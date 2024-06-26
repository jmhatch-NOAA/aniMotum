##' @title Fit Continuous-Time State-Space Models to filter Argos satellite 
##' geolocation data
##'
##' @description fits: i) a simple random walk (`rw`) ii) a correlated random walk
##' (`crw` - a random walk on velocity), or iii) a time-varying move persistence 
##' model (`mp`), all in continuous-time, to filter Argos LS, and/or KF/KS 
##' location data, GPS data, and/or generic locations with associated standard 
##' errors (e.g., processed light-level geolocation data, or high-resolution 
##' acoustic telemetry data). Location data of different types can combined in a 
##' single data frame (see details). Predicts locations at user-specified time 
##' intervals (regular or irregular).
##'
##' @param x a `data.frame`, `tibble` or `sf-tibble` of observations, depending 
##' on the tracking data type. See more in the Details section, below, and the 
##' Overview vignette \code{vignette("Overview", package = "aniMotum")}.
##' @param vmax max travel rate (m/s) to identify implausible locations
##' @param ang angles (deg) of implausible location "spikes" 
##' @param distlim lengths (m) of implausible location "spikes" 
##' @param spdf (logical) turn pre-filtering on (default; TRUE) or off
##' @param min.dt minimum allowable time difference between observations;
##' `dt <= min.dt` will be ignored by the SSM. Default is 0: all time
##' differences > 0 are allowed.
##' @param pf just pre-filter the data, do not fit the SSM (default is FALSE)
##' @param model fit a simple random walk (`rw`), correlated random walk
##' (`crw`), or a time-varying move persistence model (`mp`), all as 
##' continuous-time process models
##' @param time.step options: 1) the regular time interval, in hours, to predict to; 
##' 2) a vector of prediction times, possibly not regular, must be
##' specified as a data.frame with id and POSIXt dates; 3) NA - turns off 
##' prediction and locations are only estimated at observation times. 
##' @param emf optionally supplied data.frame of error multiplication factors for 
##' Argos location quality classes. Default behaviour is to use the factors 
##' supplied by [aniMotum::emf]
##' @param map a named list of parameters as factors that are to be fixed during
##' estimation, e.g., `list(psi = factor(NA))`
##' @param parameters a list of initial values for all model parameters and
##' unobserved states, default is to let sfilter specify these. Only play with
##' this if you know what you are doing
##' @param fit.to.subset fit the SSM to the data subset determined by 
##' [aniMotum::prefilter] (default is TRUE)
##' @param control list of control settings for the outer optimizer 
##' (see [aniMotum::ssm_control] for details)
##' @param inner.control list of control settings for the inner optimizer 
##' (see [TMB::MakeADFun] for additional details)
##' @param ... variable name arguments passed to format_data, see 
##' [aniMotum::format_data] for details 
##'
##' @details `x` is a `data.frame`, `tibble`, or `sf-tibble` with 5, 7 or 8 
##' columns (the default format), depending on the tracking data type. Argos 
##' Least-Squares and GPS data should have 5 columns in the following order: 
##' **`id`, `date`, `lc`, `lon`, `lat`**. Where `date` can be a POSIX object or text
##' string in YYYY-MM-DD HH:MM:SS format. If a text string is supplied then the
##' time zone is assumed to be `UTC`. lc (location class) can include the 
##' following values: 3, 2, 1, 0, A, B, Z, G, or GL. The latter two are for GPS 
##' locations and 'Generic Locations', respectively. Class Z values are assumed 
##' to have the same error variances as class B. By default, class `G` (GPS) 
##' locations are assumed to have error variances 10x smaller than Argos class 3 
##' variances, but unlike Argos error variances the GPS variances are the same for 
##' longitude and latitude. 
##' 
##' The [aniMotum::format_data] function can be used as a data pre-processing 
##' step or called automatically within `fit_ssm` to restructure data that is 
##' not in one of the above default formats. The minimum essential variables: 
##' **`id`, `date`, `lc`, `lon`, `lat`** must exist in the input data but they can
##' have different names and exist in a different column order. See 
##' [aniMotum::format_data] for details.
##' 
##' See [emf] for details on how to modify these assumptions. 
##' 
##' Argos Kalman Filter (or Kalman Smoother) data should have 8 columns, 
##' including the above 5 plus **`smaj`, `smin`, `eor`** that contain Argos error 
##' ellipse variables (in m for `smaj`, `smin` and deg for `eor`). 
##' 
##' Generic locations can be modelled provided each longitude and latitude 
##' (or X and Y) coordinate has a corresponding standard error. These data should 
##' have 7 columns, including the above 5 plus two extra columns, typically
##' named `x.sd`, `y.sd` that provide the standard errors for the longitude, 
##' latitude (or X, Y) coordinates. Longitude and latitude standard errors should
##' be in degrees, whereas X and Y standard errors should be in m. In either case, 
##' all `lc` values should be set to `GL` (Generic Location), the helper function
##' [format_data] will add the `lc` variable to the input data automatically. 
##' 
##' Multiple location data types can be combined in a single data frame 
##' (see the Overview vignette for examples). 
##' 
##' When data are provided as an `sf-tibble`, the user-specified projection is 
##' respected, although projected units are always transformed to km to improve 
##' SSM convergence efficiency. Otherwise, longlat data are re-projected 
##' internally to a global Mercator grid and provided as the default output. 
##' A simple `tibble`, without a geom, of `lon,lat` and `x,y` location estimates
##'  can be obtained by using [grab] with the argument `as_sf = FALSE`.
##' 
##' @return a list with components
##' * `call` the matched call
##' * `predicted` an sf tbl of predicted location states
##' * `fitted` an sf tbl of fitted locations
##' * `par` model parameter summary
##' * `data` an augmented sf tbl of the formatted input data
##' * `inits` a list of initial values
##' * `pm` the process model fit, either "rw" or "crw"
##' * `ts` time time.step in h used
##' * `opt` the object returned by the optimizer
##' * `tmb` the TMB object
##' * `rep` TMB sdreport
##' * `aic` the calculated Akaike Information Criterion
##' * `time` the processing time for sfilter
##'
##' @references
##' Jonsen ID, Patterson TA, Costa DP, et al. (2020) A continuous-time state-space
##'  model for rapid quality-control of Argos locations from animal-borne tags. 
##'  Movement Ecology 8:31 
##' 
##' Jonsen ID, McMahon CR, Patterson TA, et al. (2019) Movement responses to 
##' environment: fast inference of variation among southern elephant seals with 
##' a mixed effects model. Ecology. 100(1):e02566 
##' 
##' @examples
##' ## fit crw model to Argos LS data
##' fit <- fit_ssm(ellie, vmax = 4, model = "crw", time.step = 24, 
##' control = ssm_control(verbose = 0)) 
##' 
##' ## time series plots of fitted values and observations
##' plot(fit, what = "fitted", type = 1, ask = FALSE)
##'
##' ## 2-D tracks plots of predicted values and observations
##' plot(fit, what = "predicted", type = 2, ask = FALSE)
##'
##'
##' @export
##' @md

fit_ssm <- function(x,
                    vmax = 5,
                    ang = c(15,25),
                    distlim = c(2500,5000),
                    spdf = TRUE,
                    min.dt = 0,
                    pf = FALSE,
                    model = "crw",
                    time.step = NA,
                    emf = NULL,
                    map = NULL,
                    parameters = NULL,
                    fit.to.subset = TRUE,
                    control = ssm_control(),
                    inner.control = NULL,
                    ...
                    )
{

  dots <- list(...)
  
  stopifnot("model can only be 1 of `rw`, `crw`, or `mp`" = model %in% c("rw","crw","mp"))
  
## check args - most args handled by prefilter() & sfilter()
  if(!is.data.frame(x))  
    stop("x must be a data.frame, tibble or sf-tibble, see `?fit_ssm for details`")
  if(!is.logical(pf)) 
    stop("pf must be either FALSE (fit model) or TRUE (only run prefilter)")
  
## warnings for deprecated arguments
  if(all("verbose" %in% names(dots), !is.null(dots$verbose))) {
    warning("the `verbose` arg is deprecated & will be removed in the next version, use `control = ssm_control(verbose)` instead. See `?ssm_control for details",
            call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    control$verbose <- dots$verbose
  }
  if(all("optim" %in% names(dots), !is.null(dots$optim))) {
    warning("the `optim` arg is deprecated & will be removed in the next version, use `control = ssm_control(optim)` instead. See `?ssm_control for details",
            call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    if(dots$optim %in% c("nlminb", "optim")) control$optim <- dots$optim
    else stop("invalid optimiser specified, see ?ssm_control for options")
  }
  if(all("optMeth" %in% names(dots), !is.null(dots$optMeth))) {
    warning("the `optMeth` arg is deprecated & will be removed in the next version, use `control = ssm_control(method)` instead. See `?ssm_control for details",
            call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    if(dots$optMeth %in% c("L-BFGS-B", "BFGS", "Nelder-Mead", "CG", "SANN", "Brent"))
      control$method <- dots$optMeth
    else stop("invalid optimisation method specified, see ?ssm_control for options")
  }
  if(all("lpsi" %in% names(dots), !is.null(dots$lpsi))) {
    warning("the `lpsi` arg is deprecated & will be removed in the next version, use `control = ssm_control(lower)` instead. See `?ssm_control for details",
            call. = FALSE, immediate. = TRUE, noBreaks. = TRUE)
    control$lower <- list(l_psi = dots$lpsi)
  }
  
  ## ensure data is in expected format
  if(!inherits(x, "fG_format")) x <- format_data(x, ...) 
  
  ## apply prefilter
  fit <- lapply(split(x, x$id),
                function(xx) {
                  prefilter(x = xx,
                            vmax = vmax,
                            ang = ang,
                            distlim = distlim,
                            spdf = spdf,
                            min.dt = min.dt,
                            emf = emf)
                })
  
  ## if pf = TRUE then just prefilter & return result
  if(pf){
    fit <- try(do.call(rbind, fit))
    
    ## this shouldn't happen as invalid data should error at format_data() step
    if(inherits(fit, "try-error")) 
    stop("something went wrong trying to put prefilter-ed tracks back together")
    
  } else {
    if(control$verbose == 1)
      cat(paste0("fitting ", model, " SSM to ", length(fit), " tracks...\n"))
    if(model %in% c("crw", "rw")) {
      fit <- lapply(fit,
                    function(x) {
                      sfilter(
                        x = x,
                        model = model,
                        time.step = time.step,
                        parameters = parameters,
                        map = map,
                        fit.to.subset = fit.to.subset,
                        control = control,
                        inner.control = inner.control
                      )
                    })
      
    } else if (model == "mp") {
      fit <- lapply(fit,
                    function(x) {
                      mpfilter(
                        x = x,
                        model = model,
                        time.step = time.step,
                        parameters = parameters,
                        map = map,
                        fit.to.subset = fit.to.subset,
                        control = control,
                        inner.control = inner.control
                      )
                    })
      
    }
    
    fit <- tibble(id = names(fit), ssm = fit)
    fit <- tibble(fit,
                  converged = sapply(fit$ssm, function(x) 
                    if(length(x) == 15) {
                      x$opt$convergence == 0
                    } else if(length(x) < 15) {
                      FALSE
                    }),
                  pdHess = sapply(fit$ssm, function(x) 
                    length(x) == 15
                  ),
                  pmodel = sapply(fit$ssm, function(x) x$pm)
    )

  }

  class(fit) <- append("ssm_df", class(fit))
  return(fit)
}
