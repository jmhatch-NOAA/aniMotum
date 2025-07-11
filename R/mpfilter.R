##' @title fit the move persistence state-space model to data after passing through \code{prefilter}
##'
##' @description generates initial values for model parameters and unobserved states;
##' structures data and initial values for C++ \code{TMB} template;
##' fits state-space model; minimizes the joint log-likelihood via the selected
##' optimizer (\code{nlminb} or \code{optim}); structures and passes output
##' object to \code{fit_ssm}
##'
##' @details called by \code{fit_mpssm}, not intended for general use. \code{mpfilter} can only fit to an
##' individual track, use \code{fit_mpssm} to fit to multiple tracks (see ?fit_mpssm).
##'
##' @param x Argos data passed through prefilter()
##' @param time.step the regular time interval, in hours, to predict to.
##' Alternatively, a vector of prediction times, possibly not regular, must be
##' specified as a data.frame with id and POSIXt dates.
##' @param map a named list of parameters as factors that are to be fixed during estimation, e.g., list(psi = factor(NA))
##' @param parameters a list of initial values for all model parameters and
##' unobserved states, default is to let sfilter specify these. Only play with
##' this if you know what you are doing...
##' @param fit.to.subset fit the SSM to the data subset determined by prefilter
##' (default is TRUE)
##' @param control list of control parameters for the outer optimization (see \code{ssm_control} for details)
##' @param inner.control list of control settings for the inner optimization
##' (see ?TMB::MakeADFUN for additional details)
##' @param verbose is deprecated, use ssm_control(verbose = 1) instead, see \code{ssm_control} for details
##' @param optim is deprecated, use ssm_control(optim = "optim") instead, see \code{ssm_control} for details
##' @param optMeth is deprecated, use ssm_control(method = "L-BFGS-B") instead, see \code{ssm_control} for details
##' @param lpsi is deprecated, use ssm_control(lower = list(lpsi = -Inf)) instead, see \code{ssm_control} for details
##'
##' @importFrom TMB MakeADFun sdreport newtonOption FreeADFun
##' @importFrom stats approx cov sd predict nlminb optim na.omit filter
##' @importFrom utils flush.console
##' @importFrom tibble as_tibble
##' @importFrom sf st_crs st_coordinates st_geometry<- st_as_sf st_set_crs
##'
##'
##' @keywords internal

mpfilter <-
  function(x,
           model = model,
           time.step = 6,
           parameters = NULL,
           map = NULL,
           fit.to.subset = TRUE,
           control = ssm_control(),
           inner.control = NULL) {
    
    st <- proc.time()
    call <- match.call()
    
    ## check args
    if(!inherits(x, c("sf","tbl_df"))) stop("x must be an sf-tibble produced by `prefilter()`")
    if(!any((is.numeric(time.step) & time.step > 0) | is.na(time.step) | is.data.frame(time.step))) stop("time.step must be either: 1) a positive, non-zero value; 2) NA (to turn off predictions); or 3) a data.frame (see `?fit_ssm`)")
    if(!any(is.list(parameters) || is.null(parameters))) stop("parameters must be a named list of parameter initial values or NULL")
    if(!any(is.list(map) || is.null(map))) stop("map must be a named list of parameters to fix (turn off) in estimation or NULL")
    if(!is.logical(fit.to.subset)) stop("fit.to.subset must be TRUE (fit to prefiltered observations) or FALSE (fit to all observations)")
    if(!any(is.list(inner.control) || is.null(inner.control))) stop("inner.control must be a named list of valid newtonOptimiser control arguments or NULL")
    
    if(length(time.step) > 1 & !is.data.frame(time.step)) {
      stop("\ntime.step must be a data.frame with id's when specifying multiple prediction times")
    } else if(length(time.step) > 1 & is.data.frame(time.step)) {
      if(sum(!names(time.step) %in% c("id","date")) > 0) stop("\n time.step names must be `id` and `date`")
    }
    
    ## drop any records flagged to be ignored, if fit.to.subset is TRUE
    if(fit.to.subset) xx <- subset(x, keep)
    else xx <- x
    
    prj <- st_crs(xx)
    loc <- as.data.frame(st_coordinates(xx))
    names(loc) <- c("x","y")
    st_geometry(xx) <- NULL
    d <- cbind(xx, loc)
    d$isd <- TRUE
    
    
    if (!inherits(time.step, "data.frame") & all(!is.na(time.step))) {
      ## prediction times - assume on time.step-multiple of the hour
      tsp <- time.step * 3600
      tms <- (as.numeric(d$date) - as.numeric(d$date[1])) / tsp
      index <- floor(tms)
      
      if (time.step > 1) {
        ## time.step >= 1 h
        ## trunc so predictions start on the hour immediately prior to 1st obs
        ts <-
          data.frame(id = d$id[1],
                     date = seq(
                       trunc(d$date[1], "hour"),
                       by = tsp,
                       length.out = max(index) + 2
                     ))
      } else {
        ## time.step < 1 h
        ## trunc so predictions start on the time.step immediately prior to 1st obs
        ts1 <- trunc(d$date[1] - tsp, "mins") + tsp
        if(ts1 <= d$date[1]) {
          ts <-
            data.frame(id = d$id[1],
                       date = seq(ts1,
                                  by = tsp,
                                  length.out = max(index) + 2))
        } else {
          ts <-
            data.frame(id = d$id[1],
                       date = seq(ts1 - tsp,
                                  by = tsp,
                                  length.out = max(index) + 2))
        }
      }
      
    } else if (inherits(time.step, "data.frame") & all(!is.na(time.step))) {
      ts <- subset(time.step, id %in% unique(d$id))
      
    } else if (inherits(time.step, "data.frame") & any(is.na(time.step))) {
      stop("NA's detected in user-supplied prediction times data.frame")
    }
    
    if (all(!is.na(time.step))) {
      ## add 1 s to observation time(s) that exactly match prediction time(s)
      if (sum(d$date %in% ts$date) > 0) {
        o.times <- which(d$date %in% ts$date)
        d[o.times, "date"] <- d[o.times, "date"] + 1
      }
      
      ## merge data and prediction times
      ## add is.data flag (distinguish obs from reg states)
      d.all <- full_join(d, ts, by = c("id", "date"))
      d.all <- d.all[order(d.all$date), ]
      d.all$isd <- with(d.all, ifelse(is.na(isd), FALSE, isd))
      d.all$id <- with(d.all, ifelse(is.na(id), na.omit(unique(id))[1], id))
    } else {
      d.all <- d
    }
    
    ## calc delta times in hours for observations & interpolation points (states)
    dt <- as.numeric(difftime(d.all$date, c(as.POSIXct(NA),d.all$date[-nrow(d.all)]), 
                              units = "hours"))
    dt[1] <- 0.000001 # - 0 causes numerical issues in CRW model
    
    ## use approx & MA filter to obtain state initial values
    x.init1 <-
      approx(x = select(d, date, x),
             xout = d.all$date,
             rule = 2)$y
    x.init <-
      as.numeric(stats::filter(x.init1, rep(1, 5) / 5))
    x.na <- which(is.na(x.init))
    x.init[x.na] <- x.init1[x.na]
    
    y.init1 <-
      approx(x = select(d, date, y),
             xout = d.all$date,
             rule = 2)$y
    y.init <-
      as.numeric(stats::filter(y.init1, rep(1, 5) / 5))
    y.na <- which(is.na(y.init))
    y.init[y.na] <- y.init1[y.na]
    
    xs <- cbind(x.init, y.init)
    
    if (is.null(parameters)) {
      ## Estimate stochastic innovations
      es <- xs[-1,] - xs[-nrow(xs),]
      
      ## Estimate components of variance
      V <- cov(es)
      sigma <- sqrt(diag(V))
      rho <- V[1, 2] / prod(sigma)
      
      v <- cbind(c(0,diff(x.init)), c(0,diff(y.init))) / dt
      
      parameters <- list(
        l_sigma = log(pmax(1e-08, sigma)),
        l_rho_p = log((1 + rho) / (1 - rho)),
        X = t(xs),
        lg = rep(0, dim(xs)[1]),
        l_sigma_g = 0,
        l_psi = 0,
        l_tau = c(0, 0),
        l_rho_o = 0
      )
    } 
    
    ## start to work out which obs_mod to use for each observation
    d <- mutate(d, obs.type = factor(obs.type, 
                                        levels = c("LS","KF","GL","GPS"), 
                                        labels = c("LS","KF","GL","GPS"))
                      )
    p.obst <- table(d$obs.type) / nrow(d)
    obst <- which(table(d$obs.type) > 0)
    # favours KF when mix of few LS + many KF
    obst <- round(which(table(d$obs.type) * p.obst > 0))
    
    automap <- list(
      list(
        l_psi = factor(NA)
        ),
      list(
        l_tau = factor(c(NA, NA)),
        l_rho_o = factor(NA)
        ),
      list(
        l_tau = factor(c(NA, NA)),
        l_psi = factor(NA)
        ),
      list(
        l_psi = factor(NA)
        )
      )

    mm <- unlist(automap[obst], recursive = FALSE)
    mm[which(duplicated(names(mm)))] <- NULL
    automap <- mm
    
    if(!is.null(map)) {
      names(map) <- paste0("l_", names(map))
      map <- append(automap, map, after = 0)
    } else {
      ## ensure psi is always turned on if KF obs present in data
      ## handles cases when % LS > % KF (may be rare or non-existent)
      if(all(length(table(d.all$obs.type)) > 1, 
             "KF" %in% names(table(d.all$obs.type)))) {
        map <- automap[-which(names(automap)  == "l_psi")]
      } else {
        map <- automap
      }
    }    
    
    ## TMB - data list
    ## convert from string to integer for C++
    obs_mod <- ifelse(d.all$obs.type %in% c("LS","GPS"), 0, 
                      ifelse(d.all$obs.type == "KF", 1, 2)
                      )
    ## where is.na(obs_mod) - prediction points - set to "LS" (obs_mod = 0) so
    ##  NA's don't create an int overflow situation in C++ code. This won't matter 
    ##  as isd makes likelihood contribution goes to 0 in C++ code
    obs_mod <- ifelse(is.na(obs_mod), 0, obs_mod)
    
    data <- list(
      model_name = model,
      Y = rbind(d.all$x, d.all$y), 
      dt = dt,
      isd = as.integer(d.all$isd),
      obs_mod = as.integer(obs_mod),
      m = d.all$smin,
      M = d.all$smaj,
      c = d.all$eor,
      K = cbind(d.all$emf.x, d.all$emf.y),
      GLerr = cbind(d.all$x.sd, d.all$y.sd)
    ) 
    
    ## TMB - create objective function
    if (is.null(inner.control) | !"smartsearch" %in% names(inner.control)) {
      inner.control <- list(smartsearch = TRUE)
    }
    rnd <- c("X", "lg")
    
    obj <-
      MakeADFun(
        data,
        parameters,
        map = map,
        random = rnd,
        hessian = TRUE,
        method = control$method,
        DLL = "aniMotum",
        silent = !ifelse(control$verbose == 2, TRUE, FALSE),
        inner.control = inner.control
      )    
    
    obj$env$tracemgc <- ifelse(control$verbose == 2, TRUE, FALSE)
    
    ## add par values to trace if control$verbose = TRUE
    myfn <- function(x) {
      cat("\r", "pars:  ", round(x, 5), "     ")
      flush.console()
      obj$fn(x)
    }
    
    ## Set parameter bounds - most are -Inf, Inf
    L = c(l_sigma=c(-Inf,-Inf),
          l_rho_p=-7,
          l_sigma_g=-10,
          l_psi=-Inf,
          l_tau=c(-Inf,-Inf),
          l_rho_o=-7) ## using 2 / (1 + exp(-x)) - 1 transformation, this gives rho_o = -0.999, 0.999
    U = c(l_sigma=c(Inf,Inf),
          l_rho_p=7,
          l_sigma_g=50,
          l_psi=Inf,
          l_tau=c(Inf,Inf),
          l_rho_o=7)
    
    if(any(!is.null(control$lower))) {
      L[which(names(L) %in% names(control$lower))] <- unlist(control$lower)
    }
    if(any(!is.null(control$upper))) {
      U[which(names(U) %in% names(control$upper))] <- unlist(control$upper)
    } 
    
    names(L)[c(1:2,6:7)] <- c("l_sigma", "l_sigma", "l_tau", "l_tau")
    names(U)[c(1:2,6:7)] <- c("l_sigma", "l_sigma", "l_tau", "l_tau")
    
    # Remove inactive parameters from bounds
    L <- L[!names(L) %in% names(map)]
    U <- U[!names(U) %in% names(map)]

    ## Minimize objective function
    oldw <- getOption("warn")
    options(warn = -1)  ## turn warnings off but check if optimizer crashed & return warning at end
    
    opt <-
      switch(control$optim,
             nlminb = try(nlminb(obj$par,
                                 ifelse(control$verbose == 1, myfn, obj$fn),
                                 obj$gr,
                                 control = control$control,
                                 lower = L,
                                 upper = U
             ), silent = TRUE),
             optim = try(do.call(
               optim,
               args = list(
                 par = obj$par,
                 fn = ifelse(control$verbose == 1, myfn, obj$fn),
                 gr = obj$gr,
                 method = control$method,
                 control = control$control,
                 lower = L,
                 upper = U
               )
             ), silent = TRUE))
    cat("\n")
    
    FreeADFun(obj)
    
    ## if error then exit with limited output to aid debugging
    ## check if pdHess is FALSE at end and return warning
    rep <- try(sdreport(obj), silent = TRUE) #, skip.delta.method = !control$se
    
    options(warn = oldw) ## turn warnings back on
  
    if (!inherits(opt, "try-error") & !inherits(rep, "try-error")) {
      ## Parameters, states and the fitted values
      fxd <- summary(rep, "report")
      if(!is.null(map)) { 
        fxd <- fxd[which(!rownames(fxd) %in%
                           sapply(strsplit(names(map), "_"), function(.)
                             .[2])),]
      }
      if("sigma" %in% row.names(fxd)) {
        row.names(fxd)[which(row.names(fxd) == "sigma")] <- c("sigma_x","sigma_y")
      } 
      if("tau" %in% row.names(fxd)) {
        row.names(fxd)[which(row.names(fxd) == "tau")] <- c("tau_x","tau_y")
      }
      
      tmp <- summary(rep, "random")
      X <- tmp[rownames(tmp) %in% "X",]
      lg <- tmp[rownames(tmp) %in% "lg",]
      
      rdm <- data.frame(cbind(X[seq(1, nrow(X), by = 2), ],
                              X[seq(2, nrow(X), by = 2), ]
      ))[, c(1,3,2,4)] 
      names(rdm) <- c("x","y","x.se","y.se")
      rdm <- data.frame(rdm, logit_g = lg[,1], logit_g.se = lg[,2], g = plogis(lg[,1]))
      
      
      rdm$id <- unique(d.all$id)
      rdm$date <- d.all$date
      rdm$isd <- d.all$isd  
      rdm <- rdm[, c("id","date","x","y","x.se","y.se", "logit_g", "logit_g.se", "g", "isd")]
      
      ## coerce x,y back to sf object
      rdm <- st_as_sf(rdm, coords = c("x","y"), remove = FALSE)
      rdm <- st_set_crs(rdm, prj)
      
      # drop x,y as we have the geom
      rdm <- rdm[, c("id","date","x.se","y.se", "logit_g", "logit_g.se", "g", "isd")]
      pv <- subset(rdm, !isd)[, -8]
      #pv$gn <- with(pv, (g - min(g))  / (max(g) - min(g)))
      fv <- subset(rdm, isd)[, -8]
      #fv$gn <- with(fv, (g - min(g))  / (max(g) - min(g)))
    
      if (control$optim == "nlminb") {
        AICc <- 2 * length(opt[["par"]]) + 2 * opt[["objective"]] + 
          (2 * length(opt[["par"]])^2 + 2 * length(opt[["par"]])) / (nrow(fv) - length(opt[["par"]]))
      } else if (control$optim == "optim") {
        AICc <- 2 * length(opt[["par"]]) + 2 * opt[["value"]] + 
          (2 * length(opt[["par"]])^2 + 2 * length(opt[["par"]])) / (nrow(fv) - length(opt[["par"]]))
      }
      
      out <- list(
        call = call,
        predicted = pv,
        fitted = fv,
        par = fxd,
        data = x,
        isd = d.all$isd,
        inits = parameters,
        pm = "mp",
        ts = time.step,
        opt = opt,
        tmb = obj,
        rep = rep,
        AICc = AICc,
        optimiser = control$optim,
        time = proc.time() - st
      )
      
      if(!rep$pdHess & any(!is.na(x$smaj), !is.na(x$smin), !is.na(x$eor)) & (!"l_psi" %in% names(map))) {
      warning("Hessian was not positive-definite so some standard errors could not be calculated. Try simplifying the model by adding the following argument:
             `map = list(psi = factor(NA))`", call. = FALSE)
      } else if (!rep$pdHess) {
        warning("Hessian was not positive-definite so some standard errors could not be calculated.", call. = FALSE)
      }
    } else if (rep$pdHess & all(is.na(x$smaj), is.na(x$smin), is.na(x$eor))){
      warning("Hessian was not positive-definite so some standard errors could not be calculated.", 
              call. = FALSE)
    } else {
      
      ## if optimiser fails
      out <- list(
        call = call,
        data = x,
        inits = parameters,
        pm = model,
        ts = time.step,
        tmb = obj,
        errmsg = opt
      )
      
      if(any(!is.na(x$smaj), !is.na(x$smin), !is.na(x$eor)) & (!"l_psi" %in% names(map))) {
        warning("The optimiser failed. Try simplifying the model with the following argument: 
                                 `map = list(psi = factor(NA))`", call. = FALSE)
      } else if (all(is.na(x$smaj), is.na(x$smin), is.na(x$eor)) & (!"rho_o" %in% names(map))){
        warning("The optimiser failed. Try simplifying the model with the following argument: 
                                 `map = list(rho_o = factor(NA))`", call. = FALSE)
      } else {
        warning("The optimiser failed. You could try using a different time.step or turn off prefiltering with
                `spdf = FALSE`", call. = FALSE)
      }
    }
    
    class(out) <- append("mp_ssm", class(out))
    
    out
}