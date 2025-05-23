##' @title Coerce input data into expected `aniMotum` format
##'
##' @description format data by mapping supplied variable names to those expected by
##' `fit_ssm()`, and ensuring variables are put into the expected order. Can be 
##' run manually by user as a data pre-processing step prior to calling `fit_ssm()`
##' or can be called automatically by `fit_ssm()`. In the latter case, any custom 
##' variable names must be declared as arguments to `fit_ssm()`; see examples, below.
##'
##' @param x input data
##' @param id the name (as a quoted character string) of id variable: a unique 
##' identifier for individual (animal) track data sets.
##' @param date the name (as a quoted character string)of the date/time variable:
##' date and time (as YYYY-MM-DD HH:MM:SS) of each observation.
##' @param lc the name (as a quoted character string) of the location quality class
##' variable: Argos location quality class (values in the set: 3,2,1,0,"A","B","Z").
##' Can also include "G" for GPS data and/or "GL" for light-level geolocation (GLS)
##' and other data types.
##' @param coord the names (as quoted character strings) of the location coordinate
##' variables: defaults are c("lon","lat"), but could also be c("x","y") for planar
##'  coordinates; or if input data is an `sf` object then "geometry". If input
##'  data is an `sf` object then `coord` is set to "geometry" by default.
##' @param epar the names (as quoted character strings) of the Argos error ellipse
##' parameters: defaults are "smaj" (ellipse semi-major axis), 
##' "smin" (ellipse semi-minor axis), and "eor" (ellipse orientation). Ignored if
##' these variables are missing from the input data.
##' @param sderr the names (as quoted character strings) of provided standard 
##' errors for `lon,lat` or `x,y`: default names are `x.sd`, `y.sd`. Typically,
##' these are only provided for generic location data such as processed light-level
##' geolocations, or high-resolution acoustic detections. The argument is ignored
##' if these variables are missing from the input data.
##' @param tz the timezone the applies to the data/time variable if they are not 
##' in `tz = 'UTC'`. A list of valid timezone names can be viewed via `OlsonNames()`
##' 
##' @return a data.frame or sf-tibble of input data in expected aniMotum format. 
##' Additional columns required by `fit_ssm()`, if missing, will be added to the 
##' formatted tibble: `smaj`, `smin`, `eor`, `x.sd`, and `y.sd`.
##' 
##' @importFrom sf st_crs
##' @importFrom dplyr tibble select everything
##'
##' @examples
##' ## as a data pre-processing step
##' data(sese2_n)
##' head(sese2_n, 5)
##' d <- format_data(sese2_n, date = "time", coord = c("longitude","latitude"), 
##' tz = "America/Halifax")
##' fit <- fit_ssm(d, model = "crw", time.step = 24)
##' 
##' ## called automatically within fit_ssm()
##' fit <- fit_ssm(sese2_n, date = "time", coord = c("longitude", "latitude"), 
##' tz = "America/Halifax", model = "crw", time.step = 24)
##' @export
##' @md

format_data <- function(x,
                        id = "id",
                        date = "date",
                        lc = "lc",
                        coord = c("lon","lat"),
                        epar = c("smaj","smin","eor"),
                        sderr = c("x.sd","y.sd"),
                        tz = "UTC") {
  
  ## check that all variable names are character strings
  if(id %in% names(x)) 
    stopifnot("id must be a character string" = is.character(id))
  else {
    stop("An 'id' variable must be included in the input data\n")
  }
  stopifnot("date must be a character string" = is.character(date))
  stopifnot("lc must be a character string" = is.character(lc))
  stopifnot("coord must be a character vector with 1 or 2 elements" = all(is.character(coord)))
  stopifnot("epar must be a character vector with 3 elements" = all(is.character(epar)))
  stopifnot("sderr must be a character vector with 2 elements" = all(is.character(sderr)))
  
  ## if "lonerr", "laterr" exist in x then rename "x.sd", "y.sd"
  ## return warning that "lonerr", "laterr" are deprecated
  if(all(all(c("lonerr","laterr") %in% names(x)), 
         all(c("x.sd","y.sd") %in% sderr))) {
    warning("The error variable names `lonerr`, `laterr` are deprecated as of v 1.2 and should be changed to `x.sd`, `y.sd` to avoid this warning.
    See `?fit_ssm` for more details on `x.sd`, `y.sd` values.
    Proceeding...",
            call. = FALSE,
            immediate. = TRUE)
    err.idx <- which(names(x) %in% c("lonerr", "laterr"))
    names(x)[err.idx] <- sderr
  }
  
  ## if input is an sf data.frame (ie. a data.frame with a geometry list object) then
  ##  coerce to an sf-tibble
  if(all(!inherits(x, "sf"), "geometry" %in% names(x))) {
    if(inherits(x$geometry, "sfc")) {
      x <- st_as_sf(x)
    }
  }
  
  ## set coord = "geometry" if input data is an "sf" object
  if(inherits(x, "sf")) {
    coord <- "geometry" 
  }
  
  ## check that specified mandatory variable names are in the input data
  stopifnot("An id variable must be included in the input data; 
            see vignette('Overview', package = 'aniMotum')" = id %in% names(x))
  stopifnot("A date/time variable must be included in the input data; 
            see vignette('Overview', package = 'aniMotum')" = date %in% names(x))
  stopifnot("Coordinate variables must be included in the input data; 
            see vignette('Overview', package = 'aniMotum')" = all(coord %in% names(x)))
 
  ## if input data is an 'sf' object then check that a valid crs exists
  if(inherits(x, "sf") & is.na(st_crs(x))) {
    stop("\nCRS info is missing from input data sf object")
  }
  
  ## add lc if missing from input data
  if (!lc %in% names(x)) {
    ## Case when data are Generic Locations
    if (all(!epar %in% names(x)) & all(sderr %in% names(x))) {
      if (inherits(x, "data.frame", which = TRUE) == 1) {
        x <- data.frame(x, lc = rep("GL", nrow(x)))
      } else if (inherits(x, "tbl_df", which = TRUE) == 1) {
        x <- tibble(x, lc = "GL")
      } else if(inherits(x, "sf", which = TRUE) == 1) {
        x$lc <- rep("GL", nrow(x))
      }
      x <- x[, c(id, date, "lc", coord, sderr)]
      
      ## Case when data are GPS
    } else if (all(!epar %in% names(x)) &
               all(!sderr %in% names(x))) {
      message("Guessing that all observations are GPS locations.")
      if (inherits(x, "data.frame", which = TRUE) == 1) {
        x <- data.frame(x, lc = rep("G", nrow(x)))
      } else if (inherits(x, "tbl_df", which = TRUE) == 1) {
        x <- tibble(x, lc = rep("G", nrow(x)))
      } else if(inherits(x, "sf", which = TRUE) == 1) {
        x$lc <- rep("G", nrow(x))
      } 
      x <- x[, c(id, date, "lc", coord)]
    }
  }
  
  ## In case all format_data-added variables already exist in x
  xx <- x
  
  ## determine if there are extra variables in x
  xt.vars <- names(x)[!names(x) %in% c(id, date, lc, coord, epar, sderr)]

  if(all(!c("lon","lat") %in% coord, coord != "geometry")) {
    pos1 <- grepl("lon", coord, ignore.case = TRUE)
    pos2 <- grepl("lat", coord, ignore.case = TRUE)
    if(!any(pos1)) {
      pos1 <- grepl("x", coord, ignore.case = TRUE)
      pos2 <- grepl("y", coord, ignore.case = TRUE)
    }
    coord <- coord[c(which(pos1), which(pos2))]
  }
  
  if(all(!epar %in% names(x), !sderr %in% names(x))) {
    ## Argos LS or GPS data
    ## add expected error ellipse variables
    x$smaj <- x$smin <- x$eor <- as.double(NA)
    ## add expected sderr variables
    x$x.sd <- x$y.sd <- as.double(NA)
    xx <- x[, c(id, date, lc, coord, epar, sderr, xt.vars)]
    
    if(all(!inherits(x, "sf"), all(coord %in% c("lon","lat")))) {
      names(xx)[1:5] <- c("id","date","lc",coord)
      names(xx)[4:5] <- c("lon","lat")
    } else if(all(!inherits(x, "sf"), any(!coord %in% c("lon","lat")))) {
      names(xx)[1:5] <- c("id","date","lc","lon","lat")
    } else if(inherits(x, "sf")) {
      names(xx)[1:4] <- c("id","date","lc",coord)
    }
  } 
  if(all(epar %in% names(x), !sderr %in% names(x))) {
    ## Argos KF/KS data
    ## add expected sderr variables
    x$x.sd <- x$y.sd <- as.double(NA)
    xx <- x[, c(id, date, lc, coord, epar, sderr, xt.vars)]
    if(all(!inherits(x, "sf"), coord != "geometry")) {
      names(xx)[1:8] <- c("id","date","lc",coord,"smaj","smin","eor")
      names(xx)[4:5] <- c("lon","lat")
    } else if (inherits(x, "sf")) {
      names(xx)[1:7] <- c("id","date","lc",coord,"smaj","smin","eor")
    }
  }
  if(all(!epar %in% names(x), sderr %in% names(x))) {
    ## Generic Location data
    ## add expected error ellipse variables
    x$smaj <- x$smin <- x$eor <- as.double(NA)
    xx <- x[, c(id, date, lc, coord, epar, sderr, xt.vars)]
    if(!inherits(x, "sf")) {
      names(xx)[c(1:5, 9:10)] <- c("id","date","lc",coord,"x.sd","y.sd")
      #names(xx)[4:5] <- c("lon","lat")
    } else if (inherits(x, "sf")) {
      names(xx)[c(1:4, 8:9)] <- c("id","date","lc",coord,"x.sd","y.sd")
    }
  }

  ## in cases where user supplies id as a factor, drop any unused factor levels 
  ##    and coerce to character
  if(is.factor(xx$id)) xx$id <- droplevels(xx$id)
  xx$id <- as.character(xx$id)
  
  ## convert dates to POSIXt if not already
  if(!inherits(xx$date, "POSIXt")) {
    xx$date <- try(as.POSIXct(xx$date, tz = tz), silent = TRUE)
    if(inherits(xx$date, "try-error")) stop("dates must be in a standard format: YYYY-MM-DD HH:MM:SS")
  }
  ## order records by date
  xx <- xx[order(xx$date), ]
  
  class(xx) <- append("fG_format", class(xx))
  
  return(xx)
  
}
