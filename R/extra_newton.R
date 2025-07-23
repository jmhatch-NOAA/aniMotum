##' @title extra_newton
##'
##' @details
##' `r lifecycle::badge("experimental")`
##' 
##' @description Run extra optimization on an already fitted aniMotum object 
##' 
##' @param x a \code{ssm} fit object with class `ssm_df`
##' @param n number of extra Newton optimization steps to take
##` 
##' @return an updated \code{ssm} fit object
##' 
##' @export
##' 
##' @examples
##' # run 3 extra Newton steps
##' test_extra <- extra_newton(x = fit, n = 3)
##' 
extra_newton <- function(x, n = 1) {
  
  # ensure n > 0
  stopifnot("n must be greater than 0" = (n > 0))
  
  # ensure x is a ssm fit object
  stopifnot("x must be an `ssm_df` fit object" = inherits(x, "ssm_df"))
  
  # ensure optimizer is nlminb
  stopifnot("optimizer must be nlminb" = all(sapply(X = x$ssm, FUN = function(X) X$optimiser) == "nlminb"))
  
  # extra Newton steps
  extra_steps <- purrr::map_dfr(.x = 1:nrow(x), .f = function(.x) { 
    
    # accounting
    x_new = x[.x,]
    new_obj <- x_new$ssm[[1]]$tmb
    old_par <- new_obj$par
    new_obj$fn(old_par) # initialize TMB
    old_opt <- x_new$ssm[[1]]$opt
    Gr <- new_obj$gr(old_par)
    NLL <- old_opt$objective
    
    # check gradients for NAs
    stopifnot("NAs detected in gradient" = !any(is.na(Gr)))
    
    # run Newton steps
    for (i in 1:n) {
      g <- new_obj$gr(old_opt$par) |> 
        as.numeric()
      h <- stats::optimHess(old_opt$par, fn = new_obj$fn, gr = new_obj$gr)
      new_par <- old_opt$par - solve(h, g)
      new_objective <- new_obj$fn(new_par) 
      if (any(is.na(new_par)) | is.na(new_objective)) {
        warning(paste0("Newton step ", i, " resulted in NAs, stopping at ", i - 1))
        break
      } else {
        old_opt$par <- new_par
        old_opt$objective <- new_objective        
      }
    }
  
    # warning message
    if (old_opt$objective > NLL) message("Newton steps resulted in an increase in the NLL")
    
    # save
    x_new$ssm[[1]]$tmb <- new_obj
    x_new$ssm[[1]]$opt <- old_opt
    
    # not mapped
    not_map <- rownames(x_new$ssm[[1]]$par)
    not_map <- gsub(pattern = '_x|_y', replacement = '', x = not_map)
    
    # projection
    prj <- sf::st_crs(x_new$ssm[[1]]$fitted)
    d.all <- dplyr::bind_rows(x_new$ssm[[1]]$fitted, x_new$ssm[[1]]$predicted) |>
      dplyr::arrange(date)
      
    # re-do
    rep <- TMB::sdreport(new_obj)    
    fxd <- summary(rep, "report")
    fxd <- fxd[which(rownames(fxd) %in% not_map),]

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
    rdm$isd <- x_new$ssm[[1]]$isd  
    rdm <- rdm[, c("id","date","x","y","x.se","y.se", "logit_g", "logit_g.se", "g", "isd")]
    
    ## coerce x,y back to sf object
    rdm <- sf::st_as_sf(rdm, coords = c("x","y"), remove = FALSE)
    rdm <- sf::st_set_crs(rdm, prj)
    
    ## drop x,y as we have the geom
    rdm <- rdm[, c("id","date","x.se","y.se", "logit_g", "logit_g.se", "g", "isd")]
    pv <- subset(rdm, !isd)[, -8]
    #pv$gn <- with(pv, (g - min(g))  / (max(g) - min(g)))
    fv <- subset(rdm, isd)[, -8]
    #fv$gn <- with(fv, (g - min(g))  / (max(g) - min(g)))
    
    AICc <- 2 * length(old_opt[["par"]]) + 2 * old_opt[["objective"]] + (2 * length(old_opt[["par"]])^2 + 2 * length(old_opt[["par"]])) / (nrow(fv) - length(old_opt[["par"]]))

    # replace
    x_new$ssm[[1]]$predicted <- pv
    x_new$ssm[[1]]$fitted <- fv
    x_new$ssm[[1]]$par <- fxd
    x_new$ssm[[1]]$rep <- rep
    x_new$ssm[[1]]$AICc <- AICc  
    
    # output
    x_new

  })
  
  ## output
  extra_steps
  
}