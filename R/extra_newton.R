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
  
  # ensure convergence
  stopifnot("model(s) must have converged" = all(x$converged == TRUE))
  
  # ensure pdHess
  pdHess = sapply(X = x$ssm, FUN = function(X) X$rep$pdHess)
  id_not_pdHess = names(pdHess)[pdHess != TRUE]
  if (length(id_not_pdHess) > 0) stop(paste0('Hessian is not positive definite for id: ', paste(id_not_pdHess, collapse = ',')))
  ## stopifnot("Hessian(s) must be positive definite" = all(sapply(X = x$ssm, FUN = function(X) X$rep$pdHess) == TRUE))
  
  # ensure n > 0
  stopifnot("n must be greater than 0" = (n > 0))
  
  # ensure x is a ssm fit object
  stopifnot("x must be an `ssm_df` fit object" = inherits(x, "ssm_df"))
  
  # ensure optimizer is nlminb
  stopifnot("optimizer must be nlminb" = all(sapply(X = x$ssm, FUN = function(X) X$optimiser) == "nlminb"))
  
  # ensure model is mp or crw
  stopifnot("model(s) must be mp or crw" = all(x$pmodel %in% c('mp', 'crw')))
  
  # extra Newton steps
  extra_steps <- purrr::map_dfr(.x = 1:nrow(x), .f = function(.x) { 
    
    # accounting
    x_new <- x[.x,]
    new_obj <- x_new$ssm[[1]]$tmb
    old_opt <- x_new$ssm[[1]]$opt
    old_par <- old_opt$par
    new_obj$fn(old_par) # initialize TMB
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
    if("D" %in% row.names(fxd)) {
      row.names(fxd)[which(row.names(fxd) == "D")] <- c("D_x","D_y")
    }   
    ## separate speed+se if present
    if("sf" %in% row.names(fxd)) {
      sf <- fxd[which(row.names(fxd) == "sf"), ]
      sf[1,] <- c(NA,NA)
      if("sp" %in% row.names(fxd)) {
        sp <- fxd[which(row.names(fxd) == "sp"), ]
        sp[1,] <- c(NA,NA)
      }
      fxd <- fxd[!row.names(fxd) %in% c("sf","sp"),]
    }
    
    # different process based on model
    if (x_new$pmodel == 'mp') {
      
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
      
      ## calculate AICc
      AICc <- 2 * length(old_opt[["par"]]) + 2 * old_opt[["objective"]] + 
        (2 * length(old_opt[["par"]])^2 + 2 * length(old_opt[["par"]])) / (nrow(fv) - length(old_opt[["par"]]))
      
    } else if (x_new$pmodel == 'crw') {
      
      model = 'crw'
      switch(model,
             rw = {
               tmp <- summary(rep, "random")
               rdm <- data.frame(cbind(tmp[seq(1, nrow(tmp), by = 2), ],
                                       tmp[seq(2, nrow(tmp), by = 2), ]
               ))[, c(1,3,2,4)] 
               names(rdm) <- c("x","y","x.se","y.se")
               
               rdm$id <- unique(d.all$id)
               rdm$date <- d.all$date
               rdm$isd <- d.all$isd
               
               rdm <- rdm[, c("id","date","x","y","x.se","y.se","isd")]
             },
             crw = {
               tmp <- summary(rep, "random")
               loc <- tmp[rownames(tmp) == "mu",]
               vel <- tmp[rownames(tmp) == "v",]
               
               loc <- cbind(loc[seq(1, dim(loc)[1], by = 2),],
                            loc[seq(2, dim(loc)[1], by = 2),]) 
               loc <- as.data.frame(loc, row.names = 1:nrow(loc))
               names(loc) <- c("x", "x.se", "y", "y.se")
               
               vel <- cbind(vel[seq(1, dim(vel)[1], by = 2),],
                            vel[seq(2, dim(vel)[1], by = 2),])
               vel <- as.data.frame(vel, row.names = 1:nrow(vel))
               names(vel) <- c("u", "u.se", "v", "v.se")
               
               rdm <- cbind(loc, vel)
               rdm$id <- unique(d.all$id)
               rdm$date <- d.all$date
               rdm$isd <- x_new$ssm[[1]]$isd
               
               rdm <- rdm[, c("id", "date", "x", "y", 
                              "x.se", "y.se", "u", "v", "u.se", "v.se", "isd")]
             })
      
      ## coerce x,y back to sf object
      rdm <- st_as_sf(rdm, coords = c("x","y"), remove = FALSE)
      rdm <- st_set_crs(rdm, prj)
      
      control <- list()
      control$se <- FALSE
      message('Assuming ssm_control(se = FALSE), which is the default.')
      time.step <- x_new$ssm[[1]]$ts
      switch(model,
             rw = {
               rdm <- rdm[, c("id", "date", "x.se", "y.se", "isd")]
               ## fitted
               fv <- subset(rdm, isd)[, 1:4]
               ##predicted
               if(all(!is.na(time.step))) {
                 pv <- subset(rdm, !isd)[, 1:4]
               } else {
                 pv <- NULL
               }
             },
             crw = {
               if(control$se) {
                 sf.se <- sf[,2]
                 sf <- sf[,1]
                 if(all(!is.na(time.step))) {
                   sp.se <- sp[,2]
                   sp <- sp[,1]
                 }
               } else {
                 sf <- as.vector(new_obj$report()$sf)
                 sp <- as.vector(new_obj$report()$sp)
                 sf.se <- NA
                 sp.se <- NA
               }
               rdm <- rdm[, c("id", "date", "x.se", "y.se", 
                              "u", "v", "u.se", "v.se", "isd")]
               
               ## fitted
               fv <- subset(rdm, isd)[, -9]
               fv$s <- sf
               fv$s.se <- sf.se
               
               ## predicted
               if(all(!is.na(time.step))) {
                 pv <- subset(rdm, !isd)[, -9]
                 pv$s <- sp
                 pv$s.se <- sp.se
               } else {
                 pv <- NULL
               }
             })
      
      ## calculate AICc
      npar <- length(old_opt[["par"]])
      n <- nrow(fv)
      AICc <- 2 * npar + 2 * old_opt[["objective"]] + 2 * npar * (npar + 1)/(n - npar - 1)
      
    }
    
    ## drop sv's from fxd
    fxd <- fxd[which(row.names(fxd) != "sv"), ]
    
    # replace
    if (!is.null(pv)) x_new$ssm[[1]]$predicted <- pv
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