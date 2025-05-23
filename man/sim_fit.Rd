% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_fit.R
\name{sim_fit}
\alias{sim_fit}
\title{simulate animal tracks from a \code{ssm} fit}
\usage{
sim_fit(
  x,
  what = c("fitted", "predicted"),
  reps = 1,
  start = NULL,
  end = NULL,
  grad = NULL,
  beta = c(-300, -300),
  cpf = FALSE,
  sim_only = FALSE
)
}
\arguments{
\item{x}{a \code{ssm} fit object with class \code{ssm_df}}

\item{what}{simulate fitted (typically irregular in time) or predicted
(typically regular in time) locations}

\item{reps}{number of replicate tracks to simulate from an \code{ssm} model
fit object}

\item{start}{a 2-element vector for the simulated track start location
(lon,lat or x,y)}

\item{end}{a 2-element vector for the simulated track end location
(lon,lat or x,y)}

\item{grad}{a SpatRaster of x- and y-gradients as separate layers (see details)}

\item{beta}{a 2-element vector of parameters defining the potential function
magnitude in x- and y-directions (ignored if \code{is.null(grad)},
ie. no potential function; see details).}

\item{cpf}{logical; should simulated tracks return to their start point
(ie. a central-place forager)}

\item{sim_only}{logical, do not include \code{ssm} estimated location in output
(default is FALSE)}
}
\value{
a \code{fG_sim_fit} object containing the paths simulated from a
\code{ssm} fit object
}
\description{
simulate from the \code{rw} or \code{crw} process models to generate
either a set of x,y or lon,lat coordinates from a \code{ssm} fit with length
equal to the number of observations used in the SSM fit.
}
\details{
A potential function can be applied to the simulated paths to help
avoid locations on land (or in water), using the \code{grad} and \code{beta}
arguments. A coarse-resolution rasterStack of global x- and y-gradients of
distance to land are provided. Stronger beta parameters result in stronger
land (water) avoidance but may also introduce undesirable/unrealistic artefacts
(zig-zags) in the simulated paths. See Brillinger et al. (2012) and
\code{vignette("momentuHMM", package = "momentuHMM")} for more details on the
use of potential functions for simulating constrained animal movements.
WARNING: This application of potential functions to constrain simulated
paths is experimental, likely to change in future releases, and NOT guaranteed
to work enitrely as intended, especially if \code{cpf = TRUE}!
}
\examples{
fit <- fit_ssm(ellie, model = "crw", time.step = 24)
trs <- sim_fit(fit, what = "predicted", reps = 3)
plot(trs)

}
\references{
Brillinger DR, Preisler HK, Ager AA, Kie J (2012) The use of potential functions in modelling animal movement. In: Guttorp P., Brillinger D. (eds) Selected Works of David Brillinger. Selected Works in Probability and Statistics. Springer, New York. pp. 385-409.
}
