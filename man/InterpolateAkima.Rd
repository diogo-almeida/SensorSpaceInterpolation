% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/SensorSpaceInterpolation.R
\name{InterpolateAkima}
\alias{InterpolateAkima}
\title{Interpolate EEG/MEG data from sensors to a predefined sensor space
  using the functions in package AKIMA.}
\usage{
InterpolateAkima(potentials, x, y, xo = seq(min(x), max(x), length = 40),
  yo = seq(min(y), max(y), length = 40), grid.grain = NULL,
  extrapolate = FALSE, linear = FALSE)
}
\arguments{
\item{potentials}{A vector containing data to be interpolated.}

\item{x}{A matrix or a vector. If matrix, needs to contain the cartesian
  coordinates of the sensors. If vector, needs to contain the x-coordinates
  of the sensors.}

\item{y}{A vector. The y-coordinates of the sensors.}

\item{xo}{A vector defining the size of the grid in the x-dimension?
  Defaults to 40.}

\item{yo}{A vector defining the size of the grid in the y-dimension?
  Defaults to 40.}

\item{grid.grain}{A scalar. For the time being it is ignored (defaults to 40)}

\item{extrapolate}{A logical. Should the data be extrapolated beyond the
  convex hull created by the sensors? Defaults to FALSE. If TRUE, the
  interpolation will use the minimum spanning circle of the sensor space.}

\item{linear}{A logical. Should the function use the linear interpolation
  method? Defaults to FALSE, which uses cubic interpolation.}
}
\value{
List with interpolated potentials for each time slice.
}
\description{
Interpolate EEG/MEG data from sensors to a predefined sensor
  space.
}
\details{
Requires packages akima for linear and cubic interpolation.
}

