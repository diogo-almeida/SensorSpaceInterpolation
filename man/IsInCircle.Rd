% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/SensorSpaceInterpolation.R
\name{IsInCircle}
\alias{IsInCircle}
\title{Tests whether a point is inside the minimum spanning circle of the
  EEG/MEG sensors.}
\usage{
IsInCircle(interpolation.surface, electrode.positions, inflation = 1e-07)
}
\arguments{
\item{interpolation.surface}{A matrix containig the (x, y) locations of the
  surface to be interpolated.}

\item{electrode.positions}{A matrix containing the (x, y) locations of the
  EEG/MEG sensors.}

\item{inflation}{A scalar. This represents a small fudge factor to make sure
  that the bounding circle will include all the sensors, in case there are
  sensors that would fall exactly in the circle line.}
}
\value{
This function returns TRUE or FALSE depending on whether a point is
  inside the circumcircle.
}
\description{
Tests whether a point is inside the minimum spanning circle of
  the EEG/MEG sensors.
}
\details{
Calls function circumcircle from package tripack.
}

