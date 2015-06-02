#' @title Event-Related Potentials epoch data from Keyser and Tenke (2006).
#'   
#' @description This dataset is from one of the conditions from Keyser and Tenke
#'   (2006) paper (see reference below). It is distributed with their MATLAB
#'   package 'csd toolbox' and used in its tutorial, which can be found at 
#'   (\url{http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/tutorial.html})
#' 
#' @name erp.epoch
#' @docType data
#' @usage erp.epoch
#' @format A data frame with 3 columns and 6200 rows
#' @section Variables:
#'
#' \itemize{ 
#' \item \code{x}: x--axis coordinates for the observation.
#' \item \code{y}: y--axis coordinates for the observation.
#' \item \code{Time}: Time, in miliseconds.
#' \item \code{Electrode}: Electrode name.
#' \item \code{mV}: Grand average potential, in microvolts, from 66 healthy 
#'   adults from ``target condition for an auditory oddball task using complex 
#'   tones with a right button press.''. 
#'   (see more information at 
#'   \url{http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/tutorial.html})
#' }
#'
#' @references Kayser, J., & Tenke, C. E. (2006). Principal components analysis 
#'   of Laplacian waveforms as a generic method for identifying ERP generator 
#'   patterns: I. Evaluation with auditory oddball tasks. Clinical 
#'   Neurophysiology, 117(2), 348-368.
#'   
#' @source \url{http://psychophysiology.cpmc.columbia.edu/Software/CSDtoolbox/NR-C66-trr.dat}
#'   
"erp.epoch"