#' @title Formats the x and y locations of EEG/MEG sensors for interpolation 
#'   function.
#' 
#' @description Formats the x and y locations of EEG/MEG sensors for 
#'   interpolation function.
#' 
#' @param x A matrix or a vector. If x is a matrix, it needs to have two columns
#'   containing the (x, y) locations of the sensors, in which case y is not 
#'   required. If it is a vector, it should contain the x-coordinates of the 
#'   sensor.
#'   
#' @param y A vector containing the y-locations of the sensors. Defaults to 
#'   NULL, but will throw an error if x does not contain all the required
#'   information.
#'      
#' @return This function returns a data.frame with two fields (x, y), properly 
#'   formatted for the interpolation function.
#'   
#' @export
FormatPosition <- function(x, y = NULL) {
  #   desired input - output relationships:
  #   first argument 'x' is a [x, y] matrix
  #     -> pass x as is, y = NULL
  #   first argument 'x' is a vector of x-locations
  #     -> pass x as is, y needs to be a vector of same type and length
  #   first argument 'x' is a data.frame of [x, y] form
  #     -> coerce x to matrix, y = NULL
  #   first argument 'x' is a list of two fields {x, y}, each a vector of equal 
  #     numeric type and length
  #     -> x1 <- list$x, y1 <- list$y
  
  x1 <- NULL
  y1 <- NULL
  
  # Error checks
  if (is.null(x)) {
    stop("argument 'x' is NULL or missing!")
  }
  if (is.list(x)) {
    if (!is.null(y)) {
      warning("argument 'x' is a list or data.frame. Will ignore argument 'y'!")
    } else {
      if (length(x) < 2) {
        stop("argument 'y' is NULL, but it seems to be needed, as\n",
             "argument 'x' has only one field/column of data"
        )
      }
    }
    if (!all(c("x", "y") %in% names(x))) {
      # No fields named 'x' and 'y' on the list/data.frame
      # attempts to use first two fields
      x1 <- unlist(x[1])
      y1 <- unlist(x[2])
    } else {
      x1 <- x$x
      y1 <- x$y
    }
    if (any(!is.numeric(x1), !is.numeric(y1))) {
      stop("Problem with list 'x': \n",
           "a) It may only contain one field, or\n",
           "b) One (or both) first two columns are not numeric vectors"
      )
    }
  } else {
    # not a list, maybe a matrix?
    if (is.matrix(x)){
      if (!is.null(y)){
        warning("argument 'x' is a matrix. Will ignore argument 'y'!")
      }
      x1 <- x[, 1]
      y1 <- x[, 2]
    } else {
      # not a list and not a matrix...
      # x is a vector
      # y is necessary
      if (is.null(y)) {
        stop("argument 'y' is NULL, but it seems to be needed!")
      }
      if (any(!is.numeric(x), !is.numeric(y))) {
        stop("arguments 'x' and 'y' need to be numeric!")
      }
      if (is.matrix(y)) {
        warning("arguments 'y' is a matrix. Will use only its first column!")
        x1 <- x
        y1 <- y[, 1]
      } else {
        x1 <- x
        y1 <- y 
      }
    }
  }
  if (length(x1) != length(y1)) {
    stop("Lengh of vectors containing xy-coordinates differ!")
  }
  cbind(x = x1, y = y1)
}

#' @title Finds the bounding rectangle for the sensors in cartesian space.
#' 
#' @description Finds the bounding rectangle for the sensors in cartesian space.
#' 
#' @param sensor.space A matrix. This matrix should contain the x,y locations of
#'   the sensors.
#'   
#' @param inflation A scalar. This represents a small fudge factor to make sure
#'   that the bounding rectangle will include all the sensors, in case there are
#'   sensors that would fall exactly in the rectangle line. Defaults to 5% of 
#'   the rectangle dimensions.
#'      
#' @return This function returns a data.frame with the cartesian coordinates of
#'   the bounding rectangle for the interpolation function.
#'   
#' @export
FindBoundingRectangle <- function(sensor.space, inflation = 1.05) {
  if (!is.matrix(sensor.space)) {
    sensors <- FormatPosition(sensor.space)
  } else {
    sensors <- sensor.space
  }
  maxval.xy <- apply(apply(sensors, 2, range), 2, function(x){max(abs(x))}) * inflation
  x.range   <- c(-1 * maxval.xy[1], maxval.xy[1])
  y.range   <- c(-1 * maxval.xy[2], maxval.xy[2])
  rectangle <- rbind(x.range, y.range)
  colnames(rectangle) <- c("min", "max")
  rectangle
}

#' @title Creates a grid in cartesian space for interpolation function.
#' 
#' @description Creates a grid in cartesian space for interpolation function.
#'
#' @details This function takes a matrix of cartesian coordinates of EEG/MEG 
#'   sensors and creates a regularly-spaced grid over the sensor space bounding
#'   rectangle. A denser grid will give a very smooth-looking interpolated
#'   surface, but will take longer to compute.
#'  
#' @param sensor.space A matrix. This matrix should contain the x,y locations of
#'   the sensors.
#'   
#' @param grain A scalar. This represents that spacing between the points in the
#'   grid. Defaults to 40, to create a relatively coarse surface. Higher values
#'   will create a more fine-grained, and therefore smoother-looking surface. 
#'      
#' @return This function returns a data.frame with the cartesian coordinates of
#'   the bounding rectangle for the interpolation function.
#'   
#' @export
CreateGrid <- function(sensor.space, grain = 40L) {
  if (!is.matrix(sensor.space)) {
    sensors <- FormatPosition(sensor.space)
  } else {
    sensors <- sensor.space
  }
  
  bounding.rectangle <- FindBoundingRectangle(sensors)
  grid.xys <- apply(bounding.rectangle, 1, 
                    function(x){seq(x[1], x[2], length.out = grain)})
  grid.xys
}

#' @title Creates an interpolation surface in cartesian space for EEG/MEG 
#'   sensors.
#' 
#' @description  Creates an interpolation surface in cartesian space for EEG/MEG 
#'   sensors.
#'   
#' @param x A matrix or a vector. If x is a matrix, it needs to have two columns
#'   containig the (x, y) locations of the sensors, in which case y is not 
#'   required. If it is a vector, it should contain the x-coordinates of the 
#'   sensor.
#'   
#' @param y A vector containing the y-locations of the sensors. Defaults to 
#'   NULL, but will throw an error if x does not contain all the required
#'   information.
#'      
#' @return This function returns a data.frame with the cartesian coordinates of
#'   the surface for interpolation.
#'   
#' @export
CreateInterpolationSurface <- function(x, y = NULL) {
  if (!is.matrix(x)) {
    grid <- FormatPosition(x, y)
  } else {
    grid <- x
  }  
  unname(as.matrix(expand.grid(grid[, 1], grid[, 2])))
}

#' @title Tests whether a point is inside the convex hull created by the EEG/MEG 
#'   sensors.
#' 
#' @description Tests whether a point is inside the convex hull created by the 
#'   EEG/MEG sensors.
#'   
#' @details Calls function in.convex.hull from package tripack.
#'   
#' @param interpolation.surface A matrix containig the (x, y) locations of the 
#'   surface to be interpolated.
#'   
#' @param electrode.positions A matrix containing the (x, y) locations of the 
#'   EEG/MEG sensors.
#'      
#' @return This function returns TRUE or FALSE depending on whether a point is
#'   inside the convex hull.
#'   
#' @export
IsInConvexHull <- function(interpolation.surface, electrode.positions) {
  #if (!require("tripack")) install.packages("tripack")
  
  if (!is.matrix(electrode.positions)) {
    electrode.pos <- FormatPosition(electrode.positions)
  } else {
    electrode.pos <- electrode.positions
  }
  
  grid.query <- interpolation.surface
  x <- electrode.pos[, 1]
  y <- electrode.pos[, 2]
  electrodes.triangulation <- tripack::tri.mesh(x, y)
  tripack::in.convex.hull(electrodes.triangulation, 
                          grid.query[, 1], grid.query[, 2])
}

#' @title Tests whether a point is inside the minimum spanning circle of the 
#'   EEG/MEG sensors.
#' 
#' @description Tests whether a point is inside the minimum spanning circle of 
#'   the EEG/MEG sensors.
#'   
#' @details Calls function circumcircle from package tripack.
#'   
#' @param interpolation.surface A matrix containig the (x, y) locations of the 
#'   surface to be interpolated.
#'   
#' @param electrode.positions A matrix containing the (x, y) locations of the 
#'   EEG/MEG sensors.
#'
#' @param inflation A scalar. This represents a small fudge factor to make sure
#'   that the bounding circle will include all the sensors, in case there are
#'   sensors that would fall exactly in the circle line.
#'      
#' @return This function returns TRUE or FALSE depending on whether a point is
#'   inside the circumcircle.
#'   
#' @export
IsInCircle <- function(interpolation.surface, electrode.positions, 
                       inflation = 1e-7) {
  #if (!require("tripack")) install.packages("tripack")
  GetDistanceFromCenter <- function(point, center) {
    point - center
  }
  
  if (!is.matrix(electrode.positions)) {
    electrode.pos <- FormatPosition(electrode.positions)
  } else {
    electrode.pos <- electrode.positions
  }
  grid.query <- interpolation.surface
  x <- electrode.pos[, 1]
  y <- electrode.pos[, 2]
  
  elec.sphere <- tripack::circumcircle(x, y, num.touch = 3, plot = FALSE, 
                                       debug = FALSE)
  center.sphere <- c(elec.sphere$x, elec.sphere$y)
  squared.distances <- t(
    apply(grid.query, 1, GetDistanceFromCenter, center = center.sphere)
  )^2
  distances <- sqrt(squared.distances[, 1] + squared.distances[, 2])
  in.circle <- distances <= elec.sphere$radius + inflation
  in.circle
}

#' @title Interpolate EEG/MEG data from sensors to a predefined sensor space.
#' 
#' @description Interpolate EEG/MEG data from sensors to a predefined sensor 
#'   space.
#'   
#' @details Requires packages akima and rgcvpack for different types of data
#'   interpolation. Akima implements linear and cubic interpolation, while
#'   rgcvpack implements thin-plate spline interpolation.
#'   
#' @param potentials A vector containing data to be interpolated.
#'   
#' @param x.elec A matrix or a vector. If matrix, needs to contain the cartesian
#'   coordinates of the sensors. If vector, needs to contain the x-coordinates 
#'   of the sensors.
#'   
#' @param y.elec A vector. The y-coordinates of the sensors.
#'
#' @param grid.grain A scalar. This represents that spacing between the points in the
#'   grid. Defaults to 40, to create a relatively coarse surface. Higher values
#'   will create a more fine-grained, and therefore smoother-looking surface. 
#'      
#' @param x.vert A scalar. Defaults to NULL. Akima may use this?
#'   
#' @param y.vert A scalar. Defaults to NULL. Akima may use this?
#' 
#' @param extrapolate A logical. Should the data be extrapolated beyond the 
#'   convex hull created by the sensors? Defaults to FALSE. If TRUE, the
#'   interpolation will use the minimum spanning circle of the sensor space.
#'   
#' @param method A character string. Which interpolation method should be used?
#'   Currently accepts only three values: 'linear', for linear interpolation,
#'   'cubic' for cubic interpolation and 'tps', for thin-plate spline 
#'   interpolation.
#'   
#' @param electrode.dim A scalar. In matrix 'potentials', what dimension 
#'   represents time? Defaults to NULL, in which case the function takes a 
#'   guess.
#'
#' @param simplify A logical. Should output be simplified to save space? 
#'   Defaults to false.
#' 
#' @return List with interpolated potentials for each time slice.
#' 
#' @export
InterpolatePotentials <- function(potentials, x.elec, y.elec = NULL, 
                                  grid.grain = 40L, x.vert = NULL, 
                                  y.vert = NULL, extrapolate = FALSE, 
                                  method = c("linear", "cubic", "tps"),
                                  electrode.dim = NULL, simplify = FALSE) {
  
  #if (!require("plyr")) install.packages("plyr")
  #if (!require("akima")) install.packages("akima")
  #if (!require("rgcvpack")) install.packages("rgcvpack")
  
  if (!is.matrix(x.elec)) {
    electrode.locations <- FormatPosition(x.elec, y.elec)
  } else {
    electrode.locations <- x.elec
  }
  if (all(!is.null(x.vert), !is.null(y.vert))) {
    if (!is.matrix(x.vert)) {
      surface.vertices <- FormatPosition(x.vert, y.vert)
    } else {
      surface.vertices <- x.vert
    }
  } else {
    surface.vertices <- CreateGrid(electrode.locations, grain = grid.grain)
  }
  # The logical flag 'linear' is by default FALSE, since it is only used to
  # trigger the linear interpolation method in the Akima package
  linear <- FALSE  
  if (any(method == "linear", method == "cubic")) {
    interp.fnc <- "InterpolateAkima"
    if(method == "linear") {
      linear <- TRUE
    }
  } else {
    if (method == "tps") {
      interp.fnc <- "InterpolateTPS"
    } else {
      stop("'method' has to be either 'linear', 'cubic', or 'tps'.")
    }
  }
  x1  <- electrode.locations[, 1]
  y1  <- electrode.locations[, 2]
  z1  <- as.matrix(potentials)
  xo  <- surface.vertices[, 1]
  yo  <- surface.vertices[, 2]
  
  if (length(x1) != length(y1)) {
    stop("arguments 'x' and 'y' differ in length!")
  }
  n.electrodes <- length(x1)
  if (is.null(electrode.dim)) {
    time.dimension <- which(dim(potentials) != n.electrodes)
  } else {
    time.dimension <- electrode.dim
  }
  if (length(time.dimension) > 1) {
    warning("The data in 'potentials' has the N rows and columns,\n",
            "where N = number of electrodes!\n", 
            "Therefore, the function can not correctly guess in which dimension\n", 
            "time is represented and in which dimension the sensors are\n",
            "represented!\n",
            "In this case, the function arbitrarily picks the _rows_ dimension\n.",
            "Make sure this is the correct output!\n",
            "If not, you need to explicitly set 'electrode.dim' in the function call")
  }
  
  # Call plyr for the data structure with the user-selected interpolation
  # function.
  # cat(dim(z1), time.dimension, n.electrodes)
  out <- plyr::alply(.data = z1, .margins = time.dimension, 
               .fun = eval(parse(text = interp.fnc)), 
               x = x1, y = y1, xo = xo, yo = yo, grid.grain = grid.grain, 
               extrap = extrapolate, linear = linear)
  out
}

#' @title Interpolate EEG/MEG data from sensors to a predefined sensor space
#'   using the functions in package AKIMA.
#' 
#' @description Interpolate EEG/MEG data from sensors to a predefined sensor 
#'   space.
#'   
#' @details Requires packages akima for linear and cubic interpolation.
#'   
#' @param potentials A vector containing data to be interpolated.
#'   
#' @param x A matrix or a vector. If matrix, needs to contain the cartesian
#'   coordinates of the sensors. If vector, needs to contain the x-coordinates 
#'   of the sensors.
#'   
#' @param y A vector. The y-coordinates of the sensors.
#'
#' @param grid.grain A scalar. For the time being it is ignored (defaults to 40)
#'      
#' @param xo A vector defining the size of the grid in the x-dimension?
#'   Defaults to 40.
#'   
#' @param yo A vector defining the size of the grid in the y-dimension?
#'   Defaults to 40.
#' 
#' @param extrapolate A logical. Should the data be extrapolated beyond the 
#'   convex hull created by the sensors? Defaults to FALSE. If TRUE, the
#'   interpolation will use the minimum spanning circle of the sensor space.
#'   
#' @param linear A logical. Should the function use the linear interpolation
#'   method? Defaults to FALSE, which uses cubic interpolation.
#' 
#' @return List with interpolated potentials for each time slice.
#' 
#' @export
InterpolateAkima <- function(potentials, x, y, 
                             xo = seq(min(x), max(x), length = 40), 
                             yo = seq(min(y), max(y), length = 40),
                             grid.grain = NULL, extrapolate = FALSE, 
                             linear = FALSE) {
  # Note: 'grid.grain' is ignored in this function. It is here only to guarantee
  # interface compatibility with the InterpolateTPS function.
  
  #if (!require("akima")) install.packages("akima")
  
  if (any(length(x) != length(y), length(x) != length(potentials))) {
    stop("arguments 'x' and 'y', and 'potentials' differ in length!\n",
         paste0("length of 'x' = ", length(x), "\n"),
         paste0("length of 'y' = ", length(y), "\n"),
         paste0("length of 'potentials' = ", length(potentials), "\n")
    )
  }
  x1 <- x
  y1 <- y
  elec.loc <- data.frame(x1, y1)
  x.vertices <- xo
  y.vertices <- xo
  interpolated.potentials <- akima::interp(x = x1, y = y1, z = potentials, 
                                    xo = x.vertices, yo = y.vertices,
                                    extrap = extrapolate, linear = linear)
  interpolation.grid <- FormatPosition(x.vertices, y.vertices)
  interpolation.surface <- CreateInterpolationSurface(interpolation.grid)
  interpolated.surface <- interpolated.potentials$z
  out.zo <- interpolated.surface
  if(extrapolate) {
    mask <- IsInCircle(interpolation.surface, elec.loc)
    out.zo[!mask] <- NA
  }
  list(x = interpolated.potentials$x,
       y = interpolated.potentials$y,
       z = out.zo
  )  
}

#' @title Interpolate EEG/MEG data from sensors to a predefined sensor space
#'   using the functions in package rcgvpack.
#' 
#' @description Interpolate EEG/MEG data from sensors to a predefined sensor 
#'   space.
#'   
#' @details Requires packages rgcvpack for thin-plate spline interpolation.
#'   
#' @param potentials A vector containing data to be interpolated.
#'   
#' @param x A matrix or a vector. If matrix, needs to contain the cartesian
#'   coordinates of the sensors. If vector, needs to contain the x-coordinates 
#'   of the sensors.
#'   
#' @param y A vector. The y-coordinates of the sensors.
#'
#' @param grid.grain A scalar. For the time being it is ignored (defaults to 40)
#'      
#' @param xo A vector defining the size of the grid in the x-dimension?
#'   Defaults to 40.
#'   
#' @param yo A vector defining the size of the grid in the y-dimension?
#'   Defaults to 40.
#' 
#' @param extrapolate A logical. Should the data be extrapolated beyond the 
#'   convex hull created by the sensors? Defaults to FALSE. If TRUE, the
#'   interpolation will use the minimum spanning circle of the sensor space.
#'   
#' @param linear A logical. Should the function use the linear interpolation
#'   method? Defaults to FALSE, which uses cubic interpolation.
#' 
#' @return List with interpolated potentials for each time slice.
#' 
#' @export
InterpolateTPS <- function(potentials, x, y, xo, yo, grid.grain, 
                           extrapolate = FALSE, linear = NULL) {
  #if (!require("rgcvpack")) install.packages("rgcvpack")
  # Note: 'linear' is ignored
  #out.xo <- xo # these two vectors are going to be ignored as well.
  #out.yo <- yo
  
  elec.loc <- data.frame(x, y)
  x.vertices <- xo
  y.vertices <- xo
  interpolation.grid <- FormatPosition(x.vertices, y.vertices)
  interpolation.surface <- CreateInterpolationSurface(interpolation.grid)
  fitted.surface <- rgcvpack::fitTps(elec.loc, y = potentials)
  interpolated.surface <- rgcvpack::predict.Tps(fitted.surface, 
                                                interpolation.surface)
  out.xo <- x.vertices
  out.yo <- y.vertices
  out.zo <- matrix(interpolated.surface, nrow = length(x.vertices), 
                   byrow = FALSE)
  if(!extrapolate) {
    mask <- IsInConvexHull(interpolation.surface, elec.loc)
  } else {
    mask <- IsInCircle(interpolation.surface, elec.loc)
  }
  out.zo[!mask] <- NA
  list(x = out.xo,
       y = out.yo,
       z = out.zo
  )       
}

#' @title Simplify a list to save space
#' 
#' @description Simplify a list to save space.
#'
#' @details The output of function InterpolatePotentials is very redundant. 
#'   This function simplifies it by removing the x, y fields from every item on 
#'   the list, since they are always the same.
#'   
#' @param topography.list A list. Natural output of function 
#'   InterpolatePotentials
#' 
#' @return List with interpolated potentials for each time slice, simplified.
#' 
#' @export
SimplifyTopographyList <- function(topography.list) {
  in.list <- topography.list[[1]]
  out.list <- cbind(expand.grid(in.list$x, in.list$y), as.vector(in.list$z))
  colnames(out.list) <- names(in.list)
  out.list
}

#' @title Make a symmetric range for color bar
#' 
#' @description Make a symmetric range for color bar
#'
#' @param z.values A vector of potentials.
#' 
#' @return vector of length 2 with symmetric range.
#' 
#' @export
MakeSymmetric <- function(z.values) {
  z1 <- c(max(abs(range(z.values, na.rm = TRUE))))
  c(-z1, z1)
}

#' @title Create labels for symmetric color bar
#' 
#' @description Create labels for symmetric color bar
#'
#' @param z.values A vector of potentials.
#' 
#' @param multiplication.factor Correction factor for graphic display of units 
#'   of measurement (generally microvolts or femtotesla).
#'   
#' @param unit A string. Can be either 'mV' (for microvolts, usually for EEG 
#'   data), 'fT' (for femtotesla, usually for MEG data) or 'mVcm2' (for 
#'   microvolts over squared centimeters, usually for CSD transformed EEG data).
#' 
#' @return vector of length 3 with symmetric range around 0.
#' 
#' @export
MakeSymmetricLabels <- function(z.values, multiplication.factor = 10, 
                                unit = NULL) {
  if (is.null(unit)) {
    unit.label = unit
  }
  if (!is.null(unit)) {
    if (!any(c("mv", "fT", "mvcm2") %in% unit)) {
      stop("invalid unit!")
    } else {
      if (unit == "mv") {
        unit.label <- " mV"
      } else {
        if (unit == "fT") {
          unit.label <- " fT"
        } else {
          if (unit == "mvcm2") {
            unit.label <- " mV/cm2"
          } else{
            unit.label <- NULL
          }
        }
      }
    }
  }
  z1   <- MakeSymmetric(z.values)
  low  <- round(min(z1) * multiplication.factor)
  mid  <- paste0("0", unit.label)
  high <- round(max(z1) * multiplication.factor)
  c(low, mid, high)
}

#' @title Create range for asymmetric color bar
#' 
#' @description Create labels for symmetric color bar
#'
#' @param z.values A vector of potentials.
#' 
#' @param steps A scalar. Number of steps in color bar.
#' 
#' @return vector of length 'steps' for color bar.
#' 
#' @export
MakeSymmetricCuts <- function(z.values, steps) {
  z1   <- MakeSymmetric(z.values)
  seq(from = z1[1], to = z1[2], length = steps)
}

