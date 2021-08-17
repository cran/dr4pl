#plotting dr4pl


#' @title Make a plot of a 4PL model curve and data
#'
#' @description This function displays a dose-response curve and data. As a default,
#' the x-axis represents dose levels in log 10 scale and the y-axis represents
#' responses. The black solid line represents a dose-response curve. The blue filled
#' circles represent data points and red triangles represent outliers.
#'
#' @name plot.dr4pl
#'
#' @param x `dr4pl' object whose data and mean response function will be plotted.
#' @param type.curve Indicator of the type of a dose-response curve. "all" indicates
#' that data and a curve will be plotted while "data" indicates that only data
#' will be plotted.
#' @param text.title Character string for the title of a plot with a default set to
#'   "Dose response plot".
#' @param text.x Character string for the x-axis of the plot with a default set to
#'   "Dose".
#' @param text.y Character string for the y-axis of the plot with a default set to
#'   "Response".
#' @param indices.outlier Pass a vector indicating all indices which are outliers in
#'   the data.
#' @param breaks.x Vector of desired break points for the x-axis
#' @param breaks.y Vector of desired break points for the y-axis
#' @param ... All arguments that can normally be passed to ggplot.
#'
#' @examples
#' \dontrun{
#' dr4pl.1 <- dr4pl(Response ~ Dose, data = sample_data_1)
#'
#' plot(dr4pl.1)
#'
#' ## Able to further edit plots.
#' library(ggplot2) #needed to change color to green
#' dr4pl.1 <- dr4pl(Response ~ Dose,
#'                         data = sample_data_1,
#'                         text.title = "Sample Data Plot")
#'
#' a <- plot(dr4pl.1)
#' a + geom_point(color = "green", size = 5)
#'
#' ## Bring attention to outliers using parameter indices.outlier.
#' dr4pl.3 <- dr4pl(Response ~ Dose,
#'                  data = drc_error_3,
#'                  method.init = "Mead",
#'                  method.robust = "absolute")
#' plot(dr4pl.3, indices.outlier = c(90, 101))
#'
#' ## Change the plot title default with parameter text.title.
#' dr4pl.1 <- dr4pl::dr4pl(Response ~ Dose,
#'                         data = sample_data_1)
#' plot(dr4pl.1, text.title = "My New Dose Response plot")
#'
#' ##Change the labels of the x and y axis to your need
#' library(drc)  # Needed to load 'decontaminants' data set
#' data.hpc <- subset(decontaminants, group %in% "hpc")
#' dr4pl.hpc <- dr4pl(count~conc, data = data.hpc)
#' plot(dr4pl.hpc,
#'      text.title = "hpc Decontaminants Plot",
#'      text.x = "Concentration",
#'      text.y = "Count")
#' }
#'
#' @author Hyowon An, \email{ahwbest@gmail.com}
#' @author Justin T. Landis, \email{jtlandis314@gmail.com}
#' @author Aubrey G. Bailey, \email{aubreybailey@gmail.com}
#'
#' @export
plot.dr4pl <- function(x,
                       type.curve = "all",
                       text.title = "Dose-response plot",
                       text.x = "Dose",
                       text.y = "Response",
                       indices.outlier = NULL,
                       breaks.x = NULL,
                       breaks.y = NULL,
                       ...) {

  ### Check whether function arguments are appropriate
  if(!is.character(text.title)) {

    stop("Title text should be characters.")
  }
  if(!is.character(text.x)) {

    stop("The x-axis label text should be characters.")
  }
  if(!is.character(text.y)) {

    stop("The y-axis label text should be characters.")
  }

  ### Draw a plot
  n <- x$sample.size
  color.vec <- rep("blue", n)
  shape.vec <- rep(19, n)

  if(!is.null(indices.outlier)) {

    color.vec[indices.outlier] <- "red"
    shape.vec[indices.outlier] <- 17
  }

  a <- ggplot2::ggplot(aes(x = x$data$Dose, y = x$data$Response), data = x$data)

  if(type.curve == "all") {

    a <- a + ggplot2::stat_function(fun = MeanResponse.dr4pl_theta,
                                    args = list(theta = x$parameters),
                                    size = 1.2)
  }

  a <- a + ggplot2::geom_point(size = I(5), alpha = I(0.8), color = color.vec,
                               shape = shape.vec)

  a <- a + ggplot2::labs(title = text.title,
                         x = text.x,
                         y = text.y)

  # Set parameters for the grids
  a <- a + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 16))
  a <- a + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  a <- a + ggplot2::theme(panel.grid.major = ggplot2::element_blank())

  if(!is.null(breaks.x)) {

    a <- a + ggplot2::scale_x_log10(breaks = breaks.x)
  } else {

    a <- a + ggplot2::scale_x_log10()
  }
  if(!is.null(breaks.y)) {

    a <- a + ggplot2::scale_y_continuous(breaks = breaks.y)
  } else {

    a <- a + ggplot2:: scale_y_continuous()
  }

  a <- a + ggplot2::theme_bw()
  # Test
  # Set parameters for the titles and text / margin(top, right, bottom, left)
  a <- a + ggplot2::theme(plot.title = ggplot2::element_text(size = 20, margin = ggplot2::margin(0, 0, 10, 0)))
  a <- a + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 16, margin = ggplot2::margin(15, 0, 0, 0)))
  a <- a + ggplot2::theme(axis.title.y = ggplot2::element_text(size = 16, margin = ggplot2::margin(0, 15, 0, 0)))
  a <- a + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16))
  a <- a + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 16))

  return(a)
}
