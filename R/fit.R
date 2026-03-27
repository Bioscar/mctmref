#' Fit a Bivariate Multivariate Conditional Transformation Model
#'
#' Convenience wrapper that defines Bernstein basis functions and fits the
#' two marginal \code{mlt} models and the joint \code{mmlt} model in one call.
#'
#' @param data A \code{data.frame} with response and covariate columns.
#' @param y1 Character. Name of the first response (e.g. \code{"fpg"}).
#' @param y2 Character. Name of the second response (e.g. \code{"hba1c"}).
#' @param x  Character. Name of the continuous covariate (e.g. \code{"age"}).
#' @param bounds_y1 Numeric length-2 vector. Hard bounds for \code{y1}'s
#'   Bernstein basis. Default \code{range(data[[y1]])}.
#' @param bounds_y2 Numeric length-2 vector. Hard bounds for \code{y2}'s
#'   Bernstein basis. Default \code{range(data[[y2]])}.
#' @param bounds_x  Numeric length-2 vector. Hard bounds for the covariate
#'   basis. Default \code{range(data[[x]])}.
#' @param order_y Integer. Bernstein polynomial order for response bases.
#'   Default 6 (recommended by Klein et al. 2022).
#' @param order_x Integer. Bernstein polynomial order for covariate basis.
#'   Default 6.
#' @param support_prob Numeric length-2 vector. Quantile probabilities defining
#'   the inner support of each Bernstein basis. Default \code{c(0.05, 0.95)}.
#'
#' @return An object of class \code{mctm_fit}, a list with elements:
#' \describe{
#'   \item{m1}{Fitted marginal \code{mlt} model for \code{y1}.}
#'   \item{m2}{Fitted marginal \code{mlt} model for \code{y2}.}
#'   \item{m_full}{Fitted joint \code{mmlt} model.}
#'   \item{data}{Data used for fitting.}
#'   \item{y1, y2, x}{Variable names.}
#' }
#'
#' @seealso \code{\link{bivreg_mctm}}, \code{\link{tol_mctm}},
#'   \code{\link{mmlt_perc}}
#'
#' @references
#' Lado-Baleato et al. (2023). Multivariate reference and tolerance regions
#' based on conditional transformation models. \emph{Biometrical Journal},
#' 65, 2200229. \doi{10.1002/bimj.202200229}
#'
#' @examples
#' \dontrun{
#' library(refreg)
#' df  <- subset(refreg::aegis, dm == "no")
#' fit <- mctm_fit(df, y1 = "fpg", y2 = "hba1c", x = "age",
#'                 bounds_y1 = c(64, 183), bounds_y2 = c(4.3, 7.8),
#'                 bounds_x  = c(19, 89))
#' print(fit)
#' }
#'
#' @importFrom basefun Bernstein_basis
#' @importFrom variables numeric_var
#' @importFrom mlt ctm mlt mmlt
#' @import tram
#' @import utils
#' @export
mctm_fit <- function(data, y1, y2, x,
                     bounds_y1    = range(data[[y1]]),
                     bounds_y2    = range(data[[y2]]),
                     bounds_x     = range(data[[x]]),
                     order_y      = 6L,
                     order_x      = 6L,
                     support_prob = c(0.05, 0.95)) {

  stopifnot(is.data.frame(data))
  for (v in c(y1, y2, x))
    if (!v %in% names(data))
      stop(sprintf("Variable '%s' not found in data.", v), call. = FALSE)

  Bxlambda <- Bernstein_basis(
    numeric_var(x,
                support = stats::quantile(data[[x]], prob = support_prob),
                bounds  = bounds_x),
    order = order_x, extrapolate = TRUE)

  by1 <- Bernstein_basis(
    numeric_var(y1,
                support = stats::quantile(data[[y1]], prob = support_prob),
                bounds  = bounds_y1),
    order = order_y, ui = "increasing")

  by2 <- Bernstein_basis(
    numeric_var(y2,
                support = stats::quantile(data[[y2]], prob = support_prob),
                bounds  = bounds_y2),
    order = order_y, ui = "increasing")

  message("Fitting marginal model for ", y1, " ...")
  m1 <- mlt(ctm(response = by1, interacting = Bxlambda, todistr = "Normal"),
            data = data)

  message("Fitting marginal model for ", y2, " ...")
  m2 <- mlt(ctm(response = by2, interacting = Bxlambda, todistr = "Normal"),
            data = data)

  message("Fitting joint mmlt model ...")
  m_full <- mmlt(m1, m2, formula = Bxlambda, data = data)

  message("Done. logLik = ", round(as.numeric(logLik(m_full)), 2))

  structure(list(m1 = m1, m2 = m2, m_full = m_full,
                 data = data, y1 = y1, y2 = y2, x = x),
            class = "mctm_fit")
}

#' @export
print.mctm_fit <- function(x, ...) {
  cat("Bivariate MCTM fit\n")
  cat("  Responses :", x$y1, "and", x$y2, "\n")
  cat("  Covariate :", x$x,  "\n")
  cat("  n         :", nrow(x$data), "\n")
  cat("  logLik    :", round(as.numeric(logLik(x$m_full)), 2), "\n")
  cat("  Parameters:", length(coef(x$m_full)), "\n")
  invisible(x)
}
