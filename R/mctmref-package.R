#' mctmref: Multivariate Conditional Reference Regions via Transformation Models
#'
#' Estimates bivariate conditional reference and tolerance regions based on
#' Multivariate Conditional Transformation Models (MCTMs).
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{mctm_fit}}}{Fit a bivariate MCTM.}
#'   \item{\code{\link{mmlt_perc}}}{Marginal conditional percentile curves.}
#'   \item{\code{\link{bivreg_mctm}}}{Bivariate conditional reference region.}
#'   \item{\code{\link{tol_mctm}}}{Bivariate conditional tolerance region.}
#'   \item{\code{\link{mctm_zscores}}}{Joint Z scores and decorrelated residuals.}
#'   \item{\code{\link{mctm_mahalanobis}}}{Mahalanobis distance in Z space.}
#'   \item{\code{\link{mctm_loglik}}}{Joint log-likelihood per observation.}
#'   \item{\code{\link{mctm_loglik_grid}}}{Joint density on a response grid.}
#'   \item{\code{\link{mctm_loglik_by_age}}}{Log-likelihood profile over the covariate.}
#' }
#'
#' @references
#' Lado-Baleato, O., Cadarso-Suarez, C., Kneib, T. and Gude, F. (2023).
#' Multivariate reference and tolerance regions based on conditional
#' transformation models. \emph{Biometrical Journal}, 65, 2200229.
#' \doi{10.1002/bimj.202200229}
#'
#' @docType package
#' @name mctmref-package
#' @aliases mctmref
"_PACKAGE"
