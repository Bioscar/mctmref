#' Joint Log-Likelihood for Observed Data
#'
#' Evaluates the joint bivariate log-likelihood \eqn{\log f(y_1, y_2 \mid x)}
#' for each observation. This includes the Gaussian copula term and is
#' \strong{not} the sum of marginal log-likelihoods.
#'
#' @param m_full A fitted \code{mmlt} object.
#' @param data A \code{data.frame} with both responses and covariate.
#'   If \code{NULL}, the training data is used.
#'
#' @return A \code{data.frame} with input columns plus:
#' \describe{
#'   \item{joint_density}{Joint bivariate density value.}
#'   \item{joint_loglik}{Log of the joint density.}
#'   \item{marginal_loglik_y1}{Log marginal density for \eqn{y_1}.}
#'   \item{marginal_loglik_y2}{Log marginal density for \eqn{y_2}.}
#'   \item{copula_contribution}{Joint minus independence log-likelihood
#'     (the copula term).}
#' }
#'
#' @details
#' \code{predict(m_full, type = "density")} (no \code{margins} argument)
#' returns the joint bivariate density. The copula contribution is
#' \eqn{\log f(y_1, y_2 \mid x) - \log f_1(y_1 \mid x) - \log f_2(y_2 \mid x)}.
#'
#' @seealso \code{\link{mctm_loglik_grid}}, \code{\link{mctm_loglik_by_age}}
#'
#' @references Lado-Baleato et al. (2023). \doi{10.1002/bimj.202200229}
#'
#' @examples
#' \dontrun{
#' ll <- mctm_loglik(fit$m_full)
#' cat("Total joint log-lik:", sum(ll$joint_loglik, na.rm = TRUE), "\n")
#' cat("Copula contribution:", sum(ll$copula_contribution, na.rm = TRUE), "\n")
#' }
#'
#' @export
mctm_loglik <- function(m_full, data = NULL) {
  if (!inherits(m_full, "mmlt"))
    stop("'m_full' must be an mmlt object.", call. = FALSE)
  if (is.null(data)) data <- m_full$data
  eps  <- .Machine$double.eps
  jd   <- as.numeric(predict(m_full, newdata = data, type = "density"))
  md1  <- as.numeric(predict(m_full, newdata = data, margins = 1L, type = "density"))
  md2  <- as.numeric(predict(m_full, newdata = data, margins = 2L, type = "density"))
  lj   <- log(pmax(jd,  eps))
  lm1  <- log(pmax(md1, eps))
  lm2  <- log(pmax(md2, eps))
  cbind(data,
        joint_density       = jd,
        joint_loglik        = lj,
        marginal_loglik_y1  = lm1,
        marginal_loglik_y2  = lm2,
        copula_contribution = lj - lm1 - lm2)
}


#' Joint Log-Likelihood on a Response Grid
#'
#' Evaluates \eqn{f(y_1, y_2 \mid x)} on a grid of \eqn{(y_1, y_2)} values
#' at a fixed covariate value. Useful for density contour plots.
#'
#' @param m_full A fitted \code{mmlt} object.
#' @param x_val Numeric scalar. Covariate value.
#' @param x_name Character. Covariate variable name. Default \code{"age"}.
#' @param y1_grid Numeric vector. Grid for \eqn{y_1}. If \code{NULL}, derived
#'   from \code{mkgrid(m1)}.
#' @param y2_grid Numeric vector. Grid for \eqn{y_2}. If \code{NULL}, derived
#'   from \code{mkgrid(m2)}.
#' @param n_grid Integer. Grid size when \code{y1_grid}/\code{y2_grid} are
#'   auto-generated. Default 50.
#' @param m1 Fitted marginal \code{mlt} model for \eqn{y_1} (needed for
#'   auto grid range).
#' @param m2 Fitted marginal \code{mlt} model for \eqn{y_2}.
#'
#' @return A \code{data.frame} with \eqn{y_1}, \eqn{y_2}, the covariate,
#'   \code{joint_density}, and \code{joint_loglik}.
#'
#' @seealso \code{\link{mctm_loglik}}, \code{\link{mctm_loglik_by_age}}
#'
#' @examples
#' \dontrun{
#' g <- mctm_loglik_grid(fit$m_full, x_val = 50, m1 = fit$m1, m2 = fit$m2)
#' fpg_vals   <- sort(unique(g$fpg))
#' hba1c_vals <- sort(unique(g$hba1c))
#' contour(fpg_vals, hba1c_vals,
#'         matrix(g$joint_density, nrow = length(fpg_vals)))
#' }
#'
#' @export
mctm_loglik_grid <- function(m_full, x_val, x_name = "age",
                              y1_grid = NULL, y2_grid = NULL,
                              n_grid = 50L, m1 = NULL, m2 = NULL) {
  if (!inherits(m_full, "mmlt"))
    stop("'m_full' must be an mmlt object.", call. = FALSE)

  if ((is.null(y1_grid) || is.null(y2_grid)) && (is.null(m1) || is.null(m2)))
    stop("Supply 'm1' and 'm2' for auto grid, or provide 'y1_grid'/'y2_grid'.",
         call. = FALSE)

  y1_name <- if (!is.null(m1)) stats::variable.names(m1, "response") else m_full$names[1L]
  y2_name <- if (!is.null(m2)) stats::variable.names(m2, "response") else m_full$names[2L]

  if (is.null(y1_grid)) y1_grid <- .mctm_response_grid(m1, n = n_grid)$q
  if (is.null(y2_grid)) y2_grid <- .mctm_response_grid(m2, n = n_grid)$q

  grid <- expand.grid(stats::setNames(
    list(y1_grid, y2_grid, x_val), c(y1_name, y2_name, x_name)))
  jd   <- as.numeric(predict(m_full, newdata = grid, type = "density"))
  cbind(grid, joint_density = jd,
        joint_loglik = log(pmax(jd, .Machine$double.eps)))
}


#' Joint Log-Likelihood as a Function of the Covariate
#'
#' For a patient with known \eqn{(y_1, y_2)}, evaluates
#' \eqn{\log f(y_1, y_2 \mid x)} across a range of covariate values.
#' The maximising value is the covariate level at which the patient's
#' marker combination is most typical.
#'
#' @param m_full A fitted \code{mmlt} object.
#' @param y1_val Numeric. Observed value of the first response.
#' @param y2_val Numeric. Observed value of the second response.
#' @param x_seq Numeric vector. Covariate values to evaluate.
#'   Default \code{seq(20, 85, by = 0.5)}.
#' @param x_name Character. Covariate name. Default \code{"age"}.
#' @param m1 Fitted marginal \code{mlt} model (to get \eqn{y_1} name).
#' @param m2 Fitted marginal \code{mlt} model (to get \eqn{y_2} name).
#'
#' @return A \code{data.frame} with columns \code{x}, \code{joint_density},
#'   \code{joint_loglik}, and \code{most_likely} (logical).
#'
#' @seealso \code{\link{mctm_loglik}}, \code{\link{mctm_loglik_grid}}
#'
#' @examples
#' \dontrun{
#' prof <- mctm_loglik_by_age(fit$m_full, y1_val = 95, y2_val = 5.5,
#'                             m1 = fit$m1, m2 = fit$m2)
#' plot(prof$x, prof$joint_loglik, type = "l",
#'      xlab = "Age", ylab = "log f(fpg, hba1c | age)")
#' abline(v = prof$x[prof$most_likely], lty = 2)
#' }
#'
#' @export
mctm_loglik_by_age <- function(m_full, y1_val, y2_val,
                                x_seq  = seq(20, 85, by = 0.5),
                                x_name = "age",
                                m1 = NULL, m2 = NULL) {
  if (!inherits(m_full, "mmlt"))
    stop("'m_full' must be an mmlt object.", call. = FALSE)
  y1_name <- if (!is.null(m1)) stats::variable.names(m1, "response") else m_full$names[1L]
  y2_name <- if (!is.null(m2)) stats::variable.names(m2, "response") else m_full$names[2L]
  nd  <- stats::setNames(data.frame(y1_val, y2_val, x_seq),
                         c(y1_name, y2_name, x_name))
  jd  <- as.numeric(predict(m_full, newdata = nd, type = "density"))
  ll  <- log(pmax(jd, .Machine$double.eps))
  data.frame(x = x_seq, joint_density = jd, joint_loglik = ll,
             most_likely = ll == max(ll[is.finite(ll)]))
}
