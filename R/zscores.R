#' Extract Joint Z Scores from an MCTM
#'
#' Maps each observed \eqn{(y_1, y_2, x)} to latent normal scores
#' \eqn{(Z_1, Z_2) \sim N(0, \Sigma(x))} from the fitted joint model.
#'
#' @param m_full A fitted \code{mmlt} object.
#' @param m1 Fitted marginal \code{mlt} model for the first response.
#' @param m2 Fitted marginal \code{mlt} model for the second response.
#' @param data A \code{data.frame}. If \code{NULL}, training data from
#'   \code{m_full} is used.
#' @param remove_boundary Logical. Remove observations where \eqn{Z_1} or
#'   \eqn{Z_2} is non-finite (Bernstein basis boundary). Default \code{TRUE}.
#' @param decorrelate Logical. Also compute decorrelated residuals
#'   \eqn{U \sim N(0, I)} via \eqn{U = \mathrm{chol}(\Sigma^{-1}) Z}.
#'   Default \code{FALSE}.
#'
#' @return A \code{data.frame} with the input data columns plus:
#' \describe{
#'   \item{Z1}{Latent Z score for the first response.}
#'   \item{Z2}{Latent Z score for the second response.}
#'   \item{U1, U2}{Decorrelated residuals (only if \code{decorrelate = TRUE}).}
#' }
#'
#' @details
#' \code{Z1} and \code{Z2} have marginal distribution \eqn{\sim N(0,1)} but
#' are correlated: \eqn{\mathrm{cor}(Z_1, Z_2) \approx \rho_s(x)}. To obtain
#' independent standard normal residuals, set \code{decorrelate = TRUE}.
#'
#' Observations at the Bernstein basis boundary map to \eqn{\pm\infty}.
#' These are removed when \code{remove_boundary = TRUE}.
#'
#' @seealso \code{\link{mctm_mahalanobis}}, \code{\link{mctm_fit}}
#'
#' @references Lado-Baleato et al. (2023). \doi{10.1002/bimj.202200229}
#'
#' @examples
#' \dontrun{
#' zs <- mctm_zscores(fit$m_full, fit$m1, fit$m2, decorrelate = TRUE)
#' plot(zs$Z1, zs$Z2, main = "Z space (correlated)")
#' plot(zs$U1, zs$U2, main = "U space (independent)")
#' }
#'@importFrom stats coef logLik predict
#' @export
mctm_zscores <- function(m_full, m1, m2, data = NULL,
                          remove_boundary = TRUE,
                          decorrelate     = FALSE) {

  .mctm_check_model(m_full, m1, m2)
  if (is.null(data)) data <- .mctm_train_data(m_full, m1, m2)

  Z1 <- as.numeric(predict(m_full, newdata = data, margins = 1L, type = "trafo"))
  Z2 <- as.numeric(predict(m_full, newdata = data, margins = 2L, type = "trafo"))

  out <- cbind(data, Z1 = Z1, Z2 = Z2)

  if (remove_boundary) {
    ok <- is.finite(Z1) & is.finite(Z2)
    if (any(!ok)) message(sum(!ok), " boundary observation(s) removed.")
    out  <- out[ok, , drop = FALSE]
    Z1   <- Z1[ok];  Z2 <- Z2[ok]
    data <- data[ok, , drop = FALSE]
  }

  if (decorrelate) {
    Sigma_arr <- .mctm_get_sigma_all(m_full, data)
    Z_mat     <- cbind(Z1, Z2)
    U         <- matrix(NA_real_, nrow = nrow(data), ncol = 2L)
    for (i in seq_len(nrow(data))) {
      S      <- Sigma_arr[, , i]
      U[i, ] <- t(base::chol(base::solve(S)) %*% Z_mat[i, ])
    }
    out$U1 <- U[, 1L]
    out$U2 <- U[, 2L]
  }

  out
}


#' Mahalanobis Distance in Z Space
#'
#' Computes \eqn{D^2 = Z^\top \Sigma(x)^{-1} Z} for each observation.
#' Under the model, \eqn{D^2 \sim \chi^2_2}, so this provides a natural
#' outlier score equivalent to the reference region classification.
#'
#' @param m_full A fitted \code{mmlt} object.
#' @param m1 Fitted marginal \code{mlt} model for the first response.
#' @param m2 Fitted marginal \code{mlt} model for the second response.
#' @param data A \code{data.frame}. If \code{NULL}, training data is used.
#'
#' @return A \code{data.frame} with input columns plus:
#' \describe{
#'   \item{Z1, Z2}{Marginal Z scores.}
#'   \item{D2}{Mahalanobis distance squared.}
#'   \item{pvalue}{P-value under \eqn{\chi^2_2}.}
#'   \item{outside_95}{Logical: \code{TRUE} if \eqn{D^2 > \chi^2_2(0.95) = 5.99}.}
#' }
#'
#' @seealso \code{\link{mctm_zscores}}, \code{\link{bivreg_mctm}}
#'
#' @references Lado-Baleato et al. (2023). \doi{10.1002/bimj.202200229}
#'
#' @examples
#' \dontrun{
#' md <- mctm_mahalanobis(fit$m_full, fit$m1, fit$m2)
#' cat("Outside 95% region:", mean(md$outside_95) * 100, "%\n")
#' plot(md$age, md$D2, col = ifelse(md$outside_95, "red", "grey"))
#' abline(h = qchisq(0.95, 2), lty = 2, col = "blue")
#' }
#'
#' @export
mctm_mahalanobis <- function(m_full, m1, m2, data = NULL) {

  .mctm_check_model(m_full, m1, m2)
  zs     <- mctm_zscores(m_full, m1, m2, data = data,
                          remove_boundary = TRUE, decorrelate = FALSE)
  Z1     <- zs$Z1;  Z2 <- zs$Z2
  data_f <- zs[, !names(zs) %in% c("Z1", "Z2"), drop = FALSE]

  Sigma_arr <- .mctm_get_sigma_all(m_full, data_f)
  Z_mat     <- cbind(Z1, Z2)

  D2 <- numeric(nrow(data_f))
  for (i in seq_len(nrow(data_f))) {
    z     <- Z_mat[i, ]
    S     <- Sigma_arr[, , i]
    D2[i] <- as.numeric(t(z) %*% base::solve(S) %*% z)
  }

  pval <- stats::pchisq(D2, df = 2L, lower.tail = FALSE)
  cbind(zs, D2 = D2, pvalue = pval, outside_95 = pval < 0.05)
}
