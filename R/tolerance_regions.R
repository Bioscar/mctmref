#' Bivariate Conditional Tolerance Region from an MCTM
#'
#' Estimates a conditional tolerance region that contains \eqn{100\tau\%} of
#' the reference population with confidence \eqn{1-\alpha}, following
#' Krishnamoorthy & Mondal (2006).
#'
#' @param m_full A fitted \code{mmlt} object.
#' @param m1 Fitted marginal \code{mlt} model for the first response.
#' @param m2 Fitted marginal \code{mlt} model for the second response.
#' @param newdata A \code{data.frame} of covariate values. One region per row.
#' @param train_data Optional \code{data.frame} used to fit the model. If
#'   \code{NULL} (default), the data stored inside \code{m_full} is used.
#' @param tau Single numeric coverage level in (0, 1). Default 0.95.
#' @param alpha Uncertainty level. Default 0.05 (95\% confidence).
#' @param B Integer. Monte Carlo iterations for tolerance factor. Default 1000.
#' @param n_grid Integer. Response grid size. Default 2000.
#'
#' @return A list with one \eqn{50 \times 2} numeric matrix per row of
#'   \code{newdata}, giving the tolerance region boundary in
#'   \eqn{(y_1, y_2)} space.
#'
#' @details
#' The tolerance factor is estimated by Monte Carlo following Krishnamoorthy &
#' Mondal (2006). The tolerance ellipse is built from the empirical mean and
#' covariance of the training Z-scores, then back-transformed to response
#' space. The resulting region is always larger than the corresponding
#' reference region from \code{\link{bivreg_mctm}}.
#'
#' @seealso \code{\link{bivreg_mctm}}, \code{\link{mctm_fit}}
#'
#' @references
#' Krishnamoorthy, K. and Mondal, S. (2006). Improved tolerance factors for
#' multivariate normal distributions. \emph{Communications in Statistics:
#' Simulation and Computation}, 35(2), 461--478.
#'
#' Lado-Baleato et al. (2023). \doi{10.1002/bimj.202200229}
#'
#' @examples
#' \dontrun{
#' tol <- tol_mctm(fit$m_full, fit$m1, fit$m2,
#'                 newdata = data.frame(age = 50))
#' lines(tol[[1]], col = "navy", lty = 2, lwd = 2)
#' }
#'
#' @export
tol_mctm <- function(m_full, m1, m2, newdata, train_data = NULL,
                     tau = 0.95, alpha = 0.05, B = 1000L, n_grid = 2000L) {

  .mctm_check_model(m_full, m1, m2)
  stopifnot(length(tau) == 1L, tau > 0, tau < 1)

  gr1 <- .mctm_response_grid(m1, n = n_grid)
  gr2 <- .mctm_response_grid(m2, n = n_grid)

  td  <- .mctm_train_data(m_full, m1, m2, train_data)

  Z1 <- as.numeric(predict(m_full, newdata = td, margins = 1L, type = "trafo"))
  Z2 <- as.numeric(predict(m_full, newdata = td, margins = 2L, type = "trafo"))

  biv <- cbind(Z1, Z2)
  biv <- biv[apply(biv, 1L, function(r) all(is.finite(r))), , drop = FALSE]

  # Tolerance factor — Krishnamoorthy & Mondal (2006)
  p <- ncol(biv);  n <- nrow(biv)
  q.sq <- matrix(stats::rchisq(p * B, df = 1L), ncol = p) / n
  Lmat <- t(sapply(seq_len(B), function(i)
    base::eigen(tolerance::rwishart(n - 1L, p))$values))
  c1 <- apply((1 + q.sq) / Lmat,       1L, sum)
  c2 <- apply((1 + 2 * q.sq) / Lmat^2, 1L, sum)
  c3 <- apply((1 + 3 * q.sq) / Lmat^3, 1L, sum)
  a  <- (c2^3) / (c3^2)
  Tmc     <- (n - 1L) * (sqrt(c2 / a) * (stats::qchisq(tau, a) - a) + c1)
  tol_val <- stats::quantile(Tmc, 1 - alpha)

  # Tolerance ellipse in Z space
  mu  <- apply(biv, 2L, mean)
  sig <- stats::cov(biv)
  es  <- base::eigen(sig)
  e1  <- es$vec %*% diag(sqrt(es$val))
  th  <- seq(0, 2 * pi, length.out = 50L)
  r1  <- sqrt(tol_val)
  pol <- t(mu - (e1 %*% t(cbind(r1 * cos(th), r1 * sin(th)))))

  # Back-transform
  lapply(seq_len(nrow(newdata)), function(k) {
    nd_k <- newdata[k, , drop = FALSE]
    t1   <- .mctm_trafo_grid(m_full, gr1$q, gr1$yname, 1L, nd_k)
    t2   <- .mctm_trafo_grid(m_full, gr2$q, gr2$yname, 2L, nd_k)
    yi   <- matrix(0, nrow = nrow(pol), ncol = 2L)
    for (i in seq_len(nrow(pol))) {
      yi[i, 1L] <- .mctm_invert_trafo(pol[i, 1L], gr1$q, t1)
      yi[i, 2L] <- .mctm_invert_trafo(pol[i, 2L], gr2$q, t2)
    }
    yi
  })
}
