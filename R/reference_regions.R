#' Bivariate Conditional Reference Region from an MCTM
#'
#' Estimates the bivariate conditional reference region \eqn{R_\tau(x)} at
#' one or more covariate values, following Lado-Baleato et al. (2023).
#'
#' @param m_full A fitted \code{mmlt} object.
#' @param m1 Fitted marginal \code{mlt} model for the first response.
#' @param m2 Fitted marginal \code{mlt} model for the second response.
#' @param tau Numeric scalar or vector of coverage levels in (0, 1).
#'   Default 0.95.
#' @param newdata A \code{data.frame} of covariate values. One region is
#'   produced per row.
#' @param b Integer. Number of boundary points per region. Default 50.
#' @param n_grid Integer. Response grid size for numerical inversion.
#'   Default 2000.
#'
#' @return A named list \code{$mctm_contour} with one element per \code{tau}
#'   level (e.g. \code{"tau=0.95"}), each itself a named list with one
#'   \eqn{b \times 2} numeric matrix per row of \code{newdata}. Columns
#'   correspond to \eqn{(y_1, y_2)}.
#'
#' @details
#' The reference region (Eq. 6, Lado-Baleato et al. 2023) is:
#' \deqn{R_\tau(x) = \tilde{h}^{-1}\{Z : Z^\top \Sigma(x)^{-1} Z
#'   \leq \chi^2_2(\tau)\}}
#' where \eqn{\Sigma(x) = \code{coef(m\_full, type = "Sigma")}} is the
#' model-implied covariance at covariate value \eqn{x}.
#'
#' The ellipse boundary is generated as \eqn{Z = \mathrm{chol}(\Sigma)^\top u}
#' where \eqn{\|u\| = \sqrt{\chi^2_2(\tau)}}, then each point is back-transformed
#' by numerically inverting \code{predict(m_full, margins = j, type = "trafo")}.
#'
#' \strong{Stability:} only documented generics are used (\code{predict},
#' \code{coef}, \code{mkgrid}, \code{variable.names}). No internal slots
#' of the \code{mmlt} object are accessed directly.
#'
#' @seealso \code{\link{tol_mctm}}, \code{\link{mctm_fit}},
#'   \code{\link{mmlt_perc}}
#'
#' @references
#' Lado-Baleato et al. (2023). \doi{10.1002/bimj.202200229}
#'
#' @examples
#' \dontrun{
#' reg <- bivreg_mctm(fit$m_full, fit$m1, fit$m2,
#'                    tau = 0.95, newdata = data.frame(age = 50), b = 100)
#' plot(df$fpg, df$hba1c, col = "grey80", pch = 16)
#' lines(reg$mctm_contour[["tau=0.95"]][[1]], col = "blue", lwd = 2)
#' }
#'
#' @export
bivreg_mctm <- function(m_full, m1, m2, tau = 0.95,
                         newdata = NULL, b = 50L, n_grid = 2000L) {

  .mctm_check_model(m_full, m1, m2)
  stopifnot(is.data.frame(newdata), nrow(newdata) >= 1L)

  gr1 <- .mctm_response_grid(m1, n = n_grid)
  gr2 <- .mctm_response_grid(m2, n = n_grid)

  cont_tau <- vector("list", length(tau))

  for (kk in seq_along(tau)) {
    cont_x <- vector("list", nrow(newdata))

    for (k in seq_len(nrow(newdata))) {
      nd_k  <- newdata[k, , drop = FALSE]
      Sigma <- .mctm_get_sigma(m_full, nd_k)
      Z_ell <- .mctm_make_ellipse(Sigma, tau = tau[kk], b = b)

      t1 <- .mctm_trafo_grid(m_full, gr1$q, gr1$yname, 1L, nd_k)
      t2 <- .mctm_trafo_grid(m_full, gr2$q, gr2$yname, 2L, nd_k)

      yi <- matrix(0, nrow = b, ncol = 2L)
      for (i in seq_len(b)) {
        yi[i, 1L] <- .mctm_invert_trafo(Z_ell[i, 1L], gr1$q, t1)
        yi[i, 2L] <- .mctm_invert_trafo(Z_ell[i, 2L], gr2$q, t2)
      }
      cont_x[[k]] <- yi
    }

    names(cont_x) <- apply(newdata, 1L, function(r)
      paste0(names(newdata), "=", r, collapse = "_"))
    cont_tau[[kk]] <- cont_x
  }

  names(cont_tau) <- paste0("tau=", tau)
  list(mctm_contour = cont_tau)
}
