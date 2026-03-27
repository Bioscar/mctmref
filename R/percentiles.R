#' Marginal Conditional Percentile Curves from an MCTM
#'
#' Numerically inverts the marginal CDF \eqn{F_j(y_j \mid x)} from the joint
#' \code{mmlt} model to obtain conditional quantiles \eqn{Q_\tau(Y_j \mid x)}.
#'
#' @param m_full A fitted \code{mmlt} object.
#' @param m1 Fitted marginal \code{mlt} model for the first response.
#' @param m2 Fitted marginal \code{mlt} model for the second response.
#' @param marginal Integer, 1 or 2. Which marginal to compute.
#' @param newdata A \code{data.frame} of covariate values. Each row gives one
#'   covariate value at which percentile curves are evaluated.
#' @param tau Numeric vector of probability levels in (0, 1).
#' @param plot Logical. If \code{TRUE}, adds lines to the current plot.
#'   If \code{FALSE} (default), returns a \code{data.frame}.
#' @param n_grid Integer. Response grid size for numerical inversion.
#'   Default 10000.
#' @param ... Additional graphical arguments passed to \code{lines()} when
#'   \code{plot = TRUE}.
#'
#' @return If \code{plot = FALSE}: a \code{data.frame} with the covariate
#'   column(s) from \code{newdata} followed by one column per \code{tau} level.
#'
#' @details
#' Uses \code{predict(m_full, margins = j, type = "distribution")} on a dense
#' response grid, then finds the value closest to each requested \eqn{\tau}.
#'
#' @seealso \code{\link{mctm_fit}}, \code{\link{bivreg_mctm}}
#'
#' @references
#' Lado-Baleato et al. (2023). \doi{10.1002/bimj.202200229}
#'
#' @examples
#' \dontrun{
#' nd  <- data.frame(age = seq(20, 85, length.out = 50))
#' pct <- mmlt_perc(fit$m_full, fit$m1, fit$m2,
#'                  marginal = 1, newdata = nd,
#'                  tau = c(0.025, 0.5, 0.975))
#' plot(df$age, df$fpg, col = "grey80", pch = 16)
#' lines(pct$age, pct[["0.5"]], lwd = 2)
#' }
#' @importFrom variables mkgrid
#' @export
mmlt_perc <- function(m_full, m1, m2, marginal, newdata, tau,
                      plot = FALSE, n_grid = 10000L, ...) {

  .mctm_check_model(m_full, m1, m2)
  marginal <- as.integer(marginal)
  stopifnot(marginal %in% 1:2)

  marg_model <- if (marginal == 1L) m1 else m2
  gr         <- .mctm_response_grid(marg_model, n = n_grid)

  cond_perc <- NULL
  for (jj in seq_len(nrow(newdata))) {
    nd_grid <- cbind(
      stats::setNames(data.frame(gr$q), gr$yname),
      newdata[rep(jj, length(gr$q)), , drop = FALSE],
      row.names = NULL
    )
    d_y  <- as.numeric(predict(m_full, newdata = nd_grid,
                                margins = marginal, type = "distribution"))
    perc <- sapply(tau, function(t) gr$q[which.min(abs(d_y - t))])
    cond_perc <- rbind(cond_perc, perc)
  }

  colnames(cond_perc) <- as.character(tau)
  out <- cbind(newdata, as.data.frame(cond_perc))

  if (plot) {
    for (k in seq_along(tau))
      graphics::lines(newdata[[1L]], out[[ncol(newdata) + k]], ...)
    invisible(out)
  } else {
    return(out)
  }
}
