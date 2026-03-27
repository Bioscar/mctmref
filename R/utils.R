# Internal helpers — all tram/mlt/basefun/variables API calls are isolated here.
# If upstream packages change, only this file needs updating.

# Correct package origins:
#   Bernstein_basis  -> basefun
#   numeric_var      -> variables
#   mkgrid           -> variables
#   variable.names   -> variables
#   mlt, ctm, mmlt   -> mlt
#   predict, coef, logLik -> stats (generics dispatched on mlt/mmlt objects)

.mctm_check_model <- function(m_full, m1, m2) {
  if (!inherits(m_full, "mmlt"))
    stop("'m_full' must be an mmlt object (from mlt::mmlt()).", call. = FALSE)
  if (!inherits(m1, "mlt"))
    stop("'m1' must be an mlt object (from mlt::mlt()).", call. = FALSE)
  if (!inherits(m2, "mlt"))
    stop("'m2' must be an mlt object (from mlt::mlt()).", call. = FALSE)
}

.mctm_response_grid <- function(m_model, n = 2000L) {
  yname <- stats::variable.names(m_model, "response")
  list(q = mkgrid(m_model, n = n)[[yname]], yname = yname)
}

.mctm_get_sigma <- function(m_full, newdata_1row) {
  sigma_obj <- coef(m_full, newdata = newdata_1row, type = "Sigma")
  arr       <- as.array(sigma_obj)
  arr[, , 1]
}

.mctm_get_sigma_all <- function(m_full, newdata) {
  as.array(coef(m_full, newdata = newdata, type = "Sigma"))
}

.mctm_trafo_grid <- function(m_full, q, yname, j, nd_covar) {
  nd <- cbind(
    stats::setNames(data.frame(q), yname),
    nd_covar[rep(1L, length(q)), , drop = FALSE],
    row.names = NULL
  )
  as.numeric(predict(m_full, newdata = nd, margins = j, type = "trafo"))
}

.mctm_make_ellipse <- function(Sigma, tau, b = 50L) {
  r     <- sqrt(stats::qchisq(tau, df = 2L))
  theta <- seq(0, 2 * pi, length.out = b)
  U     <- r * cbind(cos(theta), sin(theta))
  t(t(base::chol(Sigma)) %*% t(U))
}

.mctm_invert_trafo <- function(z_target, q, t_grid) {
  q[which.min((z_target - t_grid)^2)]
}

.mctm_train_data <- function(m_full, m1, m2, train_data = NULL) {
  if (is.null(train_data)) train_data <- m_full$data
  if (is.null(train_data))
    stop("Cannot retrieve training data. Please supply 'train_data'.", call. = FALSE)
  yname1  <- stats::variable.names(m1, "response")
  yname2  <- stats::variable.names(m2, "response")
  cov1    <- setdiff(stats::variable.names(m1), yname1)
  cov2    <- setdiff(stats::variable.names(m2), yname2)
  needed  <- unique(c(yname1, yname2, cov1, cov2))
  needed  <- intersect(needed, names(train_data))
  train_data[, needed, drop = FALSE]
}
