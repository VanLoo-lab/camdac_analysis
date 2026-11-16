#' Shared helpers for ASM DMP classification and Figure 4 plotting.
#' Extracted from PhD Notes_ASM_plots_*.R and figure4/camdac_utils.R.

round2 <- function(x, digits = 0) {
  round(x, digits = digits)
}

intervalWidth <- function(lowTailPr, ICDFname, credMass, ...) {
  ICDFname(credMass + lowTailPr, ...) - ICDFname(lowTailPr, ...)
}

HDIofICDF <- function(ICDFname, credMass = 0.99, tol = 1e-4, ...) {
  incredMass <- 1.0 - credMass
  optInfo <- optimize(
    f = intervalWidth,
    interval = c(0, incredMass),
    ICDFname = ICDFname,
    credMass = credMass,
    tol = tol,
    ...
  )
  HDIlowTailPr <- optInfo$minimum
  setNames(
    ICDFname(c(HDIlowTailPr, credMass + HDIlowTailPr), ...),
    nm = c("low", "high")
  )
}

MCMC_overlap_beta <- function(M_1, UM_1, M_2, UM_2) {
  ov <- as.numeric(NA)
  if (sum(is.na(c(M_1, UM_1, M_2, UM_2))) == 0) {
    if (sum(c(M_1, UM_1)) > 0 && sum(c(M_2, UM_2)) > 0) {
      sampleVec1 <- rbeta(n = 2000, shape1 = M_1, shape2 = UM_1)
      sampleVec2 <- rbeta(n = 2000, shape1 = M_2, shape2 = UM_2)
      sortedPts_1 <- sort(sampleVec1)
      m_1 <- mean(sortedPts_1, na.rm = TRUE)
      sortedPts_2 <- sort(sampleVec2)
      m_2 <- mean(sortedPts_2, na.rm = TRUE)
      if (min(sortedPts_1) > max(sortedPts_2) || min(sortedPts_2) > max(sortedPts_1)) {
        ov <- 0
      } else if (m_1 > m_2) {
        ov <- length(which(sortedPts_1 < max(sortedPts_2))) / 2000
      } else if (m_1 < m_2) {
        ov <- length(which(sortedPts_2 < max(sortedPts_1))) / 2000
      }
    }
  }
  ov
}

#' Evan Miller beta-greater probability (normal vs bulk tumour).
h <- function(alpha_n, beta_n, alpha_b, beta_b) {
  j <- seq.int(0, round(alpha_b) - 1)
  log_vals <- (
    lbeta(alpha_n + j, beta_n + beta_b) - log(beta_b + j) -
      lbeta(1 + j, beta_b) - lbeta(alpha_n, beta_n)
  )
  1 - sum(exp(log_vals))
}

g_legend <- function(a.gplot) {
  tmp <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  tmp$grobs[[leg]]
}
