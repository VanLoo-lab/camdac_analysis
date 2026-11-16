#' @title Figure 3A–C / S7D–F: DMP calling on simulated data
#' @author Lili Larose Cadieux
#' @description Figure 3A–C: bi-allelic tumor–normal FNR/FPR and tumor–tumor
#' FNR/FPR heatmaps. S7D/E compare mono- and bi-allelic tumor–normal FNR by
#' copy number and coverage; S7F shows mono-allelic tumor–tumor heatmaps.
#'
#' Source: \code{Notes_DMP_simulation.R} (tumour–tumour + tumour–normal blocks).
#' Regenerates all six panels from
#' \code{~/rdeocd/txMeth/power_sim_results_list.RData}.
#' This follows the final layouts and keeps Figure 3A–C bi-allelic only.

x <- c("data.table", "ggplot2", "optparse", "cowplot", "Rcpp")
invisible(lapply(x, function(y) {
  suppressWarnings(
    suppressMessages(
      library(y, character.only = TRUE, quietly = TRUE,
              warn.conflicts = FALSE)
    )
  )
}))
rm(x)

resolve_path <- function(path, base_dir = NULL) {
  if (grepl("^/", path) || grepl("^~", path)) {
    return(path.expand(path))
  }
  if (is.null(base_dir)) return(path)
  file.path(base_dir, path)
}

# Evan Miller beta-vs-beta probability used in the Notes TN block. Compiling
# the row loop avoids the prohibitively slow mcmapply implementation.
Rcpp::cppFunction('
NumericVector beta_greater_prob(NumericVector an, NumericVector bn,
                                NumericVector ab, NumericVector bb,
                                IntegerVector terms) {
  int nrow = an.size();
  NumericVector out(nrow);
  for (int i = 0; i < nrow; ++i) {
    int n = terms[i];
    double total = 0.0;
    for (int j = 0; j < n; ++j) {
      double z = R::lbeta(an[i] + j, bn[i] + bb[i])
        - std::log(bb[i] + j)
        - R::lbeta(1.0 + j, bb[i])
        - R::lbeta(an[i], bn[i]);
      total += std::exp(z);
    }
    out[i] = std::max(0.0, std::min(1.0, 1.0 - total));
  }
  return out;
}')

#' Load + trim TT simulation table (Notes lines 1297–1415)
prepare_tt_df <- function(results_list) {
  df <- rbindlist(results_list, fill = TRUE)
  df <- df[!is.na(p1) & !is.na(cov_b1) & !is.na(p2) & !is.na(cov_b2)]
  # column aliases used later in Notes
  if (!"DMP_true_t1n" %in% names(df) && "DMP_true_status_t1n" %in% names(df)) {
    df[, DMP_true_t1n := DMP_true_status_t1n]
    df[, DMP_true_t2n := DMP_true_status_t2n]
  }
  if (!"DMP_pure_trimmed" %in% names(df) &&
        "DMP_pure_status_trimmed" %in% names(df)) {
    df[, DMP_pure_trimmed := DMP_pure_status_trimmed]
    df[, DMP_bulk_trimmed := DMP_bulk_status_trimmed]
  }
  df <- df[cov_n1 > 0 & cov_n2 > 0]
  df <- df[!(DMP_true_t1n == "mono-DM" & nA_1 != nB_1)]
  df <- df[!(DMP_true_t2n == "mono-DM" & nA_2 != nB_2)]

  # The Notes code orients mixed-purity pairs with the high-purity sample as
  # sample 1. It does NOT sort by CN; sorting caused the erroneous triangles.
  df[, `:=`(CN_1_new = CN_1, CN_2_new = CN_2, p1_new = p1, p2_new = p2)]
  swap <- df$p2 > 0.58 & df$p1 < 0.3
  df[swap, `:=`(
    CN_1_new = CN_2, CN_2_new = CN_1,
    p1_new = p2, p2_new = p1
  )]
  df[CN_1_new > 6, CN_1_new := 6]
  df[CN_2_new > 6, CN_2_new := 6]
  df[, labels := ifelse(
    p1_new < 0.3 & p2_new < 0.3, "low_vs_low",
    ifelse(
      p1_new > 0.58 & p2_new > 0.58, "high_vs_high",
      ifelse(
        (p1_new < 0.3 & p2_new > 0.58) | (p2_new < 0.3 & p1_new > 0.58),
        "low_vs_high", "other"
      )
    )
  )]
  df[]
}

summarise_tt_rates <- function(df, flag_rows) {
  d <- df[flag_rows]
  one <- function(status, method) {
    ans <- d[, .(
      TP = sum(get(status) == "TP", na.rm = TRUE),
      FN = sum(get(status) == "FN", na.rm = TRUE),
      TN = sum(get(status) == "TN", na.rm = TRUE),
      FP = sum(get(status) == "FP", na.rm = TRUE)
    ), by = .(CN_1_new, CN_2_new)]
    ans[, Method := method]
    ans
  }
  out <- rbindlist(list(
    one("DMP_pure_trimmed", "CAMDAC"),
    one("DMP_bulk_trimmed", "bulk")
  ))
  out[, FNR := 1 - ((TP + 1) / (TP + FN + 1))]
  out[, FPR := 1 - ((TN + 1) / (TN + FP + 1))]
  out[]
}

calculate_tt_rates <- function(df, flag_name = c("bi", "mono", "all")) {
  flag_name <- match.arg(flag_name)
  strata <- c("low_vs_low", "low_vs_high", "high_vs_high")
  rates_by_stratum <- lapply(strata, function(label) {
    d <- df[labels == label]
    if (flag_name == "bi") {
      flag <- d[, which(
        DMP_true_t1n != "mono-DM" & DMP_true_t2n != "mono-DM"
      )]
    } else if (flag_name == "mono") {
      flag <- d[, which(
        DMP_true_t1n != "bi-DM" & DMP_true_t2n != "bi-DM"
      )]
    } else {
      flag <- seq_len(nrow(d))
    }
    rates <- summarise_tt_rates(d, flag)
    rates[, `:=`(
      Purity = label,
      DM_class = ifelse(flag_name == "bi", "bi-allelic",
                        ifelse(flag_name == "mono", "mono-allelic", "all")),
      Comparison = "tumor_tumor"
    )]
    rates
  })
  rbindlist(rates_by_stratum)
}

prepare_tn_observations <- function(results_list) {
  src <- rbindlist(results_list, fill = TRUE)
  aliases <- list(
    DMP_true_t1n = "DMP_true_status_t1n",
    DMP_true_t2n = "DMP_true_status_t2n",
    DMP_pure = "DMP_pure_status",
    DMP_bulk = "DMP_bulk_status",
    DMP_pure_trimmed = "DMP_pure_status_trimmed",
    DMP_bulk_trimmed = "DMP_bulk_status_trimmed"
  )
  for (new in names(aliases)) {
    old <- aliases[[new]]
    if (!new %in% names(src) && old %in% names(src)) setnames(src, old, new)
  }

  # Reproduce the QC applied to `df` immediately before the Notes TN block.
  src <- src[
    !is.na(p1) & !is.na(cov_b1) & !is.na(p2) & !is.na(cov_b2) &
      cov_n1 > 0 & cov_n2 > 0
  ]
  src <- src[!(DMP_true_t1n == "mono-DM" & nA_1 != nB_1)]
  src <- src[!(DMP_true_t2n == "mono-DM" & nA_2 != nB_2)]
  src <- src[
    !is.na(m_t1) & !is.na(m_t2) & !is.na(m_b1) & !is.na(m_b2)
  ]
  status_columns <- c(
    "DMP_pure", "DMP_bulk", "DMP_pure_trimmed", "DMP_bulk_trimmed"
  )
  if (all(status_columns %in% names(src))) {
    src <- src[
      !is.na(DMP_pure) & !is.na(DMP_bulk) &
        !is.na(DMP_pure_trimmed) & !is.na(DMP_bulk_trimmed)
    ]
  }

  cols <- c("p", "CN", "m_n", "m_t", "m_b", "cov_n", "cov_t", "cov_b",
            "DMP_true_tn")
  extract_side <- function(side) {
    old <- c(
      paste0("p", side), paste0("CN_", side),
      paste0("m_n", side), paste0("m_t", side), paste0("m_b", side),
      paste0("cov_n", side), paste0("cov_t", side), paste0("cov_b", side),
      paste0("DMP_true_t", side, "n")
    )
    absent <- setdiff(old, names(src))
    if (length(absent)) stop("Missing TN simulation columns: ", paste(absent, collapse = ", "))
    z <- copy(src[, ..old])
    setnames(z, old, cols)
    z
  }
  dt <- rbindlist(list(extract_side("1"), extract_side("2")))
  dt <- dt[
    !is.na(m_n) & !is.na(m_t) & !is.na(m_b) &
      !is.na(cov_n) & !is.na(cov_b) & !is.na(DMP_true_tn)
  ]
  an <- round(dt$m_n * dt$cov_n) + 0.5
  bn <- round(dt$cov_n - dt$m_n * dt$cov_n) + 0.5
  ab <- round(dt$m_b * dt$cov_b) + 0.5
  bb <- round(dt$cov_b - dt$m_b * dt$cov_b) + 0.5
  # Keep R's ties-to-even rounding used by the original h() implementation.
  prob <- beta_greater_prob(an, bn, ab, bb, as.integer(round(ab)))
  significant <- prob > 0.99 | prob < 0.01
  truth <- dt$DMP_true_tn != "no"

  dt[, `:=`(
    CN = pmin(CN, 6),
    Purity = ifelse(p < 0.3, "low", "high"),
    true_dmp = truth,
    CAMDAC_call = significant & abs(m_t - m_n) > 0.2,
    bulk_call = significant & abs(m_b - m_n) > 0.2
  )]
  dt
}

summarise_tn_rates <- function(dt, flag_name = "bi",
                               group_by = c("CN", "coverage")) {
  group_by <- match.arg(group_by)
  d <- copy(dt)
  if (flag_name == "bi") d <- d[DMP_true_tn != "mono-DM"]
  if (flag_name == "mono") d <- d[DMP_true_tn != "bi-DM"]
  if (group_by == "coverage") {
    d[, Coverage_interval := cut(
      cov_b, breaks = c(0, 25, 50, 75, 100, 125, 150, Inf),
      include.lowest = TRUE
    )]
    groups <- c("Purity", "Coverage_interval")
  } else {
    groups <- c("Purity", "CN")
  }
  rate_one <- function(call_col, method) {
    z <- d[, .(
      TP = sum(true_dmp & get(call_col)),
      FN = sum(true_dmp & !get(call_col)),
      TN = sum(!true_dmp & !get(call_col)),
      FP = sum(!true_dmp & get(call_col))
    ), by = groups]
    z[, `:=`(
      Method = method,
      FNR = 1 - ((TP + 1) / (TP + FN + 1)),
      FPR = 1 - ((TN + 1) / (TN + FP + 1))
    )]
    z
  }
  rbindlist(list(
    rate_one("CAMDAC_call", "CAMDAC"),
    rate_one("bulk_call", "bulk")
  ))
}

prepare_tn_rates <- function(results_list, flag_name = "bi",
                             group_by = c("CN", "coverage")) {
  observations <- prepare_tn_observations(results_list)
  summarise_tn_rates(observations, flag_name, match.arg(group_by))
}

plot_tn <- function(rates, metric = c("FNR", "FPR")) {
  metric <- match.arg(metric)
  d <- copy(rates)
  if (metric == "FPR") d[, value := FPR * 1000] else d[, value := FNR]
  y_limits <- if (metric == "FNR") c(-0.025, 1) else c(-0.075, 3)
  ggplot(d, aes(CN, value, colour = Purity)) +
    geom_line(
      data = d[Method == "CAMDAC"],
      linewidth = 0.55, linetype = "solid",
      lineend = "butt", linejoin = "mitre"
    ) +
    geom_line(
      data = d[Method == "bulk"],
      linewidth = 0.8, linetype = "dashed",
      lineend = "butt", linejoin = "mitre"
    ) +
    scale_colour_manual(values = c(low = "grey70", high = "grey15")) +
    scale_x_continuous(
      name = "copy number", breaks = 1:6, labels = c(1:5, ">=6")
    ) +
    scale_y_continuous(
      name = if (metric == "FNR") "FNR" else expression(FPR~(10^-3)),
      breaks = if (metric == "FNR") seq(0, 1, 0.2) else 0:3,
      expand = expansion(mult = c(0, 0.03))
    ) +
    coord_cartesian(xlim = c(1, 6), ylim = y_limits) +
    theme_classic(base_size = 9) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 7)
    )
}

plot_s7_fnr <- function(rates, x_var = c("CN", "Coverage_interval"),
                        facet_dm = FALSE) {
  x_var <- match.arg(x_var)
  d <- copy(rates)
  d[, DM_class := factor(DM_class, levels = c("mono-allelic", "bi-allelic"))]
  if (x_var == "Coverage_interval") {
    d[, Coverage_interval := factor(
      as.character(Coverage_interval),
      levels = c(
        "[0,25]", "(25,50]", "(50,75]", "(75,100]",
        "(100,125]", "(125,150]", "(150,Inf]"
      )
    )]
  }
  x_aes <- if (x_var == "CN") quote(CN) else quote(Coverage_interval)
  grouping <- interaction(d$Purity, d$DM_class, drop = TRUE)
  d[, plot_group := grouping]

  p <- ggplot(
    d,
    aes(
      x = !!x_aes, y = FNR, colour = Purity,
      shape = DM_class, group = plot_group
    )
  ) +
    geom_line(
      data = d[Method == "CAMDAC"],
      aes(linetype = Method, linewidth = Method),
      lineend = "butt", linejoin = "mitre"
    ) +
    geom_line(
      data = d[Method == "bulk"],
      aes(linetype = Method, linewidth = Method),
      lineend = "butt", linejoin = "mitre"
    ) +
    geom_point(size = 1.7) +
    scale_colour_manual(
      name = NULL, values = c(low = "#f28e00", high = "#ffc04c"),
      labels = c(low = "low purity", high = "high purity")
    ) +
    scale_shape_manual(
      name = NULL, values = c("mono-allelic" = 17, "bi-allelic" = 16)
    ) +
    scale_linetype_manual(
      name = NULL, values = c(CAMDAC = "solid", bulk = "dashed")
    ) +
    scale_linewidth_manual(
      name = NULL, values = c(CAMDAC = 0.55, bulk = 0.8)
    ) +
    scale_y_continuous(
      "FNR", breaks = seq(0, 1, 0.2),
      expand = expansion(mult = c(0, 0.02))
    ) +
    coord_cartesian(ylim = c(-0.025, 1)) +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 6)
    )
  if (x_var == "CN") {
    p <- p + scale_x_continuous(
      "copy number", breaks = 1:6, labels = c(1:5, ">=6")
    )
  } else {
    p <- p + labs(x = "coverage interval") +
      theme(axis.text.x = element_text(angle = 55, hjust = 1, size = 6))
  }
  if (facet_dm) {
    p <- p + facet_wrap(
      ~DM_class, nrow = 1,
      labeller = as_labeller(c(
        "mono-allelic" = "mono-allelic",
        "bi-allelic" = "bi-allelic"
      ))
    )
  }
  p
}

plot_tt_composite <- function(rates, metric = c("FNR", "FPR")) {
  metric <- match.arg(metric)
  d <- copy(rates)
  mono_only <- "DM_class" %in% names(d) &&
    all(d$DM_class == "mono-allelic")
  if (mono_only) {
    d <- d[CN_1_new %in% c(2, 4, 6) & CN_2_new %in% c(2, 4, 6)]
    d[, `:=`(
      CN_1_plot = match(CN_1_new, c(2, 4, 6)),
      CN_2_plot = match(CN_2_new, c(2, 4, 6))
    )]
    axis_breaks <- 1:3
    axis_labels <- c("2", "4", ">=6")
  } else {
    d[, `:=`(CN_1_plot = CN_1_new, CN_2_plot = CN_2_new)]
    axis_breaks <- 1:6
    axis_labels <- c(1:5, ">=6")
  }
  d[, Purity := factor(
    Purity,
    levels = c("low_vs_low", "low_vs_high", "high_vs_high"),
    labels = c(
      "low purity sample pair", "low vs high purity sample pair",
      "high purity sample pair"
    )
  )]
  d[, Method := factor(Method, levels = c("CAMDAC", "bulk"))]
  limits <- if (metric == "FNR") c(0, 1) else c(0, 0.45)
  colours <- if (metric == "FNR") {
    c("#add8e6", "#ee7600")
  } else {
    c("#add8e6", "#ee2c2c")
  }
  ggplot(d, aes(CN_1_plot, CN_2_plot, fill = .data[[metric]])) +
    geom_tile(width = 0.92, height = 0.92) +
    facet_grid(Method ~ Purity, switch = "y") +
    scale_fill_gradient(
      name = metric, low = colours[1], high = colours[2], limits = limits,
      na.value = "white"
    ) +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels) +
    scale_y_continuous(breaks = axis_breaks, labels = axis_labels) +
    labs(x = "CpG copy number sample x", y = NULL) +
    coord_fixed() +
    theme_classic(base_size = 8) +
    theme(
      legend.position = "none",
      axis.title.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 90, face = "bold", size = 9),
      strip.text.x = element_text(size = 6.2),
      panel.spacing = grid::unit(1.5, "mm")
    )
}

plot_tt_cell <- function(rates, purity, method, metric,
                         show_y = FALSE, show_x = TRUE) {
  d <- rates[Purity == purity & Method == method]
  mono_only <- "DM_class" %in% names(d) &&
    all(d$DM_class == "mono-allelic")
  if (mono_only) {
    d <- d[CN_1_new %in% c(2, 4, 6) & CN_2_new %in% c(2, 4, 6)]
    d[, `:=`(
      CN_1_plot = match(CN_1_new, c(2, 4, 6)),
      CN_2_plot = match(CN_2_new, c(2, 4, 6))
    )]
    axis_breaks <- 1:3
    axis_labels <- c("2", "4", ">=6")
  } else {
    d[, `:=`(CN_1_plot = CN_1_new, CN_2_plot = CN_2_new)]
    axis_breaks <- 1:6
    axis_labels <- c(1:5, ">=6")
  }
  limits <- if (metric == "FNR") c(0, 1) else c(0, 0.45)
  colours <- if (metric == "FNR") {
    c("#add8e6", "#ee7600")
  } else {
    c("#add8e6", "#ee2c2c")
  }
  p <- ggplot(d, aes(CN_1_plot, CN_2_plot, fill = .data[[metric]])) +
    geom_tile(width = 0.92, height = 0.92) +
    scale_fill_gradient(
      low = colours[1], high = colours[2], limits = limits,
      na.value = "white"
    ) +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels) +
    scale_y_continuous(breaks = axis_breaks, labels = axis_labels) +
    coord_fixed() +
    theme_classic(base_size = 7) +
    theme(
      legend.position = "none",
      axis.title = element_blank(),
      axis.text.y = if (show_y) element_text(size = 6) else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      axis.text.x = if (show_x) element_text(size = 6) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank()
    )
  if (show_y) p <- p + labs(y = method)
  p
}

plot_tt_fnr_fpr <- function(rates) {
  strata <- c("low_vs_low", "low_vs_high", "high_vs_high")
  plots <- list()
  for (method in c("CAMDAC", "bulk")) {
    for (purity in strata) {
      for (metric in c("FNR", "FPR")) {
        plots[[length(plots) + 1]] <- plot_tt_cell(
          rates, purity, method, metric,
          show_y = purity == strata[1] && metric == "FNR",
          show_x = method == "bulk"
        )
      }
    }
  }
  body <- cowplot::plot_grid(plotlist = plots, ncol = 6, align = "hv")
  headings <- cowplot::ggdraw()
  heading_text <- c(
    "low purity sample pair", "low vs high purity sample pair",
    "high purity sample pair"
  )
  for (i in seq_along(heading_text)) {
    headings <- headings + cowplot::draw_label(
      heading_text[i], x = (2 * i - 1) / 6, size = 7
    )
  }
  metric_headings <- cowplot::ggdraw()
  for (i in 0:5) {
    metric_headings <- metric_headings + cowplot::draw_label(
      if (i %% 2 == 0) "FNR" else "FPR", x = (i + 0.5) / 6, size = 6
    )
  }
  x_title <- cowplot::ggdraw() + cowplot::draw_label(
    "CpG copy number sample x", size = 8
  )
  cowplot::plot_grid(
    headings, metric_headings, body, x_title,
    ncol = 1, rel_heights = c(0.10, 0.08, 1, 0.08)
  )
}

################################################################################
# CLI
################################################################################

if (sys.nframe() == 0L) {

option_list <- list(
  make_option(c("--work_dir"), type = "character",
              default = "~/data/CAMDAC_paper/figure3A-C",
              help = "Figure working directory [default %default]"),
  make_option(c("--output_dir"), type = "character",
              default = "outputs",
              help = "Output directory [default %default]"),
  make_option(c("--sim_rdata"), type = "character",
              default = "~/rdeocd/txMeth/power_sim_results_list.RData",
              help = "Cached TT simulation results [default %default]"),
  make_option(c("--panels"), type = "character",
              default = "A,B,C",
              help = "Comma-separated panels: A,B,C,D,E,F [default %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

work_dir <- path.expand(opt$work_dir)
if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
setwd(work_dir)
output_dir <- resolve_path(opt$output_dir, work_dir)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

panels <- toupper(trimws(unlist(strsplit(opt$panels, ",", fixed = TRUE))))
if (any(panels %in% c("A", "B", "C", "D", "E", "F"))) {
  rdata <- path.expand(opt$sim_rdata)
  if (!file.exists(rdata)) stop("Missing: ", rdata)
  message("Loading ", rdata)
  e <- new.env(parent = emptyenv())
  load(rdata, envir = e)
  if (!exists("results_list", envir = e)) {
    stop("results_list not in ", rdata)
  }
}

p_a <- p_b <- p_c <- p_c_fnr <- p_c_fpr <- p_d <- p_e <- p_f <- NULL
# Main Figure 3A-C is bi-allelic. S7D/E compare both classes; S7F is mono.
flag_name <- "bi"

if (any(panels %in% c("A", "B", "D", "E"))) {
  message("Calculating tumor-normal rates")
  tn_observations <- prepare_tn_observations(e$results_list)
  tn_rates <- summarise_tn_rates(tn_observations, "bi", "CN")
}

if (any(panels %in% c("A", "B"))) {
  fwrite(tn_rates, file.path(output_dir, "Fig3A-B_TN_rates.tsv"), sep = "\t")
  p_a <- plot_tn(tn_rates, "FNR")
  p_b <- plot_tn(tn_rates, "FPR")
  if ("A" %in% panels) {
    ggsave(file.path(output_dir, "Fig3A_TN_FNR_vs_copy_number.pdf"),
           p_a, width = 2.2, height = 2.3, units = "in")
  }
  if ("B" %in% panels) {
    ggsave(file.path(output_dir, "Fig3B_TN_FPR_vs_copy_number.pdf"),
           p_b, width = 2.2, height = 2.3, units = "in")
  }
}

if ("D" %in% panels) {
  s7d_rates <- rbindlist(list(
    copy(tn_rates)[, DM_class := "bi-allelic"],
    summarise_tn_rates(tn_observations, "mono", "CN")[
      , DM_class := "mono-allelic"
    ]
  ))
  fwrite(s7d_rates, file.path(output_dir, "FigS7D_TN_FNR_CN_rates.tsv"), sep = "\t")
  p_d <- plot_s7_fnr(s7d_rates, "CN", facet_dm = TRUE)
  ggsave(file.path(output_dir, "FigS7D_TN_FNR_vs_copy_number.pdf"),
         p_d, width = 7, height = 3, units = "in")
}

if ("E" %in% panels) {
  s7e_rates <- rbindlist(list(
    summarise_tn_rates(tn_observations, "bi", "coverage")[
      , DM_class := "bi-allelic"
    ],
    summarise_tn_rates(tn_observations, "mono", "coverage")[
      , DM_class := "mono-allelic"
    ]
  ))
  fwrite(
    s7e_rates,
    file.path(output_dir, "FigS7E_TN_FNR_coverage_interval_rates.tsv"),
    sep = "\t"
  )
  p_e <- plot_s7_fnr(s7e_rates, "Coverage_interval", facet_dm = TRUE)
  ggsave(file.path(output_dir, "FigS7E_TN_FNR_vs_coverage_interval.pdf"),
         p_e, width = 7, height = 3, units = "in")
}

if (any(panels %in% c("C", "F"))) {
  df <- prepare_tt_df(e$results_list)
  message("TT rows after trim: ", nrow(df))
}

if ("C" %in% panels) {
  tt_rates <- calculate_tt_rates(df, "bi")
  fwrite(tt_rates, file.path(output_dir, "Fig3C_TT_FNR_rates.tsv"), sep = "\t")
  fwrite(tt_rates, file.path(output_dir, "Fig3C_TT_FNR_FPR_rates.tsv"), sep = "\t")
  p_c_fnr <- plot_tt_composite(tt_rates, "FNR")
  p_c_fpr <- plot_tt_composite(tt_rates, "FPR")
  p_c <- plot_tt_fnr_fpr(tt_rates)
  ggsave(file.path(output_dir, "Fig3C_TT_FNR_heatmap_2x3.pdf"),
         p_c_fnr, width = 7.2, height = 3.5, units = "in")
  ggsave(file.path(output_dir, "Fig3C_TT_FPR_heatmap_2x3.pdf"),
         p_c_fpr, width = 7.2, height = 3.5, units = "in")
  ggsave(file.path(output_dir, "Fig3C_TT_FNR_FPR_heatmap_2x6.pdf"),
         p_c, width = 9, height = 3.5, units = "in")
}

if ("F" %in% panels) {
  s7f_rates <- calculate_tt_rates(df, "mono")
  fwrite(
    s7f_rates,
    file.path(output_dir, "FigS7F_TT_mono_FNR_FPR_rates.tsv"),
    sep = "\t"
  )
  p_f_fnr <- plot_tt_composite(s7f_rates, "FNR")
  p_f_fpr <- plot_tt_composite(s7f_rates, "FPR")
  p_f <- plot_tt_fnr_fpr(s7f_rates)
  ggsave(file.path(output_dir, "FigS7F_TT_mono_FNR_heatmap_2x3.pdf"),
         p_f_fnr, width = 7.2, height = 3.5, units = "in")
  ggsave(file.path(output_dir, "FigS7F_TT_mono_FPR_heatmap_2x3.pdf"),
         p_f_fpr, width = 7.2, height = 3.5, units = "in")
  ggsave(file.path(output_dir, "FigS7F_TT_mono_FNR_FPR_heatmap_2x6.pdf"),
         p_f, width = 9, height = 3.5, units = "in")
}

if (all(c("D", "E", "F") %in% panels)) {
  s7_tn <- cowplot::plot_grid(
    p_d, p_e, labels = c("D", "E"), nrow = 1,
    rel_widths = c(1, 1)
  )
  s7_tn <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(
      "Tumor vs. normal DMP (simulated)", x = 0, hjust = 0, size = 10
    ),
    s7_tn, ncol = 1, rel_heights = c(0.12, 1)
  )
  s7_tt <- cowplot::plot_grid(p_f, labels = "F", nrow = 1)
  s7_tt <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(
      "Tumor vs. tumor DMP (simulated)", x = 0, hjust = 0, size = 10
    ),
    s7_tt, ncol = 1, rel_heights = c(0.12, 1)
  )
  s7 <- cowplot::plot_grid(
    s7_tn, s7_tt, nrow = 1, rel_widths = c(3.5, 4.5)
  )
  ggsave(file.path(output_dir, "FigS7D-F_simulated_DMP.pdf"),
         s7, width = 16, height = 3.5, units = "in")
}

if (all(c("A", "B", "C") %in% panels)) {
  tn <- cowplot::plot_grid(p_a, p_b, labels = c("A", "B"), nrow = 1)
  tn <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(
      "Tumor vs. normal DMP (simulated)", x = 0, hjust = 0, size = 10
    ),
    tn, ncol = 1, rel_heights = c(0.12, 1)
  )
  tt <- cowplot::plot_grid(p_c, labels = "C", nrow = 1)
  tt <- cowplot::plot_grid(
    cowplot::ggdraw() + cowplot::draw_label(
      "Tumor vs. tumor DMP (simulated)", x = 0, hjust = 0, size = 10
    ),
    tt, ncol = 1, rel_heights = c(0.12, 1)
  )
  combined <- cowplot::plot_grid(tn, tt, nrow = 1, rel_widths = c(2, 6.2))
  ggsave(file.path(output_dir, "Fig3A-C_simulated_DMP.pdf"),
         combined, width = 13, height = 3.3, units = "in")
}

message("Done.")

} # end main
