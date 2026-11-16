#' @description Intra-tumour comparison of pure tumour (mt-mt) versus bulk (mb-mb)
#' differentially methylated positions (DMPs) for TRACERx patient CRUK0062.
#' Overlap statistics and Venn-style figures follow the logic in
#' \code{Notes_mt_mt_vs_mb_mb_DMPs_figure5FG.R}.
#'
#' @author Lili Larose Cadieux

################################################################################
# Load libraries
################################################################################

## Prefer CRUKID in defaults; resolve to Patient (LTX*) via clinical table for paths
.shared <- file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..", "_shared", "resolve_tx_ids.R")
.shared <- normalizePath(.shared, mustWork = FALSE)
if (!file.exists(.shared)) {
  .shared <- normalizePath(file.path("..", "..", "_shared", "resolve_tx_ids.R"), mustWork = FALSE)
}
if (!exists("resolve_tx_patient_id", mode = "function")) {
  if (file.exists(.shared)) source(.shared) else stop("Missing TRACERx/_shared/resolve_tx_ids.R")
}
.tx_map <- load_tx_patient_map()
x <- c(
  "DescTools", "eulerr", "ggplot2", "data.table", "limma", "dplyr", "ggforce"
)

invisible(lapply(x, function(y) {
  suppressWarnings(
    suppressMessages(
      library(y, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    )
  )
}))
rm(x)

################################################################################
# Paths and identifiers — CRUK0062 only
################################################################################

# Internal TRACERx ID and published study identifier
patient_id <- resolve_tx_patient_id("CRUK0062", .tx_map)
cruk_id <- "CRUK0062"

parent_dir <- path.expand("~/rdeocd/txMeth/")
out_dir <- path.expand("~/data/CAMDAC_paper/figure3DE/outputs")

# Patient metadata (RRBS clinical / purity / etc.)
f_nm <- file.path(
  parent_dir,
  "tx100clinicalinfo_RRBS_published_status_immune_clusters_purity_ploidy_CpG_cov.txt"
)
patient_info <- fread(f_nm, sep = "\t")

################################################################################
# Run parameters
################################################################################

reseg <- TRUE
effect_size <- 0.2

patient_id1 <- patient_id
patient_id2 <- patient_id

path1 <- file.path(parent_dir, patient_id1)
path2 <- file.path(parent_dir, patient_id2)

sample_ids1 <- unique(patient_info[Patient == patient_id1, Sample])
sample_ids2 <- unique(patient_info[Patient == patient_id2, Sample])

cols <- c(
  "chromosome", "start", "end",
  "m_t_low", "m_t_high", "m_t_corr",
  "m_b_low", "m_b_high", "m_b"
)

################################################################################
# Load and merge per-region DMP calls (all samples for CRUK0062)
################################################################################

df_list1 <- NULL
for (i in seq_along(sample_ids1)) {
  sample_id <- sample_ids1[i]
  if (reseg) {
    f_nm <- file.path(
      path1, "Methylation", sample_id, "re_segmented",
      "purified_tumour_DMP_calls.RData"
    )
  } else {
    f_nm <- file.path(
      path1, "Methylation", sample_id, "purified_tumour_DMP_calls.RData"
    )
  }
  load(f_nm)
  df_mt <- data.table(df_mt[, cols])
  colnames(df_mt)[4:length(cols)] <-
    paste(cols[4:length(cols)], patient_id1, sample_id, sep = "_")
  df_list1[[length(df_list1) + 1]] <- df_mt
  if (length(df_list1) > 1) {
    df1 <- Reduce(
      function(...) merge(..., by = c("chromosome", "start", "end"), all = TRUE),
      df_list1
    )
    df_list1 <- list(df1)
    rm(df1)
  }
  rm(df_mt)
}
df1 <- df_list1
rm(df_list1)

# Same patient: reuse merged tumour table
df2 <- df1

df <- Reduce(
  function(...) merge(..., by = c("chromosome", "start", "end"), all = TRUE),
  list(df1, df2)
)
dt <- data.table(df)
DMPs_mt <- as.data.table(df[, 1:3])
rm(df1, df2)

################################################################################
# Pure tumour (mt): corrected DMP calls per sample pair (upper triangle)
################################################################################

for (sample_id1 in sample_ids1) {
  for (sample_id2 in sample_ids2) {
    if (patient_id1 != patient_id2 ||
      (patient_id1 == patient_id2 &&
        which(sample_ids2 == sample_id2) > which(sample_ids1 == sample_id1))) {
      cat(patient_id1, "-", sample_id1, " vs. ",
        patient_id2, "-", sample_id2, "\n")
      raw <- rep(NA, nrow(DMPs_mt))
      cols1 <- paste(c("m_t_low", "m_t_high"), patient_id1, sample_id1, sep = "_")
      cols2 <- paste(c("m_t_low", "m_t_high"), patient_id2, sample_id2, sep = "_")

      if (patient_id1 == patient_id2) {
        cols1 <- paste(cols1, "x", sep = ".")
        cols2 <- paste(cols2, "y", sep = ".")
      }
      tmp <- rowSums(dt[, is.na(.SD), .SDcols = c(cols1[1], cols2[1])])
      idx <- c(1:nrow(df))[tmp == 0]
      rm(tmp)

      raw[idx] <- Overlap(df[idx, cols1], df[idx, cols2]) == 0
      if (patient_id1 != patient_id2) {
        diff <- df[, paste(c("m_t_corr"), patient_id1, sample_id1, sep = "_")] -
          df[, paste(c("m_t_corr"), patient_id2, sample_id2, sep = "_")]
      } else {
        diff <- df[, paste0("m_t_corr_", patient_id1, "_", sample_id1, ".x")] -
          df[, paste0("m_t_corr_", patient_id2, "_", sample_id2, ".y")]
      }
      corr <- abs(diff) >= effect_size & raw == TRUE

      cols_out <- c(
        paste("corrected", patient_id1, sample_id1, patient_id2, sample_id2, sep = "_")
      )
      DMPs_mt[, as.character(cols_out) := list(NA)]
      DMPs_mt[, as.character(cols_out) := as.list(corr)]
    }
  }
}
rm(raw, diff, corr)

################################################################################
# Bulk (mb): corrected DMP calls per sample pair (upper triangle)
################################################################################

DMPs_mb <- as.data.table(df[, 1:3])
for (sample_id1 in sample_ids1) {
  for (sample_id2 in sample_ids2) {
    if (patient_id1 != patient_id2 ||
      (patient_id1 == patient_id2 &&
        which(sample_ids2 == sample_id2) > which(sample_ids1 == sample_id1))) {
      cat(patient_id1, "-", sample_id1, " vs. ",
        patient_id2, "-", sample_id2, "\n")
      raw <- rep(NA, nrow(DMPs_mb))
      cols1 <- paste(c("m_b_low", "m_b_high"), patient_id1, sample_id1, sep = "_")
      cols2 <- paste(c("m_b_low", "m_b_high"), patient_id2, sample_id2, sep = "_")
      if (patient_id1 == patient_id2) {
        cols1 <- paste(cols1, "x", sep = ".")
        cols2 <- paste(cols2, "y", sep = ".")
      }

      tmp <- rowSums(dt[, is.na(.SD), .SDcols = c(cols1[1], cols2[1])])
      idx <- c(1:nrow(df))[tmp == 0]
      rm(tmp)

      raw[idx] <- Overlap(df[idx, cols1], df[idx, cols2]) == 0
      if (patient_id1 != patient_id2) {
        diff <- df[, paste(c("m_b"), patient_id1, sample_id1, sep = "_")] -
          df[, paste(c("m_b"), patient_id2, sample_id2, sep = "_")]
      } else {
        diff <- df[, paste0("m_b_", patient_id1, "_", sample_id1, ".x")] -
          df[, paste0("m_b_", patient_id2, "_", sample_id2, ".y")]
      }
      corr <- abs(diff) >= effect_size & raw == TRUE

      cols_out <- c(
        paste("corrected", patient_id1, sample_id1, patient_id2, sample_id2, sep = "_")
      )
      DMPs_mb[, as.character(cols_out) := NA]
      DMPs_mb[, as.character(cols_out) := as.list(corr)]
    }
  }
}
rm(raw, diff, corr)

################################################################################
# Union of loci with at least one pair-wise DMP call (mt and/or mb)
################################################################################

DMPs <- na.omit(cbind.data.frame(DMPs_mt, DMPs_mb))
cols_corrected <- grep("corrected_", colnames(DMPs_mt))

venn_corr <- list()
max_DMPs_corr <- max(
  DMPs_mt[, lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)),
    .SDcols = cols_corrected
  ],
  DMPs_mb[, lapply(.SD, function(x) sum(as.numeric(x), na.rm = TRUE)),
    .SDcols = cols_corrected
  ]
)

venn_overlap_list <- list()

for (sample_id1 in sample_ids1) {
  for (sample_id2 in sample_ids2) {
    if (patient_id1 != patient_id2 ||
      (patient_id1 == patient_id2 &&
        which(sample_ids2 == sample_id2) > which(sample_ids1 == sample_id1))) {
      pair_cols <- grep(
        paste0("corrected_", patient_id1, "_", sample_id1, "_",
          patient_id2, "_", sample_id2),
        colnames(DMPs_mt)
      )
      tmp <- na.omit(data.frame(
        mt = DMPs_mt[, unlist(.SD), .SDcols = pair_cols],
        mb = DMPs_mb[, unlist(.SD), .SDcols = pair_cols]
      ))
      rownames(tmp) <- seq_len(nrow(tmp))

      tmp2 <- euler(tmp, shape = "ellipse")
      venn_overlap_list[[length(venn_overlap_list) + 1]] <- data.table(
        sample_id1 = sample_id1,
        sample_id2 = sample_id2,
        perc_mt_also_mb = tmp2$original[3] / sum(tmp2$original[c(1, 3)]),
        perc_mb_also_mt = tmp2$original[3] / sum(tmp2$original[c(2, 3)]),
        mt_only = tmp2$original[1],
        mb_only = tmp2$original[2],
        both = tmp2$original[3]
      )

      f_prefix <- paste0(
        cruk_id, "_", patient_id1, "-", sample_id1, "_",
        patient_id2, "-", sample_id2
      )
      f_pdf <- paste0(f_prefix, "_mt_mt_vs_mb_mb_DMPs.pdf")
      pdf(file.path(out_dir, f_pdf), width = 3, height = 2.5)
      print(plot(tmp2, quantities = TRUE))
      dev.off()
      Sys.sleep(1)
      rm(tmp2)

      df.venn <- data.frame(
        x = c(0, NA),
        y = c(1, 1),
        labels = c("mb1/mb2\n DMPs", "mt1/mt2\n DMPs")
      )
      vdc <- vennCounts(tmp[, c(1, 2)])
      class(vdc) <- "matrix"
      df.vdc <- as.data.frame(vdc)[-1, ]
      df.vdc.ordered <- df.vdc[c(1, 3, 2), ] %>%
        dplyr::mutate(x = NA, y = c(1, 1, 1))
      rm(vdc, df.vdc)

      df.vdc.ordered$label <- rowSums(df.vdc.ordered[, 1:2])
      df.vdc.ordered$label <- ifelse(
        df.vdc.ordered[, grepl("mt", colnames(df.vdc.ordered))] == 1 &
          df.vdc.ordered$label == 1,
        "false negatives",
        ifelse(df.vdc.ordered$label == 2, "true positives", "false positives")
      )

      coordinates <- data.frame(
        label = c("mt", "shared", "mb"),
        Counts = c(
          sum(df.vdc.ordered[df.vdc.ordered$label != "false positives", "Counts"]),
          df.vdc.ordered[df.vdc.ordered$label == "true positives", "Counts"],
          sum(df.vdc.ordered[df.vdc.ordered$label != "false negatives", "Counts"])
        ),
        stringsAsFactors = FALSE
      )
      coordinates$area <- c(
        pi * (coordinates$Counts[1] / max_DMPs_corr),
        pi * (coordinates$Counts[2] / max_DMPs_corr),
        pi * (coordinates$Counts[3] / max_DMPs_corr)
      )
      coordinates$rad <- sqrt(coordinates$area / pi)
      coordinates$rad[2] <- NA
      coordinates$rad2 <- coordinates$area / pi
      coordinates$rad2[2] <- NA
      df.venn$rad <- c(coordinates$rad[3], coordinates$rad[1])

      if (coordinates$area[2] > 0.5 * coordinates$area[1]) {
        s <- seq(
          -min(coordinates$rad, na.rm = TRUE) + 1e-12,
          min(coordinates$rad, na.rm = TRUE) - 1e-12,
          0.0001
        )
        tmp2 <- pi * coordinates$rad2[1] -
          (coordinates$rad2[1] * asin(s / coordinates$rad[1])) +
          s * sqrt(coordinates$rad2[1] - (s^2)) +
          (coordinates$rad2[3] * asin(s / coordinates$rad[3])) -
          s * sqrt(coordinates$rad2[3] - (s^2)) - coordinates$area[2]
        y_coord <- s[which(tmp2 == min(tmp2[tmp2 >= 0]))]
        rm(tmp2, s)
      }
      if (coordinates$area[2] <= 0.5 * coordinates$area[1]) {
        s <- seq(
          -min(coordinates$rad, na.rm = TRUE) + 1e-12,
          min(coordinates$rad, na.rm = TRUE) - 1e-12,
          0.0001
        )
        tmp2 <- (coordinates$rad2[1] * asin(s / coordinates$rad[1])) -
          s * sqrt(coordinates$rad2[1] - (s^2)) +
          (coordinates$rad2[3] * asin(s / coordinates$rad[3])) -
          s * sqrt(coordinates$rad2[3] - (s^2)) - coordinates$area[2]
        y_coord <- s[which(tmp2 == min(tmp2[tmp2 >= 0]))]
        rm(tmp2, s)
      }

      fct2 <- function(x, x_coord, radx, rady) {
        x^2 - 2 * x_coord * x + radx - rady
      }

      x_coord <- sqrt(coordinates$rad2[3] - (y_coord^2))

      radx <- coordinates$rad2[3]
      rady <- coordinates$rad2[1]
      d <- UnirootAll(
        fct2,
        x_coord = x_coord, radx = radx, rady = rady,
        lower = -min(coordinates$rad, na.rm = TRUE) + 1e-6,
        upper = sum(coordinates$rad, na.rm = TRUE)
      )
      if (paste0(sample_id1, "_", sample_id2) %in% c(
        "R1_R5", "R1_R7", "R2_R3", "R2_R4", "R2_R7", "R3_R5",
        "R4_R5", "R3_R7", "R4_R7", "R5_R7"
      )) {
        d <- d[1]
      }
      if (length(d) == 0) {
        idx <- which.min(abs(fct2(seq(0, 1, 0.0001), x_coord, radx, rady)))
        d <- seq(0, 1, 0.0001)[idx]
        rm(idx)
      } else {
        d <- d[which.min(abs(coordinates$rad[1] + sqrt(coordinates$area / pi)[2] - d))]
      }
      df.venn$x[2] <- d
      stopifnot(is.finite(d))

      df.vdc.ordered$x <- c(
        0 - (0.5 * coordinates$rad[3]),
        x_coord,
        d + ((2 * coordinates$rad[1] - d) / 2)
      )
      df.vdc.ordered$y <- c(1 + y_coord, 1, 1 + y_coord)

      out_file <- file.path(
        out_dir,
        paste0(
          cruk_id, "_", sample_id1, "_", sample_id2,
          "_mt_mt_vs_mb_mb_DMPs_corrected.pdf"
        )
      )
      ggplot(df.venn) +
        coord_fixed() +
        theme_void() +
        theme(legend.position = "right") +
        geom_circle(
          aes(x0 = x, y0 = y, r = rad, color = labels, fill = labels),
          alpha = 0.1, size = 2
        ) +
        scale_colour_manual(
          name = "legend",
          values = c("lightgrey", "dodgerblue")
        ) +
        scale_fill_manual(
          name = "legend",
          values = c("lightgrey", "dodgerblue")
        ) +
        labs(fill = NULL) +
        annotate(
          "text",
          x = df.vdc.ordered$x,
          y = df.vdc.ordered$y,
          label = df.vdc.ordered$Counts,
          size = 5
        )
      ggsave(
        filename = out_file,
        device = "pdf",
        width = 4.5,
        height = 2.5,
        units = "in"
      )

      venn_corr[[length(venn_corr) + 1]] <- df.venn
      venn_corr[[length(venn_corr)]]$sample_id1 <- sample_id1
      venn_corr[[length(venn_corr)]]$sample_id2 <- sample_id2
    }
  }
}

rm(df.venn, df.vdc.ordered, tmp, pair_cols, coordinates, x_coord, y_coord, d)

################################################################################
# Summaries and outputs
################################################################################

df_overlaps <- rbindlist(venn_overlap_list)

cat(cruk_id, "(", patient_id, ") — median fraction of mt DMPs also mb:",
  median(df_overlaps$perc_mt_also_mb, na.rm = TRUE), "\n")
cat(cruk_id, "(", patient_id, ") — median fraction of mb DMPs also mt:",
  median(df_overlaps$perc_mb_also_mt, na.rm = TRUE), "\n")

overlap_rdata <- file.path(
  out_dir,
  paste0("intra_tumour_DMPs_", cruk_id, "_", patient_id, ".RData")
)
save(df_overlaps, file = overlap_rdata)
cat("Saved:", overlap_rdata, "\n")
rm(venn_overlap_list)

df_plot <- dplyr::bind_rows(venn_corr)
rho1 <- patient_info[Patient == patient_id1, .(Sample, rho_reseg)]
rho2 <- patient_info[Patient == patient_id2, .(Sample, rho_reseg)]
df_plot$sample_id1 <- factor(
  df_plot$sample_id1,
  levels = rho1$Sample[order(rho1$rho_reseg)],
  ordered = TRUE
)
df_plot$sample_id2 <- factor(
  df_plot$sample_id2,
  levels = rho2$Sample[order(rho2$rho_reseg)],
  ordered = TRUE
)

summary_pdf <- file.path(
  out_dir,
  paste0("summary_figure_venn_corrected_", cruk_id, "_", patient_id, ".pdf")
)
ggplot(
  data = df_plot[is.finite(df_plot$x) & !grepl("raw", df_plot$labels), ]
) +
  coord_fixed() +
  theme_void() +
  facet_wrap(
    sample_id1 ~ sample_id2,
    ncol = length(sample_ids2),
    nrow = length(sample_ids1),
    drop = FALSE
  ) +
  geom_circle(
    aes(x0 = x, y0 = y, r = rad, color = labels, fill = labels),
    alpha = 0.25,
    size = 0.5
  ) +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.line = element_blank(),
    legend.text = element_text(size = 8),
    axis.ticks = element_blank(),
    strip.background = element_blank()
  ) +
  scale_colour_manual(
    name = "legend",
    values = c(
      "mb1/mb2\n DMPs" = "lightgrey",
      "mt1/mt2\n DMPs" = "dodgerblue",
      "mb1/mb2\n raw DMPs" = "grey20",
      "mt1/mt2\n raw DMPs" = "darkblue"
    )
  ) +
  scale_fill_manual(
    name = "legend",
    values = c(
      "mb1/mb2\n DMPs" = "lightgrey",
      "mb1/mb2\n raw DMPs" = "grey20",
      "mt1/mt2\n DMPs" = "dodgerblue",
      "mt1/mt2\n raw DMPs" = "darkblue"
    )
  )
ggsave(
  filename = summary_pdf,
  device = "pdf",
  width = 5,
  height = 6.5,
  units = "in"
)
cat("Saved:", summary_pdf, "\n")

# END
