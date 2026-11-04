#' Enriched per-locus dashboard (fine-mapping) — builders + data assembly.
#'
#' Pure functions (no Shiny) so they are unit-testable against an opened bundle.
#' The UI/server that drive them live in mod_snp_assoc.R (the Fine-mapping tab).
#' Panels (stacked with cowplot, shared x): A association + LD (Phase 3),
#' B fine-mapping PIP, C genes-by-biotype. See R/mod_locus_plot.R builders.

# ---- palettes ---------------------------------------------------------------

# Distinct colour per credible set. `levels` are "CS1","CS5",… in display order.
.lp_cs_palette <- function(levels) {
  base <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02",
            "#a6761d","#666666","#1f78b4","#b2182b","#35978f","#bf812d")
  levels <- unique(levels)
  cols <- base[((seq_along(levels) - 1) %% length(base)) + 1]
  setNames(cols, levels)
}

# Biotype -> colour (gene-body track). Unknown biotypes fall back to grey.
.lp_biotype_palette <- function() {
  c(protein_coding = "#2b2b2b", lincRNA = "#1f78b4", lncRNA = "#1f78b4",
    antisense = "#35978f", processed_transcript = "#80cdc1",
    miRNA = "#33a02c", snRNA = "#6a3d9a", snoRNA = "#6a3d9a",
    pseudogene = "#9e9e9e", processed_pseudogene = "#9e9e9e",
    TR_gene = "#b15928", IG_gene = "#b15928")
}
.lp_biotype_colour <- function(bt) {
  pal <- .lp_biotype_palette()
  out <- unname(pal[as.character(bt)])
  out[is.na(out)] <- "#c0c0c0"
  out
}

.lp_is_true <- function(x) {
  if (is.logical(x)) return(x %in% TRUE)
  xc <- tolower(as.character(x))
  xc %in% c("true","t","1","yes")
}

# ---- legend grobs (extracted via a dummy plot with the real scale) ----------
.lp_legend_from <- function(g) cowplot::get_legend(
  g + ggplot2::theme(legend.position = "right",
                     legend.key.size = grid::unit(0.8, "lines"),
                     legend.margin = ggplot2::margin(2, 2, 2, 2)))

lp_ld_legend <- function(font_size = 11) {
  pal <- .lp_ld_palette()
  d <- data.frame(x = 1, y = seq_along(pal), b = factor(names(pal), levels = names(pal)))
  .lp_legend_from(ggplot2::ggplot(d, ggplot2::aes(x, y, colour = b)) + ggplot2::geom_point() +
    ggplot2::scale_colour_manual(expression(r^2 * " to lead"), values = pal) +
    ggplot2::theme_void(base_size = font_size))
}
lp_z_legend <- function(font_size = 11) {
  d <- data.frame(x = 1, y = 1, z = 0)
  .lp_legend_from(ggplot2::ggplot(d, ggplot2::aes(x, y, colour = z)) + ggplot2::geom_point() +
    ggplot2::scale_colour_gradient2("Z-score", low = "#2166ac", mid = "#f7f7f7",
                                    high = "#b2182b", midpoint = 0, limits = c(-6, 6)) +
    ggplot2::theme_void(base_size = font_size))
}
lp_prio_legend <- function(font_size = 11) {
  d <- data.frame(x = 1, y = 1, p = 0.5)
  .lp_legend_from(ggplot2::ggplot(d, ggplot2::aes(x, y, colour = p)) + ggplot2::geom_point() +
    ggplot2::scale_colour_gradient("CALDERA prob.", low = "#f2f0f7", high = "#54278f", limits = c(0, 1)) +
    ggplot2::theme_void(base_size = font_size))
}

# Label-width-aware greedy packing (ported from scripts/plot_locus_zoom.R): each
# gene's effective footprint is widened to include its label text so gene names
# in the track never collide. Assigns each gene the lowest free row.
lp_pack_rows <- function(genes, window) {
  if (is.null(genes) || nrow(genes) == 0) return(genes)
  win_w <- window[2] - window[1]
  chars_per_window <- 56; pad <- 1.4
  mid <- (genes$gstart_mb + genes$gend_mb) / 2
  lab_half <- nchar(as.character(genes$symbol)) * (win_w / chars_per_window) / 2 * pad
  genes$eff_start <- pmin(genes$gstart_mb, mid - lab_half)
  genes$eff_end   <- pmax(genes$gend_mb,   mid + lab_half)
  g <- genes[order(genes$eff_start), , drop = FALSE]
  row_end <- numeric(0); g$row <- NA_integer_
  for (i in seq_len(nrow(g))) {
    free <- which(g$eff_start[i] > row_end)
    k <- if (length(free)) free[1] else length(row_end) + 1L
    g$row[i] <- k
    row_end[k] <- g$eff_end[i]
  }
  g
}

# ---- data assembly ----------------------------------------------------------

#' Assemble everything a locus plot needs for one lead locus + L.
#' Returns NULL when there is nothing to draw; else a list with $members,
#' $baseline, $cs_summary, $genes (packed), $window (Mb), $chr, $cs_palette,
#' $all_points (for hover), and flags $has_baseline / $has_bodies.
lp_build_locus_data <- function(gd, gwas, locus, L,
                                pip_type = "marginal", pc_only = FALSE,
                                hc_methods = c("caldera","coding","nearest"),
                                caldera_hc_thr = 0.5, pad_frac = 0.05) {
  # Accept L as 1/10 or the UI radio values "L1"/"L10".
  Lval <- suppressWarnings(as.integer(gsub("[^0-9]", "", as.character(L))))
  if (length(Lval) == 0 || is.na(Lval)) Lval <- 10L
  rd <- function(block) tryCatch(gd_read(gd, gwas, block), error = function(e) NULL)
  sub_locus <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    d <- as.data.frame(df)
    keep <- as.character(d$locus) == as.character(locus)
    if ("L" %in% names(d)) keep <- keep & (suppressWarnings(as.integer(d$L)) == Lval)
    d <- d[keep, , drop = FALSE]
    if (nrow(d) == 0) NULL else d
  }

  cc <- rd("credset_coding")
  members <- sub_locus(cc$variants)
  sa <- rd("snp_assoc")
  lp <- if (!is.null(sa) && !is.null(sa$susie))
    sa$susie[[if (Lval == 1) "locus_pip_L1" else "locus_pip_L10"]] else NULL
  baseline <- NULL
  if (!is.null(lp) && nrow(lp) > 0) {
    b <- as.data.frame(lp); b <- b[as.character(b$locus) == as.character(locus), , drop = FALSE]
    if (nrow(b) > 0) baseline <- b
  }
  if (is.null(members) && is.null(baseline)) return(NULL)

  bps <- c(baseline$BP, members$BP)
  chr <- if (!is.null(members)) members$CHR[1] else baseline$CHR[1]
  lo <- min(bps, na.rm = TRUE); hi <- max(bps, na.rm = TRUE)
  pad <- max(1e4, pad_frac * (hi - lo))
  if (is.null(baseline)) {
    # No full-locus baseline (bundle predates locus_pip): the members alone give a
    # tiny span, so widen to a fine-map-sized region (~±250kb) so genes appear.
    centre <- (lo + hi) / 2; half <- max(2.5e5, (hi - lo) / 2 + pad)
    window <- c(centre - half, centre + half) / 1e6
  } else {
    window <- c(lo - pad, hi + pad) / 1e6
  }

  if (!is.null(members)) {
    yv <- if (pip_type == "alpha" && "PIP" %in% names(members))
      suppressWarnings(as.numeric(members$PIP)) else suppressWarnings(as.numeric(members$marginal_pip))
    if (all(is.na(yv))) yv <- suppressWarnings(as.numeric(members$PIP))
    members$y <- yv
    members$BP_Mb <- members$BP / 1e6
    members$pa <- if ("protein_altering" %in% names(members)) .lp_is_true(members$protein_altering) else FALSE
    members$csf <- paste0("CS", members$cs)
  }
  cs_levels <- sort(unique(c(members$csf, if (!is.null(baseline)) paste0("CS", baseline$cs[baseline$cs > 0]))))
  cs_palette <- .lp_cs_palette(cs_levels)

  # Per-CS summary: lead, size, purity (cs_avg_r2), coding share.
  cs_summary <- NULL
  st <- if (!is.null(sa) && !is.null(sa$susie)) sa$susie[[if (Lval == 1) "L1" else "L10"]] else NULL
  if (!is.null(st) && nrow(st) > 0) {
    st <- as.data.frame(st); st <- st[as.character(st$SNP) == as.character(locus), , drop = FALSE]
    if (nrow(st) > 0) cs_summary <- data.frame(
      cs = st$cs, csf = paste0("CS", st$cs),
      nsnp = if ("NSNP" %in% names(st)) st$NSNP else NA,
      purity = if ("cs_avg_r2" %in% names(st)) round(as.numeric(st$cs_avg_r2), 2) else NA,
      stringsAsFactors = FALSE)
  }
  csum <- sub_locus(cc$summary)
  if (!is.null(cs_summary) && !is.null(csum) && "frac_protein_altering" %in% names(csum)) {
    cs_summary$coding_share <- round(as.numeric(csum$frac_protein_altering[match(cs_summary$cs, csum$cs)]), 2)
  }

  # Genes overlapping the window, with HC criteria.
  gp <- tryCatch(load_gene_positions(), error = function(e) NULL)
  genes <- NULL; has_bodies <- FALSE
  if (!is.null(gp) && !is.null(gp$canonical)) {
    cg <- as.data.frame(gp$canonical)
    cg <- cg[!is.na(cg$CHR) & cg$CHR == chr, , drop = FALSE]
    has_bodies <- all(c("gene_start","gene_end") %in% names(cg))
    if (has_bodies) {
      cg$gstart_mb <- cg$gene_start / 1e6; cg$gend_mb <- cg$gene_end / 1e6
      cg <- cg[!is.na(cg$gstart_mb) & !is.na(cg$gend_mb) &
               cg$gend_mb > window[1] & cg$gstart_mb < window[2], , drop = FALSE]
    } else {
      cg$gstart_mb <- cg$BP / 1e6; cg$gend_mb <- cg$BP / 1e6
      cg <- cg[!is.na(cg$BP) & cg$gstart_mb > window[1] & cg$gstart_mb < window[2], , drop = FALSE]
    }
    if (isTRUE(pc_only)) cg <- cg[!is.na(cg$gene_biotype) & cg$gene_biotype == "protein_coding", , drop = FALSE]
    if (nrow(cg) > 0) {
      cg$body_colour <- .lp_biotype_colour(cg$gene_biotype)
      cg$hc <- ""
      cal <- rd("caldera")
      c2 <- if (!is.null(cal)) sub_locus(cal) else NULL
      # Match on Ensembl id where available (CALDERA uses GENCODE symbols like
      # CFAP410 that the older biomart can't resolve, so a symbol match misses
      # genes the biomart knows as e.g. C21orf2); fall back to canonical symbol.
      add_by_ensgid <- function(letters, ens) {
        ens <- unique(ens[!is.na(ens) & nzchar(as.character(ens))])
        if (!length(ens)) return(invisible())
        idx <- cg$ensembl_gene_id %in% ens
        cg$hc[idx] <<- paste0(cg$hc[idx], letters)
      }
      add_by_symbol <- function(letters, syms) {
        syms <- unique(harmonise_genes(as.character(syms))$symbol)
        idx <- cg$symbol %in% syms
        cg$hc[idx] <<- paste0(cg$hc[idx], letters)
      }
      if ("caldera" %in% hc_methods && !is.null(c2)) {
        hit <- suppressWarnings(as.numeric(c2$caldera)) >= caldera_hc_thr
        if ("ensgid" %in% names(c2)) add_by_ensgid("C", c2$ensgid[hit]) else add_by_symbol("C", c2$gene[hit])
      }
      if ("coding" %in% hc_methods && !is.null(members)) {
        pa <- members[members$pa & suppressWarnings(as.numeric(members$PIP)) >= 0.5, , drop = FALSE]
        add_by_symbol("M", pa$gene)
      }
      if ("nearest" %in% hc_methods && !is.null(c2) && "nearest_gene" %in% names(c2)) {
        nr <- .lp_is_true(c2$nearest_gene)
        if ("ensgid" %in% names(c2)) add_by_ensgid("N", c2$ensgid[nr]) else add_by_symbol("N", c2$gene[nr])
      }
      cg$is_hc <- nzchar(cg$hc)
      cg$gstart_mb <- pmax(cg$gstart_mb, window[1]); cg$gend_mb <- pmin(cg$gend_mb, window[2])
      genes <- lp_pack_rows(cg, window)
    }
  }

  # Panel A: association (-log10 P) over the window, coloured by r^2-to-lead.
  # Prefer the full-locus baseline (every SNP, with P) for a dense scatter; fall
  # back to the (p<1e-3-thinned) manhattan_data on bundles without locus_pip.
  assoc <- NULL
  if (!is.null(baseline) && "P" %in% names(baseline)) {
    a <- baseline; a$neglogP <- -log10(suppressWarnings(as.numeric(a$P)))
    a <- a[is.finite(a$neglogP), , drop = FALSE]
    if (nrow(a) > 0) assoc <- a[, c("SNP","CHR","BP","neglogP")]
  }
  if (is.null(assoc)) {
    gq <- rd("gwas_qc")
    md <- if (!is.null(gq)) gq$manhattan_data else NULL
    if (!is.null(md) && !is.null(md$data)) {
      a <- as.data.frame(md$data)
      a <- a[!is.na(a$CHR) & a$CHR == chr &
             a$BP / 1e6 >= window[1] & a$BP / 1e6 <= window[2], , drop = FALSE]
      if (nrow(a) > 0) assoc <- a
    }
  }
  ld2 <- if (!is.null(sa) && !is.null(sa$susie)) sa$susie$ld_to_lead else NULL
  r2map <- NULL
  if (!is.null(ld2) && nrow(ld2) > 0) {
    l2 <- as.data.frame(ld2); l2 <- l2[as.character(l2$locus) == as.character(locus), , drop = FALSE]
    if (nrow(l2) > 0) r2map <- setNames(suppressWarnings(as.numeric(l2$r2)), as.character(l2$SNP))
  }
  lead_bp <- NA_real_
  if (!is.null(assoc) && as.character(locus) %in% assoc$SNP) lead_bp <- assoc$BP[assoc$SNP == as.character(locus)][1]
  if (is.na(lead_bp) && !is.null(st) && nrow(st) > 0) lead_bp <- suppressWarnings(as.numeric(st$BP[1]))
  # Recombination-rate track for the locus (BP, rate), when packaged.
  rcb <- if (!is.null(sa) && !is.null(sa$susie)) sa$susie$recomb else NULL
  recomb <- NULL
  if (!is.null(rcb) && nrow(rcb) > 0) {
    r <- as.data.frame(rcb); r <- r[as.character(r$locus) == as.character(locus), c("BP","rate")]
    if (nrow(r) >= 2) recomb <- r[order(r$BP), ]
  }

  # Points frame for hover (baseline + members).
  all_points <- NULL
  if (!is.null(members))
    all_points <- data.frame(BP_Mb = members$BP_Mb, y = members$y, SNP = members$SNP,
                             cs = members$csf, pip = round(members$y, 4),
                             consequence = if ("consequence" %in% names(members)) members$consequence else NA,
                             gene = if ("gene" %in% names(members)) members$gene else NA,
                             stringsAsFactors = FALSE)

  list(members = members, baseline = baseline, cs_summary = cs_summary,
       genes = genes, window = window, chr = chr, cs_palette = cs_palette,
       all_points = all_points, pip_type = pip_type,
       assoc = assoc, r2map = r2map, lead_bp = lead_bp, recomb = recomb,
       has_baseline = !is.null(baseline), has_bodies = has_bodies,
       has_assoc = !is.null(assoc), has_ld = !is.null(r2map), locus = locus, L = Lval)
}

# ---- credible-set shaded bands (shared across panels) -----------------------

.lp_cs_bands <- function(ld) {
  m <- ld$members
  if (is.null(m)) return(NULL)
  sp <- do.call(rbind, lapply(split(m, m$csf), function(d) data.frame(
    csf = d$csf[1], xmin = min(d$BP_Mb), xmax = max(d$BP_Mb), stringsAsFactors = FALSE)))
  sp$col <- ld$cs_palette[sp$csf]
  sp
}

# ---- Panel A: association (-log10 P) coloured by LD r^2 to lead --------------

# LocusZoom-style r^2 palette (matches scripts/plot_locus_zoom.R).
.lp_ld_palette <- function() c("<0.2"="#357EBD","0.2-0.4"="#46B8DA","0.4-0.6"="#5CB85C",
                               "0.6-0.8"="#F0AD4E","0.8-1"="#D43F3A")
.lp_index_colour <- "#9632B8"

build_locus_assoc_panel <- function(ld, font_size = 14, point_size = 1.3,
                                    theme_fn = ggplot2::theme_bw) {
  a <- ld$assoc
  if (is.null(a) || nrow(a) == 0) return(NULL)
  a$BP_Mb <- a$BP / 1e6
  a$neglogP <- if ("neglogP" %in% names(a)) suppressWarnings(as.numeric(a$neglogP)) else -log10(suppressWarnings(as.numeric(a$P)))
  r2 <- if (!is.null(ld$r2map)) unname(ld$r2map[as.character(a$SNP)]) else rep(NA_real_, nrow(a))
  pal <- .lp_ld_palette()
  bins <- cut(r2, breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf),
              labels = names(pal), right = FALSE)
  a$col <- unname(pal[as.character(bins)]); a$col[is.na(a$col)] <- "grey70"
  bands <- .lp_cs_bands(ld)
  p <- ggplot2::ggplot()
  if (!is.null(bands)) p <- p + ggplot2::geom_rect(
    data = bands, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = bands$col, alpha = 0.12, inherit.aes = FALSE)
  ymax <- max(c(a$neglogP, -log10(5e-9)), na.rm = TRUE) + 0.5
  sf <- ymax / 100
  rc <- ld$recomb
  p <- p +
    ggplot2::geom_point(data = a, ggplot2::aes(x = BP_Mb, y = neglogP), colour = a$col, size = point_size) +
    ggplot2::geom_hline(yintercept = -log10(5e-8), colour = "red", linetype = "solid", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(1e-5), colour = "blue", linetype = "dashed", linewidth = 0.3)
  # Recombination-rate line (secondary axis) drawn ON TOP of the dense scatter so
  # it stays visible; scaled to the -log10(P) range, LocusZoom-style.
  if (!is.null(rc) && nrow(rc) >= 2) p <- p + ggplot2::geom_line(
    data = rc, ggplot2::aes(x = BP / 1e6, y = rate * sf),
    colour = "#2B7BBA", linewidth = 0.8, alpha = 0.85)
  lead_row <- a[a$SNP == as.character(ld$locus), , drop = FALSE]
  if (nrow(lead_row) == 0 && !is.na(ld$lead_bp)) lead_row <- a[which.min(abs(a$BP - ld$lead_bp)), , drop = FALSE]
  if (nrow(lead_row) > 0) p <- p + ggplot2::geom_point(
    data = lead_row, ggplot2::aes(x = BP_Mb, y = neglogP),
    shape = 23, size = point_size * 3.2, fill = .lp_index_colour, colour = "black")
  p +
    ggplot2::scale_y_continuous(
      name = expression(-log[10](italic(P))),
      sec.axis = ggplot2::sec_axis(~ . / sf, name = "Recomb. (cM/Mb)",
                                   breaks = c(0, 20, 40, 60, 80, 100))) +
    ggplot2::coord_cartesian(xlim = ld$window, ylim = c(0, ymax)) +
    ggplot2::labs(x = NULL) +
    theme_fn(base_size = font_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(4, 8, 2, 8))
}

# ---- Panel B: fine-mapping PIP ---------------------------------------------

build_locus_finemap_panel <- function(ld, font_size = 14, point_size = 1.6,
                                      theme_fn = ggplot2::theme_bw) {
  if (is.null(ld) || (is.null(ld$members) && is.null(ld$baseline))) return(NULL)
  bands <- .lp_cs_bands(ld)
  p <- ggplot2::ggplot()
  if (!is.null(bands)) p <- p + ggplot2::geom_rect(
    data = bands, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = bands$col, alpha = 0.12, inherit.aes = FALSE)
  if (!is.null(ld$baseline)) {
    b <- ld$baseline; b <- b[suppressWarnings(as.numeric(b$pip)) > 1e-4, , drop = FALSE]
    if (nrow(b) > 0) p <- p + ggplot2::geom_point(
      data = b, ggplot2::aes(x = BP / 1e6, y = pip), colour = "grey70",
      size = point_size * 0.7, alpha = 0.5)
  }
  if (!is.null(ld$members)) {
    m <- ld$members; m$col <- ld$cs_palette[m$csf]
    p <- p +
      ggplot2::geom_point(data = m, ggplot2::aes(x = BP_Mb, y = y, colour = col, shape = pa),
                          size = point_size + 0.6) +
      ggplot2::scale_colour_identity() +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none")
    lab <- do.call(rbind, lapply(split(m, m$csf), function(d) d[which.max(d$y), , drop = FALSE]))
    lab$txt <- paste0(lab$csf)
    pa_lab <- m[m$pa, , drop = FALSE]
    if (nrow(pa_lab) > 0) pa_lab$txt <- paste0(ifelse(!is.na(pa_lab$gene) & nzchar(pa_lab$gene),
                                                      paste0(pa_lab$consequence, " in ", pa_lab$gene), pa_lab$consequence))
    labs_all <- rbind(lab[, c("BP_Mb","y","txt")], if (nrow(pa_lab)) pa_lab[, c("BP_Mb","y","txt")] else NULL)
    p <- p + ggrepel::geom_text_repel(data = labs_all, ggplot2::aes(x = BP_Mb, y = y, label = txt),
                                      size = font_size * 0.26, min.segment.length = 0,
                                      max.overlaps = Inf, box.padding = 0.4, segment.size = 0.3)
  }
  if (!is.null(ld$cs_summary)) {
    txt <- apply(ld$cs_summary, 1, function(r) {
      bits <- c(paste0(r[["csf"]], ": ", r[["nsnp"]], " SNPs"),
                if (!is.na(r[["purity"]])) paste0("purity ", r[["purity"]]),
                if ("coding_share" %in% names(r) && !is.na(r[["coding_share"]])) paste0("coding share ", r[["coding_share"]]))
      paste(bits, collapse = " | ")
    })
    lab_df <- data.frame(x = ld$window[1], y = 0.98 - 0.07 * (seq_along(txt) - 1), txt = txt,
                         col = ld$cs_palette[ld$cs_summary$csf], stringsAsFactors = FALSE)
    p <- p + ggplot2::geom_text(data = lab_df, ggplot2::aes(x = x, y = y, label = txt),
                                colour = lab_df$col, hjust = 0, size = font_size * 0.24)
  }
  p +
    ggplot2::scale_y_continuous(limits = c(-0.02, 1.02), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
    ggplot2::coord_cartesian(xlim = ld$window) +
    ggplot2::labs(x = NULL, y = "PIP") +
    theme_fn(base_size = font_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(4, 8, 2, 8))
}

# ---- Panel C: genes by biotype, HC shading + badges -------------------------

build_locus_gene_panel <- function(ld, font_size = 14, theme_fn = ggplot2::theme_bw) {
  g <- ld$genes
  if (is.null(g) || nrow(g) == 0) return(NULL)
  nrows <- max(g$row)
  w <- ld$window
  p <- ggplot2::ggplot()
  bands <- .lp_cs_bands(ld)
  if (!is.null(bands)) p <- p + ggplot2::geom_rect(
    data = bands, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = bands$col, alpha = 0.12, inherit.aes = FALSE)
  # High-confidence highlight box behind the gene body (filled + padded so even a
  # small gene is clearly flagged).
  hc <- g[g$is_hc, , drop = FALSE]
  if (nrow(hc) > 0) {
    xpad <- (w[2] - w[1]) * 0.004
    hc$bx1 <- hc$gstart_mb - xpad; hc$bx2 <- hc$gend_mb + xpad
    p <- p + ggplot2::geom_rect(
      data = hc, ggplot2::aes(xmin = bx1, xmax = bx2, ymin = row - 0.33, ymax = row + 0.33),
      fill = "#fdecc9", colour = "#e39a2e", linewidth = 1.0, inherit.aes = FALSE)
  }
  if (isTRUE(ld$has_bodies)) {
    # Gene body as a filled bar coloured by biotype, with a strand-direction arrow.
    p <- p + ggplot2::geom_rect(
      data = g, ggplot2::aes(xmin = gstart_mb, xmax = gend_mb, ymin = row - 0.12, ymax = row + 0.12,
                             fill = body_colour), colour = NA) +
      ggplot2::scale_fill_identity()
    ga <- g[!is.na(g$strand), , drop = FALSE]
    if (nrow(ga) > 0) {
      ext <- (w[2] - w[1]) * 0.012
      ga$ax  <- ifelse(ga$strand == 1, ga$gend_mb, ga$gstart_mb)
      ga$axe <- ifelse(ga$strand == 1, pmin(ga$gend_mb + ext, w[2]), pmax(ga$gstart_mb - ext, w[1]))
      p <- p + ggplot2::geom_segment(
        data = ga, ggplot2::aes(x = ax, xend = axe, y = row, yend = row),
        arrow = grid::arrow(length = grid::unit(0.06, "cm"), type = "closed"),
        colour = "grey25", linewidth = 0.3)
    }
  } else {
    p <- p + ggplot2::geom_point(
      data = g, ggplot2::aes(x = gstart_mb, y = row, colour = body_colour), size = 2.5) +
      ggplot2::scale_colour_identity()
  }
  # Label every gene — the label-width-aware packing keeps names from colliding.
  # Label font ~ axis tick-label size (base_size * 0.8, converted pt -> mm).
  p <- p + ggplot2::geom_text(
    data = g, ggplot2::aes(x = (gstart_mb + gend_mb) / 2, y = row + 0.44,
                           label = symbol, fontface = ifelse(is_hc, "bold", "plain")),
    size = font_size * 0.28, colour = "grey15", vjust = 0)
  badges <- g[g$is_hc, , drop = FALSE]
  if (nrow(badges) > 0) p <- p + ggplot2::geom_text(
    data = badges, ggplot2::aes(x = (gstart_mb + gend_mb) / 2, y = row - 0.44, label = hc),
    size = font_size * 0.24, vjust = 1, colour = "#b5651d", fontface = "bold")
  p +
    ggplot2::scale_y_continuous(limits = c(0.3, nrows + 0.9), breaks = NULL) +
    ggplot2::coord_cartesian(xlim = w) +
    ggplot2::labs(x = paste0("Chromosome ", ld$chr, " position (Mb)"), y = NULL) +
    theme_fn(base_size = font_size) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(2, 8, 4, 8))
}

# ---- Panel D: gene x method association grid (positional) -------------------

.lp_z_hex <- function(z, lim = 6) {
  z <- suppressWarnings(as.numeric(z)); z[is.na(z)] <- 0
  zc <- pmax(-lim, pmin(lim, z)) / (2 * lim) + 0.5
  m <- grDevices::colorRamp(c("#2166ac", "#f7f7f7", "#b2182b"))(zc)
  grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
}
.lp_prio_hex <- function(p) {
  p <- suppressWarnings(as.numeric(p)); p[is.na(p)] <- 0; p <- pmax(0, pmin(1, p))
  m <- grDevices::colorRamp(c("#f2f0f7", "#54278f"))(p)
  grDevices::rgb(m[, 1], m[, 2], m[, 3], maxColorValue = 255)
}

#' build_mol_assoc_data for the locus window, with the method-facet Group label
#' and positional gene order attached. Shared by the grid builder + its controls.
.lp_grid_data <- function(gd, gwas, cf, ld) {
  dat <- tryCatch(build_mol_assoc_data(gd, gwas, cf), error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  dat <- as.data.frame(dat)
  pos <- tryCatch(resolve_feature_positions(unique(as.character(dat$ID)), gd, gwas, cf), error = function(e) NULL)
  if (is.null(pos) || nrow(pos) == 0) return(NULL)
  pos <- pos[!is.na(pos$CHR) & pos$CHR == ld$chr &
             pos$BP / 1e6 >= ld$window[1] & pos$BP / 1e6 <= ld$window[2], , drop = FALSE]
  if (nrow(pos) == 0) return(NULL)
  dat <- dat[as.character(dat$ID) %in% pos$ID, , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)
  dat$Group <- make_group_labels(dat$Method, dat$Type)
  attr(dat, "gene_order") <- pos$ID[order(pos$BP)]
  dat
}

# Methods (Groups) + expression/splicing + protein panels in the window, for the
# grid's selectors (mirrors the Molecular Associations summary controls).
lp_grid_options <- function(gd, gwas, cf, ld) {
  dat <- .lp_grid_data(gd, gwas, cf, ld)
  if (is.null(dat)) return(list(groups = character(0), expr = character(0), protein = character(0)))
  list(groups  = intersect(get_group_order(), unique(dat$Group)),
       expr    = sort(unique(as.character(dat$Panel[dat$Type %in% c("Expr.", "Splice")]))),
       protein = sort(unique(as.character(dat$Panel[dat$Type == "Protein"]))))
}

#' Positional gene x method association grid, FACETED by method (like the
#' Molecular Associations summary), filterable by method + panel. Returns NULL
#' if nothing in the window after filtering.
build_locus_assoc_grid <- function(gd, gwas, cf, ld, methods = NULL, panels = NULL,
                                   font_size = 14, theme_fn = ggplot2::theme_bw) {
  dat <- .lp_grid_data(gd, gwas, cf, ld)
  if (is.null(dat)) return(NULL)
  gene_order <- attr(dat, "gene_order")
  # `panels` filters only the tissue panels (TWAS/SMR/PWAS); method-level rows
  # (SuSiE/MAGMA/Nearest/CALDERA) are kept whenever their method is selected.
  tissue <- unique(as.character(dat$Panel[dat$Type %in% c("Expr.", "Splice", "Protein")]))
  if (!is.null(methods)) dat <- dat[dat$Group %in% methods, , drop = FALSE]
  if (!is.null(panels))  dat <- dat[!(as.character(dat$Panel) %in% tissue) |
                                     as.character(dat$Panel) %in% panels, , drop = FALSE]
  if (nrow(dat) == 0) return(NULL)
  dat$ID <- factor(as.character(dat$ID), levels = gene_order)
  dat$Group <- factor(dat$Group, levels = intersect(get_group_order(), unique(dat$Group)))

  has_prio <- "Prio" %in% names(dat)
  cls <- ifelse(grepl("^CALDERA", dat$Method), "caldera",
         ifelse(grepl("^SuSiE", dat$Method) | dat$Method == "Nearest\nGene", "presence", "assoc"))
  dat$fill <- "#dddddd"
  dat$fill[cls == "assoc"] <- .lp_z_hex(dat$Z[cls == "assoc"])
  dat$fill[cls == "presence"] <- "#2ca25f"
  if (has_prio) dat$fill[cls == "caldera"] <- .lp_prio_hex(dat$Prio[cls == "caldera"])

  p <- ggplot2::ggplot(dat, ggplot2::aes(x = ID, y = Panel)) +
    ggplot2::geom_point(colour = dat$fill, size = 3.3)
  sig <- dat[dat$Sig %in% TRUE, , drop = FALSE]
  if (nrow(sig)) p <- p + ggplot2::geom_point(data = sig, shape = 1, size = 4.4, colour = "black", stroke = 0.7)
  clc <- dat[dat$Coloc %in% TRUE, , drop = FALSE]
  if (nrow(clc)) p <- p + ggplot2::geom_point(data = clc, shape = 0, size = 5.0, colour = "black", stroke = 0.7)
  p <- p +
    ggplot2::facet_grid(rows = ggplot2::vars(Group), scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_y_discrete(position = "right") +   # panel/tissue labels on the right
    ggplot2::labs(x = NULL, y = NULL) +
    theme_fn(base_size = font_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   strip.placement = "outside",
                   strip.text.y.left = ggplot2::element_text(angle = 0, size = font_size * 0.7),
                   panel.grid.major = ggplot2::element_line(colour = "grey92"),
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.spacing = ggplot2::unit(2, "pt"))
  # Z-score + CALDERA-probability legends to the right of the grid.
  legs <- cowplot::plot_grid(lp_z_legend(font_size), if (has_prio) lp_prio_legend(font_size) else NULL,
                             ncol = 1, align = "v")
  cowplot::plot_grid(p, legs, ncol = 2, rel_widths = c(1, 0.16))
}

# Facet-row + gene counts (for dynamic plot height/width), tissue-aware filter.
lp_grid_dims <- function(gd, gwas, cf, ld, methods = NULL, panels = NULL) {
  dat <- .lp_grid_data(gd, gwas, cf, ld)
  if (is.null(dat)) return(list(nrows = 0L, ngenes = 0L))
  tissue <- unique(as.character(dat$Panel[dat$Type %in% c("Expr.", "Splice", "Protein")]))
  if (!is.null(methods)) dat <- dat[dat$Group %in% methods, , drop = FALSE]
  if (!is.null(panels))  dat <- dat[!(as.character(dat$Panel) %in% tissue) |
                                     as.character(dat$Panel) %in% panels, , drop = FALSE]
  if (nrow(dat) == 0) return(list(nrows = 0L, ngenes = 0L))
  list(nrows = length(unique(paste(dat$Group, dat$Panel))),
       ngenes = length(unique(as.character(dat$ID))))
}

# ---- Assembly ---------------------------------------------------------------

#' Stack the available panels with cowplot (shared x). `panels` is an ordered
#' list of ggplots (NULLs dropped). Returns a cowplot gtable (printable) or NULL.
assemble_locus_plot <- function(panels, rel_heights = NULL) {
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (length(panels) == 0) return(NULL)
  if (length(panels) == 1) return(panels[[1]])
  if (is.null(rel_heights)) rel_heights <- rep(1, length(panels))
  cowplot::plot_grid(plotlist = panels, ncol = 1, align = "v", axis = "lr",
                     rel_heights = rel_heights[seq_along(panels)])
}
