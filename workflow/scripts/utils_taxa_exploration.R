###############################################################################
# PERMANOVA
###############################################################################

###############################################################################
# PCOA
###############################################################################

plot_beta_div <- function(physeq, dist, group, show_centroids = TRUE, ellipse = TRUE,
                          theme = c("pubr", "prism", "classic"), title = NULL, pairwise_adonis = FALSE,
                          add_pval = FALSE, physeq_object = TRUE) {
  set.seed(123)

  #' @title Plot beta diversity metrics
  #' @description Create a PCoA with centroids using phyloseq.
  #' @param physeq Phyloseq object, containing otu table, metadata and taxa table.
  #' @param dist Distances. Can be obtained using `get_dist()` function. Distance between samples
  #' @param group Variable to group the sample in the PCoA plot (grouping variable)
  #' @param method_ord Ordination method. Default is "PCoA".
  #'  Currently supported method options are: c("DCA", "CCA", "RDA", "CAP", "DPCoA", "NMDS", "MDS", "PCoA")
  #' @param show_centroids If True, plot centroids
  #' @param ellipse If True, show group ellipses
  #' @param theme ggplot theme
  #' @param title plot title
  #' @param pcoa_type The plot type. Default is "samples".
  #' The currently supported options are c("samples", "sites", "species", "taxa", "biplot", "split", "scree")
  #' @param xaxis The principal component to plot in the x axis. Default is "PCoA1"
  #' @param yaxis The principal component to plot in the y axis. Default is "PCoA2"
  #' @return PCoA ggplot object

  #
  ## Settings
  #

  set.seed(123)
  df_comp <- NULL

  # Metadata
  if (physeq_object == TRUE) {
    metadata <- sample_data(physeq)
  } else {
    metadata <- physeq
  }

  theme <- match.arg(theme)

  #
  ## filter distance and metadata to keep only matching sample ids
  #
  ids <- intersect(rownames(metadata), rownames(as.matrix(dist)))

  metadata <- metadata[ids, ]
  dist <- dist(as.matrix(dist)[ids, ids])

  adonis_res <- adonis2(dist ~ metadata[[group]],
    permutations = 10000,
    na.action = na.omit,
    by = "terms"
  )

  #
  ## Compute PCOA
  #

  # Compute distance between samples
  mod <- vegan::betadisper(dist, metadata[[group]])
  d <- permutest(mod)
  centroids <- data.frame(group = rownames(mod$centroids), data.frame(mod$centroids))
  vectors <- data.frame(group = mod$group, data.frame(mod$vectors))

  #
  ## Plot PCoA
  #
  p <- ggplot() +
    geom_point(
      data = vectors, alpha = 0.8,
      aes(color = group, x = PCoA1, y = PCoA2)
    ) +
    guides(color = guide_legend(title = group))

  # theme
  p <- switch(theme,
    "prism" = p + ggprism::theme_prism(),
    "pubr" = p + theme_pubr(),
    "classic" = p + theme_classic()
  )

  if (show_centroids) {
    p <- p + geom_point(data = centroids, aes(color = group, x = PCoA1, y = PCoA2, size = 5)) +
      guides(size = FALSE)
  }

  if (ellipse) {
    p <- p + stat_ellipse(
      data = vectors, mapping = aes(color = group, fill = group, x = PCoA1, y = PCoA2),
      type = "norm", level = 0.9, geom = "polygon", alpha = 0.04
    ) +
      guides(fill = guide_legend(title = group))
  }

  if (!is.null(title)) {
    p <- p + ggtitle(title)
  }



  if (pairwise_adonis) {
    # all comparisons
    comp <- combn(levels(as.factor(metadata[[group]])), 2)

    df_comp <- setNames(
      data.frame(matrix(ncol = 4, nrow = 0)),
      c("group1", "group2", "adonis_pval", "permutest_pval")
    )

    for (i in 1:ncol(comp)) {
      pair <- comp[, i]


      m <- metadata[metadata[[group]] %in% pair, ]

      # filter distance and metadata to keep only matching sample ids

      ids <- intersect(rownames(m), rownames(as.matrix(dist)))

      m <- m[ids, ]
      dist_filt <- dist(as.matrix(dist)[ids, ids])

      adonis_pair_res <- adonis2(dist_filt ~ m[[group]],
        na.action = na.omit,
        by = "terms"
      )
      pval <- adonis_pair_res$`Pr(>F)`[1]

      mod_pair <- betadisper(dist_filt, m[[group]])
      d_pair <- permutest(mod_pair)
      pempval <- d_pair$tab$`Pr(>F)`[1]

      df_comp[nrow(df_comp) + 1, ] <- c(pair, pval, pempval)
    }

    df_comp$adonis_p.adj <- p.adjust(df_comp$adonis_pval, method = "BH")

    df_comp <- df_comp %>%
      mutate(adonis_pval = as.numeric(adonis_pval)) %>%
      mutate(adonis_p.signif = case_when(
        adonis_pval <= 0.001 ~ "***",
        adonis_pval > 0.001 & adonis_pval <= 0.01 ~ "**",
        adonis_pval > 0.01 & adonis_pval <= 0.05 ~ "*",
        adonis_pval > 0.05 & adonis_pval <= 0.1 ~ ".",
        adonis_pval > 0.1 ~ "ns"
      )) %>%
      mutate(permutest_pval = as.numeric(permutest_pval)) %>%
      mutate(permutest_p.signif = case_when(
        permutest_pval <= 0.001 ~ "***",
        permutest_pval > 0.001 & permutest_pval <= 0.01 ~ "**",
        permutest_pval > 0.01 & permutest_pval <= 0.05 ~ "*",
        permutest_pval > 0.05 & permutest_pval <= 0.1 ~ ".",
        permutest_pval > 0.1 ~ "ns"
      )) %>%
      mutate(adonis_p.adj = as.numeric(adonis_p.adj)) %>%
      mutate(adonis_p.adj.signif = case_when(
        adonis_p.adj <= 0.001 ~ "***",
        adonis_p.adj > 0.001 & adonis_p.adj <= 0.01 ~ "**",
        adonis_p.adj > 0.01 & adonis_p.adj <= 0.05 ~ "*",
        adonis_p.adj > 0.05 & adonis_p.adj <= 0.1 ~ ".",
        adonis_p.adj > 0.1 ~ "ns"
      ))
  }


  if (add_pval) {
    pval_text <- paste0("PERMANOVA, R2=", round(adonis_res$R2[1], 2), ", p = ", round(adonis_res$`Pr(>F)`[1], 3))

    p <- p + annotate(
      geom = "text",
      label = pval_text,
      x = Inf, y = -Inf, hjust = 1, vjust = -1
    )

    if (pairwise_adonis) {
      p <- plot_grid(p,
        plot_grid(
          tableGrob(get_pairwise_matrix(df_comp, values_from = "adonis_p.adj"), rows = NULL),
          tableGrob(get_pairwise_matrix(df_comp, values_from = "adonis_p.adj.signif"), rows = NULL),
          # tableGrob(get_pairwise_matrix(df_comp, values_from = 'permutest_p.signif'), rows = NULL),
          nrow = 2,
          labels = c("adonis p.adj", "adonis p.adj.signif")
        ),
        ncol = 2
      )
    }
  }

  p <- p +
    theme(
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      plot.title = element_text(size = 12),
      axis.line = element_line(size = 0.6),
      axis.ticks = element_line(size = 0.6),
      plot.margin = unit(c(5.5, 12, 5.5, 5.5), "pt")
    )

  return(list("plot" = p, "adonis" = adonis_res, "permutest" = d, "pairwise" = df_comp))
}

plot_alpha_divers <- function(phy_obj, group_var, group_var2 = NULL, y = "value",
                              alpha_meas = c("Chao1", "Shannon", "Simpson"),
                              paired = TRUE,
                              add_pval = TRUE, plot_stats = TRUE,
                              title = "alpha-diversity", plot_paired = TRUE,
                              save = TRUE, outfile = NULL, wrap_by = "~ diversity_index",
                              color_pal = c("#C25E7B", "#72939E")) {
  suppressPackageStartupMessages(library("rstatix"))
  
  # Compute alpha-diversity measures
  p <- plot_richness(phy_obj, group_var, measures = alpha_meas)

  # Compute wilcoxon test
  if (is.null(group_var2)) {
    formula_test <- as.formula(paste0("value ~ ", group_var))
    pw_wilcoxon <- p$data %>%
      rename(diversity_index = variable) %>%
      group_by(diversity_index) %>%
      wilcox_test(formula_test, paired = paired) %>%
      add_xy_position(x = group_var, scales = "free")
    group_var2 <- group_var
  } else {
    formula_test <- as.formula(paste0("value ~ ", group_var2))
    pw_wilcoxon <- p$data %>%
      rename(diversity_index = variable) %>%
      group_by(across(all_of(c("diversity_index", group_var)))) %>%
      wilcox_test(formula_test, paired = paired) %>%
      add_xy_position(x = group_var, scales = "free")
  }

  alpha_divers <- p$data %>%
    rename(diversity_index = variable) %>%
    ggplot(aes(x = get(group_var), y = value, colour = get(group_var2))) +
    geom_boxplot(alpha = 0.5, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(seed = 3922), alpha = 0.5) +
    ggprism::theme_prism(
      base_fontface = "bold",
      base_line_size = 1,
      base_size = 28
    ) +
    theme(
      text = element_text(size = 15),
      legend.position = "top",
      legend.title = element_text()
    ) +
    ylab("Alpha-diversity Index") +
    facet_wrap(as.formula(wrap_by), scales = "free") +
    ggtitle(title) +
    scale_color_manual(values = color_pal) +
    guides(fill = guide_legend(title = group_var2))

  factor_levels <- factor(p$data[[group_var2]]) %>% levels()
  # Filter non-significant values
  if (length(factor_levels) > 2) {
    label_adj <- "p.adj.signif"
    pw_wilcoxon <- pw_wilcoxon %>%
      filter(`p.adj` < 0.05) %>%
      add_xy_position(x = group_var, scales = "free",
          dodge = 0.8, step.increase = 0.035)
  } else {
    label_adj <- "p"
    pw_wilcoxon <- pw_wilcoxon %>%
      filter(p < 0.05) %>%
      add_xy_position(x = group_var, scales = "free")
  }

  if(plot_stats == TRUE){
    af <- alpha_divers +
      ggpubr::stat_pvalue_manual(pw_wilcoxon,
        label = label_adj,
        hide.ns = T, tip.length = 0.005,
        step.group.by = "diversity_index"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
      theme(legend.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=11),
            plot.title = element_text(size=13))
  }else{
    af <- alpha_divers +
      theme(legend.title = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size=11),
            plot.title = element_text(size=13))
  }


  if (save == TRUE) {
    ggsave(filename = outfile, plot = af,
            height = 12, width = 16, units = "in",
            device = "png", scale = 0.8)
  }
  pw_wilcoxon <- apply(pw_wilcoxon,2,as.character) %>% data.frame()
  return(list("plot" = af, "tab" = pw_wilcoxon))
}

plot_alpha_divers_paired <- function(phy_obj, group_var, group_var2 = NULL, y = "value",
                              alpha_meas = c("Chao1", "Shannon", "Simpson"),
                              paired = TRUE, add_to_sample_data = TRUE,
                              add_pval = TRUE, plot_stats = TRUE,
                              title = "alpha-diversity", plot_paired = TRUE,
                              save = TRUE, outfile = NULL, wrap_by = "~ diversity_index",
                              color_pal = c("#C25E7B", "#72939E")) {

  suppressPackageStartupMessages(library("rstatix"))

  # Add alpha-diversity measure to the metadata
  p <- phyloseq::plot_richness(phy_obj, group_var, measures = alpha_meas)

  if (add_to_sample_data == TRUE) {
    for (index in alpha_meas) {
      d <- p$data %>% filter(variable == index)
      phyloseq::sample_data(phy_obj)[[index]] <- d$value
    }
  }

  # Compute wilcoxon test
  variable <- "diversity_index"
  
  # Plot
  formula_wrap <- as.formula(paste0("~ variable +", group_var))
  aplot_paired <- p$data %>%
    ggpaired(
      x = group_var2, y = "value", id = "PatientID", color = group_var2,
      line.color = "gray", label = "", font.label = list(size = 8),
      repel = TRUE,
      line.size = 0.5, palette = "jco"
    ) +
    ggprism::theme_prism() +
    # ylab(paste(measure,'Index'))+
    facet_wrap(formula_wrap, scales = "free", labeller = label_both) +
    ggtitle(title) +
    scale_color_manual(values = color_pal)

  print(aplot_paired)

  if (save == TRUE) {
    filename_pair <- paste0(outdir, prefix, "_paired.png")
    ggsave(filename = filename_pair, plot = aplot_paired,
          height = 12, width = 16, units = "in",
          device = "png", scale = 0.8)
  }
  
}
###############################################################################
# Relative abundance
###############################################################################
rel_ab <- function(x) {
  ##
  # Compute relative abundance
  ##
  if (sum(x) == 0) {
    return(x)
  } else {
    return(100 * x / sum(x))
  }
}

get_representative <- function(d, cumsum_threshold = 90) {
  # Input: data frame with raw counts

  d$Abundance <- rowSums(d)
  d$Rel_abundance <- (d$Abundance / sum(d$Abundance)) * 100
  d <- d %>% dplyr::arrange(desc(Rel_abundance))
  d$Cum_sum <- cumsum(d$Rel_abundance)

  d <- d %>%
    filter(Cum_sum <= cumsum_threshold) %>%
    select(Rel_abundance, Cum_sum)

  return(d)
}

get_lowest_level <- function(
    x, label_plot = FALSE,
    clades = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
  #' @title get lowest level
  #' @description Get the lowest taxonomical label available
  #' @param x vector with clades names
  #' @return named str
  #' @export

  # Record order
  names(x) <- clades

  # Remove empty strings
  x_filt <- x[nzchar(x)]

  # Lowest level
  lowest_level <- names(x_filt[length(x_filt)])
  lowest_label <- x_filt[length(x_filt)]

  if (label_plot == TRUE) {
    label <- paste0(lowest_level[1], ":", lowest_label)
    return(label)
  } else {
    return(lowest_level)
  }
}



get_npalette <- function(n) {
  #' @title get colors
  #' @description Obtain a color palette with x distinguishable colors
  #' @return vector with x colors
  #' @export
  suppressPackageStartupMessages(library(Polychrome))

  P <- createPalette(n, c("#ff0000", "#00ff00", "#0000ff"))
  return(P)
}

get_stacked_data <- function(phy_object, otu_names, level = 7) {
  #' @title Data stacked barplots
  #' @description Obtain dataframe from phyloseq object to plot stacked ggbarplots
  #' @param phy_object phyloseq object
  #' @param otu_names vector OTU names (complete lineage) to be labeled in the plot at level_str.
  #'    OTU names not present in the otu_names vector will be label as 'Others' in the plot.
  #' @param level_str str Taxonomical level. Eg. 'Genus', 'Species'
  #' @return data.frame
  #' @export
  clades_sub <- colnames(phyloseq::tax_table(phy_object))
  level_str <- clades_sub[as.numeric(level)]

  d <- psmelt(tax_glom(phy_object, taxrank = level_str, NArm = FALSE))

  d$lowest_label <- apply(d[clades_sub], 1, get_lowest_level,
    label_plot = TRUE, clades = clades_sub)

  d$otu_labels <- ifelse(d$OTU %in% otu_names,
    ifelse(d[[level_str]] != "", d[[level_str]], d$lowest_label),
    "Others"
  )

  # Set order
  otu <- unique(d[d$OTU %in% otu_names, c("OTU", level_str, "lowest_label")])
  rownames(otu) <- otu$OTU
  otu <- otu[otu_names, ]

  otu_order <- ifelse(otu[[level_str]] != "",
                      otu[[level_str]], otu[["lowest_label"]])
  # print(otu_order)

  d$otu_labels <- factor(d$otu_labels, levels = c(otu_order, "Others"))

  return(d)
}

get_stacked_plot <- function(
    d, title = NULL, color_pallete = NULL, ncol_legend = 4,
    legend_position = "none", circular = FALSE, rotate = FALSE,
    level_str = "Species", x_axis = "Sample", weave_factor = NULL,
    facet = NULL, facet_ncol = NULL) {

  suppressPackageStartupMessages(library("ggh4x"))

  if(!is.null(weave_factor)){
    begin <- d %>%
      arrange(c(.data[[weave_factor]])) %>%
      ggplot(aes(x = ggh4x::weave_factors(.data[[x_axis]], .data[[weave_factor]]), 
        y = Abundance, fill = otu_labels)) +
      guides(x = "axis_nested")
  }else{
    begin <- d %>%
      ggplot(aes(x = .data[[x_axis]], y = Abundance, fill = otu_labels))
  }
  stacked_plot <- begin +
    geom_bar(stat = "identity", position = position_stack(reverse = F)) +
    scale_fill_manual(level_str, values = color_pallete) +
    scale_y_continuous(limits=c(0,100), expand = c(0,0)) +
    theme_classic() +
    theme(
      legend.position = legend_position,
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      # legend.spacing.x = unit(0.8, 'cm'),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.background = element_rect(fill = "transparent", color = NA),
      axis.text.y = element_text(size = 11),
      # upper, right, bottom, left
      plot.margin = unit(c(8.5, 30.5, 8.5, 20.5), "pt"),
      axis.ticks = element_line(linewidth = 0.25),
      axis.title.y = element_blank()
      # axis.title.x = element_blank()
      # axis.line = element_blank()
    ) +
    guides(fill = guide_legend(ncol = ncol_legend)) +
    ggtitle(title) +
    xlab(x_axis)

  if (circular == TRUE) {
    stacked_plot <- stacked_plot + coord_polar(start = 0)
  } else {
    stacked_plot <- stacked_plot +
      ylab("Relative Abundance (%)")
  }


  if (!is.null(facet)){
    stacked_plot <- stacked_plot + facet_grid(cols = vars(.data[[facet]])) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11)
      )
  }

  if (rotate == TRUE){
    stacked_plot <- stacked_plot + coord_flip()
  }
  
  if (x_axis == "Sample") {
    stacked_plot <- stacked_plot +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 0.5, size = 4),
        ggh4x.axis.nestline.x = element_line(size = 1.2),
        ggh4x.axis.nesttext.x = element_text(size = 12, angle = 0, hjust = 0.5),
      )
  }

  return(stacked_plot)
}

get_boxplots <- function(
    phy_object, taxa_sub, split_by, level_str,
    get_relative = TRUE) {
  # filter taxa
  phy_pb_top <- prune_taxa(rownames(tax_table(phy_object)) %in% taxa_sub, phy_object)

  if (get_relative == TRUE) {
    phy_pb_top <- transform_sample_counts(phy_pb_top, rel_ab)
  }

  formula_wrap <- as.formula(paste0("~ ", level_str))

  # Boxplot
  boxplots <- list()
  ps <- psmelt(phy_pb_top)
  for (i in 1:length(taxa_sub)) {
    otu <- taxa_sub[i]
    bp <- ps %>%
      filter(OTU == otu) %>%
      ggplot(aes(x = get(split_by), y = Abundance, color = get(split_by))) +
      geom_boxplot(outlier.shape = NA) +
      geom_point(position = position_jitterdodge(seed = 3922), alpha = 0.5) +
      labs(x = "", y = "Relative Abundance\n") +
      ggprism::theme_prism() +
      theme(
        plot.title = element_text(size = 9),
        text = element_text(size = 18),
        strip.text = element_text(face = "bold"),
        legend.position = "top"
      ) +
      facet_wrap(formula_wrap, scales = "free", ncol = 3) +
      ggtitle(otu)

    boxplots[[i]] <- bp
  }


  return(boxplots)
}
