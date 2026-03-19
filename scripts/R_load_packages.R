# Load packages
library(ggplot2)
library(dplyr)
library(plyr)
library(tidyr)
library(stringr)
library(scales)
library(vegan)
library(knitr)
library(readr)
library(kableExtra)
library(randomcoloR)
library(GUniFrac)
library(ggdendro)
library(cowplot)


# Requires special installations
library(phyloseq) # BiocManager
#library(metagenomeSeq) # BiocManager

#library(btools) # github
library(pairwiseAdonis) # github
library(metagMisc) # gitHub

# Custom
# source some scripts
# source("Final_R_analysis/16S/changeSILVATaxaNames.R")
# source("Final_R_analysis/16S/R_analysis/MergeLowAbund.R")
# source("Final_R_analysis/16S/uw_unifrac.R")
# source("Final_R_analysis/16S/w_unifrac.R")

# Not installed
#library(randomForest)
#library(microbiome) # needed?
changeSILVAtaxa_w_species <- function(x) {
  # remove the D__ etc...
  tax.clean <- data.frame(row.names = row.names(x),
                          Domain = str_replace(x[,1], "d__",""),
                          Kingdom = str_replace(x[,2], "k__",""),
                          Phylum = str_replace(x[,3], "p__",""),
                          Class = str_replace(x[,4], "c__",""),
                          Order = str_replace(x[,5], "o__",""),
                          Family = str_replace(x[,6], "f__",""),
                          Genus = str_replace(x[,7], "g__",""),
                          Species = str_replace(x[,8], "s__",""),
                          stringsAsFactors = FALSE)
}



# Define a function to run pairwise.wilcox.test for each unique factor
nested_pairwise_wilcox <- function(data, measure_col, group_col, nesting_factor) {
  
  # Split the data by the unique values of the nesting factor
  split_data <- split(data, data[[nesting_factor]])
  
  # Apply the pairwise Wilcoxon test to each subset and store the results
  results <- lapply(split_data, function(subset_data) {
    test_result <- pairwise.wilcox.test(subset_data[[measure_col]], subset_data[[group_col]], p.adjust.method = "BH")
    return(test_result)
  })
  
  return(results)
}


perform_pairwise_adonis_blocking <- function(dist_object, data, grouping_var, blocking_var) {
  # Get unique levels of the blocking factor
  blocking_levels <- unique(data[[blocking_var]])
  
  # Initialize a list to store the pairwise comparison results
  pairwise_results <- list()
  
  # Loop through each level of the blocking factor and perform pairwise comparisons
  for (level in blocking_levels) {
    # Obtain the indices corresponding to the current level of the blocking factor
    indices <- which(data[[blocking_var]] == level)
    
    # Obtain the names of observations for the current level
    obs_names <- rownames(data)[indices]
    
    # Subset the data for the current level
    subset_data <- data[indices, ]
    
    # Convert dist object to a square matrix
    square_matrix <- as.matrix(dist_object)
    
    # Subset square matrix
    subset_matrix <- square_matrix[obs_names, obs_names]
    
    # Convert matrix back to dist
    subset_dist_object <- as.dist(subset_matrix)
    
    # Construct the formula for pairwise comparison
    formula <- paste("subset_dist_object ~", grouping_var)
    
    # Perform pairwise comparison for the subset
    pairwise_result <- pairwise.adonis2(
      as.formula(formula),  # subset dist and specify the grouping variable
      data = subset_data,  
      nperm = 9999,
      p.adjust.methods = "BH" # Benjamini-Hochberg correction
    )
    
    # Store the pairwise result
    pairwise_results[[as.character(level)]] <- pairwise_result
  }
  
  # Return the pairwise comparison results
  return(pairwise_results)
}

perform_pairwise_PERMDISP_blocking <- function(dist_object, data, grouping_var, blocking_var) {
  # Get unique levels of the blocking factor
  blocking_levels <- unique(data[[blocking_var]])
  
  # Initialize a list to store the pairwise comparison results
  permdisp_results <- list()
  
  # Loop through each level of the blocking factor and perform pairwise comparisons
  for (level in blocking_levels) {
    # Obtain the indices corresponding to the current level of the blocking factor
    indices <- which(data[[blocking_var]] == level)
    
    # Obtain the names of observations for the current level
    obs_names <- rownames(data)[indices]
    
    # Subset the data for the current level
    subset_data <- data[indices, ]
    
    # Convert dist object to a square matrix
    square_matrix <- as.matrix(dist_object)
    
    # Subset square matrix
    subset_matrix <- square_matrix[obs_names, obs_names]
    
    # Convert matrix back to dist
    subset_dist_object <- as.dist(subset_matrix)
    
    # Perform betadisper for the subset
    disper_result <- betadisper(
      subset_dist_object, subset_data[[grouping_var]]  # subset dist and specify the grouping variable
    )
    
    #Permute p-values for betadisper
    permdisp_result <- permutest(
      disper_result,
      permutations = 9999,
      pairwise=T
    )
    
    # Store the pairwise result
    permdisp_results[[as.character(level)]] <- permdisp_result
  }
  
  # Return the pairwise comparison results
  return(permdisp_results)
}

## Function for count plot with sig difference letters

plot_wilcox_letters <- function(df, yvar, ylab, palette,
                                test_var      = "Group",
                                alpha         = 0.05,
                                letter_nudge  = 0.05,
                                jitter_width  = 0.08,
                                letter_size   = 10) {
  
  df <- dplyr::filter(df, !is.na(.data[[yvar]]))
  y_sym   <- rlang::ensym(yvar)
  grp_sym <- rlang::sym(test_var)
  
  df[[test_var]] <- droplevels(factor(df[[test_var]]))
  groups <- levels(df[[test_var]])
  k      <- length(groups)
  if (k < 2) stop("Need ≥ 2 groups in ", test_var, ".")
  
  ## 1) pairwise Wilcoxon (BH)
  pw <- pairwise.wilcox.test(df[[yvar]], df[[test_var]],
                             p.adjust.method = "BH")
  
  ## 2) full symmetric p-value matrix
  p_full <- matrix(1, k, k, dimnames = list(groups, groups))
  rn <- rownames(pw$p.value); cn <- colnames(pw$p.value)
  for (r in rn) for (c in cn) p_full[r, c] <- p_full[c, r] <- pw$p.value[r, c]
  diag(p_full) <- 1
  
  print(pw)
  ## 3) fix undefined cells only
  nonfin <- !is.finite(p_full)
  if (any(nonfin)) {
    warning("Pairwise Wilcoxon had undefined p-values; setting those cells to 1 (no difference).")
    p_full[nonfin] <- 1
  }
  
  ## 4) compact-letter display
  letters <- multcompView::multcompLetters(p_full, threshold = alpha)$Letters
  letters <- letters[groups]  # keep panel order
  
  ## 5) y-positions
  letter_df <- tibble::tibble(!!grp_sym := groups, Letters = letters) |>
    dplyr::left_join(
      df |>
        dplyr::group_by(!!grp_sym) |>
        dplyr::summarise(ypos = max(.data[[yvar]], na.rm = TRUE) * (1 + letter_nudge),
                         .groups = "drop"),
      by = test_var
    )
  print(letter_df)
  ## 6) plot
  ggplot2::ggplot(
    df,
    ggplot2::aes(x = !!grp_sym, y = !!y_sym, fill = !!grp_sym, colour = !!grp_sym)
  ) +
    ggplot2::geom_boxplot(alpha = .30, linewidth = 1.25, outlier.shape = NA) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(data = letter_df,
                       ggplot2::aes(y = ypos, label = Letters),
                       size = letter_size, vjust = 0, colour = "black") +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_colour_manual(values = palette) +
    ggplot2::labs(y = ylab) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::scale_x_discrete(expand = c(0, 1)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.margin  = ggplot2::margin(1, 1, 1, 1, "lines"),
      panel.border = ggplot2::element_rect(colour = "black", linewidth = 2),
      legend.position = "none",
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x  = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 24),
      axis.text.y  = ggplot2::element_text(size = 16, colour = "black"),
      axis.ticks.y = ggplot2::element_line(colour = "black", linewidth = 1)
    )
}





plot_nested_wilcox_letters <- function(df,
                                       yvar,
                                       ylab,
                                       palette,
                                       group_var      = "SampleGroup",   # x-axis
                                       test_var       = "Type",          # dodge / letters
                                       compare_scope  = c("within", "global"),
                                       alpha          = 0.05,
                                       dodge_w        = 0.8,
                                       letter_nudge   = 0.08,
                                       jitter_w       = 0.15,
                                       letter_size    = 5,
                                       add_group_lines = TRUE) {
  
  compare_scope <- match.arg(compare_scope)
  
  ## ───── tidy symbols & factors ─────────────────────────────────────
  y_sym    <- rlang::ensym(yvar)
  grp_sym  <- rlang::sym(group_var)
  type_sym <- rlang::sym(test_var)
  
  df <- dplyr::filter(df, !is.na(.data[[yvar]]))
  df[[group_var]] <- droplevels(factor(df[[group_var]]))
  df[[test_var]]  <- droplevels(factor(df[[test_var]]))
  
  dodge <- ggplot2::position_dodge(width = dodge_w)
  
  ## ───── helper to build symmetric p-value matrix ──────────────────
  build_pmat <- function(pv, lv) {
    k  <- length(lv)
    pm <- matrix(1, k, k, dimnames = list(lv, lv))
    rn <- rownames(pv); cn <- colnames(pv)
    for (r in rn) {
      for (c in cn) {
        pm[r, c] <- pm[c, r] <- pv[r, c]
      }
    }
    diag(pm) <- 1
    pm[!is.finite(pm)] <- 1
    pm
  }
  
  ## ───── compute compact-letter display(s) ─────────────────────────
  if (compare_scope == "global") {
    
    df$combo <- interaction(df[[group_var]], df[[test_var]], sep = " | ")
    lv_combo <- levels(df$combo)
    
    pw <- pairwise.wilcox.test(df[[yvar]], df$combo, p.adjust.method = "BH")
    print(pw)
    
    pm <- build_pmat(pw$p.value, lv_combo)
    letters_all <- multcompView::multcompLetters(pm, threshold = alpha)$Letters
    
    letter_df <- tidyr::separate(
      tibble::tibble(combo = names(letters_all),
                     Letters = unname(letters_all)),
      combo, into = c(group_var, test_var), sep = " \\| ", remove = FALSE
    ) |>
      dplyr::left_join(
        df |>
          dplyr::group_by(!!grp_sym, !!type_sym) |>
          dplyr::summarise(ypos = max(.data[[yvar]]) * (1 + letter_nudge),
                           .groups = "drop"),
        by = c(group_var, test_var)
      )
    
  } else {  ## compare_scope == "within"
    
    letter_df <- purrr::map_dfr(
      levels(df[[group_var]]),
      function(g) {
        sub <- dplyr::filter(df, .data[[group_var]] == g)
        
        ## run Wilcoxon only if ≥2 levels:
        if (dplyr::n_distinct(sub[[test_var]]) < 2) {
          tibble::tibble(!!grp_sym  := g,
                         !!type_sym := unique(sub[[test_var]]),
                         Letters    = "",
                         ypos       = max(sub[[yvar]]) * (1 + letter_nudge))
        } else {
          pw <- pairwise.wilcox.test(sub[[yvar]], sub[[test_var]],
                                     p.adjust.method = "BH")
          cat("\n>>> Wilcoxon for", g, "\n")
          print(pw)
          
          lv <- levels(sub[[test_var]])
          pm <- build_pmat(pw$p.value, lv)
          let <- multcompView::multcompLetters(pm, threshold = alpha)$Letters[lv]
          
          tibble::tibble(
            !!grp_sym  := g,
            !!type_sym := lv,
            Letters    = unname(let),
            ypos       = max(sub[[yvar]]) * (1 + letter_nudge)
          )
        }
      })
  }
  
  ## ───── build plot ────────────────────────────────────────────────
  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = !!grp_sym,
                 y = !!y_sym,
                 fill = !!type_sym,
                 colour = !!type_sym)
  ) +
    ggplot2::geom_boxplot(position = dodge,
                          alpha = .25,
                          width  = 0.7,
                          linewidth = 1.1,
                          outlier.shape = NA) +
    ggplot2::geom_point(position = dodge,
                        size = 2,
                        alpha = 0.7) +
    ggplot2::geom_text(
      data = letter_df,
      ggplot2::aes(x = !!grp_sym,
                   y = ypos,
                   label = Letters,
                   group = !!type_sym),
      position = dodge,
      size = letter_size,
      vjust = 0,
      colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = palette, name = test_var) +
    ggplot2::scale_colour_manual(values = palette, name = test_var) +
    ggplot2::labs(y = ylab, x = NULL) +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  ## optional vertical separators
  if (add_group_lines && dplyr::n_distinct(df[[group_var]]) > 1) {
    vpos <- seq(1.5,
                length(levels(df[[group_var]])) - 0.5,
                by = 1)
    p <- p + ggplot2::geom_vline(xintercept = vpos,
                                 colour = "grey70",
                                 linetype = "dashed")
  }
  
  p
}


