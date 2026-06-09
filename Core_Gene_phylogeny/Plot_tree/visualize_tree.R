#!/usr/bin/env Rscript
# =============================================================================
# visualize_tree.R
# ----------------
# Step 4 of the genome pipeline. Visualizes phylogenetic trees produced by
# snp_phylogeny.py (Parsnp) or core_gene_phylogeny.py (RAxML-NG) using ggtree.
#
# Usage:
#   Rscript visualize_tree.R --tree <tree_file> --output <prefix>
#                            [--metadata <file>] [--type snp|core]
#                            [--title "My Tree"]
#
# Arguments:
#   --tree        Path to Newick tree file.
#                 SNP tree:       snp_phylogeny/parsnp.tree
#                 Core gene tree: core_phylogeny/raxml/tree.raxml.support
#   --output      Output file prefix (default: tree_plot).
#                 Produces <prefix>.pdf and <prefix>.png
#   --metadata    Optional metadata file (CSV or TSV).
#                 First column: strain/tip names.
#                 Additional columns: grouping variables (e.g. host, location).
#                 If provided, tips are colored by the first grouping column.
#   --type        Tree type: 'snp' or 'core' (affects default title). Default: auto-detect.
#   --title       Custom plot title.
#   --no-bootstrap  Hide bootstrap/support values on nodes.
#   --width       Plot width in inches (default: 10)
#   --height      Plot height in inches (default: 8)
#
# Output:
#   <prefix>.pdf    Publication-quality vector PDF
#   <prefix>.png    High-resolution PNG (300 dpi)
#
# Required R packages (install once):
#   install.packages("BiocManager")
#   BiocManager::install("ggtree")
#   BiocManager::install("treeio")
#   install.packages(c("ggplot2", "dplyr", "optparse", "RColorBrewer", "ape"))
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
})

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

option_list <- list(
  make_option(c("--tree"), type = "character", default = NULL,
              help = "Path to Newick tree file [required]"),
  make_option(c("--output"), type = "character", default = "tree_plot",
              help = "Output file prefix [default: tree_plot]"),
  make_option(c("--metadata"), type = "character", default = NULL,
              help = "Optional metadata file (CSV or TSV). First column = tip names."),
  make_option(c("--type"), type = "character", default = "auto",
              help = "Tree type: 'snp' or 'core' [default: auto-detect from filename]"),
  make_option(c("--title"), type = "character", default = NULL,
              help = "Custom plot title"),
  make_option(c("--no-bootstrap"), action = "store_true", default = FALSE,
              dest = "no_bootstrap",
              help = "Hide bootstrap/support values on nodes"),
  make_option(c("--width"), type = "double", default = 10,
              help = "Plot width in inches [default: 10]"),
  make_option(c("--height"), type = "double", default = 8,
              help = "Plot height in inches [default: 8]")
)

parser <- OptionParser(
  usage = "Rscript visualize_tree.R --tree <file> --output <prefix> [options]",
  option_list = option_list,
  description = paste(
    "\nVisualize phylogenetic trees from the genome pipeline.",
    "\nProduces both PDF and PNG output.",
    "\nIf --metadata is provided, tips are colored by the first grouping column."
  )
)

opt <- parse_args(parser)

# Print help if no tree provided
if (is.null(opt$tree)) {
  print_help(parser)
  quit(status = 0)
}

# ---------------------------------------------------------------------------
# Check and load required packages
# ---------------------------------------------------------------------------

required_pkgs <- c("ggtree", "treeio", "ggplot2", "dplyr", "ape", "RColorBrewer")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  cat("\n[ERROR] Missing required R packages:", paste(missing_pkgs, collapse = ", "), "\n")
  cat("        Install with:\n")
  cat("            install.packages('BiocManager')\n")
  cat("            BiocManager::install(c('ggtree', 'treeio'))\n")
  cat("            install.packages(c('ggplot2', 'dplyr', 'ape', 'RColorBrewer', 'optparse'))\n\n")
  quit(status = 1)
}

suppressPackageStartupMessages({
  library(ggtree)
  library(treeio)
  library(ggplot2)
  library(dplyr)
  library(ape)
  library(RColorBrewer)
})

# ---------------------------------------------------------------------------
# Validate inputs
# ---------------------------------------------------------------------------

if (!file.exists(opt$tree)) {
  cat("\n[ERROR] Tree file not found:", opt$tree, "\n")
  cat("        Run snp_phylogeny.py or core_gene_phylogeny.py first.\n\n")
  quit(status = 1)
}

# ---------------------------------------------------------------------------
# Auto-detect tree type from filename
# ---------------------------------------------------------------------------

tree_type <- opt$type
if (tree_type == "auto") {
  if (grepl("parsnp", basename(opt$tree), ignore.case = TRUE)) {
    tree_type <- "snp"
  } else if (grepl("raxml|core", basename(opt$tree), ignore.case = TRUE)) {
    tree_type <- "core"
  } else {
    tree_type <- "unknown"
  }
}

# ---------------------------------------------------------------------------
# Load tree
# ---------------------------------------------------------------------------

cat("[1/4] Loading tree:", opt$tree, "\n")

tree <- tryCatch({
  # treeio::read.newick handles bootstrap values in node labels
  read.newick(opt$tree, node.label = "support")
}, error = function(e) {
  tryCatch({
    read.tree(opt$tree)
  }, error = function(e2) {
    cat("[ERROR] Could not read tree file:", conditionMessage(e2), "\n")
    quit(status = 1)
  })
})

n_tips <- length(tree$tip.label)
cat("        Tips:", n_tips, "\n")

# ---------------------------------------------------------------------------
# Load and match metadata (optional)
# ---------------------------------------------------------------------------

metadata <- NULL
group_col <- NULL

if (!is.null(opt$metadata)) {
  cat("[2/4] Loading metadata:", opt$metadata, "\n")

  if (!file.exists(opt$metadata)) {
    cat("      [WARN] Metadata file not found, proceeding without it.\n")
  } else {
    # Auto-detect delimiter
    first_line <- readLines(opt$metadata, n = 1)
    delim <- if (grepl("\t", first_line)) "\t" else ","

    metadata <- read.delim(opt$metadata, sep = delim, stringsAsFactors = FALSE,
                           check.names = FALSE)

    # Standardize: first column is the tip name column
    colnames(metadata)[1] <- "label"
    group_col <- if (ncol(metadata) >= 2) colnames(metadata)[2] else NULL

    # Report matching
    matched <- sum(metadata$label %in% tree$tip.label)
    cat("        Columns:", paste(colnames(metadata), collapse = ", "), "\n")
    cat("        Matched", matched, "of", nrow(metadata), "rows to tree tips\n")

    if (matched == 0) {
      cat("      [WARN] No metadata labels matched tree tip names.\n")
      cat("             Tree tip example:", tree$tip.label[1], "\n")
      cat("             Metadata label example:", metadata$label[1], "\n")
      metadata <- NULL
      group_col <- NULL
    }
  }
} else {
  cat("[2/4] No metadata provided — plotting tree with bootstrap values only\n")
}

# ---------------------------------------------------------------------------
# Build title
# ---------------------------------------------------------------------------

plot_title <- opt$title
if (is.null(plot_title)) {
  if (tree_type == "snp") {
    plot_title <- "Whole-genome SNP Phylogeny (Parsnp)"
  } else if (tree_type == "core") {
    plot_title <- "Core Gene Phylogeny (Roary + RAxML-NG)"
  } else {
    plot_title <- "Phylogenetic Tree"
  }
}

# ---------------------------------------------------------------------------
# Build ggtree plot
# ---------------------------------------------------------------------------

cat("[3/4] Building plot\n")

# Color palette — expand if more than 8 groups
make_palette <- function(n) {
  if (n <= 8) {
    RColorBrewer::brewer.pal(max(3, n), "Set2")[1:n]
  } else if (n <= 12) {
    RColorBrewer::brewer.pal(n, "Set3")
  } else {
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n)
  }
}

# Base tree — rectangular cladogram layout
p <- ggtree(tree, layout = "rectangular", branch.length = "branch.length") +
  theme_tree2() +
  theme(
    plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.title  = element_text(size = 11, face = "bold"),
    legend.text   = element_text(size = 10),
    legend.position = "right",
    axis.text.x   = element_text(size = 8)
  )

# Add bootstrap/support values on internal nodes (if present and not suppressed)
if (!opt$no_bootstrap) {
  has_support <- !is.null(tree$node.label) &&
                 any(nchar(trimws(tree$node.label)) > 0) &&
                 any(!is.na(suppressWarnings(as.numeric(tree$node.label))))

  if (has_support) {
    # Convert labels to numeric, round, only show >= 50
    p <- p + geom_nodelab(
      aes(label = ifelse(
        !is.na(suppressWarnings(as.numeric(label))) &
        suppressWarnings(as.numeric(label)) >= 50,
        round(as.numeric(label)),
        ""
      )),
      size = 2.8,
      hjust = 1.2,
      vjust = -0.4,
      color = "grey30"
    )
  }
}

# Add tip labels and optional group coloring
if (!is.null(metadata) && !is.null(group_col)) {
  # Merge metadata onto tree
  p <- p %<+% metadata

  groups <- unique(metadata[[group_col]])
  groups <- sort(groups[!is.na(groups)])
  pal <- make_palette(length(groups))
  names(pal) <- groups

  p <- p +
    geom_tiplab(aes(color = .data[[group_col]]),
                size = 3, offset = 0.001, align = FALSE) +
    scale_color_manual(values = pal, name = group_col, na.value = "grey50") +
    labs(
      title    = plot_title,
      subtitle = paste0(n_tips, " genomes  |  colored by: ", group_col)
    )
} else {
  p <- p +
    geom_tiplab(size = 3, offset = 0.001, color = "black") +
    labs(
      title    = plot_title,
      subtitle = paste0(n_tips, " genomes")
    )
}

# Add scale bar
p <- p + geom_treescale(fontsize = 3, linesize = 0.5, offset = 0.5)

# Expand x axis so tip labels don't get clipped
p <- p + xlim(0, max(p$data$x, na.rm = TRUE) * 1.4)

# ---------------------------------------------------------------------------
# Save output
# ---------------------------------------------------------------------------

cat("[4/4] Saving output\n")

out_pdf <- paste0(opt$output, ".pdf")
out_png <- paste0(opt$output, ".png")

# Create output directory if needed
out_dir <- dirname(opt$output)
if (nchar(out_dir) > 0 && out_dir != ".") {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

ggsave(out_pdf, plot = p, width = opt$width, height = opt$height,
       device = "pdf", useDingbats = FALSE)
cat("        PDF:", out_pdf, "\n")

ggsave(out_png, plot = p, width = opt$width, height = opt$height,
       dpi = 300, device = "png")
cat("        PNG:", out_png, "\n")

cat("\n====================================================\n")
cat("Done.\n")
cat("Open", out_pdf, "for publication-quality vector output.\n\n")
