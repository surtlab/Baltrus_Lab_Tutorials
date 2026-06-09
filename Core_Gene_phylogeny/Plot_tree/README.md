# Tree Visualization — `visualize_tree.R`

Step 4 of the genome pipeline. Visualizes phylogenetic trees produced by `snp_phylogeny.py` (Parsnp) or `core_gene_phylogeny.py` (RAxML-NG) using [ggtree](https://bioconductor.org/packages/release/bioc/html/ggtree.html). Produces publication-quality PDF and PNG output.

---

## Requirements

R packages — install once:

```r
install.packages("BiocManager")
BiocManager::install(c("ggtree", "treeio"))
install.packages(c("ggplot2", "dplyr", "ape", "RColorBrewer", "optparse"))
```

Or via conda (recommended — already included in `environment.yml`):

```bash
mamba env update -f environment.yml
```

---

## Usage

```bash
Rscript visualize_tree.R --tree <tree_file> --output <prefix> [options]
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `--tree` | *(required)* | Path to Newick tree file |
| `--output` | `tree_plot` | Output file prefix — produces `<prefix>.pdf` and `<prefix>.png` |
| `--type` | auto-detect | Tree type: `snp` or `core` (affects default title) |
| `--title` | auto | Custom plot title |
| `--no-bootstrap` | off | Hide bootstrap/support values on nodes |
| `--metadata` | — | Optional metadata file for tip coloring (see below) |
| `--width` | `10` | Plot width in inches |
| `--height` | `8` | Plot height in inches |

---

## Examples

**SNP tree (Parsnp):**
```bash
Rscript visualize_tree.R \
  --tree snp_phylogeny/parsnp.tree \
  --output snp_phylogeny/snp_tree \
  --type snp \
  --title "Pantoea SNP Phylogeny"
```

**Core gene tree (RAxML-NG):**
```bash
Rscript visualize_tree.R \
  --tree core_phylogeny/raxml/tree.raxml.support \
  --output core_phylogeny/core_tree \
  --type core \
  --title "Pantoea Core Gene Phylogeny"
```

**With explicit paths (recommended):**
```bash
Rscript visualize_tree.R \
  --tree /Users/surt_lab/Desktop/genome_download_and_phylogeny/SNPtrees/parsnp.tree \
  --output /Users/surt_lab/Desktop/genome_download_and_phylogeny/SNPtrees/snp_tree \
  --type snp

Rscript visualize_tree.R \
  --tree /Users/surt_lab/Desktop/genome_download_and_phylogeny/core_gene_phylogeny/raxml/tree.raxml.support \
  --output /Users/surt_lab/Desktop/genome_download_and_phylogeny/core_gene_phylogeny/core_tree \
  --type core
```

**Without bootstrap values:**
```bash
Rscript visualize_tree.R \
  --tree core_phylogeny/raxml/tree.raxml.support \
  --output core_phylogeny/core_tree \
  --no-bootstrap
```

**Custom dimensions (wider for many taxa):**
```bash
Rscript visualize_tree.R \
  --tree core_phylogeny/raxml/tree.raxml.support \
  --output core_phylogeny/core_tree \
  --width 14 --height 12
```

---

## Output

Both files are written to the directory of `--output`:

| File | Description |
|------|-------------|
| `<prefix>.pdf` | Vector PDF — best for publication, scalable to any size |
| `<prefix>.png` | 300 dpi PNG — good for presentations and reports |

---

## Bootstrap values

Bootstrap support values are shown automatically on internal nodes when present in the tree file (RAxML-NG `.support` trees). Only values ≥50 are displayed to reduce clutter. Use `--no-bootstrap` to hide them entirely.

The Parsnp tree (`parsnp.tree`) does not contain bootstrap values — the script detects this and skips the bootstrap label step automatically.

---

## Metadata / tip coloring

If you have a metadata file with strain groupings (host, location, treatment, etc.), tips can be colored by group:

```bash
Rscript visualize_tree.R \
  --tree core_phylogeny/raxml/tree.raxml.support \
  --output core_phylogeny/core_tree_colored \
  --metadata strains_metadata.csv
```

**Metadata file format** — CSV or TSV, auto-detected:

```
strain,host
PANS_1_5_NZ_JABDYS010000005.1,tomato
PNA_03_1_NZ_JABDZT010000003.1,pepper
DBL1720_ICMP_10132_NZ_JARNMU010000002.1,bean
```

The first column must match the tip labels in the tree exactly. The second column is used for coloring. Additional columns are loaded but not currently plotted.

> **Note:** Tip labels in the tree have sanitized names (parentheses and special characters replaced with underscores). Make sure your metadata strain names match the sanitized versions.

---

## Rooting the tree

The script plots trees as-is (unrooted layout in a rectangular format). To root the tree before plotting, you can do it in R directly:

**Root on an outgroup:**
```r
library(ape)
library(treeio)

tree <- read.newick("core_phylogeny/raxml/tree.raxml.support", node.label = "support")
tree <- root(tree, outgroup = "DBL1721_ATCC_35400_NZ_JARNMT010000001.1", resolve.root = TRUE)
write.tree(tree, "core_phylogeny/raxml/tree.raxml.support.rooted")
```

Then plot the rooted tree:
```bash
Rscript visualize_tree.R \
  --tree core_phylogeny/raxml/tree.raxml.support.rooted \
  --output core_phylogeny/core_tree_rooted \
  --type core
```

**Midpoint root:**
```r
library(phangorn)
tree <- midpoint(tree)
write.tree(tree, "core_phylogeny/raxml/tree.raxml.support.midpoint")
```

---

## Troubleshooting

**`Error in library(ggtree)`**
Package not installed. Run:
```r
BiocManager::install("ggtree")
```

**Tip labels cut off**
Increase `--width`:
```bash
Rscript visualize_tree.R --tree ... --output ... --width 14
```

**Bootstrap values not showing**
The tree file may not contain support values. Check with:
```bash
head -1 your_tree.tree
```
If node labels are absent or all zeros, use `--no-bootstrap`.

**Metadata tips not matching**
Tree tip labels have sanitized names — parentheses and spaces are replaced with underscores. Check an example tip label:
```bash
head -1 snp_phylogeny/parsnp.tree | grep -o '[A-Za-z0-9_\.]*' | head -5
```
