# Core Gene Phylogeny — `core_gene_phylogeny.py`

Step 3b of the genome pipeline. Builds a core-gene phylogenetic tree from Bakta GFF3 annotations using a three-stage pipeline:

1. **[Roary](https://sanger-pathogens.github.io/Roary/)** — pan-genome analysis and core gene alignment
2. **ModelTest-NG** *(optional)* — substitution model selection
3. **[RAxML-NG](https://github.com/amkozlov/raxml-ng)** — maximum likelihood tree inference with bootstrapping

---

## How it works

Roary identifies orthologous genes across all input genomes, defines the core genome (genes present in ≥99% of strains), and produces a concatenated multi-FASTA alignment of all core genes. RAxML-NG then uses this alignment to infer a maximum likelihood phylogenetic tree with bootstrap support values.

This approach is more robust than SNP-based methods for divergent datasets, since it uses shared gene content rather than whole-genome alignment.

---

## Requirements

- `phylo_tree` conda environment (see main [README.md](README.md) for setup)
- Completed Bakta annotations from `annotate_genomes.py`
- The `perl-file-find-rule` module installed via cpanm (required by Roary):

```bash
cpanm File::Find::Rule
perl -e "use File::Find::Rule; print 'OK\n'"  # verify
```

---

## Usage

```bash
python core_gene_phylogeny.py [annotations_dir] [--output core_phylogeny] [options]
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `annotations_dir` | `./annotations` | Root folder of Bakta annotations from `annotate_genomes.py` |
| `--output` | `./core_phylogeny` | Output folder for all results |
| `--threads` | `4` | CPU threads for Roary and RAxML-NG |
| `--identity` | `95` | Minimum BLASTP % identity for Roary gene clustering |
| `--bootstrap` | `100` | Number of RAxML-NG bootstrap replicates |
| `--alignment` | — | Path to existing alignment — skips Roary entirely |
| `--model` | `GTR+G` | Substitution model for RAxML-NG |
| `--run-modeltest` | off | Auto-select best model with ModelTest-NG (see notes) |
| `--force` | off | Re-run even if output already exists |

---

## Examples

**Minimal — auto-detects `./annotations`:**
```bash
python core_gene_phylogeny.py
```

**Explicit paths and more threads:**
```bash
python core_gene_phylogeny.py annotations/ \
    --output core_phylogeny/ \
    --threads 8 \
    --bootstrap 200
```

**More diverse genomes — lower identity threshold:**
```bash
python core_gene_phylogeny.py --identity 80
```

**Skip Roary — use an existing alignment:**

Roary is the slowest step (~30-90 min for 25 genomes). If it has already run, the script auto-detects the alignment and skips to RAxML-NG automatically. You can also pass the path explicitly:

```bash
python core_gene_phylogeny.py \
    --alignment core_phylogeny/roary/roary_out/core_gene_alignment.aln
```

**Specify a substitution model:**
```bash
python core_gene_phylogeny.py --model GTR+I+G4
```

---

## Output

```
core_phylogeny/
    roary/
        roary_out/
            core_gene_alignment.aln     Core gene alignment ← input to RAxML-NG
            summary_statistics.txt      Pan-genome summary (core/shell/cloud counts)
            gene_presence_absence.csv   Gene presence/absence matrix
        roary_error.log                 Roary log
    raxml/
        tree.raxml.support              Best ML tree with bootstrap values ← main output
        tree.raxml.bestTree             Best ML tree without bootstrap values
        core_gene_alignment_fixed.aln   Alignment after length validation/fixing
        tree.raxml.log                  RAxML-NG full log
```

---

## Understanding the pan-genome summary

The `summary_statistics.txt` file from Roary breaks the pan-genome into four categories:

| Category | Definition | Notes |
|----------|-----------|-------|
| **Core** | Present in 99–100% of strains | Used for the phylogeny |
| **Soft core** | Present in 95–99% of strains | Nearly universal |
| **Shell** | Present in 15–95% of strains | Variable, strain-specific |
| **Cloud** | Present in <15% of strains | Rare / HGT / unique genes |

A high number of cloud genes relative to total genes indicates a diverse dataset with significant accessory genome variation.

---

## Substitution model

The default model is `GTR+G` (General Time Reversible + Gamma rate variation), which is the standard for bacterial phylogenomics and appropriate for most datasets.

**ModelTest-NG** can auto-select the best-fit model, but:
- It is slow on large alignments (>1 Mb can take hours)
- It is not compatible with Apple Silicon (osx-arm64)

To use it on supported hardware:
```bash
python core_gene_phylogeny.py --run-modeltest
```

Other commonly used models you can specify with `--model`:
- `GTR+G` — standard, appropriate for most datasets
- `GTR+I+G4` — adds invariant sites, good for alignments with conserved regions
- `HKY+G` — simpler model, faster, less flexible

---

## Visualizing the tree

**FigTree** — free desktop app, good for quick interactive viewing:
<https://github.com/rambaut/figtree/releases>

**iTOL (Interactive Tree of Life)** — web-based, publication-quality:
<https://itol.embl.de>

**Pipeline visualization script:**
```bash
Rscript visualize_tree.R \
    --tree core_phylogeny/raxml/tree.raxml.support \
    --output core_phylogeny/core_tree
```

With metadata for tip coloring:
```bash
Rscript visualize_tree.R \
    --tree core_phylogeny/raxml/tree.raxml.support \
    --output core_phylogeny/core_tree \
    --metadata strains_metadata.csv
```

---

## Notes on divergent outgroups

If your dataset includes a distantly related strain as an outgroup (e.g. a different species), it will appear with a very long branch in the unrooted tree. This is biologically correct. When visualizing, root the tree on the outgroup:

```r
# In R, after loading the tree:
tree <- root(tree, outgroup = "outgroup_strain_name", resolve.root = TRUE)
```

---

## Runtime expectations

Roary is by far the slowest step. Expected times for 25 bacterial genomes:

| Step | Time |
|------|------|
| Roary (protein extraction + BLAST + MCL) | 30–90 min |
| Roary (MAFFT core alignment) | 5–15 min |
| RAxML-NG (100 bootstraps, GTR+G) | 10–30 min |

The script automatically skips Roary on re-runs if an alignment already exists.

---

## Troubleshooting

**`Can't locate File/Find/Rule.pm`**
The Roary Perl dependency is missing. Install with cpanm:
```bash
cpanm File::Find::Rule
```

**`Roary finished but core gene alignment not found`**
Roary created a timestamped output directory. The script searches recursively — re-running will find it automatically. You can also pass the path directly with `--alignment`.

**`FASTA file does not contain equal size sequences`**
One or more sequences in the alignment are shorter than the others (missing data). The script automatically detects this, excludes sequences with <99% alignment coverage, and fixes minor length differences before running RAxML-NG.

**`modeltest-ng: illegal hardware instruction`**
ModelTest-NG is not compatible with Apple Silicon. Use the default `GTR+G` model (already the default) or specify a model with `--model`.

**Roary keeps creating timestamped directories**
Delete the `roary/` folder before re-running Roary:
```bash
rm -rf core_phylogeny/roary/
```

**`Already exists` — want to re-run:**
```bash
python core_gene_phylogeny.py --force
```
