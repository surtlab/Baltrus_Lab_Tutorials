# SNP Phylogeny — `snp_phylogeny.py`

Step 3a of the genome pipeline. Builds a whole-genome SNP-based phylogenetic tree using [Parsnp](https://github.com/marbl/parsnp), which performs core-genome alignment directly from assembled FASTA files. No raw reads required.

---

## How it works

Parsnp aligns all input genomes to a reference, identifies single nucleotide polymorphisms (SNPs) across the core genome, and builds a maximum likelihood phylogenetic tree. The output tree is in Newick format, compatible with any tree viewer.

---

## Requirements

- `phylo_tree` conda environment (see main [README.md](README.md) for setup)
- Assembled genome FASTA files from `download_genomes.py`

Install Parsnp if not already in your environment:
```bash
conda install -c bioconda parsnp
```

---

## Usage

```bash
python snp_phylogeny.py [genomes_dir] [--output snp_phylogeny] [options]
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `genomes_dir` | `./genomes` | Folder of `.fasta` files from `download_genomes.py` |
| `--output` | `./snp_phylogeny` | Output folder for results |
| `--reference` | auto-select | Strain name or FASTA stem to use as the reference genome |
| `--threads` | `4` | CPU threads |
| `--no-partition` | off | Disable Parsnp partitioning — recommended for small datasets (<20 genomes) |
| `--force` | off | Re-run even if output already exists |

---

## Examples

**Minimal — auto-detects `./genomes`, Parsnp picks reference:**
```bash
python snp_phylogeny.py
```

**Explicit paths:**
```bash
python snp_phylogeny.py genomes/ --output snp_phylogeny/ --threads 8
```

**Specify a reference genome:**
```bash
python snp_phylogeny.py --reference E_coli_K12
```
The reference name is matched against the FASTA filename stem — partial matches work (e.g. `E_coli` will match `E_coli_K12_NC000913.fasta`).

**Small dataset (<20 genomes):**
```bash
python snp_phylogeny.py --no-partition
```

---

## Output

```
snp_phylogeny/
    parsnp.tree      Newick phylogenetic tree ← main output
    parsnp.xmfa      Core genome alignment (multi-FASTA)
    parsnp.vcf       SNP calls in VCF format
    parsnp.ggr       Harvest archive (visualize with Gingr)
    parsnp_error.log Full error log if the run fails
```

---

## Visualizing the tree

**FigTree** — free desktop app, good for quick interactive viewing:
<https://github.com/rambaut/figtree/releases>

**iTOL (Interactive Tree of Life)** — web-based, publication-quality:
<https://itol.embl.de>

**Pipeline visualization script:**
```bash
Rscript visualize_tree.R --tree snp_phylogeny/parsnp.tree --output snp_phylogeny/snp_tree
```

---

## Notes on reference selection

Parsnp's default auto-selection (`-r !`) picks a reference internally from the input genomes. The reference **is included** in the output tree — it's not treated differently from other genomes in the final result.

If you have an outgroup or a well-assembled type strain, specifying it as the reference with `--reference` can improve alignment quality.

---

## Notes on divergent genomes

Parsnp works best for closely related strains (same species, <5% divergence). If your dataset includes distantly related genomes:

- Parsnp's ANI filter may exclude them by default — the script uses `--curated` to force inclusion of all genomes regardless of distance.
- Very divergent strains (e.g. a different species used as an outgroup) will appear with long branch lengths in the tree. This is expected and biologically meaningful — they can be used to root the tree in your visualization script.
- For more divergent datasets, the core gene approach (`core_gene_phylogeny.py`) is generally more appropriate.

---

## Troubleshooting

**`parsnp not found`**
Parsnp is not installed or the `phylo_tree` environment is not active:
```bash
mamba activate phylo_tree
conda install -c bioconda parsnp
```

**`Less than 2 input sequences provided`**
Parsnp rejected some input files. The script copies FASTAs to a clean temp directory and sanitizes filenames/headers before running — if this error persists, check that your FASTA files are valid and non-empty.

**`Too divergent for parsnp` / very small core genome**
Your genomes are too different for whole-genome SNP alignment. Switch to the core gene approach:
```bash
python core_gene_phylogeny.py
```

**`Already exists` — want to re-run:**
```bash
python snp_phylogeny.py --force
```
