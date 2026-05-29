# Genome Download & Phylogeny Pipeline

A modular pipeline for downloading bacterial genomes from NCBI, annotating them with Bakta, and (steps coming soon) building a phylogenetic tree.

---

## Pipeline Overview

| Step | Script | Description |
|------|--------|-------------|
| 1 | `download_genomes.py` | Download complete genome assemblies from NCBI |
| 2 | `annotate_genomes.py` | Annotate genomes with Bakta |
| 3 | *(coming soon)* | Phylogenetic tree inference |

---

## Requirements

- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Mambaforge](https://github.com/conda-forge/miniforge)
- Internet connection (for genome downloads and database setup)
- ~2 GB disk space minimum (light Bakta database)

> **Windows users:** Bakta is not supported on Windows natively. Use [WSL2](https://docs.microsoft.com/en-us/windows/wsl/install) and run everything from the WSL2 terminal.

---

## Setup

### 1. Clone or download the pipeline scripts

Place all scripts in the same working directory, e.g.:

```
genome_download_and_phylogeny/
    download_genomes.py
    annotate_genomes.py
    setup_env.py
    environment.yml
    README.md
```

### 2. Create the conda environment

The setup script auto-detects your platform and installs the correct package versions:

```bash
python setup_env.py
```

This handles platform-specific version pinning automatically:

| Platform | Notes |
|----------|-------|
| Apple Silicon (osx-arm64) | Installs bakta 1.11.x + pyhmmer 0.10.x |
| Intel Mac (osx-64) | Installs latest available bakta |
| Linux x86-64 | Installs latest available bakta |
| Windows | Exits with WSL2 instructions |

To update or remove the environment later:

```bash
python setup_env.py --update   # update packages
python setup_env.py --remove   # remove environment
```

### 3. Activate the environment

```bash
mamba activate phylo_pipeline
# or
conda activate phylo_pipeline
```

---

## Step 1 — Download Genomes

### Input file format

A plain text file with two columns (tab, comma, or whitespace separated):

```
# Column 1: strain name    Column 2: GenBank/RefSeq accession
E_coli_K12                 NC_000913
Salmonella_LT2             AE006468
Listeria_EGD               CP012591.1
```

- Lines starting with `#` and blank lines are ignored
- Accessions for **individual contigs** are automatically resolved to the full parent assembly — you don't need to look up assembly accessions manually

### Usage

```bash
python download_genomes.py <input_file> [output_dir] --email your@email.com
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `input_file` | *(required)* | Two-column strain/accession file |
| `output_dir` | `./genomes` | Folder for downloaded FASTA files |
| `--email` | *(required)* | Your email — required by NCBI policy |
| `--delay` | `0.4` | Seconds between NCBI requests |
| `--api-key` | — | NCBI API key (raises rate limit to 10 req/sec) |

### Example

```bash
python download_genomes.py strains.txt --email you@institution.edu
```

Or with explicit output directory:

```bash
python download_genomes.py strains.txt ./genomes --email you@institution.edu
```

### Output

```
genomes/
    StrainName_Accession.fasta    # one file per strain (full assembly)
    download_log.txt              # summary of successes and failures
```

### Notes

- Already-downloaded genomes are skipped automatically (safe to re-run after interruption)
- Empty/partial files from interrupted downloads are detected and re-downloaded
- Prefers RefSeq (GCF) assemblies over GenBank (GCA) when both exist

---

## Step 2 — Annotate Genomes

Runs [Bakta](https://github.com/oschwengers/bakta) on every `.fasta` file in the genomes folder.

### Bakta database

The database is **not** included with the conda install and must be downloaded separately. The script handles this automatically on first run — if no database is found it downloads the light database (~1.3 GB) to `~/.bakta_db/db-light`.

To download manually:

```bash
# Light database (~1.3 GB compressed, good for most use cases)
bakta_db download --output ~/bakta_db --type light

# Full database (~32 GB, better annotation for novel/understudied taxa)
bakta_db download --output ~/bakta_db --type full
```

Set the path permanently so you never need to pass `--db`:

```bash
# Add to ~/.zshrc or ~/.bash_profile
export BAKTA_DB=~/bakta_db/db-light
```

### AMRFinderPlus

After downloading the database, initialise the AMRFinderPlus sub-database (required once):

```bash
amrfinder_update --force_update --database ~/bakta_db/db-light/amrfinderplus-db
```

If this fails, the pipeline will continue with AMR detection skipped — all other annotation features work normally.

### Usage

```bash
python annotate_genomes.py [genomes_dir] [--db /path/to/bakta/db]
```

**Arguments:**

| Argument | Default | Description |
|----------|---------|-------------|
| `genomes_dir` | `./genomes` | Folder of `.fasta` files from step 1 |
| `--db` | auto-detect | Path to Bakta database |
| `--output` | `./annotations` | Root folder for annotation results |
| `--threads` | `4` | CPU threads per Bakta job |
| `--force` | off | Re-annotate genomes that already have output |
| `--no-skip-amr` | — | Enable AMR detection (requires AMRFinderPlus db) |

### Example

```bash
# Minimal — database auto-detected or downloaded
python annotate_genomes.py

# With explicit paths
python annotate_genomes.py genomes/ \
    --db ~/bakta_db/db-light \
    --output annotations/ \
    --threads 8
```

### Output

```
annotations/
    StrainName_Accession/
        StrainName_Accession.gbff    # GenBank flat file
        StrainName_Accession.gff3    # GFF3 annotation
        StrainName_Accession.faa     # Protein FASTA (CDS amino acid sequences)
        StrainName_Accession.tsv     # Annotation summary table
        StrainName_Accession.json    # Machine-readable annotation
        ...                          # Other Bakta outputs
    annotation_log.txt               # Summary of successes and failures
```

### Notes

- Already-annotated genomes are skipped (checks for non-empty `.gbff`, `.gff3`, and `.faa`)
- Use `--force` to re-annotate everything
- If a run fails partway through, simply re-run — completed genomes will be skipped

---

## Running the Full Pipeline

```bash
# Activate the environment
mamba activate phylo_pipeline

# Step 1 — download genomes
python download_genomes.py strains.txt --email you@institution.edu

# Step 2 — annotate
python annotate_genomes.py --db ~/bakta_db/db-light --threads 8
```

---

## Troubleshooting
**`No NCBI account`**
Sign up for an account at [NCBI Signup](https://www.ncbi.nlm.nih.gov/myncbi/)

**`ModuleNotFoundError: No module named 'Bio'`**
The conda environment is not activated. Run `mamba activate phylo_pipeline` first.

**`bakta: error: unrecognized arguments`**
Your bakta version may be incompatible with the arguments used. Re-create the environment with `python setup_env.py --remove` then `python setup_env.py`.

**`ERROR: wrong database version detected`**
The installed bakta version requires a newer database. Download a fresh database:
```bash
bakta_db download --output ~/bakta_db_new --type light
python annotate_genomes.py --db ~/bakta_db_new/db-light
```

**`AttributeError: 'str' object has no attribute 'decode'`**
Bakta and pyhmmer version mismatch. Re-create the environment:
```bash
python setup_env.py --remove
python setup_env.py
```

**AMRFinderPlus errors**
AMR gene detection is skipped by default and does not affect core annotation. To enable it, set up the AMRFinderPlus database:
```bash
amrfinder_update --force_update --database ~/bakta_db/db-light/amrfinderplus-db
python annotate_genomes.py --no-skip-amr
```
