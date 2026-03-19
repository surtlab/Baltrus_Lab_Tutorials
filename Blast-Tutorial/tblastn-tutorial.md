# BLAST Tutorial: Building a Database and Running tblastn

This tutorial walks through downloading example data from the Baltrus Lab GitHub repository, building a local BLAST nucleotide database, and searching it with a protein query using `tblastn`.

**About the example data:**
The nucleotide database is built from a draft genome of a *Pantoea ananatis* strain. The protein query contains sequences of the jplate and sheath proteins, which surround the tail fiber of the tailocin locus. This search is designed to identify genomic regions in *P. ananatis* that encode proteins similar to these tailocin structural components.

**What you will need:**
- BLAST+ installed ([see the BLAST installation guide](https://github.com/surtlab/Baltrus_Lab_Tutorials/tree/main/Blast-Tutorial))
- `curl` or `wget` installed ([see the example data download guide](https://github.com/surtlab/Baltrus_Lab_Tutorials/tree/main/Blast-Tutorial/Example-Data))
- A terminal (macOS) or Command Prompt / PowerShell (Windows)

---

## Step 1: Set up a working directory

Create a folder to keep your data and results organized:

```bash
mkdir blast-tutorial
cd blast-tutorial
```

---

## Step 2: Download the example data

Download the genome file (nucleotide database subject):

```bash
# macOS/Windows with curl
curl -O https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/JABEBR01.fna

# alternatively with wget
wget https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/JABEBR01.fna
```

Download the protein query file:

```bash
# macOS/Windows with curl
curl -O https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/tail_fiber_bracket_ATCC35400.prot.fa

# alternatively with wget
wget https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/tail_fiber_bracket_ATCC35400.prot.fa
```

Confirm both files downloaded correctly:

```bash
# macOS
ls -lh

# Windows
dir
```

You should see both `JABEBR01.fna` and `tail_fiber_bracket_ATCC35400.prot.fa` listed with non-zero file sizes.

---

## Step 3: Build a BLAST nucleotide database

Use `makeblastdb` to format the genome file as a searchable BLAST database:

```bash
makeblastdb -in JABEBR01.fna -dbtype nucl -out JABEBR01_db
```

**What each flag does:**

| Flag | Description |
|---|---|
| `-in JABEBR01.fna` | The input file to build the database from |
| `-dbtype nucl` | Specifies a nucleotide database (use `prot` for protein) |
| `-out JABEBR01_db` | The name prefix for the output database files |

After running this you will see several new files in your folder (`JABEBR01_db.nhr`, `JABEBR01_db.nin`, `JABEBR01_db.nsq`, etc.) — these are the database index files that BLAST uses internally.

---

## Step 4: Run tblastn

`tblastn` searches a **protein query** against a **nucleotide database** by translating the database on the fly in all six reading frames. This is useful for finding genomic regions that encode proteins similar to your query, even when no annotation exists.

```bash
tblastn -query tail_fiber_bracket_ATCC35400.prot.fa -db JABEBR01_db -out tblastn_results.txt -outfmt 6 -evalue 1e-5
```

**What each flag does:**

| Flag | Description |
|---|---|
| `-query` | Your protein query file |
| `-db` | The BLAST database you built in Step 3 |
| `-out` | Output file to write results to |
| `-outfmt 6` | Tabular output format (easy to read and parse) |
| `-evalue 1e-5` | E-value cutoff — only report hits with an E-value below this threshold |

---

## Step 5: View your results

```bash
# macOS
cat tblastn_results.txt

# Windows
type tblastn_results.txt
```

With `-outfmt 6` each line is a hit with 12 tab-separated columns:

| Column | Description |
|---|---|
| 1 | Query sequence ID |
| 2 | Subject sequence ID |
| 3 | % identity |
| 4 | Alignment length |
| 5 | Mismatches |
| 6 | Gap opens |
| 7 | Query start |
| 8 | Query end |
| 9 | Subject start |
| 10 | Subject end |
| 11 | E-value |
| 12 | Bit score |

A lower E-value indicates a more significant hit. Hits with E-values below `1e-10` are generally considered reliable.

---

## Expected results

Using the example data with `-outfmt 6` and `-evalue 1e-5`, your `tblastn_results.txt` should contain two hits:

```
sheath_ATCC35400        NZ_JABEBR010000003.1    100.000 389     0       0       389     292810  291644  0.0     747
baseplateJ_ATCC35400    NZ_JABEBR010000003.1    99.338  302     2       0       302     299408  298503  0.0     572
```

**Interpreting these results:**

The sheath protein hit shows 100% identity over the full 389 amino acid query length, and the baseplateJ protein shows 99.3% identity over 302 amino acids with only 2 mismatches — both extremely strong matches. Both hits are on contig `NZ_JABEBR010000003.1` and both have E-values of 0.0 and high bit scores, indicating highly significant matches.

Notice that for both hits the subject start coordinate is higher than the subject end coordinate (e.g. 292810 → 291644). This tells you the hits are on the **minus strand**, meaning the tailocin locus is encoded on the reverse complement strand of this contig.

The fact that both proteins hit the same contig in close proximity is consistent with what we would expect biologically — the tailocin locus is a gene cluster and its structural components should be encoded near each other in the genome.

---

## Adjusting search sensitivity

You can relax or tighten the search by changing the E-value cutoff:

```bash
# More stringent — only very strong hits
tblastn -query tail_fiber_bracket_ATCC35400.prot.fa -db JABEBR01_db -out tblastn_results_strict.txt -outfmt 6 -evalue 1e-10

# More permissive — catches weaker/more distant hits
tblastn -query tail_fiber_bracket_ATCC35400.prot.fa -db JABEBR01_db -out tblastn_results_permissive.txt -outfmt 6 -evalue 1e-3
```

---

## Troubleshooting

**"makeblastdb: command not found"**
BLAST is not in your PATH. If you installed via conda, make sure your environment is activated: `conda activate blast`

**"No hits found" or empty results file**
Try relaxing the E-value cutoff to `1e-3` and rerun. If still no hits, double-check that the database built correctly by confirming the `.nhr`, `.nin`, and `.nsq` files exist in your working directory.

**Results file is empty but no error was shown**
Confirm you are using the Raw GitHub URLs (starting with `raw.githubusercontent.com`) when downloading — downloading the regular GitHub page URL will give you an HTML file instead of sequence data.
