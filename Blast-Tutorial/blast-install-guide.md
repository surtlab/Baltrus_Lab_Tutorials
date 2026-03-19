# Installing NCBI BLAST on Windows and macOS

BLAST (Basic Local Alignment Search Tool) is a suite of command-line programs developed by NCBI for comparing nucleotide and protein sequences against databases. This guide covers downloading and installing the **BLAST+ command-line tools** on Windows and macOS.

---

## Windows

### Step 1: Download the installer

1. Go to the NCBI BLAST download page:
   [https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
2. Download the file ending in `win64.exe`, for example:
   `ncbi-blast-2.16.0+-win64.exe`

### Step 2: Run the installer

1. Double-click the downloaded `.exe` file
2. Accept the license agreement
3. Choose an install location (default is `C:\Program Files\NCBI\blast-2.16.0+`)
4. Click **Install**, then **Finish**

### Step 3: Add BLAST to your PATH

To run BLAST from any terminal window, add the `bin` folder to your system PATH:

1. Open the Start menu and search for **Environment Variables**
2. Click **Edit the system environment variables**
3. Click **Environment Variables** at the bottom
4. Under **System variables**, find and select **Path**, then click **Edit**
5. Click **New** and add the path to the BLAST `bin` folder, e.g.:
   ```
   C:\Program Files\NCBI\blast-2.16.0+\bin
   ```
6. Click **OK** on all dialogs to save

### Step 4: Verify the installation

Open a new **Command Prompt** or **PowerShell** window and run:

```bash
blastn -version
```

You should see output like `blastn: 2.16.0+`.

---

## macOS

### Option A: Download directly from NCBI (recommended)

#### Step 1: Download the package

1. Go to the NCBI BLAST download page:
   [https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
2. Download the correct file for your Mac:
   - **Apple Silicon (M1/M2/M3):** `ncbi-blast-2.16.0+-aarch64-apple-macosx.tar.gz`
   - **Intel Mac:** `ncbi-blast-2.16.0+-x64-macosx.tar.gz`

Not sure which chip you have? Click the Apple menu → **About This Mac** and look under Chip or Processor.

#### Step 2: Extract and install

Open Terminal and run:

```bash
# Move to your Downloads folder
cd ~/Downloads

# Extract the archive (update filename to match your download)
tar -xzf ncbi-blast-2.16.0+-x64-macosx.tar.gz

# Move to a standard location
sudo mv ncbi-blast-2.16.0+ /usr/local/ncbi-blast
```

#### Step 3: Add BLAST to your PATH

Add the BLAST `bin` directory to your shell profile so it's available in every terminal session:

```bash
# For zsh (default on macOS Catalina and later)
echo 'export PATH="/usr/local/ncbi-blast/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc

# For bash (older macOS or if you use bash)
echo 'export PATH="/usr/local/ncbi-blast/bin:$PATH"' >> ~/.bash_profile
source ~/.bash_profile
```

### Option B: Install via Homebrew

If you have [Homebrew](https://brew.sh) installed, this is a one-liner:

```bash
brew install blast
```

Homebrew handles the PATH automatically.

### Step 4: Verify the installation

```bash
blastn -version
```

You should see output like `blastn: 2.16.0+`.

---

## Optional: Install via conda

If you already have conda (Miniconda or Anaconda) installed, this is the quickest way to get BLAST running on either Windows or macOS — no manual PATH setup required.

BLAST is available on the **bioconda** channel, which is a conda channel maintained specifically for bioinformatics tools.

### Step 1: Add the required channels (one-time setup)

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### Step 2: Install BLAST

It is strongly recommended to install BLAST into a dedicated environment rather than your base environment:

```bash
# Create a new environment and install BLAST in one step
conda create -n blast -c bioconda blast

# Activate the environment
conda activate blast
```

Or if you prefer to install into an existing environment:

```bash
conda install -c bioconda blast
```

### Step 3: Verify the installation

```bash
blastn -version
```

You should see output like `blastn: 2.16.0+`.

### Notes

- conda handles all PATH configuration automatically within the active environment — no manual edits to `.zshrc` or system variables needed
- To use BLAST in future sessions, just run `conda activate blast` first
- To update BLAST later: `conda update -c bioconda blast`

---

## Running your first BLAST search

Once installed, here is a basic example searching a nucleotide query against a local database:

```bash
# Run a nucleotide BLAST search
blastn -query my_sequence.fasta -db nt -out results.txt -outfmt 6 -evalue 1e-5

# Run a protein BLAST search
blastp -query my_protein.fasta -db nr -out results.txt -outfmt 6

# Run a remote search against NCBI servers (no local database needed)
blastn -query my_sequence.fasta -db nt -remote -out results.txt
```

**Common output formats (`-outfmt`):**

| Value | Format |
|---|---|
| `0` | Default pairwise alignment |
| `6` | Tabular (easy to parse) |
| `7` | Tabular with comments |
| `10` | CSV |

---

## BLAST programs quick reference

| Program | Query | Database | Use case |
|---|---|---|---|
| `blastn` | Nucleotide | Nucleotide | DNA vs DNA |
| `blastp` | Protein | Protein | Protein vs protein |
| `blastx` | Nucleotide | Protein | Translated DNA vs protein |
| `tblastn` | Protein | Nucleotide | Protein vs translated DNA |
| `tblastx` | Nucleotide | Nucleotide | Translated DNA vs translated DNA |

---

## Troubleshooting

**"command not found" after installation**
The BLAST `bin` directory is likely not in your PATH. Re-check Step 3 for your platform and make sure you opened a new terminal window after editing your profile.

**macOS security warning ("cannot be opened because the developer cannot be verified")**
Go to **System Settings → Privacy & Security** and click **Allow Anyway** next to the blocked item. This is a Gatekeeper warning for unsigned binaries and is expected for NCBI tools.

**Remote searches timing out**
NCBI rate-limits remote BLAST queries. For repeated or large searches, download a local database instead.
