# Downloading Example Data

This folder contains example data files for use with Baltrus Lab tutorials. This guide covers how to download individual files directly to your machine from the terminal using `wget` or `curl` — no git required.

---

## Before you start

You need either `wget` or `curl` installed. Availability depends on your operating system:

| | curl | wget |
|---|---|---|
| macOS | Pre-installed | Not included — install via conda or Homebrew |
| Windows 10/11 | Pre-installed (PowerShell and Command Prompt) | Not included — install via conda |

**Check if you have curl:**
```bash
curl --version
```

**Check if you have wget:**
```bash
wget --version
```

**Install wget via conda (works on both macOS and Windows):**
```bash
conda install -c conda-forge wget
```

**Install wget via Homebrew (macOS only):**
```bash
brew install wget
```

### Windows-specific notes

- Open **Command Prompt** or **PowerShell** to run the commands in this guide — both work fine
- `curl` is available natively on Windows 10 version 1803 and later. If you are on an older version of Windows, install it via conda: `conda install -c conda-forge curl`
- On Windows, use `dir` instead of `ls` to list files in a folder, and use `cd %HOMEPATH%` instead of `cd ~` to navigate to your home folder
- File paths on Windows use backslashes (`\`) but the download commands themselves work the same as on macOS

---

## How to download a file

### Step 1: Get the file URL from GitHub

1. Navigate to the file you want in this repository on GitHub
2. Click on the filename to open it
3. Click the **Raw** button in the top right of the file viewer
4. Copy the URL from your browser — it will look something like:
   ```
   https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/example.fasta
   ```

### Step 2: Download with curl

```bash
curl -O https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/example.fasta
```

The `-O` flag saves the file using its original filename. To save it with a different name use `-o`:

```bash
curl -o my_example.fasta https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/example.fasta
```

### Step 2 (alternative): Download with wget

```bash
wget https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/example.fasta
```

wget saves the file using its original filename by default. To save with a different name use `-O`:

```bash
wget -O my_example.fasta https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/example.fasta
```

---

## Downloading multiple files at once

If you need several files, you can run multiple commands in sequence. For example:

```bash
curl -O https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/example.fasta
curl -O https://raw.githubusercontent.com/surtlab/Baltrus_Lab_Tutorials/main/Blast-Tutorial/Example-Data/example_db.fasta
```

---

## Tips

- By default, files download to whatever directory you are currently in. Use `cd` to navigate to your desired folder before downloading.
- To check where you are: `pwd` (macOS) or `cd` with no arguments (Windows)
- To check that the file downloaded correctly: `ls -lh` (macOS) or `dir` (Windows) — this shows file sizes so you can confirm the file isn't empty

---

## Troubleshooting

**"curl: command not found" or "wget: command not found"**
See the installation options at the top of this guide.

**Downloaded file is empty or very small**
The URL may be incorrect. Make sure you are using the **Raw** URL from GitHub (starting with `raw.githubusercontent.com`) and not the regular GitHub page URL.

**"Permission denied" when saving the file**
You may not have write permission in the current directory. Navigate to your home folder first:

macOS:
```bash
cd ~
```
Windows:
```bash
cd %HOMEPATH%
```
Then re-run the download command.
