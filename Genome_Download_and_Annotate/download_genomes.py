#!/usr/bin/env python3
"""
download_genomes.py
-------------------
Reads a two-column text file (strain name + GenBank/RefSeq accession) and
downloads the COMPLETE genome assembly for each entry as a nucleotide FASTA.

The accession can be for a single contig, a chromosome, or a whole genome —
the script resolves it to the full assembly via the NCBI Assembly database
and downloads the complete *_genomic.fna.gz, then decompresses it.

Usage:
    python download_genomes.py <input_file> [output_dir] [--email your@email.com]

Arguments:
    input_file   Two-column file (tab, comma, or whitespace separated).
                 Column 1: strain name  |  Column 2: accession number
    output_dir   Folder for output FASTA files. Default: ./genomes
    --email      Email for NCBI Entrez (required by NCBI policy).
    --delay      Seconds between NCBI requests (default 0.4; max ~3/sec without API key).
    --api-key    Optional NCBI API key (raises rate limit to 10 req/sec).

Output:
    One .fasta per strain: <output_dir>/<strain>_<accession>.fasta
    A summary log:         <output_dir>/download_log.txt
"""

import argparse
import gzip
import os
import shutil
import sys
import time
import urllib.request
from pathlib import Path

try:
    from Bio import Entrez
except ImportError:
    print(
        "\n[ERROR] Missing required dependency: biopython\n"
        "\n"
        "  Install it with:\n"
        "      pip install biopython\n"
        "  or, if using conda:\n"
        "      conda install -c conda-forge biopython\n"
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_input_file(filepath: str) -> list[tuple[str, str]]:
    """Parse the two-column input file, auto-detecting delimiter."""
    entries = []
    with open(filepath, "r") as fh:
        for line_num, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "\t" in line:
                parts = line.split("\t")
            elif "," in line:
                parts = line.split(",")
            else:
                parts = line.split()
            if len(parts) < 2:
                print(f"  [WARNING] Line {line_num} skipped (< 2 columns): {line!r}")
                continue
            entries.append((parts[0].strip(), parts[1].strip()))
    return entries


def sanitize_name(name: str) -> str:
    """Replace characters unsafe in filenames."""
    for ch in r'/\:*?"<>| ':
        name = name.replace(ch, "_")
    return name


# ---------------------------------------------------------------------------
# NCBI helpers
# ---------------------------------------------------------------------------

def accession_to_assembly_ftp(accession: str) -> tuple[str | None, str | None]:
    """
    Given any nucleotide accession (contig, chromosome, or whole genome),
    resolve it to the NCBI Assembly FTP base path.

    Returns (ftp_base_url, assembly_name) or (None, None) on failure.

    Strategy:
      1. esearch the accession in nuccore to get its internal UID.
      2. elink nuccore UID → assembly UID.
      3. esummary the assembly UID to get the FTP path.
      Prefer RefSeq (GCF) over GenBank (GCA).
    """
    # Step 1 – nuccore UID
    handle = Entrez.esearch(db="nuccore", term=accession, retmax=1)
    result = Entrez.read(handle)
    handle.close()

    if not result["IdList"]:
        return None, None
    nuc_uid = result["IdList"][0]

    # Step 2 – linked assembly UID
    handle = Entrez.elink(dbfrom="nuccore", db="assembly", id=nuc_uid)
    link_result = Entrez.read(handle)
    handle.close()

    asm_uid = None
    for linkset in link_result:
        for db_link in linkset.get("LinkSetDb", []):
            if db_link.get("LinkName") == "nuccore_assembly":
                links = db_link.get("Link", [])
                if links:
                    asm_uid = links[0]["Id"]
                    break
        if asm_uid:
            break

    if not asm_uid:
        return None, None

    # Step 3 – assembly summary → FTP path
    handle = Entrez.esummary(db="assembly", id=asm_uid, report="full")
    summary = Entrez.read(handle, validate=False)
    handle.close()

    doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
    asm_name = str(doc.get("AssemblyName", ""))

    # Prefer RefSeq; fall back to GenBank
    ftp = str(doc.get("FtpPath_RefSeq", ""))
    if not ftp or ftp == "na":
        ftp = str(doc.get("FtpPath_GenBank", ""))
    if not ftp or ftp == "na":
        return None, asm_name

    return ftp, asm_name


def download_and_decompress(ftp_base: str, dest_fasta: Path) -> None:
    """
    Download <ftp_base>/<asm>_genomic.fna.gz and decompress to dest_fasta.
    """
    asm_name = ftp_base.rstrip("/").split("/")[-1]
    gz_url = f"{ftp_base}/{asm_name}_genomic.fna.gz"
    gz_tmp = dest_fasta.with_suffix(".fna.gz")

    # Download
    urllib.request.urlretrieve(gz_url, gz_tmp)

    # Decompress
    with gzip.open(gz_tmp, "rb") as f_in, open(dest_fasta, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    gz_tmp.unlink()


# ---------------------------------------------------------------------------
# Main per-entry logic
# ---------------------------------------------------------------------------

def process_entry(strain: str, accession: str, output_dir: Path) -> tuple[bool, str]:
    """
    Resolve accession -> full assembly FTP -> download FASTA.
    Returns (success, message).
    """
    dest = output_dir / f"{sanitize_name(strain)}_{sanitize_name(accession)}.fasta"

    if dest.exists():
        if dest.stat().st_size > 0:
            return True, f"Already exists, skipped: {dest.name}"
        else:
            # Previous run left an empty/partial file — remove and re-download
            dest.unlink()
            print(f"    [WARN] Found empty file from previous run, re-downloading...")

    # Resolve to assembly FTP
    ftp_base, asm_name = accession_to_assembly_ftp(accession)

    if not ftp_base:
        # Fallback: if no assembly link exists (e.g. old/direct accessions),
        # download just the single record as FASTA
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        if not fasta_data.strip() or not fasta_data.startswith(">"):
            return False, f"No assembly found and direct fetch failed for: {accession}"
        dest.write_text(fasta_data)
        seq_count = fasta_data.count(">")
        return True, f"No assembly link; fetched single record ({seq_count} seq) → {dest.name}"

    # Download full assembly
    download_and_decompress(ftp_base, dest)

    # Count sequences in downloaded file
    seq_count = sum(1 for line in open(dest) if line.startswith(">"))
    label = asm_name if asm_name else "assembly"
    return True, f"Downloaded {label} ({seq_count} sequences) → {dest.name}"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Download complete genome assemblies from NCBI for a list of accessions.\n\n"
            "INPUT FILE FORMAT:\n"
            "  Two columns separated by a tab, comma, or whitespace.\n"
            "  Column 1: strain name   Column 2: GenBank/RefSeq accession\n"
            "  Lines beginning with '#' and blank lines are ignored.\n\n"
            "  Example:\n"
            "    E_coli_K12      NC_000913\n"
            "    Salmonella_LT2  AE006468\n"
            "    Listeria_EGD    CP012591.1\n\n"
            "NOTES:\n"
            "  Accessions for individual contigs are automatically resolved to the\n"
            "  full parent assembly. Output is one .fasta per strain.\n"
            "  Requires: pip install biopython"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input_file", help="Two-column file: strain_name<delim>accession")
    parser.add_argument("output_dir", nargs="?", default="genomes",
                        help="Output folder (default: ./genomes)")
    parser.add_argument("--email", default=None,
                        help="Email for NCBI Entrez — required by NCBI policy")
    parser.add_argument("--delay", type=float, default=0.4,
                        help="Seconds between NCBI requests (default: 0.4; max ~3/sec without API key)")
    parser.add_argument("--api-key", default=None,
                        help="NCBI API key (optional; raises rate limit to 10 req/sec)")

    # Print full help when called with no arguments
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    # Email is required — warn loudly if missing
    if not args.email:
        parser.print_usage()
        print("\n[ERROR] --email is required by NCBI policy.")
        print("        Example: --email your@institution.edu")
        print("        Run with no arguments to see full help.\n")
        sys.exit(1)

    Entrez.email = args.email
    if args.api_key:
        Entrez.api_key = args.api_key

    input_path = Path(args.input_file)
    if not input_path.exists():
        parser.print_usage()
        print(f"\n[ERROR] Input file not found: {input_path}")
        print("        Run with no arguments to see full help.\n")
        sys.exit(1)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input file : {input_path}")
    print(f"Output dir : {output_dir.resolve()}")
    print(f"NCBI email : {args.email}")
    print()

    entries = parse_input_file(str(input_path))
    if not entries:
        sys.exit("[ERROR] No valid entries found in the input file.")

    # Scan output dir for already-downloaded genomes
    existing = {f.name for f in output_dir.glob("*.fasta") if f.stat().st_size > 0}
    expected_names = {
        f"{sanitize_name(strain)}_{sanitize_name(acc)}.fasta"
        for strain, acc in entries
    }
    already_done = expected_names & existing
    to_download = len(entries) - len(already_done)

    print(f"Found {len(entries)} entr{'y' if len(entries) == 1 else 'ies'} in input file.")
    if already_done:
        print(f"  {len(already_done)} already downloaded in '{output_dir}', will skip.")
    print(f"  {to_download} to download.\n")

    log_lines = []
    success_count = 0
    fail_count = 0

    for i, (strain, accession) in enumerate(entries, 1):
        print(f"[{i}/{len(entries)}] {strain} | {accession}")
        try:
            ok, msg = process_entry(strain, accession, output_dir)
        except Exception as exc:
            ok, msg = False, f"Unexpected error: {exc}"

        status = "OK  " if ok else "FAIL"
        print(f"    [{status}] {msg}")
        log_lines.append(f"[{status}] {strain}\t{accession}\t{msg}")

        if ok:
            success_count += 1
        else:
            fail_count += 1

        if i < len(entries):
            time.sleep(args.delay)

    log_path = output_dir / "download_log.txt"
    with open(log_path, "w") as lf:
        lf.write("\n".join(log_lines) + "\n")

    print()
    print("=" * 50)
    print(f"Done. {success_count} succeeded, {fail_count} failed.")
    print(f"Log written to: {log_path}")


if __name__ == "__main__":
    main()
