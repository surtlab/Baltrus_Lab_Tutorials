#!/usr/bin/env python3
"""
snp_phylogeny.py
----------------
Step 3a of the genome pipeline. Builds a whole-genome SNP-based phylogeny
using Parsnp for core-genome alignment and tree inference.

Parsnp aligns all input genomes, identifies core-genome SNPs, and outputs
a phylogenetic tree. No raw reads required — assembled FASTA files only.

Usage:
    python snp_phylogeny.py [genomes_dir] [--output snp_phylogeny]
                            [--reference strain_name] [--threads 4]

Arguments:
    genomes_dir     Folder of .fasta files from download_genomes.py.
                    Default: ./genomes
    --output        Output folder for Parsnp results. Default: ./snp_phylogeny
    --reference     Strain name (or FASTA filename stem) to use as reference.
                    Default: Parsnp auto-selects from input genomes.
    --threads       CPU threads. Default: 4
    --no-partition  Disable Parsnp's partitioning (use for small genome sets,
                    <20 genomes). Default: partitioning is on.

Output:
    snp_phylogeny/
        parsnp.tree         Newick phylogenetic tree
        parsnp.xmfa         Core genome alignment
        parsnp.vcf          SNP calls (VCF format)
        parsnp.ggr          Harvest archive (for Gingr visualization)
"""

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Characters illegal in RAxML taxon names
RAXML_ILLEGAL = set("():;,[]'") | {chr(9), chr(13), chr(32)}


def sanitize_taxon_name(name: str) -> str:
    """Replace characters illegal in RAxML/Newick taxon names."""
    result = ""
    for ch in name:
        result += "_" if ch in RAXML_ILLEGAL else ch
    # Collapse multiple consecutive underscores
    while "__" in result:
        result = result.replace("__", "_")
    return result.strip("_")


def find_tool(name: str) -> bool:
    """Check if a tool is on PATH."""
    return shutil.which(name) is not None


def check_dependencies() -> list[str]:
    """Return list of any missing required tools."""
    missing = []
    for tool in ["parsnp"]:
        if not find_tool(tool):
            missing.append(tool)
    return missing


def find_reference(genomes_dir: Path, ref_name: str | None) -> Path | None:
    """
    Locate a reference FASTA.
    If ref_name given, find the matching file in genomes_dir.
    Otherwise return None (let Parsnp auto-select).
    """
    if not ref_name:
        return None

    fastas = list(genomes_dir.glob("*.fasta"))
    # Match by stem or by partial name
    for f in fastas:
        if ref_name == f.stem or ref_name in f.stem:
            return f

    print(f"  [WARN] Reference '{ref_name}' not found in {genomes_dir}.")
    print(f"         Available genomes:")
    for f in sorted(fastas):
        print(f"           {f.stem}")
    print("         Parsnp will auto-select a reference.\n")
    return None


def already_done(output_dir: Path) -> bool:
    """Check if Parsnp has already produced a tree."""
    tree = output_dir / "parsnp.tree"
    return tree.exists() and tree.stat().st_size > 0


def run_parsnp(
    genomes_dir: Path,
    output_dir: Path,
    reference: Path | None,
    threads: int,
    no_partition: bool,
) -> tuple[bool, str]:
    """
    Run Parsnp on all FASTA files in genomes_dir.
    Copies FASTAs to a system temp directory (outside output_dir) with
    sanitized filenames and headers so RAxML never sees illegal characters.
    Returns (success, message).
    """
    import shutil as _shutil
    import tempfile as _tempfile

    fastas = list(genomes_dir.glob("*.fasta"))
    if not fastas:
        return False, f"No .fasta files found in {genomes_dir}"

    output_dir.mkdir(parents=True, exist_ok=True)

    # Use a system temp dir OUTSIDE output_dir — Parsnp scans its own output
    # dir and gets confused by any extra files or subdirs inside it.
    tmp_dir = Path(_tempfile.mkdtemp(prefix="parsnp_input_"))

    try:
        sanitized = {}  # original path -> sanitized copy path
        for f in fastas:
            clean_stem = sanitize_taxon_name(f.stem)
            dest = tmp_dir / (clean_stem + ".fasta")
            # Rewrite FASTA with sanitized sequence headers too
            with open(f) as fin, open(dest, "w") as fout:
                for line in fin:
                    if line.startswith(">"):
                        header_id = line[1:].split()[0]
                        fout.write(f">{sanitize_taxon_name(header_id)}\n")
                    else:
                        fout.write(line)
            sanitized[f] = dest

        san_reference = sanitized.get(reference) if reference is not None else None

        # When using -r ! (auto-select), ALL genomes go in the temp dir.
        # Parsnp picks one as reference internally and includes it in the tree.
        # When a specific reference is given, keep it in the dir too — Parsnp
        # needs it there even when also passing it via -r.
        # Bottom line: always put every genome in the temp dir.
        # Parsnp handles reference deduplication internally.

        # Pass the temp dir as -d — Parsnp scans it for .fasta files
        cmd = [
            "parsnp",
            "-d", str(tmp_dir),
            "-o", str(output_dir),
            "-p", str(threads),
            "--force-overwrite",
        ]

        if san_reference:
            cmd += ["-r", str(san_reference)]
        else:
            cmd += ["-r", "!"]

        if no_partition:
            cmd.append("--no-partition")

        print(f"  Running: {' '.join(cmd)}\n")

        result = subprocess.run(cmd, text=True, capture_output=True)

        if result.returncode != 0:
            err_log = output_dir / "parsnp_error.log"
            with open(err_log, "w") as ef:
                ef.write("=== STDOUT ===\n")
                ef.write(result.stdout or "(empty)\n")
                ef.write("\n=== STDERR ===\n")
                ef.write(result.stderr or "(empty)\n")
            err_lines = (result.stderr or result.stdout or "").strip().splitlines()
            err_tail = "\n        ".join(err_lines[-15:]) if err_lines else "(no output)"
            return False, (
                f"Parsnp exited with code {result.returncode}:\n"
                f"        {err_tail}\n"
                f"        Full log: {err_log}"
            )

        if not already_done(output_dir):
            return False, "Parsnp finished but parsnp.tree not found."

        vcf = output_dir / "parsnp.vcf"
        snp_count = ""
        if vcf.exists():
            with open(vcf) as f:
                n = sum(1 for line in f if not line.startswith("#"))
            snp_count = f" ({n:,} SNPs)"

        return True, f"Tree built{snp_count} \u2192 {output_dir / 'parsnp.tree'}"

    except FileNotFoundError:
        return False, (
            "parsnp not found. Install with:\n"
            "        conda install -c bioconda parsnp"
        )
    except Exception as exc:
        return False, f"Unexpected error: {exc}"

    finally:
        _shutil.rmtree(tmp_dir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(
        description="Build a whole-genome SNP phylogeny with Parsnp (pipeline step 3a).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "EXAMPLES:\n"
            "  # Minimal — uses ./genomes, auto-selects reference:\n"
            "  python snp_phylogeny.py\n\n"
            "  # Explicit reference genome:\n"
            "  python snp_phylogeny.py --reference E_coli_K12\n\n"
            "  # Custom paths and threads:\n"
            "  python snp_phylogeny.py genomes/ --output snp_phylogeny/ --threads 8\n\n"
            "OUTPUT FILES:\n"
            "  parsnp.tree    Newick tree — open in FigTree, iTOL, or any tree viewer\n"
            "  parsnp.xmfa    Core genome alignment\n"
            "  parsnp.vcf     SNP calls\n"
            "  parsnp.ggr     Harvest archive (visualize with Gingr)\n"
        ),
    )
    parser.add_argument(
        "genomes_dir", nargs="?", default="genomes",
        help="Folder of .fasta files from download_genomes.py (default: ./genomes)",
    )
    parser.add_argument(
        "--output", default="snp_phylogeny",
        help="Output folder for Parsnp results (default: ./snp_phylogeny)",
    )
    parser.add_argument(
        "--reference", default=None,
        help="Strain name or FASTA stem to use as reference (default: Parsnp auto-selects)",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="CPU threads (default: 4)",
    )
    parser.add_argument(
        "--no-partition", action="store_true",
        help="Disable Parsnp partitioning — use for small sets (<20 genomes)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run even if output already exists",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    genomes_dir = Path(args.genomes_dir)
    output_dir = Path(args.output)

    # Validate input
    if not genomes_dir.exists():
        parser.print_usage()
        print(f"\n[ERROR] Genomes folder not found: {genomes_dir}")
        print("        Run download_genomes.py first, or pass the correct path.\n")
        sys.exit(1)

    fastas = list(genomes_dir.glob("*.fasta"))
    if not fastas:
        print(f"[ERROR] No .fasta files found in: {genomes_dir}")
        sys.exit(1)

    if len(fastas) < 3:
        print(f"[WARN] Only {len(fastas)} genome(s) found — Parsnp works best with ≥3 genomes.")

    # Check dependencies
    missing = check_dependencies()
    if missing:
        print(
            f"\n[ERROR] Missing required tools: {', '.join(missing)}\n"
            "        Install with:\n"
            "            conda install -c bioconda parsnp\n"
        )
        sys.exit(1)

    # Check if already done
    if not args.force and already_done(output_dir):
        print(f"[INFO] Output already exists at {output_dir}/parsnp.tree")
        print("       Use --force to re-run.")
        sys.exit(0)

    # Resolve reference
    reference = find_reference(genomes_dir, args.reference)

    print(f"\nGenomes dir  : {genomes_dir.resolve()} ({len(fastas)} genomes)")
    print(f"Output dir   : {output_dir.resolve()}")
    print(f"Reference    : {reference.stem if reference else 'auto-select'}")
    print(f"Threads      : {args.threads}")
    print(f"Partitioning : {'off' if args.no_partition else 'on'}")
    print()

    start = time.time()
    ok, msg = run_parsnp(
        genomes_dir=genomes_dir,
        output_dir=output_dir,
        reference=reference,
        threads=args.threads,
        no_partition=args.no_partition,
    )
    elapsed = time.time() - start

    print()
    print("=" * 50)
    if ok:
        print(f"[OK  ] {msg}")
        print(f"       Time: {elapsed:.0f}s")
        print(f"\nView tree in FigTree or upload parsnp.tree to https://itol.embl.de")
    else:
        print(f"[FAIL] {msg}")
        sys.exit(1)


if __name__ == "__main__":
    main()
