#!/usr/bin/env python3
"""
core_gene_phylogeny.py
----------------------
Step 3b of the genome pipeline. Builds a core-gene phylogeny from Bakta
GFF3 annotations using:

    Roary       → pan-genome analysis, core gene alignment
    ModelTest-NG → best-fit substitution model selection
    RAxML-NG    → maximum likelihood tree inference

Usage:
    python core_gene_phylogeny.py [annotations_dir] [--output core_phylogeny]
                                  [--threads 4] [--identity 95]
                                  [--bootstrap 100]

Arguments:
    annotations_dir   Root folder of Bakta annotations from annotate_genomes.py.
                      Default: ./annotations
    --output          Output folder. Default: ./core_phylogeny
    --threads         CPU threads for Roary and RAxML-NG. Default: 4
    --identity        Minimum BLASTP % identity for Roary clustering. Default: 95
    --bootstrap       Number of bootstrap replicates for RAxML-NG. Default: 100
    --force           Re-run even if output already exists.

Output:
    core_phylogeny/
        roary/                          Roary pan-genome output
            core_gene_alignment.aln     Core gene alignment (FASTA)
            summary_statistics.txt      Pan-genome summary
            gene_presence_absence.csv   Presence/absence matrix
        modeltest/
            modeltest.out               Model selection results
            modeltest.log               Full ModelTest-NG log
        raxml/
            tree.raxml.bestTree         Best ML tree (Newick)
            tree.raxml.support          Tree with bootstrap support
            tree.raxml.log              RAxML-NG log
"""

import argparse
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# Characters illegal in shell commands / Roary internal processing
ILLEGAL_CHARS = set('():;,[]') | {chr(34), chr(39), chr(9), chr(13), chr(32)}


def sanitize_name(name: str) -> str:
    """Replace characters that cause shell/Roary issues."""
    result = ""
    for ch in name:
        result += "_" if ch in ILLEGAL_CHARS else ch
    while "__" in result:
        result = result.replace("__", "_")
    return result.strip("_")


def find_tool(name: str) -> bool:
    return shutil.which(name) is not None


def check_dependencies() -> list[str]:
    missing = []
    for tool in ["roary", "modeltest-ng", "raxml-ng"]:
        if not find_tool(tool):
            missing.append(tool)
    return missing


def collect_gff_files(annotations_dir: Path) -> list[Path]:
    """
    Find all .gff3 files produced by Bakta.
    Bakta puts them at: annotations/<strain>/<strain>.gff3
    """
    gffs = sorted(annotations_dir.rglob("*.gff3"))
    # Exclude any empty files
    return [g for g in gffs if g.stat().st_size > 0]


def already_done(output_dir: Path) -> bool:
    """Check if a final support tree exists."""
    support = output_dir / "raxml" / "tree.raxml.support"
    best = output_dir / "raxml" / "tree.raxml.bestTree"
    return (support.exists() and support.stat().st_size > 0) or \
           (best.exists() and best.stat().st_size > 0)


# ---------------------------------------------------------------------------
# Step 1: Roary
# ---------------------------------------------------------------------------

def run_roary(
    gff_files: list[Path],
    output_dir: Path,
    threads: int,
    identity: int,
) -> tuple[bool, Path | None]:
    """
    Run Roary to compute pan-genome and core gene alignment.
    Copies GFF files to a temp directory with sanitized names (.gff extension)
    so parentheses and other special characters don't break Roary's shell calls.
    Returns (success, path_to_alignment).
    """
    import shutil as _shutil
    import tempfile as _tempfile

    output_dir.mkdir(parents=True, exist_ok=True)

    # Copy GFF3 files to a temp dir with:
    # 1. Sanitized filenames (no parentheses or special chars)
    # 2. .gff extension (Roary works more reliably with .gff than .gff3)
    tmp_dir = Path(_tempfile.mkdtemp(prefix="roary_input_"))
    name_map = {}  # sanitized stem -> original stem (for output remapping)

    try:
        for gff in gff_files:
            clean = sanitize_name(gff.stem)
            dest = tmp_dir / (clean + ".gff")
            _shutil.copy2(gff, dest)
            name_map[clean] = gff.stem

        # Remove existing roary_out dir so Roary doesn't create a timestamped copy
        import shutil as _shutil2
        roary_out = output_dir / "roary_out"
        if roary_out.exists():
            _shutil2.rmtree(roary_out)
        roary_out.mkdir()

        cmd = [
            "roary",
            "-f", str(roary_out),
            "-e",           # multiFASTA alignment of core genes
            "--mafft",      # use MAFFT (faster than PRANK)
            "-p", str(threads),
            "-i", str(identity),
            "-v",
            str(tmp_dir / "*.gff"),  # glob — roary expands this
        ]

        # Roary wants a glob or list of files; pass each explicitly
        cmd = [
            "roary",
            "-f", str(roary_out),
            "-e",
            "--mafft",
            "-p", str(threads),
            "-i", str(identity),
            "-v",
        ] + [str(f) for f in sorted(tmp_dir.glob("*.gff"))]

        print(f"  Running Roary on {len(gff_files)} GFF files (sanitized names)...")
        print(f"  Command: roary -f {roary_out} -e --mafft -p {threads} -i {identity} *.gff\n")

        # Set PERL5LIB explicitly so Roary's internal scripts find all modules
        # regardless of how the conda environment was activated
        import os as _os
        import shutil as _sh
        perl_bin = _sh.which("perl") or ""
        conda_prefix = str(Path(perl_bin).parent.parent) if perl_bin else ""
        perl_env = dict(_os.environ)
        if conda_prefix:
            perl5lib = ":".join([
                f"{conda_prefix}/lib/perl5/site_perl/5.22.0/darwin-thread-multi-2level",
                f"{conda_prefix}/lib/perl5/site_perl/5.22.0",
                f"{conda_prefix}/lib/perl5/5.22.0/darwin-thread-multi-2level",
                f"{conda_prefix}/lib/perl5/5.22.0",
            ])
            existing = perl_env.get("PERL5LIB", "")
            perl_env["PERL5LIB"] = f"{perl5lib}:{existing}" if existing else perl5lib

        result = subprocess.run(cmd, text=True, capture_output=True, env=perl_env)

        # Always write log
        err_log = output_dir / "roary_error.log"
        with open(err_log, "w") as ef:
            ef.write("=== STDOUT ===\n")
            ef.write(result.stdout or "(empty)\n")
            ef.write("\n=== STDERR ===\n")
            ef.write(result.stderr or "(empty)\n")

        if result.returncode != 0:
            err_lines = (result.stderr or result.stdout or "").strip().splitlines()
            err_tail = "\n        ".join(err_lines[-15:]) if err_lines else "(no output)"
            print(f"  [FAIL] Roary error:\n        {err_tail}")
            return False, None

        # Find the core gene alignment — Roary often creates a timestamped subdir
        # Search recursively under roary_out and also the parent output_dir
        search_roots = [roary_out, output_dir]
        for root in search_roots:
            if not root.exists():
                continue
            # Check root itself and all immediate subdirs
            candidates = [root] + [d for d in root.iterdir() if d.is_dir()]
            for search_dir in sorted(candidates, key=lambda p: p.stat().st_mtime, reverse=True):
                for name in ["core_gene_alignment.aln", "core_alignment.aln"]:
                    aln = search_dir / name
                    if aln.exists() and aln.stat().st_size > 0:
                        summary = search_dir / "summary_statistics.txt"
                        if summary.exists():
                            print(summary.read_text())
                        return True, aln

        # Last resort: recursive search anywhere under output_dir
        for aln in sorted(output_dir.rglob("core_gene_alignment.aln")) + sorted(output_dir.rglob("core_alignment.aln")):
            if aln.stat().st_size > 0:
                return True, aln

        print("  [FAIL] Roary finished but core gene alignment not found.")
        print(f"  Contents of {roary_out}:")
        if roary_out.exists():
            for f in sorted(roary_out.rglob("*.aln"))[:10]:
                print(f"    {f}")
        return False, None

    except FileNotFoundError:
        print("  [FAIL] roary not found. Install: conda install -c bioconda roary")
        return False, None
    except Exception as exc:
        print(f"  [FAIL] Unexpected error: {exc}")
        return False, None
    finally:
        _shutil.rmtree(tmp_dir, ignore_errors=True)


# ---------------------------------------------------------------------------
# Step 2: ModelTest-NG
# ---------------------------------------------------------------------------

def run_modeltest(
    alignment: Path,
    output_dir: Path,
    threads: int,
) -> tuple[bool, str]:
    """
    Run ModelTest-NG to select the best substitution model.
    Returns (success, model_string).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(output_dir / "modeltest")

    cmd = [
        "modeltest-ng",
        "--input", str(alignment),
        "--output", prefix,
        "--threads", str(threads),
        "--datatype", "nt",
        "--input-format", "fasta",
        "--force",
    ]

    print(f"  Running ModelTest-NG to select best substitution model...")

    try:
        result = subprocess.run(cmd, text=True, capture_output=True)

        # Save full log regardless
        with open(f"{prefix}.log", "w") as lf:
            lf.write(result.stdout or "")
            lf.write(result.stderr or "")

        if result.returncode != 0:
            err_lines = (result.stderr or result.stdout or "").strip().splitlines()
            err_tail = "\n        ".join(err_lines[-10:]) if err_lines else "(no output)"
            print(f"  [WARN] ModelTest-NG failed:\n        {err_tail}")
            print("  [WARN] Falling back to GTR+G model.\n")
            return False, "GTR+G"

        # ModelTest-NG writes results to <prefix>.out file, not stdout
        # Try reading from the output file first, fall back to stdout
        out_text = result.stdout or ""
        out_file = Path(f"{prefix}.out")
        if out_file.exists():
            out_text = out_file.read_text()

        # Parse best model from output (BIC criterion)
        model = parse_best_model(out_text)
        if not model:
            # Also check the log file
            log_file = Path(f"{prefix}.log")
            if log_file.exists():
                model = parse_best_model(log_file.read_text())

        if not model:
            print("  [WARN] Could not parse best model from ModelTest-NG output.")
            print("  [WARN] Falling back to GTR+G.\n")
            return False, "GTR+G"

        print(f"  Best model (BIC): {model}\n")
        return True, model

    except FileNotFoundError:
        print("  [WARN] modeltest-ng not found — falling back to GTR+G.\n")
        return False, "GTR+G"
    except Exception as exc:
        print(f"  [WARN] ModelTest-NG error: {exc} — falling back to GTR+G.\n")
        return False, "GTR+G"


def parse_best_model(output: str) -> str | None:
    """
    Parse the best-fit model (BIC) from ModelTest-NG stdout.
    Looks for lines like: 'Best model according to BIC: GTR+I+G4'
    """
    for line in output.splitlines():
        if "best model according to bic" in line.lower():
            parts = line.strip().split()
            if parts:
                return parts[-1]
    # Also try AIC as fallback
    for line in output.splitlines():
        if "best model according to aic" in line.lower():
            parts = line.strip().split()
            if parts:
                return parts[-1]
    return None


# ---------------------------------------------------------------------------
# Alignment validation
# ---------------------------------------------------------------------------

def fix_alignment(alignment: Path, output_dir: Path, length_threshold: float = 0.99) -> Path:
    """
    Check that all sequences in the alignment are the same length.
    Sequences shorter than length_threshold * max_length are excluded
    (they have too much missing data to be reliable in the tree).
    Sequences only slightly shorter are padded with gap characters.
    Returns path to the (possibly fixed) alignment.
    """
    records = {}
    current = None
    with open(alignment) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                current = line[1:].split()[0]
                records[current] = []
            elif current:
                records[current].append(line)

    seqs = {k: "".join(v) for k, v in records.items()}
    lengths = {k: len(v) for k, v in seqs.items()}
    unique_lengths = set(lengths.values())

    if len(unique_lengths) == 1:
        print(f"  Alignment OK: {len(seqs)} sequences, {list(unique_lengths)[0]} bp each")
        return alignment

    max_len = max(unique_lengths)
    min_len = min(unique_lengths)
    threshold = int(max_len * length_threshold)

    print(f"  [WARN] Unequal sequence lengths detected ({min_len} - {max_len} bp)")

    excluded = {k for k, l in lengths.items() if l < threshold}
    if excluded:
        print(f"  Excluding {len(excluded)} sequence(s) with <{length_threshold*100:.0f}% alignment coverage:")
        for name in sorted(excluded):
            pct = lengths[name] / max_len * 100
            print(f"    {name}: {lengths[name]} bp ({pct:.1f}%)")

    kept = {k: v for k, v in seqs.items() if k not in excluded}
    # Pad any remaining minor length differences
    kept_max = max(len(v) for v in kept.values())
    padded_count = sum(1 for v in kept.values() if len(v) < kept_max)
    if padded_count:
        print(f"  Padding {padded_count} sequence(s) with minor length differences...")

    output_dir.mkdir(parents=True, exist_ok=True)
    fixed = output_dir / "core_gene_alignment_fixed.aln"
    with open(fixed, "w") as f:
        for name, seq in kept.items():
            padded = seq + "-" * (kept_max - len(seq))
            f.write(f">{name}\n{padded}\n")

    print(f"  Fixed alignment: {len(kept)} sequences, {kept_max} bp")
    print(f"  Written to: {fixed.name}")
    return fixed


# ---------------------------------------------------------------------------
# Step 3: RAxML-NG
# ---------------------------------------------------------------------------

def run_raxml(
    alignment: Path,
    model: str,
    output_dir: Path,
    threads: int,
    bootstrap: int,
) -> tuple[bool, str]:
    """
    Run RAxML-NG for ML tree inference with bootstrapping.
    Returns (success, message).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = str(output_dir / "tree")

    cmd = [
        "raxml-ng",
        "--all",                    # ML search + bootstrapping in one run
        "--msa", str(alignment),
        "--msa-format", "FASTA",    # Roary outputs FASTA format alignment
        "--data-type", "DNA",
        "--model", model,
        "--prefix", prefix,
        "--threads", str(threads),
        "--bs-trees", str(bootstrap),
        "--force", "model_lh_impr",  # suppress convergence warnings on short alignments
    ]

    print(f"  Running RAxML-NG (model: {model}, bootstrap: {bootstrap})...")

    try:
        result = subprocess.run(cmd, text=True, capture_output=True)

        # Save log
        with open(f"{prefix}.raxml.log", "w") as lf:
            lf.write(result.stdout or "")
            lf.write(result.stderr or "")

        if result.returncode != 0:
            err_lines = (result.stderr or result.stdout or "").strip().splitlines()
            err_tail = "\n        ".join(err_lines[-15:]) if err_lines else "(no output)"
            return False, (
                f"RAxML-NG exited with code {result.returncode}:\n"
                f"        {err_tail}\n"
                f"        Full log: {prefix}.raxml.log"
            )

        # Find output tree
        support = Path(f"{prefix}.raxml.support")
        best = Path(f"{prefix}.raxml.bestTree")
        tree_path = support if support.exists() else best

        if not tree_path.exists():
            return False, "RAxML-NG finished but no tree file found."

        return True, f"Tree with bootstrap support → {tree_path}"

    except FileNotFoundError:
        return False, (
            "raxml-ng not found. Install with:\n"
            "        conda install -c bioconda raxml-ng"
        )
    except Exception as exc:
        return False, f"Unexpected error: {exc}"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Build a core-gene phylogeny with Roary + ModelTest-NG + RAxML-NG (step 3b).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "EXAMPLES:\n"
            "  # Minimal — uses ./annotations automatically:\n"
            "  python core_gene_phylogeny.py\n\n"
            "  # Explicit paths and settings:\n"
            "  python core_gene_phylogeny.py annotations/ \\\n"
            "      --output core_phylogeny/ --threads 8 --bootstrap 200\n\n"
            "  # Lower identity threshold for more diverse genomes:\n"
            "  python core_gene_phylogeny.py --identity 80\n\n"
            "OUTPUT FILES:\n"
            "  raxml/tree.raxml.support    Best ML tree with bootstrap values\n"
            "  raxml/tree.raxml.bestTree   Best ML tree (no bootstrap)\n"
            "  roary/core_gene_alignment.aln  Core gene alignment\n"
            "  roary/gene_presence_absence.csv  Pan-genome matrix\n"
        ),
    )
    parser.add_argument(
        "annotations_dir", nargs="?", default="annotations",
        help="Root folder of Bakta annotations from annotate_genomes.py (default: ./annotations)",
    )
    parser.add_argument(
        "--output", default="core_phylogeny",
        help="Output folder (default: ./core_phylogeny)",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="CPU threads for Roary and RAxML-NG (default: 4)",
    )
    parser.add_argument(
        "--identity", type=int, default=95,
        help="Minimum BLASTP %% identity for Roary gene clustering (default: 95)",
    )
    parser.add_argument(
        "--bootstrap", type=int, default=100,
        help="Number of bootstrap replicates for RAxML-NG (default: 100)",
    )
    parser.add_argument(
        "--alignment", default=None,
        help="Path to existing core gene alignment — skips Roary entirely",
    )
    parser.add_argument(
        "--model", default="GTR+G",
        help="Substitution model for RAxML-NG (default: GTR+G)",
    )
    parser.add_argument(
        "--run-modeltest", action="store_true",
        dest="run_modeltest",
        help="Run ModelTest-NG to auto-select best substitution model (slow, may not work on Apple Silicon)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-run even if output already exists",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    annotations_dir = Path(args.annotations_dir)
    output_dir = Path(args.output)
    roary_dir = output_dir / "roary"
    modeltest_dir = output_dir / "modeltest"
    raxml_dir = output_dir / "raxml"

    # Validate input
    if not annotations_dir.exists():
        parser.print_usage()
        print(f"\n[ERROR] Annotations folder not found: {annotations_dir}")
        print("        Run annotate_genomes.py first, or pass the correct path.\n")
        sys.exit(1)

    gff_files = collect_gff_files(annotations_dir)
    if not gff_files:
        print(f"[ERROR] No .gff3 files found under: {annotations_dir}")
        print("        Make sure annotate_genomes.py has completed successfully.")
        sys.exit(1)

    if len(gff_files) < 3:
        print(f"[WARN] Only {len(gff_files)} GFF3 file(s) found — Roary works best with ≥3 genomes.")

    # Check dependencies
    missing = check_dependencies()
    if missing:
        print(
            f"\n[ERROR] Missing required tools: {', '.join(missing)}\n"
            "        Install with:\n"
            "            conda install -c bioconda roary raxml-ng modeltest-ng\n"
        )
        sys.exit(1)

    # Check if already done
    if not args.force and already_done(output_dir):
        print(f"[INFO] Output already exists at {output_dir}/raxml/")
        print("       Use --force to re-run.")
        sys.exit(0)

    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"\nAnnotations  : {annotations_dir.resolve()} ({len(gff_files)} genomes)")
    print(f"Output dir   : {output_dir.resolve()}")
    print(f"Threads      : {args.threads}")
    print(f"Identity     : {args.identity}%")
    print(f"Bootstrap    : {args.bootstrap}")
    print()

    total_start = time.time()

    # -------------------------
    # Step 1: Roary (or use existing alignment)
    # -------------------------
    alignment = None

    # --alignment flag: use a pre-existing alignment, skip Roary entirely
    if hasattr(args, 'alignment') and args.alignment:
        alignment = Path(args.alignment)
        if not alignment.exists():
            print(f"[ERROR] --alignment file not found: {alignment}")
            sys.exit(1)
        print(f"[1/3] Using provided alignment: {alignment}\n")

    # Auto-detect: search output dir for alignment from a previous Roary run
    if alignment is None:
        found = list(output_dir.rglob("core_gene_alignment.aln")) + list(output_dir.rglob("core_alignment.aln"))
        found = [f for f in found if f.stat().st_size > 0]
        if found:
            alignment = sorted(found, key=lambda f: f.stat().st_mtime, reverse=True)[0]
            print(f"[1/3] Found existing Roary alignment, skipping Roary:")
            print(f"      {alignment}\n")

    # Run Roary if no alignment found
    if alignment is None:
        print("[1/3] Running Roary (pan-genome + core gene alignment)...")
        roary_start = time.time()
        roary_ok, alignment = run_roary(
            gff_files=gff_files,
            output_dir=roary_dir,
            threads=args.threads,
            identity=args.identity,
        )
        roary_elapsed = time.time() - roary_start

        if not roary_ok or alignment is None:
            print(f"[FAIL] Roary failed. Check {roary_dir}/roary_error.log")
            sys.exit(1)

        print(f"[OK  ] Core gene alignment -> {alignment}  ({roary_elapsed:.0f}s)\n")

    # -------------------------
    # Step 2: Model selection
    # -------------------------
    if args.run_modeltest:
        print("[2/3] Running ModelTest-NG (model selection)...")
        model_start = time.time()
        _, model = run_modeltest(
            alignment=alignment,
            output_dir=modeltest_dir,
            threads=args.threads,
        )
        model_elapsed = time.time() - model_start
        print(f"[OK  ] Model selected: {model}  ({model_elapsed:.0f}s)\n")
    else:
        model = args.model
        print(f"[2/3] Using substitution model: {model}")
        print("      (Use --run-modeltest to auto-select, or --model to specify)\n")

    # -------------------------
    # Step 3: RAxML-NG
    # -------------------------
    print("[3/3] Running RAxML-NG (ML tree + bootstrapping)...")
    print("  Validating alignment...")
    alignment = fix_alignment(alignment, raxml_dir)
    raxml_start = time.time()
    raxml_ok, msg = run_raxml(
        alignment=alignment,
        model=model,
        output_dir=raxml_dir,
        threads=args.threads,
        bootstrap=args.bootstrap,
    )
    raxml_elapsed = time.time() - raxml_start

    total_elapsed = time.time() - total_start

    print()
    print("=" * 50)
    if raxml_ok:
        print(f"[OK  ] {msg}")
        print(f"       Total time: {total_elapsed:.0f}s")
        print(f"\nView tree in FigTree or upload to https://itol.embl.de")
        print(f"Pan-genome matrix: {roary_dir}/gene_presence_absence.csv")
    else:
        print(f"[FAIL] {msg}")
        sys.exit(1)


if __name__ == "__main__":
    main()
