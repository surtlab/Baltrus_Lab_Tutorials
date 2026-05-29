#!/usr/bin/env python3
"""
annotate_genomes.py
-------------------
Step 2 of the genome pipeline. Runs Bakta to annotate every .fasta file
in the genomes/ folder produced by download_genomes.py.

Outputs per genome (in annotations/<strain_accession>/):
    <prefix>.gbff   GenBank flat file
    <prefix>.gff3   GFF3 annotation
    <prefix>.faa    Protein FASTA (CDS amino acid sequences)
    + all other default Bakta outputs (TSV, JSON, etc.)

Usage:
    python annotate_genomes.py [genomes_dir] [--db /path/to/bakta/db]
                               [--output annotations] [--threads 4]
                               [--conda-env bakta] [--force]

Arguments:
    genomes_dir     Folder of .fasta files from download_genomes.py.
                    Default: ./genomes
    --db            Path to Bakta database directory. Can also be set via
                    the BAKTA_DB environment variable. If neither is provided
                    and no database is auto-detected, the light database will
                    be downloaded automatically to ~/.bakta_db/db-light.
    --output        Root folder for annotation results. Default: ./annotations
    --threads       CPU threads per Bakta job. Default: 4
    --conda-env     Name of the conda/mamba environment that has Bakta.
                    Default: bakta
    --force         Re-annotate genomes that already have an output folder.

Output layout:
    annotations/
        <strain>_<accession>/
            <strain>_<accession>.gbff
            <strain>_<accession>.gff3
            <strain>_<accession>.faa
            ... (other Bakta outputs)
        annotation_log.txt
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path


# ---------------------------------------------------------------------------
# Bakta location
# ---------------------------------------------------------------------------

def find_bakta(conda_env: str) -> str:
    """
    Return the command prefix needed to invoke bakta.
    Priority:
      1. bakta already on PATH (env already activated)
      2. conda run -n <env> bakta
    Returns a prefix string (may be empty) or None if not found.
    """
    if shutil.which("bakta"):
        return ""  # already on PATH

    for conda_bin in ["conda", "mamba", "micromamba"]:
        if shutil.which(conda_bin):
            return f"{conda_bin} run -n {conda_env}"

    return None


def check_bakta(cmd_prefix: str) -> bool:
    """Verify bakta is callable and print its version."""
    if cmd_prefix is None:
        return False
    cmd = f"{cmd_prefix} bakta --version".strip()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode == 0:
        print(f"  Bakta found: {result.stdout.strip()}")
        return True
    return False


def find_bakta_python(cmd_prefix: str) -> str | None:
    """
    Return the path to the Python interpreter in the same env as bakta.
    Used to locate the bakta package directory for db auto-detection.
    """
    if cmd_prefix:
        # Using conda run — ask that env for its python
        cmd = f"{cmd_prefix} python3 -c \"import sys; print(sys.executable)\""
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
        if r.returncode == 0:
            return r.stdout.strip()
    else:
        # bakta is on PATH — its python is in the same bin directory
        bakta_bin = shutil.which("bakta")
        if bakta_bin:
            env_bin = Path(bakta_bin).parent
            for py in ["python3", "python"]:
                p = env_bin / py
                if p.exists():
                    return str(p)
        # Fall back to current interpreter
        return sys.executable
    return None


# ---------------------------------------------------------------------------
# Database resolution
# ---------------------------------------------------------------------------

DEFAULT_DB_DIR = Path.home() / ".bakta_db"


def db_is_valid(db_path: Path) -> bool:
    """Check that a path looks like a valid bakta database (has version.json)."""
    return db_path.is_dir() and (db_path / "version.json").exists()


def find_db_in_bakta_package(cmd_prefix: str) -> Path | None:
    """
    Look for a db/ or db-light/ folder inside bakta's own package directory.
    This is where the db lives if it was placed alongside the conda install.
    """
    py = find_bakta_python(cmd_prefix)
    if not py:
        return None

    r = subprocess.run(
        [py, "-c", "import bakta, os; print(os.path.dirname(bakta.__file__))"],
        capture_output=True, text=True
    )
    if r.returncode != 0:
        return None

    pkg_dir = Path(r.stdout.strip())
    for candidate in ["db", "db-light"]:
        p = pkg_dir / candidate
        if db_is_valid(p):
            return p
    return None


def update_amrfinder(db_path: Path, cmd_prefix: str) -> bool:
    """
    Run amrfinder_update to initialise AMRFinderPlus's internal database.
    This is required once after every fresh bakta_db download.
    Returns True on success.
    """
    amr_db = db_path / "amrfinderplus-db"
    print(f"  Initialising AMRFinderPlus database at {amr_db} ...")
    print("  (This downloads ~500 MB and is required once after install)\n")
    cmd = f"{cmd_prefix} amrfinder_update --force_update --database {amr_db}".strip()
    result = subprocess.run(cmd, shell=True, text=True)
    if result.returncode != 0:
        print(
            "  [WARN] AMRFinderPlus update failed — AMR gene detection will be skipped by bakta.\n"
            "         To fix manually, run:\n"
            f"             amrfinder_update --force_update --database {amr_db}\n"
        )
        return False
    print("  AMRFinderPlus database ready.\n")
    return True


def download_light_db(target_dir: Path, cmd_prefix: str) -> Path | None:
    """
    Download the Bakta light database to target_dir using bakta_db,
    then initialise the AMRFinderPlus sub-database.
    Returns the path to the downloaded db, or None on failure.
    """
    target_dir.mkdir(parents=True, exist_ok=True)
    db_dest = target_dir / "db-light"

    if db_is_valid(db_dest):
        print(f"  Light database already exists at {db_dest}, skipping download.")
        # Still check AMRFinder is initialised (may have failed on a previous run)
        amr_db = db_dest / "amrfinderplus-db"
        if not any(amr_db.glob("*/AMR.fa")) if amr_db.exists() else True:
            update_amrfinder(db_dest, cmd_prefix)
        return db_dest

    print(f"  Downloading Bakta light database to {target_dir} ...")
    print("  (This is ~1.3 GB and may take a few minutes)\n")

    cmd = f"{cmd_prefix} bakta_db download --output {target_dir} --type light".strip()
    result = subprocess.run(cmd, shell=True, text=True)

    # bakta_db exits non-zero when only the AMRFinder step fails —
    # check whether the core db was actually written before giving up
    resolved = None
    if db_is_valid(db_dest):
        resolved = db_dest
    else:
        for p in sorted(target_dir.iterdir()):
            if p.is_dir() and db_is_valid(p):
                resolved = p
                break

    if resolved is None:
        print(f"\n[ERROR] Database download failed (exit code {result.returncode}).")
        return None

    # Run amrfinder_update separately to finish what bakta_db may have left incomplete
    update_amrfinder(resolved, cmd_prefix)
    return resolved


def resolve_db(args_db: str | None, cmd_prefix: str) -> str | None:
    """
    Determine the database path to use, in priority order:

      1. --db argument
      2. BAKTA_DB environment variable
      3. db/ or db-light/ next to the bakta package (conda install default)
      4. Already-downloaded light db at ~/.bakta_db/db-light
      5. Auto-download the light db to ~/.bakta_db/db-light

    Returns the resolved path string, or None to let bakta try its own logic
    (shouldn't happen after step 5, but kept as a safety valve).
    """
    # 1. Explicit --db
    if args_db:
        p = Path(args_db).expanduser().resolve()
        if db_is_valid(p):
            print(f"  Using database from --db: {p}")
            return str(p)
        else:
            print(f"[ERROR] --db path does not look like a valid Bakta database: {p}")
            print("        Expected a directory containing version.json")
            sys.exit(1)

    # 2. BAKTA_DB environment variable
    env_db = os.environ.get("BAKTA_DB")
    if env_db:
        p = Path(env_db).expanduser().resolve()
        if db_is_valid(p):
            print(f"  Using database from BAKTA_DB: {p}")
            return str(p)
        else:
            print(f"  [WARN] BAKTA_DB is set but path is not a valid database: {p}")

    # 3. Check bakta's package directory
    print("  Checking for database next to bakta package install...")
    pkg_db = find_db_in_bakta_package(cmd_prefix)
    if pkg_db:
        print(f"  Found database in bakta package dir: {pkg_db}")
        return str(pkg_db)

    # 4. Previously downloaded light db
    cached = DEFAULT_DB_DIR / "db-light"
    if db_is_valid(cached):
        print(f"  Found previously downloaded light database: {cached}")
        return str(cached)

    # 5. Auto-download light db
    print("\n  No Bakta database found. Automatically downloading the light database.")
    print(f"  Install location: {DEFAULT_DB_DIR}")
    print("  (To use a full database instead, cancel and run:")
    print("       bakta_db download --output ~/bakta_db --type full")
    print("   then pass --db ~/bakta_db/db)\n")

    result = download_light_db(DEFAULT_DB_DIR, cmd_prefix)
    if result:
        print(f"\n  Database ready: {result}\n")
        return str(result)

    print(
        "\n[ERROR] Could not auto-download the database.\n"
        "        Please download it manually and pass the path with --db:\n"
        "            bakta_db download --output ~/bakta_db --type light\n"
        "            python annotate_genomes.py --db ~/bakta_db/db-light\n"
    )
    sys.exit(1)


# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

def already_annotated(output_subdir: Path) -> bool:
    """
    Return True if this genome has a complete annotation:
    requires a non-empty .gbff, .gff3, and .faa to all be present.
    """
    if not output_subdir.exists():
        return False
    required_suffixes = {".gbff", ".gff3", ".faa"}
    found = {
        f.suffix for f in output_subdir.iterdir()
        if f.suffix in required_suffixes and f.stat().st_size > 0
    }
    return required_suffixes == found


def run_bakta(
    fasta: Path,
    prefix: str,
    output_subdir: Path,
    db_path: str,
    threads: int,
    cmd_prefix: str,
    skip_amr: bool = True,
) -> tuple[bool, str]:
    """
    Run Bakta on a single genome FASTA.
    Returns (success, message).
    """
    output_subdir.mkdir(parents=True, exist_ok=True)

    cmd_parts = []
    if cmd_prefix:
        cmd_parts.extend(cmd_prefix.split())

    cmd_parts += [
        "bakta",
        "--output", str(output_subdir),
        "--prefix", prefix,
        "--threads", str(threads),
        "--force",
        "--db", db_path,
    ]
    cmd_parts.append(str(fasta))

    try:
        result = subprocess.run(cmd_parts, capture_output=True, text=True)

        if result.returncode != 0:
            err_log = output_subdir / "bakta_error.log"
            with open(err_log, "w") as ef:
                ef.write("=== STDOUT ===\n")
                ef.write(result.stdout or "(empty)\n")
                ef.write("\n=== STDERR ===\n")
                ef.write(result.stderr or "(empty)\n")
            err_lines = (result.stderr or result.stdout or "").strip().splitlines()
            err_tail = "\n        ".join(err_lines[-15:]) if err_lines else "(no output)"
            return False, (
                f"Bakta exited with code {result.returncode}:\n"
                f"        {err_tail}\n"
                f"        Full error log: {err_log}"
            )

        if not already_annotated(output_subdir):
            return False, "Bakta finished but expected output files (.gbff/.gff3/.faa) not found."

        return True, f"Annotation complete -> {output_subdir.name}/"

    except FileNotFoundError:
        return False, "bakta executable not found. Check --conda-env or activate your environment."
    except Exception as exc:
        return False, f"Unexpected error: {exc}"


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Annotate genomes with Bakta (pipeline step 2).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "EXAMPLES:\n"
            "  # Minimal — bakta env active, db auto-detected or downloaded:\n"
            "  python annotate_genomes.py\n\n"
            "  # Explicit db path:\n"
            "  python annotate_genomes.py genomes/ --db ~/bakta_db/db-light\n\n"
            "  # Custom conda env name and thread count:\n"
            "  python annotate_genomes.py --conda-env my_bakta_env --threads 8\n\n"
            "DATABASE NOTES:\n"
            "  If no database is found via --db or BAKTA_DB, the script will\n"
            "  automatically download the light database (~1.3 GB) to:\n"
            f"      {DEFAULT_DB_DIR / 'db-light'}\n\n"
            "  To use the full database instead (~32 GB, better for novel taxa):\n"
            "      bakta_db download --output ~/bakta_db --type full\n"
            "      python annotate_genomes.py --db ~/bakta_db/db\n"
        ),
    )
    parser.add_argument(
        "genomes_dir", nargs="?", default="genomes",
        help="Folder of .fasta files from download_genomes.py (default: ./genomes)",
    )
    parser.add_argument(
        "--db", default=None,
        help="Path to Bakta database. Auto-detected or downloaded if not provided.",
    )
    parser.add_argument(
        "--output", default="annotations",
        help="Root output folder for annotations (default: ./annotations)",
    )
    parser.add_argument(
        "--threads", type=int, default=4,
        help="CPU threads per Bakta run (default: 4)",
    )
    parser.add_argument(
        "--conda-env", default="bakta",
        help="Conda/mamba environment name containing Bakta (default: bakta)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="Re-annotate genomes that already have complete output folders",
    )
    parser.add_argument(
        "--skip-amr", action="store_true", default=True,
        help="Skip AMR gene detection (default: on; use --no-skip-amr to enable)",
    )
    parser.add_argument(
        "--no-skip-amr", dest="skip_amr", action="store_false",
        help="Enable AMR gene detection (requires AMRFinderPlus database to be set up)",
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()

    genomes_dir = Path(args.genomes_dir)
    if not genomes_dir.exists():
        parser.print_usage()
        print(f"\n[ERROR] Genomes folder not found: {genomes_dir}")
        print("        Run download_genomes.py first, or pass the correct path.\n")
        sys.exit(1)

    fastas = sorted(genomes_dir.glob("*.fasta"))
    if not fastas:
        print(f"[ERROR] No .fasta files found in: {genomes_dir}")
        sys.exit(1)

    output_root = Path(args.output)
    output_root.mkdir(parents=True, exist_ok=True)

    # Locate bakta binary
    print("\nLocating Bakta...")
    cmd_prefix = find_bakta(args.conda_env)
    if not check_bakta(cmd_prefix):
        print(
            "\n[ERROR] Could not find Bakta.\n"
            f"        Expected conda environment: '{args.conda_env}'\n"
            "        Activate it with:  conda activate <env>\n"
            "        Or pass a different name with --conda-env <name>\n"
        )
        sys.exit(1)

    # Resolve database — auto-downloads light db if nothing found
    print("\nResolving Bakta database...")
    db_path = resolve_db(args.db, cmd_prefix)

    print(f"\nGenomes dir  : {genomes_dir.resolve()}")
    print(f"Output dir   : {output_root.resolve()}")
    print(f"Database     : {db_path}")
    print(f"Threads      : {args.threads}")
    print()

    # Tally already-done vs to-do
    skipped_names = []
    todo = []
    for fasta in fastas:
        prefix = fasta.stem
        subdir = output_root / prefix
        if not args.force and already_annotated(subdir):
            skipped_names.append(prefix)
        else:
            todo.append(fasta)

    print(f"Found {len(fastas)} genome(s).")
    if skipped_names:
        print(f"  {len(skipped_names)} already annotated, will skip.")
    print(f"  {len(todo)} to annotate.\n")

    if not todo:
        print("Nothing to do. Use --force to re-annotate existing results.")
        sys.exit(0)

    log_lines = []
    success_count = 0
    fail_count = 0

    for i, fasta in enumerate(todo, 1):
        prefix = fasta.stem
        subdir = output_root / prefix
        print(f"[{i}/{len(todo)}] {prefix}")

        start = time.time()
        ok, msg = run_bakta(
            fasta=fasta,
            prefix=prefix,
            output_subdir=subdir,
            db_path=db_path,
            threads=args.threads,
            cmd_prefix=cmd_prefix,
            skip_amr=args.skip_amr,
        )
        elapsed = time.time() - start

        status = "OK  " if ok else "FAIL"
        detail = f"{msg}  ({elapsed:.0f}s)"
        print(f"    [{status}] {detail}")
        log_lines.append(f"[{status}] {prefix}\t{detail}")

        if ok:
            success_count += 1
        else:
            fail_count += 1

    log_path = output_root / "annotation_log.txt"
    with open(log_path, "w") as lf:
        lf.write("\n".join(log_lines) + "\n")

    print()
    print("=" * 50)
    print(f"Done. {success_count} succeeded, {fail_count} failed.")
    if skipped_names:
        print(f"       {len(skipped_names)} skipped (already annotated).")
    print(f"Log written to: {log_path}")


if __name__ == "__main__":
    main()
