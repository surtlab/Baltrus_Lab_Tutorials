#!/usr/bin/env python3
"""
setup_env.py
------------
Creates the phylo_pipeline conda environment with the correct package
versions for your platform.

Usage:
    python setup_env.py           # create the environment
    python setup_env.py --update  # update an existing environment
    python setup_env.py --remove  # remove the environment

Platforms handled:
    osx-arm64   Apple Silicon Mac (M1/M2/M3/M4)
                bakta 1.11.x + pyhmmer 0.10.x
                (blast >=2.17 not yet on osx-arm64; 1.11.x uses blast 2.16)

    osx-64      Intel Mac
                bakta >=1.9 (latest available)

    linux-64    Linux x86-64
                bakta >=1.9 (latest available)

    win-*       Windows
                bakta is not supported on Windows natively.
                WSL2 (Windows Subsystem for Linux) is required.
"""

import argparse
import platform
import shutil
import subprocess
import sys
from pathlib import Path

ENV_NAME = "phylo_pipeline"
BASE_CHANNELS = ["-c", "bioconda", "-c", "conda-forge"]
BASE_DEPS = [
    "python>=3.10",
    "biopython>=1.81",
]

# Platform-specific bakta/pyhmmer pins
# osx-arm64: blast >=2.17 not available, so cap bakta at 1.11.x
# and pin pyhmmer to 0.10.x (what bakta 1.11 expects)
PLATFORM_DEPS = {
    "osx-arm64": [
        "bakta>=1.11,<1.12",
        "pyhmmer=0.10.*",
    ],
    "osx-64": [
        "bakta>=1.9",
    ],
    "linux-64": [
        "bakta>=1.9",
    ],
    "linux-aarch64": [
        "bakta>=1.9",
    ],
}


def detect_platform() -> str:
    """Detect the conda platform string for this machine."""
    system = platform.system().lower()
    machine = platform.machine().lower()

    if system == "windows":
        return "win-64" if "64" in machine else "win-32"
    elif system == "darwin":
        # Apple Silicon
        if machine in ("arm64", "aarch64"):
            return "osx-arm64"
        return "osx-64"
    elif system == "linux":
        if machine in ("arm64", "aarch64"):
            return "linux-aarch64"
        return "linux-64"
    return "unknown"


def find_conda() -> str | None:
    """Return the path to mamba, micromamba, or conda (in that order)."""
    for cmd in ["mamba", "micromamba", "conda"]:
        if shutil.which(cmd):
            return cmd
    return None


def env_exists(conda: str) -> bool:
    result = subprocess.run(
        [conda, "env", "list"], capture_output=True, text=True
    )
    return ENV_NAME in result.stdout


def create_env(conda: str, extra_deps: list[str]):
    cmd = (
        [conda, "create", "-n", ENV_NAME, "--yes"]
        + BASE_CHANNELS
        + [d for dep in BASE_DEPS + extra_deps for d in [dep]]
    )
    # conda create takes packages as positional args
    cmd = [conda, "create", "-n", ENV_NAME, "--yes"] + BASE_CHANNELS + BASE_DEPS + extra_deps
    print(f"\nRunning: {' '.join(cmd)}\n")
    result = subprocess.run(cmd)
    return result.returncode == 0


def update_env(conda: str, extra_deps: list[str]):
    cmd = [conda, "install", "-n", ENV_NAME, "--yes"] + BASE_CHANNELS + BASE_DEPS + extra_deps
    print(f"\nRunning: {' '.join(cmd)}\n")
    result = subprocess.run(cmd)
    return result.returncode == 0


def remove_env(conda: str):
    cmd = [conda, "env", "remove", "-n", ENV_NAME, "--yes"]
    print(f"\nRunning: {' '.join(cmd)}\n")
    result = subprocess.run(cmd)
    return result.returncode == 0


def main():
    parser = argparse.ArgumentParser(
        description=f"Set up the {ENV_NAME} conda environment for your platform.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--update", action="store_true", help="Update packages in existing environment")
    parser.add_argument("--remove", action="store_true", help="Remove the environment")
    args = parser.parse_args()

    # Detect platform
    plat = detect_platform()
    print(f"Detected platform : {plat}")

    # Windows check
    if plat.startswith("win"):
        print(
            "\n[ERROR] Bakta is not supported on Windows natively.\n"
            "        Please use WSL2 (Windows Subsystem for Linux) and re-run\n"
            "        this script from inside the WSL2 terminal.\n"
            "        Setup guide: https://docs.microsoft.com/en-us/windows/wsl/install\n"
        )
        sys.exit(1)

    if plat == "unknown":
        print("[WARN] Could not detect platform, defaulting to linux-64 package set.")
        plat = "linux-64"

    # Get platform-specific deps
    extra_deps = PLATFORM_DEPS.get(plat, PLATFORM_DEPS["linux-64"])
    if plat == "osx-arm64":
        print(
            "  Apple Silicon detected — using bakta 1.11.x + pyhmmer 0.10.x\n"
            "  (blast >=2.17 is not yet available for osx-arm64)\n"
        )
    else:
        print(f"  Using latest available bakta for {plat}\n")

    # Find conda/mamba
    conda = find_conda()
    if not conda:
        print(
            "[ERROR] Could not find conda, mamba, or micromamba.\n"
            "        Please install Miniconda or Mambaforge first:\n"
            "        https://github.com/conda-forge/miniforge\n"
        )
        sys.exit(1)
    print(f"Using installer    : {conda}")

    # Remove
    if args.remove:
        if not env_exists(conda):
            print(f"\n[INFO] Environment '{ENV_NAME}' does not exist, nothing to remove.")
            sys.exit(0)
        print(f"\nRemoving environment '{ENV_NAME}'...")
        if remove_env(conda):
            print(f"Environment '{ENV_NAME}' removed.")
        else:
            print("[ERROR] Removal failed.")
            sys.exit(1)
        sys.exit(0)

    # Update
    if args.update:
        if not env_exists(conda):
            print(f"\n[INFO] Environment '{ENV_NAME}' not found — creating instead of updating.")
        else:
            print(f"\nUpdating environment '{ENV_NAME}'...")
            if update_env(conda, extra_deps):
                print(f"\nEnvironment '{ENV_NAME}' updated successfully.")
                print(f"Activate with: {conda} activate {ENV_NAME}")
            else:
                print("\n[ERROR] Update failed.")
                sys.exit(1)
            sys.exit(0)

    # Create
    if env_exists(conda):
        print(
            f"\n[INFO] Environment '{ENV_NAME}' already exists.\n"
            f"       Use --update to update packages, or --remove to start fresh.\n"
        )
        sys.exit(0)

    print(f"\nCreating environment '{ENV_NAME}'...")
    if create_env(conda, extra_deps):
        print(f"\nEnvironment '{ENV_NAME}' created successfully!")
        print(f"\nActivate with:\n    {conda} activate {ENV_NAME}\n")
    else:
        print("\n[ERROR] Environment creation failed.")
        sys.exit(1)


if __name__ == "__main__":
    main()
