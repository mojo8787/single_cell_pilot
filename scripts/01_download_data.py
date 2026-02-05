#!/usr/bin/env python3
"""
Download GSE231935 (M3-Seq bacterial single-cell RNA-seq) from NCBI GEO.
Processed data: CSV files in TAR archive (~42 Mb).
"""

import os
import tarfile
from pathlib import Path

import requests

# Configuration
GEO_ACCESSION = "GSE231935"
GEO_URL = f"https://www.ncbi.nlm.nih.gov/geo/download/?acc={GEO_ACCESSION}&format=file"
FTP_FALLBACK = f"ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE231nnn/{GEO_ACCESSION}/suppl/{GEO_ACCESSION}_RAW.tar"

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
TAR_PATH = DATA_DIR / f"{GEO_ACCESSION}_RAW.tar"


def download_geo_data() -> Path:
    """Download GSE231935 processed data from GEO. Returns path to TAR file."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    if TAR_PATH.exists():
        print(f"Data already exists at {TAR_PATH}. Skipping download.")
        return TAR_PATH

    print(f"Downloading {GEO_ACCESSION} from GEO...")
    try:
        response = requests.get(GEO_URL, allow_redirects=True, timeout=120)
        response.raise_for_status()

        # GEO may return HTML with a redirect link; check content type
        content_type = response.headers.get("Content-Type", "")
        if "text/html" in content_type:
            # Try FTP fallback
            print("HTTP returned HTML, trying FTP fallback...")
            response = requests.get(FTP_FALLBACK, allow_redirects=True, timeout=120)
            response.raise_for_status()

        with open(TAR_PATH, "wb") as f:
            f.write(response.content)

        print(f"Downloaded {TAR_PATH.stat().st_size / 1e6:.1f} Mb to {TAR_PATH}")
        return TAR_PATH

    except requests.RequestException as e:
        raise RuntimeError(f"Download failed: {e}") from e


def extract_tar(tar_path: Path) -> Path:
    """Extract TAR archive to data directory."""
    extract_dir = DATA_DIR / GEO_ACCESSION
    if extract_dir.exists() and any(extract_dir.iterdir()):
        print(f"Data already extracted to {extract_dir}. Skipping extraction.")
        return extract_dir

    print(f"Extracting {tar_path}...")
    with tarfile.open(tar_path, "r:*") as tf:
        tf.extractall(DATA_DIR)

    # GEO TAR often extracts to a subfolder; list what we got
    extracted = list(DATA_DIR.iterdir())
    for item in extracted:
        if item.is_dir() and item.name != GEO_ACCESSION:
            # Might be GSE231935_RAW or similar
            pass
        print(f"  {item.name}")

    # Find the directory with extracted files
    for item in DATA_DIR.iterdir():
        if item.is_dir() and item.name.startswith(GEO_ACCESSION):
            return item
        if item.is_dir():
            return item

    return DATA_DIR


def main():
    tar_path = download_geo_data()
    extract_dir = extract_tar(tar_path)
    print(f"Done. Data in {extract_dir}")


if __name__ == "__main__":
    main()
