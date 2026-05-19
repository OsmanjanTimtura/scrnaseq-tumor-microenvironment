"""Data loading and persistence helpers.

The default :data:`DATA_DIR` is the repo's ``data/`` folder, but both functions
accept absolute paths so the working subsample can live outside the repo
(typical on Windows, where users keep multi-GB scientific data on a separate
drive). See ``data/README.md`` for the recommended layout.
"""
from __future__ import annotations

from pathlib import Path
from typing import Optional

import anndata as ad


# Default repo-relative data dir. Subsampling output and processed h5ads can
# live here; the raw full atlas (~13 GB) usually lives elsewhere.
DATA_DIR = Path(__file__).resolve().parent.parent / "data"

# Canonical filenames used across notebooks 01–04.
DEFAULT_SUBSAMPLE_NAME = "luca_subsample_100k.h5ad"
DEFAULT_FILTERED_NAME = "processed/luca_filtered.h5ad"


def _resolve(filename: str | Path, base: Path) -> Path:
    """If ``filename`` is absolute, return it as-is; otherwise join with ``base``."""
    p = Path(filename)
    return p if p.is_absolute() else base / p


def load_adata(
    filename: str | Path,
    data_dir: Optional[Path] = None,
) -> ad.AnnData:
    """Load an h5ad file.

    Parameters
    ----------
    filename
        Either a bare filename (e.g. ``"luca_subsample_100k.h5ad"``) that is
        joined with ``data_dir``, or an absolute path that is used as-is.
    data_dir
        Override the default repo ``data/`` directory. Ignored if ``filename``
        is absolute.

    Raises
    ------
    FileNotFoundError
        Points the user at ``data/README.md`` for download/subsample steps.
    """
    base = data_dir or DATA_DIR
    path = _resolve(filename, base)
    if not path.exists():
        raise FileNotFoundError(
            f"Expected file not found: {path}\n"
            f"See data/README.md for download + subsample instructions."
        )
    return ad.read_h5ad(path)


def save_adata(
    adata: ad.AnnData,
    filename: str | Path,
    data_dir: Optional[Path] = None,
    compression: str | None = "gzip",
) -> Path:
    """Save an AnnData object and return the output path.

    Creates parent directories as needed. Defaults to gzip compression so the
    repo's ``data/processed/`` artifacts stay small enough to inspect quickly.
    """
    base = data_dir or DATA_DIR
    path = _resolve(filename, base)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, compression=compression)
    return path
