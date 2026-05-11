from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .core import SpectralFilter
from .profiles import expand_profile_to_coefficients


def save_spectral_filter(path, spectral_filter: SpectralFilter) -> None:
    """Save a reusable diagonal spectral filter to a compressed NumPy ``.npz`` file."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    metadata = json.dumps(spectral_filter.to_metadata(), sort_keys=True)
    np.savez_compressed(
        target,
        profile=np.asarray(spectral_filter.profile, dtype=np.float64),
        coefficient_weights=spectral_filter.coefficient_weights,
        metadata=np.array(metadata),
    )


def load_spectral_filter(path) -> SpectralFilter:
    """Load a diagonal spectral filter saved by :func:`save_spectral_filter`."""
    with np.load(Path(path), allow_pickle=False) as data:
        profile = np.asarray(data["profile"], dtype=np.float64)
        metadata = json.loads(str(data["metadata"].item()))
    return SpectralFilter(
        profile,
        kind=str(metadata.get("kind", "custom")),
        cutoff=metadata.get("cutoff"),
        selectivity=metadata.get("selectivity"),
        low_pass_exponent=metadata.get("low_pass_exponent"),
        gaussian_exponent=metadata.get("gaussian_exponent"),
    )


def save_filter_matrix(path, profile) -> None:
    """Save diagonal coefficient weights implied by a wave-number profile.

    This is a lightweight diagonal matrix format, not the FULLPOS stretched
    geometry LFI matrix format produced by ``CPFPFILTER``.
    """
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    profile_arr = np.asarray(profile, dtype=np.float64)
    np.savez_compressed(
        target,
        profile=profile_arr,
        coefficient_weights=expand_profile_to_coefficients(profile_arr),
        metadata=np.array(
            json.dumps(
                {
                    "version": 1,
                    "kind": "diagonal_profile_matrix",
                    "ntrunc": int(profile_arr.size - 1),
                    "layout": "ectrans_global_real_imag_by_total_wave_number",
                },
                sort_keys=True,
            )
        ),
    )


def load_filter_matrix(path) -> dict:
    """Load a diagonal profile matrix saved by :func:`save_filter_matrix`."""
    with np.load(Path(path), allow_pickle=False) as data:
        return {
            "profile": np.asarray(data["profile"], dtype=np.float64),
            "coefficient_weights": np.asarray(data["coefficient_weights"], dtype=np.float64),
            "metadata": json.loads(str(data["metadata"].item())),
        }
