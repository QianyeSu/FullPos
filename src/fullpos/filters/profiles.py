from __future__ import annotations

import numpy as np


def gaussian_filter_profile(
    ntrunc: int,
    *,
    cutoff: int | None = None,
    exponent: float = 1.0,
) -> np.ndarray:
    """Return the FULLPOS ``FPFILTER`` Gaussian profile by total wave number."""
    ntrunc = _validate_ntrunc(ntrunc)
    cutoff = ntrunc if cutoff is None else _validate_cutoff(cutoff, ntrunc)
    scale = min(ntrunc, cutoff)
    if scale <= 0:
        return np.ones(ntrunc + 1, dtype=np.float64)
    n = np.arange(ntrunc + 1, dtype=np.float64)
    return np.exp(-0.5 * float(exponent) * n**2 / float(scale * scale))


def low_pass_filter_profile(
    ntrunc: int,
    *,
    cutoff: int,
    selectivity: float = 1.0,
    exponent: float = 1.0,
) -> np.ndarray:
    """Return the FULLPOS ``FPFILTER`` low-pass/THX profile by total wave number."""
    ntrunc = _validate_ntrunc(ntrunc)
    cutoff = _validate_cutoff(cutoff, ntrunc)
    n = np.arange(ntrunc + 1, dtype=np.float64)
    z = np.exp(-float(exponent) * float(selectivity))
    return 0.5 * (1.0 - np.tanh(z * (n - float(cutoff))))


def fpfilter_profile(
    ntrunc: int,
    *,
    kind: str = "gaussian",
    cutoff: int | None = None,
    selectivity: float = 1.0,
    low_pass_exponent: float = 1.0,
    gaussian_exponent: float = 1.0,
) -> np.ndarray:
    """Build a 1D FULLPOS-style profile indexed by total wave number."""
    normalized = _normalize_kind(kind)
    if normalized == "gaussian":
        return gaussian_filter_profile(ntrunc, cutoff=cutoff, exponent=gaussian_exponent)
    if cutoff is None:
        raise ValueError("cutoff is required for low-pass spectral filters")
    return low_pass_filter_profile(
        ntrunc,
        cutoff=cutoff,
        selectivity=selectivity,
        exponent=low_pass_exponent,
    )


def expand_profile_to_coefficients(profile: np.ndarray) -> np.ndarray:
    """Expand total-wave-number weights to ECTRANS real/imag coefficient weights."""
    weights = np.asarray(profile, dtype=np.float64)
    if weights.ndim != 1 or weights.size == 0:
        raise ValueError("profile must be a non-empty 1D array")
    ntrunc = weights.size - 1
    coeff_weights = np.empty((ntrunc + 1) * (ntrunc + 2), dtype=np.float64)
    pos = 0
    for m in range(ntrunc + 1):
        for n in range(m, ntrunc + 1):
            coeff_weights[pos] = weights[n]
            coeff_weights[pos + 1] = weights[n]
            pos += 2
    return coeff_weights


def infer_ntrunc_from_nspec2(nspec2: int) -> int:
    """Infer triangular truncation from an ECTRANS real-valued coefficient count."""
    nspec2 = int(nspec2)
    if nspec2 <= 0:
        raise ValueError("coefficient count must be positive")
    ntrunc = int((np.sqrt(1 + 4 * nspec2) - 3) / 2)
    if (ntrunc + 1) * (ntrunc + 2) != nspec2:
        raise ValueError("coefficient count is not a triangular spectral layout")
    return ntrunc


def _normalize_kind(kind: str) -> str:
    normalized = str(kind).lower().replace("-", "_")
    aliases = {
        "classic": "gaussian",
        "gauss": "gaussian",
        "lowpass": "low_pass",
        "thx": "low_pass",
    }
    normalized = aliases.get(normalized, normalized)
    if normalized not in {"gaussian", "low_pass"}:
        raise ValueError("filter kind must be 'gaussian' or 'low_pass'")
    return normalized


def _validate_ntrunc(ntrunc: int) -> int:
    value = int(ntrunc)
    if value < 0:
        raise ValueError("ntrunc must be non-negative")
    return value


def _validate_cutoff(cutoff: int, ntrunc: int) -> int:
    value = int(cutoff)
    if value < 0 or value > ntrunc:
        raise ValueError(f"cutoff must be between 0 and {ntrunc}, got {cutoff}")
    return value
