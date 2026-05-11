from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..interpolation import DEFAULT_CHUNK_SIZE
from ..spectral import spectral_fit, spectral_synthesis
from .profiles import expand_profile_to_coefficients, fpfilter_profile, infer_ntrunc_from_nspec2


@dataclass(frozen=True)
class SpectralFilter:
    """Reusable diagonal spectral filter in ECTRANS coefficient space."""

    profile: np.ndarray
    kind: str = "custom"
    cutoff: int | None = None
    selectivity: float | None = None
    low_pass_exponent: float | None = None
    gaussian_exponent: float | None = None

    def __post_init__(self) -> None:
        """Validate and freeze the filter profile as float64."""
        profile = np.asarray(self.profile, dtype=np.float64)
        if profile.ndim != 1 or profile.size == 0:
            raise ValueError("profile must be a non-empty 1D array")
        object.__setattr__(self, "profile", profile.copy())

    @classmethod
    def gaussian(
        cls,
        ntrunc: int,
        *,
        cutoff: int | None = None,
        exponent: float = 1.0,
    ) -> "SpectralFilter":
        """Create a reusable FULLPOS-style Gaussian spectral filter."""
        return cls(
            fpfilter_profile(
                ntrunc,
                kind="gaussian",
                cutoff=cutoff,
                gaussian_exponent=exponent,
            ),
            kind="gaussian",
            cutoff=cutoff,
            gaussian_exponent=float(exponent),
        )

    @classmethod
    def low_pass(
        cls,
        ntrunc: int,
        *,
        cutoff: int,
        selectivity: float = 1.0,
        exponent: float = 1.0,
    ) -> "SpectralFilter":
        """Create a reusable FULLPOS-style low-pass/THX spectral filter."""
        return cls(
            fpfilter_profile(
                ntrunc,
                kind="low_pass",
                cutoff=cutoff,
                selectivity=selectivity,
                low_pass_exponent=exponent,
            ),
            kind="low_pass",
            cutoff=int(cutoff),
            selectivity=float(selectivity),
            low_pass_exponent=float(exponent),
        )

    @classmethod
    def from_profile(cls, profile, *, kind: str = "custom") -> "SpectralFilter":
        """Create a filter from explicit total-wave-number weights."""
        return cls(np.asarray(profile, dtype=np.float64), kind=str(kind))

    @property
    def ntrunc(self) -> int:
        """Triangular truncation implied by the profile length."""
        return int(self.profile.size - 1)

    @property
    def nspec2(self) -> int:
        """Number of real-valued ECTRANS coefficients for this truncation."""
        return (self.ntrunc + 1) * (self.ntrunc + 2)

    @property
    def coefficient_weights(self) -> np.ndarray:
        """Profile expanded to the ECTRANS real/imag coefficient layout."""
        return expand_profile_to_coefficients(self.profile)

    def apply_coefficients(self, coefficients, *, axis: int = -1) -> np.ndarray:
        """Apply this diagonal filter to spectral coefficients."""
        arr = np.asarray(coefficients)
        if arr.ndim == 0:
            raise ValueError("coefficients must have at least one dimension")
        coeff_axis = int(axis)
        if coeff_axis < 0:
            coeff_axis += arr.ndim
        if coeff_axis < 0 or coeff_axis >= arr.ndim:
            raise ValueError("axis is out of bounds")
        if arr.shape[coeff_axis] != self.nspec2:
            actual = infer_ntrunc_from_nspec2(arr.shape[coeff_axis])
            raise ValueError(f"filter is T{self.ntrunc}, but coefficients are T{actual}")
        shape = [1] * arr.ndim
        shape[coeff_axis] = self.nspec2
        return arr * self.coefficient_weights.reshape(shape)

    def apply(
        self,
        values,
        *,
        grid,
        axis=-1,
        chunk_size=DEFAULT_CHUNK_SIZE,
        missing_policy="error",
    ) -> np.ndarray:
        """Fit grid values, apply this filter in spectral space, and synthesize back."""
        coeffs = spectral_fit(
            values,
            grid=grid,
            ntrunc=self.ntrunc,
            axis=axis,
            chunk_size=chunk_size,
            missing_policy=missing_policy,
        )
        return spectral_synthesis(
            self.apply_coefficients(coeffs),
            grid=grid,
            ntrunc=self.ntrunc,
            axis=-1,
            chunk_size=chunk_size,
        )

    def to_metadata(self) -> dict:
        """Return serializable metadata describing this filter."""
        return {
            "version": 1,
            "kind": self.kind,
            "ntrunc": self.ntrunc,
            "cutoff": self.cutoff,
            "selectivity": self.selectivity,
            "low_pass_exponent": self.low_pass_exponent,
            "gaussian_exponent": self.gaussian_exponent,
            "layout": "ectrans_global_real_imag_by_total_wave_number",
        }
