from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping

import numpy as np


@dataclass(frozen=True)
class GaussianGrid:
    """Description of a supported ECMWF Gaussian grid."""

    name: str
    n: int
    kind: str
    nlat: int
    work_nlon: int
    pl: tuple[int, ...] | None = None

    @property
    def is_reduced(self) -> bool:
        """Whether this grid is stored as a reduced packed Gaussian grid."""
        return self.pl is not None

    @property
    def size(self) -> int:
        """Total number of horizontal grid points."""
        if self.pl is None:
            return self.nlat * self.work_nlon
        return int(sum(self.pl))


def octahedral_pl(n: int) -> np.ndarray:
    """Return ECMWF octahedral reduced Gaussian row lengths for ``O<n>``."""
    if n <= 0:
        raise ValueError("octahedral resolution must be positive")
    north = 20 + 4 * np.arange(2 * n // 2, dtype=np.int64)
    return np.concatenate([north, north[::-1]])


def classic_reduced_grid_from_pl(n: int, pl) -> GaussianGrid:
    """Build a classic reduced Gaussian grid from GRIB row lengths."""
    if n <= 0:
        raise ValueError("grid resolution must be positive")
    arr = np.asarray(pl, dtype=np.int64)
    if arr.ndim != 1:
        arr = arr.reshape(-1)
    if arr.size == 0 or arr.size % 2:
        raise ValueError("GRIB_pl must be a non-empty 1D array with 2*N rows")
    if np.any(arr <= 0):
        raise ValueError("GRIB_pl values must be positive")
    inferred_n = int(arr.size // 2)
    if int(n) != inferred_n:
        raise ValueError(f"GRIB_N={n} is inconsistent with GRIB_pl length {arr.size}")
    if not np.array_equal(arr, arr[::-1]):
        raise ValueError("classic reduced Gaussian grid must be symmetric north/south")
    return GaussianGrid(
        name=f"N{int(n)}",
        n=int(n),
        kind="classic_reduced",
        nlat=2 * int(n),
        work_nlon=int(arr.max()),
        pl=tuple(int(v) for v in arr),
    )


def gaussian_latitudes(nlat: int) -> np.ndarray:
    """Return Gaussian latitude coordinates in north-to-south order."""
    if nlat <= 0:
        raise ValueError("nlat must be positive")
    nodes, _weights = np.polynomial.legendre.leggauss(nlat)
    return np.rad2deg(np.arcsin(nodes))[::-1]


def regular_longitudes(nlon: int) -> np.ndarray:
    """Return regular longitude coordinates in degrees east."""
    if nlon <= 0:
        raise ValueError("nlon must be positive")
    return np.arange(nlon, dtype=np.float64) * (360.0 / nlon)


def parse_grid(name: str) -> GaussianGrid:
    """Parse a grid name such as ``O320`` or ``F320``."""
    if not isinstance(name, str) or not name.strip():
        raise ValueError("grid name must be a non-empty string")

    normalized = name.strip().upper()
    prefix = normalized[0]
    if prefix not in {"O", "F", "N"}:
        raise ValueError("Only O* octahedral and F* regular Gaussian grids are supported")

    try:
        n = int(normalized[1:])
    except ValueError as exc:
        raise ValueError(f"invalid Gaussian grid name: {name!r}") from exc
    if n <= 0:
        raise ValueError(f"grid resolution must be positive: {name!r}")

    nlat = 2 * n
    if prefix == "O":
        pl = octahedral_pl(n)
        return GaussianGrid(
            name=normalized,
            n=n,
            kind="octahedral",
            nlat=nlat,
            work_nlon=int(pl.max()),
            pl=tuple(int(v) for v in pl),
        )

    regular_name = f"F{n}" if prefix == "F" else normalized
    return GaussianGrid(
        name=regular_name,
        n=n,
        kind="regular",
        nlat=nlat,
        work_nlon=4 * n,
        pl=None,
    )


def infer_grid_from_attrs(attrs: Mapping | None) -> str | GaussianGrid:
    """Infer a supported grid from GRIB-style xarray attributes."""
    if not attrs:
        raise ValueError("source_grid is required when it cannot be inferred from metadata")

    grid_type = _string_attr(attrs, "GRIB_gridType", "gridType")
    n = _int_attr(attrs, "GRIB_N", "N")
    pl = _pl_attr(attrs, "GRIB_pl", "pl")

    if grid_type == "reduced_gg":
        if pl is None:
            raise ValueError(
                "reduced Gaussian source grid requires GRIB_pl metadata for safe inference"
            )
        return _infer_reduced_grid_from_pl(pl, n)

    if grid_type == "regular_gg":
        if n is None:
            raise ValueError("regular Gaussian source grid requires GRIB_N metadata")
        _validate_regular_point_count(attrs, n)
        return f"F{n}"

    if pl is not None:
        return _infer_reduced_grid_from_pl(pl, n)

    if n is not None and grid_type in {None, ""}:
        return f"N{n}"

    if grid_type:
        raise ValueError(f"unsupported GRIB_gridType for spectral regridding: {grid_type!r}")
    raise ValueError("source_grid is required when it cannot be inferred from metadata")


def infer_grid_name_from_attrs(attrs: Mapping | None) -> str:
    """Infer a supported grid name from GRIB-style xarray attributes."""
    grid = infer_grid_from_attrs(attrs)
    return grid.name if isinstance(grid, GaussianGrid) else grid


def infer_grid_name_from_shape(sizes: Mapping[str, int], dims: tuple[str, ...]) -> str:
    """Infer a regular Gaussian ``F`` grid from latitude/longitude dimensions."""
    if "latitude" not in dims or "longitude" not in dims:
        raise ValueError("source_grid is required when it cannot be inferred from metadata")
    nlat = int(sizes["latitude"])
    nlon = int(sizes["longitude"])
    if nlat <= 0 or nlat % 2:
        raise ValueError("regular Gaussian latitude dimension must be a positive even size")
    n = nlat // 2
    if nlon != 4 * n:
        raise ValueError(
            "regular Gaussian longitude dimension must be 4*N for automatic inference"
        )
    return f"F{n}"


def _infer_octahedral_n_from_pl(pl: np.ndarray, n: int | None) -> int:
    if pl.ndim != 1 or pl.size == 0 or pl.size % 2:
        raise ValueError("GRIB_pl must be a non-empty 1D array with 2*N rows")
    inferred_n = int(pl.size // 2)
    if n is not None and int(n) != inferred_n:
        raise ValueError(f"GRIB_N={n} is inconsistent with GRIB_pl length {pl.size}")
    expected = octahedral_pl(inferred_n)
    if not np.array_equal(pl.astype(np.int64), expected):
        raise ValueError(
            "reduced Gaussian grid is not an ECMWF octahedral O-grid; "
            "only O* reduced and F* regular Gaussian grids are currently supported"
        )
    return inferred_n


def _infer_reduced_grid_from_pl(pl: np.ndarray, n: int | None) -> str | GaussianGrid:
    if pl.ndim != 1 or pl.size == 0 or pl.size % 2:
        raise ValueError("GRIB_pl must be a non-empty 1D array with 2*N rows")
    if np.any(pl <= 0):
        raise ValueError("GRIB_pl values must be positive")
    inferred_n = int(pl.size // 2) if n is None else int(n)
    if pl.size != 2 * inferred_n:
        raise ValueError(f"GRIB_N={n} is inconsistent with GRIB_pl length {pl.size}")
    if np.array_equal(pl.astype(np.int64), octahedral_pl(inferred_n)):
        return f"O{inferred_n}"
    return classic_reduced_grid_from_pl(inferred_n, pl)


def _validate_regular_point_count(attrs: Mapping, n: int) -> None:
    number_of_points = _int_attr(attrs, "GRIB_numberOfPoints", "numberOfPoints")
    if number_of_points is None:
        return
    expected = 8 * int(n) * int(n)
    if int(number_of_points) != expected:
        raise ValueError(
            f"regular Gaussian point count {number_of_points} is inconsistent with F{n} "
            f"(expected {expected})"
        )


def _string_attr(attrs: Mapping, *names: str) -> str | None:
    value = _first_attr(attrs, *names)
    if value is None:
        return None
    return str(value).strip()


def _int_attr(attrs: Mapping, *names: str) -> int | None:
    value = _first_attr(attrs, *names)
    if value is None:
        return None
    arr = np.asarray(value)
    if arr.size != 1:
        raise ValueError(f"{names[0]} must be scalar")
    return int(arr.reshape(-1)[0])


def _pl_attr(attrs: Mapping, *names: str) -> np.ndarray | None:
    value = _first_attr(attrs, *names)
    if value is None:
        return None
    arr = np.asarray(value, dtype=np.int64)
    if arr.ndim != 1:
        arr = arr.reshape(-1)
    return arr


def _first_attr(attrs: Mapping, *names: str):
    for name in names:
        if name in attrs:
            return attrs[name]
    return None
