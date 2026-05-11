from __future__ import annotations

import numpy as np


def periodic_interp_row(row: np.ndarray, target_nlon: int) -> np.ndarray:
    """Linearly interpolate one periodic longitude row to ``target_nlon`` points."""
    row = np.asarray(row, dtype=np.float32)
    source_nlon = row.size
    if source_nlon == target_nlon:
        return row.copy()
    source_x = np.arange(source_nlon + 1, dtype=np.float64) * (360.0 / source_nlon)
    target_x = np.arange(target_nlon, dtype=np.float64) * (360.0 / target_nlon)
    periodic_values = np.concatenate([row, row[:1]])
    return np.interp(target_x, source_x, periodic_values).astype(np.float32)


def unpack_reduced(values: np.ndarray, pl: np.ndarray, target_nlon: int) -> np.ndarray:
    """Unpack a reduced Gaussian packed field to a regular 2D display array."""
    values = np.asarray(values, dtype=np.float32)
    if values.ndim != 1:
        raise ValueError("packed reduced field must be 1D")
    if values.size != int(pl.sum()):
        raise ValueError(f"packed field has {values.size} points, expected {int(pl.sum())}")

    out = np.empty((pl.size, target_nlon), dtype=np.float32)
    offset = 0
    for j, row_nlon in enumerate(pl.astype(int)):
        row = values[offset : offset + row_nlon]
        out[j, :] = periodic_interp_row(row, target_nlon)
        offset += row_nlon
    return out


def pack_reduced(values: np.ndarray, pl: np.ndarray) -> np.ndarray:
    """Pack a regular 2D field into reduced Gaussian row lengths."""
    values = np.asarray(values, dtype=np.float32)
    if values.ndim != 2:
        raise ValueError("regular field to pack must be 2D")
    if values.shape[0] != pl.size:
        raise ValueError(f"regular field has {values.shape[0]} rows, expected {pl.size}")

    rows = [periodic_interp_row(values[j, :], int(row_nlon)) for j, row_nlon in enumerate(pl)]
    return np.concatenate(rows).astype(np.float32)
