import numpy as np
import pytest

from fullpos.grids import (
    GaussianGrid,
    classic_reduced_grid_from_pl,
    infer_grid_from_attrs,
    infer_grid_name_from_attrs,
    octahedral_pl,
    parse_grid,
)


def test_parse_o320() -> None:
    grid = parse_grid("O320")
    assert grid.nlat == 640
    assert grid.work_nlon == 1296
    assert grid.size == 421120
    assert grid.is_reduced


def test_parse_f320() -> None:
    grid = parse_grid("F320")
    assert grid.nlat == 640
    assert grid.work_nlon == 1280
    assert grid.size == 819200
    assert not grid.is_reduced


def test_parse_n320_alias_for_regular_gaussian() -> None:
    grid = parse_grid("N320")
    assert grid.name == "N320"
    assert grid.nlat == 640
    assert grid.work_nlon == 1280
    assert grid.size == 819200
    assert not grid.is_reduced


def test_octahedral_max_row_uses_ecmwf_4n_plus_16_rule() -> None:
    assert parse_grid("O96").work_nlon == 400
    assert parse_grid("O160").work_nlon == 656
    assert parse_grid("O320").work_nlon == 1296
    assert parse_grid("F160").work_nlon == 640


def test_octahedral_pl_o320_matches_expected_edges() -> None:
    pl = octahedral_pl(320)
    assert pl.size == 640
    assert int(pl.sum()) == 421120
    assert pl[:5].tolist() == [20, 24, 28, 32, 36]
    assert pl[-5:].tolist() == [36, 32, 28, 24, 20]
    assert np.array_equal(pl, pl[::-1])


def test_classic_reduced_grid_from_pl_builds_gaussian_grid() -> None:
    pl = np.array([10, 12, 14, 16, 16, 14, 12, 10], dtype=np.int64)
    grid = classic_reduced_grid_from_pl(4, pl)

    assert isinstance(grid, GaussianGrid)
    assert grid.name == "N4"
    assert grid.kind == "classic_reduced"
    assert grid.nlat == 8
    assert grid.work_nlon == 16
    assert grid.size == int(pl.sum())
    assert grid.is_reduced


def test_infer_reduced_gaussian_from_grib_attrs() -> None:
    attrs = {
        "GRIB_gridType": "reduced_gg",
        "GRIB_N": 4,
        "GRIB_pl": octahedral_pl(4),
        "GRIB_numberOfPoints": int(octahedral_pl(4).sum()),
    }
    assert infer_grid_name_from_attrs(attrs) == "O4"


def test_infer_classic_reduced_gaussian_from_grib_attrs() -> None:
    pl = np.array([10, 12, 14, 16, 16, 14, 12, 10], dtype=np.int64)
    attrs = {
        "GRIB_gridType": "reduced_gg",
        "GRIB_N": 4,
        "GRIB_pl": pl,
        "GRIB_numberOfPoints": int(pl.sum()),
    }
    grid = infer_grid_from_attrs(attrs)

    assert isinstance(grid, GaussianGrid)
    assert grid.name == "N4"
    assert grid.kind == "classic_reduced"
    assert infer_grid_name_from_attrs(attrs) == "N4"


def test_infer_regular_gaussian_from_grib_attrs() -> None:
    attrs = {
        "GRIB_gridType": "regular_gg",
        "GRIB_N": 4,
        "GRIB_numberOfPoints": 128,
    }
    assert infer_grid_name_from_attrs(attrs) == "F4"


def test_infers_classic_reduced_gaussian_from_non_octahedral_attrs() -> None:
    attrs = {
        "GRIB_gridType": "reduced_gg",
        "GRIB_N": 4,
        "GRIB_pl": np.full(8, 16),
    }
    assert infer_grid_name_from_attrs(attrs) == "N4"
