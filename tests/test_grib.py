from __future__ import annotations

import sys
import types

import numpy as np
import pytest
import xarray as xr

from fullpos import to_grib


class _FakeEccodes:
    def __init__(self, messages):
        self.messages = [dict(message) for message in messages]
        self.positions = {}
        self.writes = []
        self.releases = []

    def module(self):
        module = types.ModuleType("eccodes")
        module.codes_grib_new_from_file = self.codes_grib_new_from_file
        module.codes_get = self.codes_get
        module.codes_clone = self.codes_clone
        module.codes_set_values = self.codes_set_values
        module.codes_write = self.codes_write
        module.codes_release = self.codes_release
        return module

    def codes_grib_new_from_file(self, handle):
        index = self.positions.get(handle, 0)
        if index >= len(self.messages):
            return None
        self.positions[handle] = index + 1
        return {
            "keys": dict(self.messages[index]),
            "source_index": index,
            "values": None,
            "clone": False,
        }

    def codes_get(self, gid, key):
        if key not in gid["keys"]:
            raise KeyError(key)
        return gid["keys"][key]

    def codes_clone(self, gid):
        return {
            "keys": dict(gid["keys"]),
            "source_index": gid["source_index"],
            "values": None,
            "clone": True,
        }

    def codes_set_values(self, gid, values):
        gid["values"] = np.asarray(values).copy()

    def codes_write(self, gid, handle):
        self.writes.append(
            {
                "shortName": gid["keys"]["shortName"],
                "source_index": gid["source_index"],
                "values": gid["values"].copy(),
            }
        )
        handle.write(b"GRIB")

    def codes_release(self, gid):
        self.releases.append((gid["source_index"], gid["clone"]))


def _install_fake_eccodes(monkeypatch, messages):
    fake = _FakeEccodes(messages)
    monkeypatch.setitem(sys.modules, "eccodes", fake.module())
    return fake


def test_to_grib_dataarray_clones_template_and_sets_values(tmp_path, monkeypatch) -> None:
    fake = _install_fake_eccodes(
        monkeypatch,
        [{"shortName": "t", "numberOfPoints": 4}],
    )
    template = tmp_path / "template.grib2"
    output = tmp_path / "out.grib2"
    template.write_bytes(b"TEMPLATE")
    field = xr.DataArray(
        np.array([[1.0, 2.0], [3.0, 4.0]], dtype=np.float32),
        dims=("latitude", "longitude"),
        name="t",
    )

    result = to_grib(field, output, template=template)

    assert result == output
    assert output.read_bytes() == b"GRIB"
    assert len(fake.writes) == 1
    assert fake.writes[0]["shortName"] == "t"
    assert fake.writes[0]["source_index"] == 0
    np.testing.assert_array_equal(fake.writes[0]["values"], [1.0, 2.0, 3.0, 4.0])


def test_to_grib_dataset_writes_matching_template_messages_in_order(
    tmp_path, monkeypatch
) -> None:
    fake = _install_fake_eccodes(
        monkeypatch,
        [
            {"shortName": "t", "numberOfPoints": 2},
            {"shortName": "q", "numberOfPoints": 2},
            {"shortName": "t", "numberOfPoints": 2},
            {"shortName": "u", "numberOfPoints": 2},
        ],
    )
    template = tmp_path / "template.grib2"
    output = tmp_path / "out.grib2"
    template.write_bytes(b"TEMPLATE")
    dataset = xr.Dataset(
        {
            "t": (("time", "values"), np.array([[1.0, 2.0], [3.0, 4.0]])),
            "u": (("values",), np.array([5.0, 6.0])),
        }
    )

    to_grib(dataset, output, template=template)

    assert output.read_bytes() == b"GRIBGRIBGRIB"
    assert [write["shortName"] for write in fake.writes] == ["t", "t", "u"]
    assert [write["source_index"] for write in fake.writes] == [0, 2, 3]
    np.testing.assert_array_equal(fake.writes[0]["values"], [1.0, 2.0])
    np.testing.assert_array_equal(fake.writes[1]["values"], [3.0, 4.0])
    np.testing.assert_array_equal(fake.writes[2]["values"], [5.0, 6.0])


def test_to_grib_rejects_template_count_mismatch_by_default(tmp_path, monkeypatch) -> None:
    _install_fake_eccodes(
        monkeypatch,
        [
            {"shortName": "t", "numberOfPoints": 2},
            {"shortName": "t", "numberOfPoints": 2},
        ],
    )
    template = tmp_path / "template.grib2"
    output = tmp_path / "out.grib2"
    template.write_bytes(b"TEMPLATE")
    field = xr.DataArray(np.array([1.0, 2.0]), dims=("values",), name="t")

    with pytest.raises(ValueError, match="has 1 field"):
        to_grib(field, output, template=template)


def test_to_grib_strict_false_writes_leading_template_subset(tmp_path, monkeypatch) -> None:
    fake = _install_fake_eccodes(
        monkeypatch,
        [
            {"shortName": "t", "numberOfPoints": 2},
            {"shortName": "t", "numberOfPoints": 2},
        ],
    )
    template = tmp_path / "template.grib2"
    output = tmp_path / "out.grib2"
    template.write_bytes(b"TEMPLATE")
    field = xr.DataArray(np.array([1.0, 2.0]), dims=("values",), name="t")

    to_grib(field, output, template=template, strict=False)

    assert output.read_bytes() == b"GRIB"
    assert [write["source_index"] for write in fake.writes] == [0]


def test_to_grib_append_mode_preserves_existing_bytes(tmp_path, monkeypatch) -> None:
    _install_fake_eccodes(
        monkeypatch,
        [{"shortName": "t", "numberOfPoints": 2}],
    )
    template = tmp_path / "template.grib2"
    output = tmp_path / "out.grib2"
    template.write_bytes(b"TEMPLATE")
    output.write_bytes(b"OLD")
    field = xr.DataArray(np.array([1.0, 2.0]), dims=("values",), name="t")

    to_grib(field, output, template=template, append=True)

    assert output.read_bytes() == b"OLDGRIB"


def test_to_grib_requires_dataset_variables_to_exist(tmp_path, monkeypatch) -> None:
    _install_fake_eccodes(
        monkeypatch,
        [{"shortName": "t", "numberOfPoints": 2}],
    )
    template = tmp_path / "template.grib2"
    output = tmp_path / "out.grib2"
    template.write_bytes(b"TEMPLATE")
    dataset = xr.Dataset({"t": (("values",), np.array([1.0, 2.0]))})

    with pytest.raises(KeyError, match="variables not found"):
        to_grib(dataset, output, template=template, variables=["q"])
