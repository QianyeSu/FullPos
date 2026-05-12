from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

from .errors import FullposBackendError

__all__ = ["to_grib"]


@dataclass(frozen=True)
class _TemplateMessage:
    index: int
    short_name: str
    number_of_points: int


@dataclass(frozen=True)
class _VariableRecord:
    short_name: str
    fields: np.ndarray


def to_grib(
    obj,
    path,
    *,
    template,
    variables=None,
    append: bool = False,
    strict: bool = True,
) -> Path:
    """Write xarray data to GRIB by cloning messages from a template file.

    GRIB is a message-based format whose grid, level, time, packing, and
    discipline metadata are not safely reconstructable from xarray metadata
    alone. This writer therefore requires ``template=`` and preserves each
    selected template message as-is except for its data values.

    Parameters
    ----------
    obj:
        Input ``xarray.DataArray`` or ``xarray.Dataset``.
    path:
        Output GRIB file path.
    template:
        Existing GRIB/GRIB2 file whose messages are cloned.
    variables:
        Optional Dataset variable names to write. DataArray input ignores this
        argument and writes the single array.
    append:
        Append to an existing GRIB file instead of replacing it.
    strict:
        When true, each written variable must have exactly the same number of
        fields as matching template messages. When false, the first matching
        template messages are used and extra template messages are skipped.

    Returns
    -------
    pathlib.Path
        The output path.
    """
    eccodes = _import_eccodes()
    output_path = Path(path)
    template_path = Path(template)
    if not template_path.exists():
        raise FileNotFoundError(f"template GRIB file not found: {template_path}")

    template_messages = _read_template_messages(eccodes, template_path)
    if not template_messages:
        raise ValueError(f"template GRIB file contains no messages: {template_path}")

    records = _prepare_records(
        obj,
        template_messages=template_messages,
        variables=variables,
        strict=strict,
    )
    if not records:
        raise ValueError("no variables were selected for GRIB output")

    queues = {record.short_name: deque(record.fields) for record in records}
    output_path.parent.mkdir(parents=True, exist_ok=True)
    mode = "ab" if append else "wb"

    with template_path.open("rb") as template_handle, output_path.open(mode) as output_handle:
        for message in template_messages:
            queue = queues.get(message.short_name)
            if not queue:
                _release_next_template_message(eccodes, template_handle)
                continue
            gid = eccodes.codes_grib_new_from_file(template_handle)
            if gid is None:
                raise ValueError("template GRIB ended unexpectedly while writing")
            clone = None
            try:
                clone = eccodes.codes_clone(gid)
                values = np.asarray(queue.popleft(), dtype=np.float64)
                eccodes.codes_set_values(clone, values)
                eccodes.codes_write(clone, output_handle)
            finally:
                if clone is not None:
                    eccodes.codes_release(clone)
                eccodes.codes_release(gid)

    remaining = {name: len(queue) for name, queue in queues.items() if queue}
    if remaining:
        raise ValueError(f"template GRIB did not provide enough messages: {remaining}")
    return output_path


def _import_eccodes():
    try:
        import eccodes
    except ImportError as exc:  # pragma: no cover - depends on local environment
        raise FullposBackendError(
            "to_grib requires the ecCodes Python bindings. Install the optional "
            "GRIB dependencies with `pip install fullpos[grib]` or install "
            "`eccodes` in the active environment."
        ) from exc
    return eccodes


def _read_template_messages(eccodes, template_path: Path) -> list[_TemplateMessage]:
    messages: list[_TemplateMessage] = []
    with template_path.open("rb") as handle:
        while True:
            gid = eccodes.codes_grib_new_from_file(handle)
            if gid is None:
                break
            try:
                short_name = _get_required_key(eccodes, gid, "shortName", len(messages))
                number_of_points = int(
                    _get_required_key(eccodes, gid, "numberOfPoints", len(messages))
                )
                messages.append(
                    _TemplateMessage(
                        index=len(messages),
                        short_name=str(short_name),
                        number_of_points=number_of_points,
                    )
                )
            finally:
                eccodes.codes_release(gid)
    return messages


def _get_required_key(eccodes, gid, key: str, index: int):
    try:
        value = eccodes.codes_get(gid, key)
    except Exception as exc:
        raise ValueError(f"template message {index} is missing GRIB key {key!r}") from exc
    if value is None:
        raise ValueError(f"template message {index} has empty GRIB key {key!r}")
    return value


def _release_next_template_message(eccodes, template_handle) -> None:
    gid = eccodes.codes_grib_new_from_file(template_handle)
    if gid is not None:
        eccodes.codes_release(gid)


def _prepare_records(
    obj,
    *,
    template_messages: list[_TemplateMessage],
    variables,
    strict: bool,
) -> list[_VariableRecord]:
    if isinstance(obj, xr.DataArray):
        if variables is not None:
            raise ValueError("variables is only supported for xarray Dataset input")
        return [
            _prepare_data_array_record(
                obj,
                template_messages=template_messages,
                explicit_short_name=None,
                strict=strict,
            )
        ]

    if isinstance(obj, xr.Dataset):
        selected = list(obj.data_vars) if variables is None else [str(name) for name in variables]
        missing = [name for name in selected if name not in obj.data_vars]
        if missing:
            raise KeyError(f"variables not found in dataset: {missing}")
        records = [
            _prepare_data_array_record(
                obj[name],
                template_messages=template_messages,
                explicit_short_name=name,
                strict=strict,
            )
            for name in selected
        ]
        short_names = [record.short_name for record in records]
        duplicates = [name for name, count in Counter(short_names).items() if count > 1]
        if duplicates:
            raise ValueError(f"multiple variables map to the same GRIB shortName: {duplicates}")
        return records

    raise TypeError("obj must be an xarray DataArray or Dataset")


def _prepare_data_array_record(
    data_array: xr.DataArray,
    *,
    template_messages: list[_TemplateMessage],
    explicit_short_name: str | None,
    strict: bool,
) -> _VariableRecord:
    short_name = _infer_short_name(data_array, template_messages, explicit_short_name)
    matches = [message for message in template_messages if message.short_name == short_name]
    if not matches:
        raise ValueError(f"template GRIB contains no messages for shortName={short_name!r}")

    point_counts = {message.number_of_points for message in matches}
    if len(point_counts) != 1:
        raise ValueError(
            f"template messages for shortName={short_name!r} use mixed numberOfPoints"
        )
    number_of_points = point_counts.pop()
    fields = _data_array_to_fields(data_array, number_of_points)

    if strict and fields.shape[0] != len(matches):
        raise ValueError(
            f"variable {short_name!r} has {fields.shape[0]} field(s), but template "
            f"contains {len(matches)} matching message(s); use strict=False to write "
            "a leading subset"
        )
    if fields.shape[0] > len(matches):
        raise ValueError(
            f"variable {short_name!r} has {fields.shape[0]} field(s), but template "
            f"only contains {len(matches)} matching message(s)"
        )
    return _VariableRecord(short_name=short_name, fields=fields)


def _infer_short_name(
    data_array: xr.DataArray,
    template_messages: list[_TemplateMessage],
    explicit_short_name: str | None,
) -> str:
    for key in ("GRIB_shortName", "shortName"):
        value = data_array.attrs.get(key)
        if value not in (None, ""):
            return str(value)
    if data_array.name not in (None, ""):
        return str(data_array.name)
    if explicit_short_name not in (None, ""):
        return str(explicit_short_name)
    short_names = {message.short_name for message in template_messages}
    if len(short_names) == 1:
        return next(iter(short_names))
    raise ValueError(
        "cannot infer GRIB shortName for unnamed DataArray; set the DataArray "
        "name or GRIB_shortName attr"
    )


def _data_array_to_fields(data_array: xr.DataArray, number_of_points: int) -> np.ndarray:
    horizontal_dims = _infer_horizontal_dims(data_array, number_of_points)
    leading_dims = tuple(dim for dim in data_array.dims if dim not in horizontal_dims)
    ordered = data_array.transpose(*leading_dims, *horizontal_dims)
    values = np.asarray(ordered.values)
    horizontal_shape = values.shape[-len(horizontal_dims) :]
    horizontal_size = int(np.prod(horizontal_shape, dtype=np.int64))
    if horizontal_size != number_of_points:
        raise ValueError(
            f"horizontal shape {horizontal_shape} has {horizontal_size} point(s), "
            f"but template message requires {number_of_points}"
        )
    return np.ascontiguousarray(values.reshape((-1, number_of_points)))


def _infer_horizontal_dims(data_array: xr.DataArray, number_of_points: int) -> tuple[str, ...]:
    dims = tuple(data_array.dims)
    sizes = data_array.sizes
    if "values" in sizes and int(sizes["values"]) == number_of_points:
        return ("values",)
    lat_dim = _find_dim(dims, ("latitude", "lat"))
    lon_dim = _find_dim(dims, ("longitude", "lon"))
    if lat_dim is not None and lon_dim is not None:
        if int(sizes[lat_dim]) * int(sizes[lon_dim]) == number_of_points:
            return (lat_dim, lon_dim)
    if dims and int(sizes[dims[-1]]) == number_of_points:
        return (dims[-1],)
    if len(dims) >= 2 and int(sizes[dims[-2]]) * int(sizes[dims[-1]]) == number_of_points:
        return (dims[-2], dims[-1])
    raise ValueError(
        f"cannot find horizontal dimension(s) matching template numberOfPoints="
        f"{number_of_points}"
    )


def _find_dim(dims: tuple[str, ...], names: tuple[str, ...]) -> str | None:
    for name in names:
        if name in dims:
            return name
    return None
