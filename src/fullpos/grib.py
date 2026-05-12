from __future__ import annotations

from collections import Counter, deque
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import xarray as xr

from .errors import FullposBackendError

__all__ = ["to_grib"]

_PACKING_ALIASES = {
    "template": None,
    "simple": "grid_simple",
    "ccsds": "grid_ccsds",
    "aec": "grid_ccsds",
    "second_order": "grid_second_order",
    "jpeg": "grid_jpeg",
    "png": "grid_png",
    "ieee": "grid_ieee",
}


@dataclass(frozen=True)
class _TemplateMessage:
    index: int
    short_name: str
    number_of_points: int
    param_id: int | None = None


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
    packing_type: str | None = None,
    bits_per_value: int | None = None,
    edition: int | None = None,
    centre=None,
    sub_centre: int | None = None,
    generating_process_identifier: int | None = None,
    data_date=None,
    data_time=None,
    step_type: str | None = None,
    step_range=None,
    forecast_time: int | None = None,
    param_id: int | None = None,
    short_name: str | None = None,
    type_of_level: str | None = None,
    level=None,
    key_overrides: Mapping[str, object] | None = None,
) -> Path:
    """Write xarray data to GRIB by cloning messages from a template file.

    GRIB is a message-based format whose grid, level, time, packing, and
    discipline metadata are not safely reconstructable from xarray metadata
    alone. This writer therefore requires ``template=`` and preserves each
    selected template message as-is except for its data values.

    The template must already match the desired output geometry and message
    layout. For example, write an ``O480`` result with an ``O480`` template, not
    with the original ``O320`` input file. Run ``regrid`` or
    ``vertical_interpolate`` first, then pass the processed xarray object here.

    Parameters
    ----------
    obj:
        Input ``xarray.DataArray`` or ``xarray.Dataset``. Supported horizontal
        layouts are packed reduced-grid data with a ``values`` dimension, or
        two-dimensional regular data with ``latitude``/``longitude`` or
        ``lat``/``lon`` dimensions. Leading dimensions such as ``time`` and
        ``hybrid`` are flattened in xarray order and written as one GRIB
        message per horizontal field.
    path:
        Output GRIB file path.
    template:
        Existing GRIB/GRIB2 file whose messages are cloned. Template messages
        are matched by GRIB ``shortName``. Grid, level, time, product
        definition, bitmap, and packing metadata are inherited from the
        template unless an explicit packing override is supplied.
    variables:
        Optional Dataset variable names to write. DataArray input ignores this
        argument and writes the single array. Dataset variable names are used as
        GRIB ``shortName`` unless the variable has ``GRIB_shortName`` or
        ``shortName`` in ``attrs``.
    append:
        Append to an existing GRIB file instead of replacing it.
    strict:
        When true, each written variable must have exactly the same number of
        fields as matching template messages. When false, the first matching
        template messages are used and extra template messages are skipped.
    packing_type:
        Optional ecCodes ``packingType`` override. ``None`` or ``"template"``
        preserves the template packing. Common aliases include ``"simple"``
        for ``"grid_simple"`` and ``"ccsds"``/``"aec"`` for
        ``"grid_ccsds"``.
    bits_per_value:
        Optional ecCodes ``bitsPerValue`` override. Leave unset to preserve
        the template precision/packing settings.
    edition:
        Optional GRIB edition override, normally ``1`` or ``2``. By default the
        cloned template edition is preserved. Prefer a GRIB2 template over
        converting GRIB1 messages at write time.
    centre, sub_centre:
        Optional producing centre overrides, mapped to ecCodes ``centre`` and
        ``subCentre``.
    generating_process_identifier:
        Optional provenance override mapped to
        ``generatingProcessIdentifier``. This can be used to mark fields as
        produced by a local FullPos workflow.
    data_date, data_time:
        Optional reference date/time overrides mapped to ``dataDate`` and
        ``dataTime``. ecCodes commonly represents these as ``YYYYMMDD`` and
        ``HHMM`` integers.
    step_type, step_range, forecast_time:
        Optional forecast-step overrides mapped to ``stepType``, ``stepRange``,
        and ``forecastTime``. ``step_range`` and ``forecast_time`` are mutually
        exclusive because they can describe the same GRIB step in different
        ways.
    param_id, short_name:
        Optional parameter identity overrides mapped to ``paramId`` and
        ``shortName``. These are mutually exclusive because ecCodes derives one
        from the other and conflicting values can create invalid metadata.
        These overrides are also used to select matching template messages, and
        can only be used with a single DataArray or one selected Dataset
        variable.
    type_of_level, level:
        Optional level overrides mapped to ``typeOfLevel`` and ``level``.
        Prefer inheriting these from the template unless the new message layout
        is known to match.
    key_overrides:
        Low-level ecCodes key/value overrides applied after the named
        arguments. This is intended for advanced metadata keys only; do not use
        it to change geometry keys such as ``numberOfPoints`` or ``pl`` unless
        the values and template message are already consistent.

    Returns
    -------
    pathlib.Path
        The output path.

    Examples
    --------
    Write a regridded DataArray with template packing preserved:

    >>> from fullpos import regrid, to_grib
    >>> out = regrid(ds["t"].isel(time=0), target_grid="O480")
    >>> output_path = to_grib(out, "t_o480.grib2", template="template_o480_t.grib2")

    Write selected Dataset variables. The template must contain matching
    ``shortName`` messages for ``t``, ``u``, and ``v``:

    >>> output_path = to_grib(
    ...     ds_out,
    ...     "tuv_o480.grib2",
    ...     template="template_o480_tuv.grib2",
    ...     variables=["t", "u", "v"],
    ... )

    Request GRIB2 CCSDS/AEC packing through ecCodes:

    >>> output_path = to_grib(
    ...     out,
    ...     "t_o480_ccsds.grib2",
    ...     template="template_o480_t.grib2",
    ...     packing_type="ccsds",
    ...     bits_per_value=16,
    ... )

    Override common metadata keys while preserving the rest of the template:

    >>> output_path = to_grib(
    ...     out,
    ...     "t_o480_20250102.grib2",
    ...     template="template_o480_t.grib2",
    ...     centre="ecmf",
    ...     generating_process_identifier=255,
    ...     data_date=20250102,
    ...     data_time=600,
    ...     step_type="instant",
    ...     forecast_time=0,
    ... )

    Use low-level ecCodes keys only when you know the target message metadata:

    >>> output_path = to_grib(
    ...     out,
    ...     "t_o480_custom.grib2",
    ...     template="template_o480_t.grib2",
    ...     key_overrides={"localDefinitionNumber": 1},
    ... )

    Use ``strict=False`` only when intentionally writing a leading subset of
    matching template messages, for example 10 levels into the first 10
    ``shortName="t"`` template messages:

    >>> output_path = to_grib(
    ...     t_ml1_10,
    ...     "t_ml1_10.grib2",
    ...     template="template_ml137.grib2",
    ...     strict=False,
    ... )

    Notes
    -----
    ``packing_type="ccsds"`` maps to ecCodes ``grid_ccsds`` and requires a
    GRIB2 template plus an ecCodes build with libaec support. If ecCodes cannot
    encode the requested packing, the underlying ecCodes error is raised.
    """
    eccodes = _import_eccodes()
    output_path = Path(path)
    template_path = Path(template)
    normalized_packing_type = _normalize_packing_type(packing_type)
    normalized_bits_per_value = _normalize_bits_per_value(bits_per_value)
    message_key_overrides = _build_message_key_overrides(
        edition=edition,
        centre=centre,
        sub_centre=sub_centre,
        generating_process_identifier=generating_process_identifier,
        data_date=data_date,
        data_time=data_time,
        step_type=step_type,
        step_range=step_range,
        forecast_time=forecast_time,
        param_id=param_id,
        short_name=short_name,
        type_of_level=type_of_level,
        level=level,
        key_overrides=key_overrides,
    )
    if not template_path.exists():
        raise FileNotFoundError(f"template GRIB file not found: {template_path}")

    template_messages = _read_template_messages(eccodes, template_path)
    if not template_messages:
        raise ValueError(f"template GRIB file contains no messages: {template_path}")

    records = _prepare_records(
        obj,
        template_messages=template_messages,
        variables=variables,
        output_param_id=param_id,
        output_short_name=short_name,
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
                if normalized_packing_type is not None:
                    eccodes.codes_set(clone, "packingType", normalized_packing_type)
                if normalized_bits_per_value is not None:
                    eccodes.codes_set(clone, "bitsPerValue", normalized_bits_per_value)
                for key, value in message_key_overrides:
                    eccodes.codes_set(clone, key, value)
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


def _normalize_packing_type(packing_type: str | None) -> str | None:
    if packing_type is None:
        return None
    normalized = str(packing_type).strip().lower()
    if not normalized:
        return None
    if normalized in _PACKING_ALIASES:
        return _PACKING_ALIASES[normalized]
    if normalized.startswith("grid_"):
        return normalized
    raise ValueError(
        "packing_type must be None, 'template', 'simple', 'ccsds', or an "
        "ecCodes grid_* packingType"
    )


def _normalize_bits_per_value(bits_per_value: int | None) -> int | None:
    if bits_per_value is None:
        return None
    value = int(bits_per_value)
    if value < 0:
        raise ValueError("bits_per_value must be non-negative")
    return value


def _build_message_key_overrides(
    *,
    edition,
    centre,
    sub_centre,
    generating_process_identifier,
    data_date,
    data_time,
    step_type,
    step_range,
    forecast_time,
    param_id,
    short_name,
    type_of_level,
    level,
    key_overrides,
) -> list[tuple[str, object]]:
    if param_id is not None and short_name is not None:
        raise ValueError("param_id and short_name are mutually exclusive")
    if step_range is not None and forecast_time is not None:
        raise ValueError("step_range and forecast_time are mutually exclusive")

    items: list[tuple[str, object]] = []
    _add_override(items, "edition", _normalize_edition(edition))
    _add_override(items, "centre", centre)
    _add_override(items, "subCentre", _normalize_int("sub_centre", sub_centre))
    _add_override(
        items,
        "generatingProcessIdentifier",
        _normalize_int("generating_process_identifier", generating_process_identifier),
    )
    _add_override(items, "dataDate", _normalize_int_like("data_date", data_date))
    _add_override(items, "dataTime", _normalize_int_like("data_time", data_time))
    _add_override(items, "stepType", step_type)
    _add_override(items, "stepRange", step_range)
    _add_override(items, "forecastTime", _normalize_int("forecast_time", forecast_time))
    _add_override(items, "paramId", _normalize_int("param_id", param_id))
    _add_override(items, "shortName", short_name)
    _add_override(items, "typeOfLevel", type_of_level)
    _add_override(items, "level", level)

    if key_overrides is not None:
        if not isinstance(key_overrides, Mapping):
            raise TypeError("key_overrides must be a mapping of ecCodes key names to values")
        for key, value in key_overrides.items():
            if value is not None:
                items.append((str(key), value))
    return items


def _add_override(items: list[tuple[str, object]], key: str, value) -> None:
    if value is not None:
        items.append((key, value))


def _normalize_edition(edition: int | None) -> int | None:
    if edition is None:
        return None
    value = int(edition)
    if value not in {1, 2}:
        raise ValueError("edition must be 1 or 2")
    return value


def _normalize_int(name: str, value) -> int | None:
    if value is None:
        return None
    return int(value)


def _normalize_int_like(name: str, value) -> int | None:
    if value is None:
        return None
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            raise ValueError(f"{name} must not be empty")
        return int(stripped)
    return int(value)


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
                        param_id=_get_optional_int_key(eccodes, gid, "paramId"),
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


def _get_optional_int_key(eccodes, gid, key: str) -> int | None:
    try:
        value = eccodes.codes_get(gid, key)
    except Exception:
        return None
    if value is None:
        return None
    return int(value)


def _release_next_template_message(eccodes, template_handle) -> None:
    gid = eccodes.codes_grib_new_from_file(template_handle)
    if gid is not None:
        eccodes.codes_release(gid)


def _prepare_records(
    obj,
    *,
    template_messages: list[_TemplateMessage],
    variables,
    output_param_id: int | None,
    output_short_name: str | None,
    strict: bool,
) -> list[_VariableRecord]:
    if isinstance(obj, xr.DataArray):
        if variables is not None:
            raise ValueError("variables is only supported for xarray Dataset input")
        return [
            _prepare_data_array_record(
                obj,
                template_messages=template_messages,
                explicit_param_id=output_param_id,
                explicit_short_name=output_short_name,
                strict=strict,
            )
        ]

    if isinstance(obj, xr.Dataset):
        selected = list(obj.data_vars) if variables is None else [str(name) for name in variables]
        missing = [name for name in selected if name not in obj.data_vars]
        if missing:
            raise KeyError(f"variables not found in dataset: {missing}")
        if (output_param_id is not None or output_short_name is not None) and len(selected) != 1:
            raise ValueError(
                "param_id and short_name can only be used with one selected Dataset variable"
            )
        records = [
            _prepare_data_array_record(
                obj[name],
                template_messages=template_messages,
                explicit_param_id=output_param_id,
                explicit_short_name=output_short_name or name,
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
    explicit_param_id: int | None,
    explicit_short_name: str | None,
    strict: bool,
) -> _VariableRecord:
    short_name, matches = _match_template_messages(
        data_array,
        template_messages=template_messages,
        explicit_param_id=explicit_param_id,
        explicit_short_name=explicit_short_name,
    )
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


def _match_template_messages(
    data_array: xr.DataArray,
    *,
    template_messages: list[_TemplateMessage],
    explicit_param_id: int | None,
    explicit_short_name: str | None,
) -> tuple[str, list[_TemplateMessage]]:
    if explicit_param_id is not None:
        matches = [message for message in template_messages if message.param_id == explicit_param_id]
        if not matches:
            raise ValueError(f"template GRIB contains no messages for paramId={explicit_param_id!r}")
        short_names = {message.short_name for message in matches}
        if len(short_names) != 1:
            raise ValueError(
                f"template messages for paramId={explicit_param_id!r} use mixed shortName values"
            )
        return next(iter(short_names)), matches

    short_name = _infer_short_name(data_array, template_messages, explicit_short_name)
    matches = [message for message in template_messages if message.short_name == short_name]
    return short_name, matches


def _infer_short_name(
    data_array: xr.DataArray,
    template_messages: list[_TemplateMessage],
    explicit_short_name: str | None,
) -> str:
    if explicit_short_name not in (None, ""):
        return str(explicit_short_name)
    for key in ("GRIB_shortName", "shortName"):
        value = data_array.attrs.get(key)
        if value not in (None, ""):
            return str(value)
    if data_array.name not in (None, ""):
        return str(data_array.name)
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
