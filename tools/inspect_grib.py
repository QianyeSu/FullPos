from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path

import eccodes
import numpy as np


DEFAULT_PATH = Path(
    r"L:\ERA5_test\era5_reanalysis_model_level_20250102_packing_CCSDS_O320.grib2"
)


def inspect(path: Path = DEFAULT_PATH) -> None:
    keys = [
        "shortName",
        "typeOfLevel",
        "level",
        "dataDate",
        "dataTime",
        "gridType",
        "N",
        "Ny",
        "numberOfPoints",
        "packingType",
    ]
    groups: Counter[tuple[object, ...]] = Counter()
    levels: defaultdict[tuple[object, object], set[object]] = defaultdict(set)
    first_pl = None
    count = 0

    with path.open("rb") as handle:
        while True:
            gid = eccodes.codes_grib_new_from_file(handle)
            if gid is None:
                break
            count += 1
            row = {}
            for key in keys:
                try:
                    row[key] = eccodes.codes_get(gid, key)
                except Exception:
                    row[key] = None
            groups[
                (
                    row["shortName"],
                    row["typeOfLevel"],
                    row["gridType"],
                    row["packingType"],
                )
            ] += 1
            levels[(row["shortName"], row["typeOfLevel"])].add(row["level"])
            if first_pl is None:
                first_pl = np.asarray(eccodes.codes_get_array(gid, "pl"), dtype=np.int64)
            eccodes.codes_release(gid)

    print(f"path: {path}")
    print(f"messages: {count}")
    print("groups:")
    for key, value in groups.most_common():
        print(f"  {value}: {key}")
    print("levels:")
    for key, values in sorted(levels.items()):
        values_sorted = sorted(values)
        print(f"  {key}: {len(values_sorted)} levels, {values_sorted[:3]} ... {values_sorted[-3:]}")
    if first_pl is not None:
        print("first message pl:")
        print(f"  rows: {first_pl.size}")
        print(f"  points: {int(first_pl.sum())}")
        print(f"  min/max row length: {int(first_pl.min())}/{int(first_pl.max())}")
        print(f"  symmetric: {bool(np.array_equal(first_pl, first_pl[::-1]))}")
        print(f"  first 20: {first_pl[:20].tolist()}")
        print(f"  last 20: {first_pl[-20:].tolist()}")


if __name__ == "__main__":
    inspect()
