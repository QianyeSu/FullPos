from __future__ import annotations

import argparse
import json
from pathlib import Path


DEFAULT_REPORT = Path("validation_roundtrip_o320_o480_o320_t_137.json")


def main() -> None:
    args = _parse_args()
    report = json.loads(args.report.read_text(encoding="utf-8"))
    failures = _check_report(
        report,
        max_global_relative_rmse=args.max_global_relative_rmse,
        max_level_relative_rmse=args.max_level_relative_rmse,
        require_finite=not args.allow_nonfinite,
    )

    if failures:
        print("Regression check failed:")
        for failure in failures:
            print(f"  - {failure}")
        raise SystemExit(1)

    print("Regression check passed:")
    print(f"  global relative_rmse={report['relative_rmse']:.6g}")
    if "by_level" in report:
        worst = max(report["by_level"], key=lambda row: row["relative_rmse"])
        print(
            "  worst level "
            f"level={worst['level']} "
            f"relative_rmse={worst['relative_rmse']:.6g} "
            f"rmse={worst['rmse']:.6g} "
            f"max_abs={worst['max_abs']:.6g}"
        )


def _check_report(
    report: dict,
    *,
    max_global_relative_rmse: float,
    max_level_relative_rmse: float,
    require_finite: bool,
) -> list[str]:
    failures = []
    if require_finite and not report.get("finite", False):
        failures.append("global finite flag is false")

    global_relative_rmse = float(report["relative_rmse"])
    if global_relative_rmse > max_global_relative_rmse:
        failures.append(
            f"global relative_rmse {global_relative_rmse:.6g} exceeds "
            f"{max_global_relative_rmse:.6g}"
        )

    by_level = report.get("by_level")
    if by_level is None:
        failures.append("report does not contain by_level metrics")
        return failures

    worst = max(by_level, key=lambda row: row["relative_rmse"])
    worst_relative_rmse = float(worst["relative_rmse"])
    if worst_relative_rmse > max_level_relative_rmse:
        failures.append(
            f"level {worst['level']} relative_rmse {worst_relative_rmse:.6g} exceeds "
            f"{max_level_relative_rmse:.6g}"
        )
    if require_finite:
        bad_levels = [row["level"] for row in by_level if not row.get("finite", False)]
        if bad_levels:
            failures.append(f"non-finite levels: {bad_levels}")
    return failures


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check saved FullPos validation metrics.")
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    parser.add_argument("--max-global-relative-rmse", type=float, default=0.01)
    parser.add_argument("--max-level-relative-rmse", type=float, default=0.03)
    parser.add_argument("--allow-nonfinite", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    main()
