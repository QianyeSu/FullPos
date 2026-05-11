from __future__ import annotations

import argparse
import json

from fullpos import doctor


def main() -> None:
    args = _parse_args()
    report = doctor()
    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
    else:
        _print_report(report)
    raise SystemExit(0 if report["ok"] else 1)


def _print_report(report: dict) -> None:
    info = report["backend_info"]
    print("FullPos doctor")
    print(f"  backend: {info['backend']}")
    print(f"  python: {info['python']}")
    print(f"  platform: {info['platform']}")
    print(f"  numpy: {info['numpy']}")
    print(f"  native_module: {info['native_module']}")
    print(f"  native_prefix: {info['native_prefix']}")
    print(f"  native_runtime_platform: {info['native_runtime_platform']}")
    print(f"  native_runtime_dir: {info['native_runtime_dir']}")
    if info["external_runtime_libraries"]:
        print("  external_runtime_libraries:")
        for name, path in info["external_runtime_libraries"].items():
            print(f"    {name}: {path}")
    print("checks:")
    for check in report["checks"]:
        status = "PASS" if check["ok"] else "FAIL"
        print(f"  [{status}] {check['name']}: {check['detail']}")


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check FullPos native backend runtime health.")
    parser.add_argument("--json", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    main()
