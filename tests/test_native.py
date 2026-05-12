from pathlib import Path

import fullpos.native as native_module
from fullpos.native import native_library_names_for_platform, native_prefix, native_runtime_dir


def test_native_prefix_uses_environment_override(monkeypatch) -> None:
    monkeypatch.setenv("FULLPOS_NATIVE_PREFIX", r"Z:\custom-fullpos")
    assert native_prefix() == Path(r"Z:\custom-fullpos")
    assert str(native_runtime_dir()).startswith(r"Z:\custom-fullpos")


def test_native_prefix_prefers_bundled_runtime(monkeypatch, tmp_path) -> None:
    package_dir = tmp_path / "fullpos"
    runtime_dir = package_dir / "_native" / ("bin" if native_module.sys.platform == "win32" else "lib")
    runtime_dir.mkdir(parents=True)
    for name in native_library_names_for_platform():
        (runtime_dir / name).write_text("", encoding="utf-8")

    monkeypatch.delenv("FULLPOS_NATIVE_PREFIX", raising=False)
    monkeypatch.setattr(native_module, "__file__", str(package_dir / "native.py"))

    assert native_prefix() == package_dir / "_native"
    assert native_runtime_dir() == runtime_dir
