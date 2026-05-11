from pathlib import Path

from fullpos.native import native_prefix, native_runtime_dir


def test_native_prefix_uses_environment_override(monkeypatch) -> None:
    monkeypatch.setenv("FULLPOS_NATIVE_PREFIX", r"Z:\custom-fullpos")
    assert native_prefix() == Path(r"Z:\custom-fullpos")
    assert str(native_runtime_dir()).startswith(r"Z:\custom-fullpos")
