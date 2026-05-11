from fullpos import backend_info, doctor


def test_backend_info_reports_native_status() -> None:
    info = backend_info()
    assert "python" in info
    assert "platform" in info
    assert "numpy" in info
    assert info["backend"] == "native"
    assert "native_available" in info
    assert info["native_module"] == "fullpos._ectrans"
    assert "native_runtime_platform" in info
    assert "native_prefix" in info
    assert "native_prefix_exists" in info
    assert "native_runtime_dir" in info
    assert "required_native_libraries" in info
    assert "required_native_libraries_present" in info
    assert "external_runtime_libraries" in info
    assert "external_runtime_libraries_present" in info


def test_doctor_runs_native_smoke_check() -> None:
    report = doctor()
    assert "checks" in report
    names = {check["name"] for check in report["checks"]}
    assert "native prefix exists" in names
    assert "native module importable" in names
    assert "required native libraries present" in names
    assert "external runtime libraries present" in names
    assert "native smoke regrid N4 -> N4" in names
    assert report["ok"]
