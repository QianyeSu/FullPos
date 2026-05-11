from tools.check_regression import _check_report


def test_check_report_accepts_current_thresholds() -> None:
    report = {
        "finite": True,
        "relative_rmse": 0.004,
        "by_level": [
            {"level": 1, "finite": True, "relative_rmse": 0.001},
            {"level": 137, "finite": True, "relative_rmse": 0.024},
        ],
    }
    failures = _check_report(
        report,
        max_global_relative_rmse=0.01,
        max_level_relative_rmse=0.03,
        require_finite=True,
    )
    assert failures == []


def test_check_report_rejects_bad_level() -> None:
    report = {
        "finite": True,
        "relative_rmse": 0.004,
        "by_level": [
            {"level": 137, "finite": True, "relative_rmse": 0.04},
        ],
    }
    failures = _check_report(
        report,
        max_global_relative_rmse=0.01,
        max_level_relative_rmse=0.03,
        require_finite=True,
    )
    assert failures
    assert "level 137" in failures[0]
