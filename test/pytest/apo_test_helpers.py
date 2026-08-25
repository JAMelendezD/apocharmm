# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from __future__ import annotations

from collections.abc import Callable, Sequence
import math
import os
from pathlib import Path
from typing import TypeVar

import apocharmm as apo

T = TypeVar("T")
ExceptionT = TypeVar("ExceptionT", bound=BaseException)

_APO_STATUS_NAMES: dict[int, str] = {
    apo.APO_STATUS_OK: "APO_STATUS_OK",
    apo.APO_STATUS_INVALID_ARGUMENT: "APO_STATUS_INVALID_ARGUMENT",
    apo.APO_STATUS_RUNTIME_ERROR: "APO_STATUS_RUNTIME_ERROR",
    apo.APO_STATUS_CUDA_ERROR: "APO_STATUS_CUDA_ERROR",
    apo.APO_STATUS_NOT_INITIALIZED: "APO_STATUS_NOT_INITIALIZED",
    apo.APO_STATUS_NOT_IMPLEMENTED: "APO_STATUS_NOT_IMPLEMENTED",
}


def get_repo_root() -> Path:
    repo_root: str | None = os.environ.get("APOCHARMM_REPO_ROOT")
    if repo_root is not None:
        return Path(repo_root)

    return Path(__file__).resolve().parents[2]


def get_data_dir() -> Path:
    data_path: str | None = os.environ.get("APOCHARMM_TEST_DATA_DIR")
    if data_path is not None:
        return Path(data_path)

    return get_repo_root() / "test/data"


def get_toppar_dir() -> Path:
    toppar_path: str | None = os.environ.get("APOCHARMM_TOPPAR_DIR")
    if toppar_path is not None:
        return Path(toppar_path)

    return get_repo_root() / "toppar"


def write_text_file(path: Path, contents: str) -> None:
    path.write_text(contents, encoding="utf-8")
    return


def require_file(path: Path) -> str:
    if not path.is_file():
        raise FileNotFoundError(f"required test file does not exist: {path}")

    return str(path)


def assert_file_created(path: Path) -> None:
    if not path.is_file():
        raise AssertionError(f"expected subscriber output file does not exist: {path}")

    if path.stat().st_size <= 0:
        raise AssertionError(f"expected subscriber output file is empty: {path}")

    return


def remove_if_exists(path: Path) -> None:
    if path.exists():
        path.unlink()
    return


def assert_equal(label: str, observed: T, expected: T) -> None:
    if observed != expected:
        raise AssertionError(f"{label}: expected {expected!r}, observed {observed!r}")

    return


def assert_close(
    label: str, observed: float, expected: float, tolerance: float
) -> None:
    if not math.isfinite(observed):
        raise AssertionError(f"{label}: observed value is not finite: {observed}")

    if abs(observed - expected) > tolerance:
        raise AssertionError(
            f"{label}: expected {expected:.12g}, obersvered {observed:.12g}, tolerance {tolerance:.3g}"
        )

    return


def assert_sequence_close(
    label: str, observed: Sequence[float], expected: Sequence[float], tolerance: float
) -> None:
    assert_equal(f"{label} length", len(observed), len(expected))

    for i in range(len(observed)):
        assert_close(f"{label}[{i}]", float(observed[i]), float(expected[i]), tolerance)

    return


def assert_nested_sequence_close(
    label: str,
    observed: Sequence[Sequence[float]],
    expected: Sequence[Sequence[float]],
    tolerance: float,
) -> None:
    assert_equal(f"{label} outer length", len(observed), len(expected))

    for i in range(len(expected)):
        assert_sequence_close(f"{label}[{i}]", observed[i], expected[i], tolerance)

    return


def assert_finite_scalar(label: str, value: float) -> None:
    if not math.isfinite(value):
        raise AssertionError(f"{label} is not finite: {value}")

    return


def assert_finite_sequence(label: str, values: Sequence[float]) -> None:
    for i, value in enumerate(values):
        assert_finite_scalar(f"{label}[{i}]", float(value))

    return


def assert_finite_nested_sequence(
    label: str, values_values: Sequence[Sequence[float]]
) -> None:
    for i, values in enumerate(values_values):
        assert_finite_sequence(f"{label}[{i}]", values)


def assert_finite_nonnegative_scalar(label: str, value: float) -> None:
    assert_finite_scalar(label, value)

    if value < 0.0:
        raise AssertionError(f"{label} is negative: {value}")

    return


def assert_finite_temperature(label: str, temperature: float) -> None:
    assert_finite_nonnegative_scalar(f"{label} temperature", temperature)
    return


def expect_exception(
    label: str, exception_type: type[ExceptionT], action: Callable[[], object]
) -> ExceptionT:
    try:
        action()
    except exception_type as exc:
        print(f"Expected exception observed: {label}")
        return exc
    except Exception as exc:
        raise AssertionError(
            f"{label}: expected {exception_type.__name__}, observed {type(exc).__name__}: {exc}"
        ) from exc

    raise AssertionError(
        f"{label}: expected {exception_type.__name__}, but no exception was raised"
    )


def expect_apo_error(
    label: str,
    action: Callable[[], object],
    expected_status: int,
    diagnostic_substring: str,
    expected_context: str | None = None,
) -> apo.ApoCharmmError:
    error: apo.ApoCharmmError = expect_exception(label, apo.ApoCharmmError, action)

    assert_equal(f"{label} status", error.status, expected_status)

    expected_status_name: str | None = _APO_STATUS_NAMES.get(expected_status)
    if expected_status_name is None:
        raise AssertionError(
            f"{label}: no expected status name is defined for apo_status {expected_status}"
        )

    assert_equal(f"{label} status name", error.status_name, expected_status_name)

    if expected_context is not None:
        assert_equal(f"{label} operation context", error.context, expected_context)
        assert_equal(
            f"{label} rendered context occurrence count",
            error.message.count(expected_context),
            1,
        )

    if diagnostic_substring not in error.native_diagnostic:
        raise AssertionError(
            f"{label}: expected native diagnostic to contain {diagnostic_substring!r}, observed {error.native_diagnostic!r}"
        )

    if "ERROR:" in error.native_diagnostic:
        raise AssertionError(
            f"{label}: native diagnostic contains a manual ERROR: prefix"
        )

    if error.native_diagnostic.endswith("\n"):
        raise AssertionError(f"{label}: native diagnostic has a trailing newline")

    return error


def expect_invalid_argument(
    label: str,
    action: Callable[[], object],
    diagnostic_substring: str,
    expected_context: str | None = None,
) -> apo.ApoCharmmError:
    return expect_apo_error(
        label,
        action,
        apo.APO_STATUS_INVALID_ARGUMENT,
        diagnostic_substring,
        expected_context=expected_context,
    )
