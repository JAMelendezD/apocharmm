# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from __future__ import annotations

from pathlib import Path
import sys

import apocharmm as apo

from python_api_test_helpers import (
    assert_equal,
    expect_exception,
    expect_invalid_argument,
    get_repo_root,
    remove_if_exists,
)

TIME_STEP: float = 0.001


def check_report_frequency(path: Path) -> None:
    print("Checking Subscriber report-frequency validation...")

    remove_if_exists(path)
    subscriber = apo.RestartSubscriber(path)

    try:
        assert_equal(
            "default Subscriber report frequency", subscriber.getReportFrequency(), 1000
        )

        subscriber.setReportFrequency(17)
        assert_equal(
            "updated Subscriber report frequency", subscriber.getReportFrequency(), 17
        )

        zero_error = expect_invalid_argument(
            "Subscriber rejects zero report frequency",
            lambda: subscriber.setReportFrequency(0),
            "Subscriber report frequency must be positive; observed 0",
        )
        assert_equal(
            "Subscriber zero report frequency context",
            zero_error.context,
            "Subscriber.setReportFrequency(report_frequency)",
        )
        assert_equal(
            "Subscriber report frequency after zero rejection",
            subscriber.getReportFrequency(),
            17,
        )

        negative_error = expect_invalid_argument(
            "Subscriber rejects negative report frequency",
            lambda: subscriber.setReportFrequency(-1),
            "Subscriber report frequency must be positive; observed -1",
        )
        assert_equal(
            "Subscriber negative report frequency context",
            negative_error.context,
            "Subscriber.setReportFrequency(report_frequency)",
        )
        assert_equal(
            "Subscriber report frequency after negative rejection",
            subscriber.getReportFrequency(),
            17,
        )

        expect_exception(
            "Subscriber rejects non-int report frequency",
            TypeError,
            lambda: subscriber.setReportFrequency(1.5),
        )
        expect_exception(
            "Subscriber rejects too-large report frequency",
            ValueError,
            lambda: subscriber.setReportFrequency(2**31),
        )
        expect_exception(
            "Subscriber rejects too-small report frequency",
            ValueError,
            lambda: subscriber.setReportFrequency(-(2**31) - 1),
        )
    finally:
        subscriber.close()
        remove_if_exists(path)

    return


def check_output_open_failure() -> None:
    print("Checking Subscriber output-open failure...")

    error = expect_exception(
        "Subscriber reports output-open failure",
        apo.ApoCharmmError,
        lambda: apo.RestartSubscriber("."),
    )

    assert_equal(
        "Subscriber output-open status", error.status, apo.APO_STATUS_RUNTIME_ERROR
    )
    assert_equal(
        "Subscriber output-open status name",
        error.status_name,
        "APO_STATUS_RUNTIME_ERROR",
    )
    assert_equal(
        "Subscriber output-open context",
        error.context,
        "RestartSubscriber construction",
    )

    expected_diagnostic: str = "Failed to open subscriber output file for writing: ."
    if expected_diagnostic not in error.native_diagnostic:
        raise AssertionError(
            "Subscriber output-open diagnostic: expected "
            f"{expected_diagnostic!r}, observed {error.native_diagnostic!r}"
        )

    if "ERROR:" in error.native_diagnostic:
        raise AssertionError(
            "Subscriber output-open diagnostic contains a manual ERROR: prefix"
        )

    if error.native_diagnostic.endswith("\n"):
        raise AssertionError("Subscriber output-open diagnostic has a trailing newline")

    return


def check_integrator_attachment(path: Path) -> None:
    print("Checking Subscriber integrator attachment errors...")

    remove_if_exists(path)
    subscriber = apo.RestartSubscriber(path)
    integrator = apo.CudaLangevinThermostatIntegrator(TIME_STEP)

    try:
        missing_error = expect_invalid_argument(
            "CudaIntegrator rejects unsubscribing a missing Subscriber",
            lambda: integrator.unsubscribe(subscriber),
            f'Subscriber not found (file "{path}")',
        )
        assert_equal(
            "missing Subscriber unsubscribe context",
            missing_error.context,
            "CudaIntegrator.unsubscribe(subscriber)",
        )

        integrator.subscribe(subscriber)

        duplicate_error = expect_invalid_argument(
            "CudaIntegrator rejects duplicate Subscriber attachment",
            lambda: integrator.subscribe(subscriber),
            "Subscriber is already subscribed to this CudaIntegrator",
        )
        assert_equal(
            "duplicate Subscriber attachment context",
            duplicate_error.context,
            "CudaIntegrator.subscribe(subscriber)",
        )

        integrator.unsubscribe(subscriber)

        second_missing_error = expect_invalid_argument(
            "CudaIntegrator rejects repeated Subscriber removal",
            lambda: integrator.unsubscribe(subscriber),
            f'Subscriber not found (file "{path}")',
        )
        assert_equal(
            "repeated Subscriber removal context",
            second_missing_error.context,
            "CudaIntegrator.unsubscribe(subscriber)",
        )
    finally:
        integrator.close()
        subscriber.close()
        remove_if_exists(path)

    return


def main(argc: int, argv: list[str]) -> int:
    output_dir: Path = get_repo_root() / "test/pytest"
    output_dir.mkdir(parents=True, exist_ok=True)

    check_report_frequency(output_dir / "tmp_python_api_subscriber.rst")
    check_output_open_failure()
    check_integrator_attachment(output_dir / "tmp_python_api_subscriber_attachment.rst")

    print("\033[32m" + "PASS: Subscriber Python API tests completed." + "\033[0m")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
