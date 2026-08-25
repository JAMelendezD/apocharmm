# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from __future__ import annotations

import sys

import apocharmm as apo

import apo_test_helpers as apo_test

TIME_STEP: float = 0.002


def check_base_cuda_integrator_has_no_valid_handle() -> None:
    print("Checking base CudaIntegrator invalid-handle behavior...")

    integrator = apo.CudaIntegrator()

    apo_test.expect_exception(
        "base CudaIntegrator.integrator_handle rejects an empty handle",
        RuntimeError,
        lambda: integrator.integrator_handle,
    )
    apo_test.expect_exception(
        "base CudaIntegrator.setTimeStep rejects an empty handle",
        RuntimeError,
        lambda: integrator.setTimeStep(TIME_STEP),
    )
    apo_test.expect_exception(
        "base CudaIntegrator.propagate rejects an empty handle",
        RuntimeError,
        lambda: integrator.propagate(1),
    )
    apo_test.expect_exception(
        "base CudaIntegrator.initializeFromRestartFile rejects an empty handle",
        RuntimeError,
        lambda: integrator.initializeFromRestartFile("restart.rst"),
    )

    integrator.close()

    apo_test.expect_exception(
        "closed base CudaIntegrator still rejects an empty handle",
        RuntimeError,
        lambda: integrator.setTimeStep(TIME_STEP),
    )

    return


def check_python_side_validation() -> None:
    print("Checking base CudaIntegrator Python-side validation...")

    integrator = apo.CudaIntegrator()

    apo_test.expect_exception(
        "CudaIntegrator.setCharmmContext rejects non-CharmmContext",
        TypeError,
        lambda: integrator.setCharmmContext(object()),
    )
    apo_test.expect_exception(
        "CudaIntegrator.subscribe rejects non-Subscriber",
        TypeError,
        lambda: integrator.subscribe(object()),
    )
    apo_test.expect_exception(
        "CudaIntegrator.propagate rejects step count below C int range",
        ValueError,
        lambda: integrator.propagate(-(2**31) - 1),
    )
    apo_test.expect_exception(
        "CudaIntegrator.propagate rejects step count outside C int range",
        ValueError,
        lambda: integrator.propagate(2**31),
    )

    return


def main(argc: int, argv: list[str]) -> int:
    check_base_cuda_integrator_has_no_valid_handle()
    check_python_side_validation()

    print(
        "\033[32m"
        + "PASS: base CudaIntegrator Python API validation tests completed."
        + "\033[0m"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
