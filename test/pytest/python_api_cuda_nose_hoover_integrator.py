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
    require_file,
    assert_close,
    assert_nested_sequence_close,
    assert_finite_temperature,
    expect_exception,
)

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
RANDOM_SEED: int = 314159
TEMPERATURE: float = 300.0
TIME_STEP: float = 0.001
REFERENCE_TEMPERATURE: float = 315.0
NOSE_HOOVER_PISTON_MASS: float = 500.0
NUM_STEPS: int = 5
TOLERANCE: float = 1.0e-10
DETERMINISTIC_TOLERANCE: float = 0.0


def create_context(repo_root: Path) -> apo.CharmmContext:
    prm_path: str = require_file(repo_root / "test/data/toppar_water_ions.str")
    psf_path: str = require_file(repo_root / "test/data/nacl_pair.psf")
    crd_path: str = require_file(repo_root / "test/data/nacl_pair.cor")

    prm = apo.CharmmParameters(prm_path)
    psf = apo.CharmmPsf(psf_path)
    crd = apo.CharmmCrd(crd_path)

    fm = apo.ForceManager(psf, prm)
    fm.setBoxDimensions(BOX_DIMENSIONS)

    ctx = apo.CharmmContext(fm)
    ctx.setCoordinates(crd)
    ctx.useHolonomicConstraints(False)
    ctx.setRandomSeedForVelocities(RANDOM_SEED)
    ctx.assignVelocitiesAtTemperature(TEMPERATURE)

    return ctx


def create_integrator() -> apo.CudaNoseHooverIntegrator:
    integrator = apo.CudaNoseHooverIntegrator(TIME_STEP)
    integrator.setTimeStep(TIME_STEP)
    integrator.setReferenceTemperature(REFERENCE_TEMPERATURE)
    integrator.setNoseHooverPistonMass(NOSE_HOOVER_PISTON_MASS)
    integrator.useOldTemperature(False)
    return integrator


def check_setters_and_getters() -> None:
    print("Checking CudaNoseHooverIntegrator setters/getters...")

    integrator = create_integrator()

    assert_close(
        "CudaNoseHooverIntegrator.getReferenceTemperature",
        integrator.getReferenceTemperature(),
        REFERENCE_TEMPERATURE,
        TOLERANCE,
    )
    assert_close(
        "CudaNoseHooverIntegrator.getNoseHooverPistonMass",
        integrator.getNoseHooverPistonMass(),
        NOSE_HOOVER_PISTON_MASS,
        TOLERANCE,
    )

    integrator.resetAverageTemperature()
    assert_close(
        "CudaNoseHooverIntegrator.getAverageTemperature after reset",
        integrator.getAverageTemperature(),
        0.0,
        TOLERANCE,
    )

    integrator.close()

    return


def check_validation(repo_root: Path) -> None:
    print("Checking CudaNoseHooverIntegrator validation...")

    integrator = create_integrator()

    expect_exception(
        "CudaNoseHooverIntegrator.setCharmmContext rejects non-CharmmContext",
        TypeError,
        lambda: integrator.setCharmmContext(object()),
    )
    expect_exception(
        "CudaNoseHooverIntegrator.propagate rejects negative step count",
        ValueError,
        lambda: integrator.propagate(-1),
    )
    expect_exception(
        "CudaNoseHooverIntegrator.propagate rejects missing CharmmContext",
        apo.ApoCharmmError,
        lambda: integrator.propagate(1),
    )
    expect_exception(
        "CudaNoseHooverIntegrator.initializeFromRestartFile rejects missing CharmmContext",
        apo.ApoCharmmError,
        lambda: integrator.initializeFromRestartFile("missing.rst"),
    )

    ctx = create_context(repo_root)
    integrator.setCharmmContext(ctx)

    expect_exception(
        "CudaNoseHooverIntegrator.setCharmmContext rejects second context",
        apo.ApoCharmmError,
        lambda: integrator.setCharmmContext(ctx),
    )
    expect_exception(
        "CudaNoseHooverIntegrator.initializeFromRestartFile rejects missing restart file",
        apo.ApoCharmmError,
        lambda: integrator.initializeFromRestartFile("missing.rst"),
    )

    integrator.close()

    return


def check_short_propagation(repo_root: Path) -> None:
    print("Checking CudaNoseHooverIntegrator short propagation...")

    ctx = create_context(repo_root)
    integrator = create_integrator()

    integrator.setCharmmContext(ctx)
    integrator.resetAverageTemperature()

    initial_coordinates_charges: list[list[float]] = ctx.getCoordinatesCharges()

    integrator.propagate(NUM_STEPS)

    final_coordinates_charges: list[list[float]] = ctx.getCoordinatesCharges()

    assert_finite_temperature(
        "CharmmContext.computeTemperature", ctx.computeTemperature()
    )
    assert_finite_temperature(
        "CudaNoseHooverIntegrator.getInstantaneousTemperature",
        integrator.getInstantaneousTemperature(),
    )
    assert_finite_temperature(
        "CudaNoseHooverIntegrator.getAverageTemperature",
        integrator.getAverageTemperature(),
    )

    if final_coordinates_charges[0][:3] == initial_coordinates_charges[0][:3]:
        raise AssertionError("short propagation did not move atom 0 coordinates")

    integrator.close()

    return


def check_deterministic_trajectories(repo_root: Path) -> None:
    print("Checking deterministic Nose-Hoover trajectories...")

    ctx1 = create_context(repo_root)
    ctx2 = create_context(repo_root)

    assert_nested_sequence_close(
        "initial coordinates/charges",
        ctx1.getCoordinatesCharges(),
        ctx2.getCoordinatesCharges(),
        DETERMINISTIC_TOLERANCE,
    )

    integrator1 = create_integrator()
    integrator2 = create_integrator()

    integrator1.setCharmmContext(ctx1)
    integrator2.setCharmmContext(ctx2)

    assert_close(
        "deterministic initial reference temperature",
        integrator1.getReferenceTemperature(),
        integrator2.getReferenceTemperature(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_close(
        "deterministic initial piston mass",
        integrator1.getNoseHooverPistonMass(),
        integrator2.getNoseHooverPistonMass(),
        DETERMINISTIC_TOLERANCE,
    )

    initial_coordinates_charges1: list[list[float]] = ctx1.getCoordinatesCharges()
    initial_coordinates_charges2: list[list[float]] = ctx2.getCoordinatesCharges()

    integrator1.propagate(NUM_STEPS)
    integrator2.propagate(NUM_STEPS)

    final_coordinates_charges1: list[list[float]] = ctx1.getCoordinatesCharges()
    final_coordinates_charges2: list[list[float]] = ctx2.getCoordinatesCharges()

    assert_nested_sequence_close(
        "deterministic final coordinates/charges",
        final_coordinates_charges1,
        final_coordinates_charges2,
        DETERMINISTIC_TOLERANCE,
    )
    assert_close(
        "deterministic final average temperature",
        integrator1.getAverageTemperature(),
        integrator2.getAverageTemperature(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_close(
        "deterministic final instantaneous temperature",
        integrator1.getInstantaneousTemperature(),
        integrator2.getInstantaneousTemperature(),
        DETERMINISTIC_TOLERANCE,
    )

    if final_coordinates_charges1[0][:3] == initial_coordinates_charges1[0][:3]:
        raise AssertionError("deterministic trajectory 1 did not move atom 0")
    if final_coordinates_charges2[0][:3] == initial_coordinates_charges2[0][:3]:
        raise AssertionError("deterministic trajectory 2 did not move atom 0")

    integrator1.close()
    integrator2.close()

    return


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")

    check_setters_and_getters()
    check_validation(repo_root)
    check_short_propagation(repo_root)
    check_deterministic_trajectories(repo_root)

    print(
        "\033[32m"
        + "PASS: CudaNoseHooverIntegrator Python API tests completed."
        + "\033[0m"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
