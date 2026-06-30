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
    get_data_path,
    require_file,
    assert_equal,
    assert_close,
    assert_finite_scalar,
    assert_finite_temperature,
    assert_sequence_close,
    assert_finite_sequence,
    assert_finite_nested_sequence,
    assert_nested_sequence_close,
    expect_exception,
)

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
RANDOM_SEED: int = 314159
TEMPERATURE: float = 300.0
TIME_STEP: float = 0.001
REFERENCE_TEMPERATURE: float = 315.0
NOSE_HOOVER_PISTON_MASS: float = 500.0
LANGEVIN_PISTON_MASS: list[float] = [5.0]
LANGEVIN_PISTON_FRICTION: float = 5.0
LANGEVIN_PISTON_SEED: int = 987654321
REFERENCE_PRESSURE_TENSOR: list[list[float]] = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
]
NUM_STEPS: int = 3
TOLERANCE: float = 1.0e-10
DETERMINISTIC_TOLERANCE: float = 1.0e-14


def create_context() -> apo.CharmmContext:
    prm_path: str = require_file(get_data_path() / "toppar_water_ions.str")
    psf_path: str = require_file(get_data_path() / "nacl_pair.psf")
    crd_path: str = require_file(get_data_path() / "nacl_pair.cor")

    prm = apo.CharmmParameters(prm_path)
    psf = apo.CharmmPsf(psf_path)
    crd = apo.CharmmCrd(crd_path)

    ctx = apo.CharmmContext(psf, prm)
    ctx.setBoxDimensions(BOX_DIMENSIONS)
    ctx.setCoordinates(crd)
    ctx.useHolonomicConstraints(False)
    ctx.setRandomSeedForVelocities(RANDOM_SEED)
    ctx.assignVelocitiesAtTemperature(TEMPERATURE)

    return ctx


def create_integrator() -> apo.CudaLangevinPistonIntegrator:
    integrator = apo.CudaLangevinPistonIntegrator(TIME_STEP)
    integrator.useNoseHooverThermostat(True)
    integrator.setReferenceTemperature(REFERENCE_TEMPERATURE)
    integrator.setNoseHooverPistonMass(NOSE_HOOVER_PISTON_MASS)
    integrator.setReferencePressure(REFERENCE_PRESSURE_TENSOR)
    integrator.setCrystalType(apo.CrystalType.CUBIC)
    integrator.setLangevinPistonMass(LANGEVIN_PISTON_MASS)
    integrator.setLangevinPistonFrictionSeed(LANGEVIN_PISTON_SEED)
    integrator.setLangevinPistonFriction(LANGEVIN_PISTON_FRICTION)
    return integrator


def check_setters_and_getters() -> None:
    print("Checking CudaLangevinPistonIntegrator setters/getters...")

    integrator = create_integrator()

    assert_close(
        "CudaLangevinPistonIntegrator.getReferenceTemperature",
        integrator.getReferenceTemperature(),
        REFERENCE_TEMPERATURE,
        TOLERANCE,
    )
    assert_close(
        "CudaLangevinPistonIntegrator.getNoseHooverPistonMass",
        integrator.getNoseHooverPistonMass(),
        NOSE_HOOVER_PISTON_MASS,
        TOLERANCE,
    )
    assert_nested_sequence_close(
        "CudaLangevinPistonIntegrator.getReferencePressureTensor",
        integrator.getReferencePressureTensor(),
        REFERENCE_PRESSURE_TENSOR,
        TOLERANCE,
    )
    assert_equal(
        "CudaLangevinPistonIntegrator.getCrystalType",
        integrator.getCrystalType(),
        apo.CrystalType.CUBIC,
    )
    assert_sequence_close(
        "CudaLangevinPistonIntegrator.getLangevinPistonMass",
        integrator.getLangevinPistonMass(),
        LANGEVIN_PISTON_MASS,
        TOLERANCE,
    )
    assert_equal(
        "CudaLangevinPistonIntegrator.getLangevinPistonFrictionSeed",
        integrator.getLangevinPistonFrictionSeed(),
        LANGEVIN_PISTON_SEED,
    )

    integrator.resetAverages()

    assert_close(
        "CudaLangevinPistonIntegrator.getAverageTemperature after reset",
        integrator.getAverageTemperature(),
        0.0,
        TOLERANCE,
    )
    assert_nested_sequence_close(
        "CudaLangevinPistonIntegrator.getAveragePressureTensor after reset",
        integrator.getAveragePressureTensor(),
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
        TOLERANCE,
    )
    assert_close(
        "CudaLangevinPistonIntegrator.getAveragePressureScalar after reset",
        integrator.getAveragePressureScalar(),
        0.0,
        TOLERANCE,
    )

    integrator.close()

    return


def check_validation() -> None:
    print("Checking CudaLagnevinPistonIntegrator validation...")

    integrator = apo.CudaLangevinPistonIntegrator(TIME_STEP)

    expect_exception(
        "setReferencePressure rejects wrong row length",
        ValueError,
        lambda: integrator.setReferencePressure(
            [[1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        ),
    )
    expect_exception(
        "setReferencePressure rejects wrong flattened length",
        apo.ApoCharmmError,
        lambda: integrator.setReferencePressure([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
    )
    expect_exception(
        "setReferencePressure rejects non-finite entries",
        apo.ApoCharmmError,
        lambda: integrator.setReferencePressure(
            [[1.0, 0.0, 0.0], [0.0, float("inf"), 0.0], [0.0, 0.0, 1.0]]
        ),
    )
    expect_exception(
        "setCrystalType rejects invalid crystal type",
        ValueError,
        lambda: integrator.setCrystalType(999),
    )
    expect_exception(
        "setLangevinPistonMass rejects missing crystal type",
        apo.ApoCharmmError,
        lambda: integrator.setLangevinPistonMass([1.0]),
    )
    expect_exception(
        "setLangevinPistonFriction rejects missing crystal type",
        apo.ApoCharmmError,
        lambda: integrator.setLangevinPistonFriction(1.0),
    )
    expect_exception(
        "setLangevinPistonFrictionSeed rejects negative seed",
        ValueError,
        lambda: integrator.setLangevinPistonFrictionSeed(-1),
    )

    integrator.setCrystalType(apo.CrystalType.CUBIC)

    expect_exception(
        "setLangevinPistonMass rejects wrong length for cubic crystal",
        apo.ApoCharmmError,
        lambda: integrator.setLangevinPistonMass([1.0, 2.0]),
    )
    expect_exception(
        "setLangevinPistonMass rejects negative mass",
        apo.ApoCharmmError,
        lambda: integrator.setLangevinPistonMass([-1.0]),
    )

    missing_crystal_integrator = apo.CudaLangevinPistonIntegrator(TIME_STEP)
    ctx = create_context()
    expect_exception(
        "setCharmmContext rejects missing crystal type",
        apo.ApoCharmmError,
        lambda: missing_crystal_integrator.setCharmmContext(ctx),
    )

    integrator.close()
    missing_crystal_integrator.close()

    return


def check_short_propagation() -> None:
    print("Checking CudaLangevinPistonIntegrator short propagation...")

    ctx = create_context()
    integrator = create_integrator()

    integrator.setCharmmContext(ctx)
    integrator.resetAverages()

    initial_coordinates_charges: list[list[float]] = ctx.getCoordinatesCharges()

    integrator.propagate(NUM_STEPS)

    final_coordinates_charges: list[list[float]] = ctx.getCoordinatesCharges()

    assert_finite_temperature(
        "CharmmContext.computeTemperature", ctx.computeTemperature()
    )
    assert_finite_temperature(
        "CudaLangevinPistonIntegrator.getInstantaneousTemperature",
        integrator.getInstantaneousTemperature(),
    )
    assert_finite_temperature(
        "CudaLangevinPistonIntegrator.getAverageTemperature",
        integrator.getAverageTemperature(),
    )

    assert_finite_nested_sequence(
        "CudaLangevinPistonIntegrator.getInstantaneousPressureTensor",
        integrator.getInstantaneousPressureTensor(),
    )
    assert_finite_scalar(
        "CudaLangevinPistonIntegrator.getInstantaneousPressureScalar",
        integrator.getInstantaneousPressureScalar(),
    )
    assert_finite_nested_sequence(
        "CudaLangevinPistonIntegrator.getAveragePressureTensor",
        integrator.getAveragePressureTensor(),
    )
    assert_finite_scalar(
        "CudaLangevinPistonIntegrator.getAveragePressureScalar",
        integrator.getAveragePressureScalar(),
    )

    if final_coordinates_charges[0][:3] == initial_coordinates_charges[0][:3]:
        raise AssertionError("short propagation did not move atom 0 coordinates")

    integrator.close()

    return


def check_deterministic_trajectories() -> None:
    print("Checking deterministic Langevin piston trajectories...")

    ctx1 = create_context()
    ctx2 = create_context()

    assert_nested_sequence_close(
        "initial coordinates/charges",
        ctx1.getCoordinatesCharges(),
        ctx2.getCoordinatesCharges(),
        DETERMINISTIC_TOLERANCE,
    )

    integrator1 = create_integrator()
    integrator1.setCharmmContext(ctx1)

    integrator2 = create_integrator()
    integrator2.setCharmmContext(ctx2)

    assert_close(
        "deterministic initial reference temperature",
        integrator1.getReferenceTemperature(),
        integrator2.getReferenceTemperature(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_nested_sequence_close(
        "deterministic initial reference pressure tensor",
        integrator1.getReferencePressureTensor(),
        integrator2.getReferencePressureTensor(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_sequence_close(
        "deterministic initial Langevin piston mass",
        integrator1.getLangevinPistonMass(),
        integrator2.getLangevinPistonMass(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_equal(
        "deterministic initial Langevin piston seed",
        integrator1.getLangevinPistonFrictionSeed(),
        integrator2.getLangevinPistonFrictionSeed(),
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
        "deterministic final instantaneous temperature",
        integrator1.getInstantaneousTemperature(),
        integrator2.getInstantaneousTemperature(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_close(
        "deterministic final average temperature",
        integrator1.getAverageTemperature(),
        integrator2.getAverageTemperature(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_nested_sequence_close(
        "determinstic final instantaneous pressure tensor",
        integrator1.getInstantaneousPressureTensor(),
        integrator2.getInstantaneousPressureTensor(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_close(
        "determinstic final instantaneous pressure scalar",
        integrator1.getInstantaneousPressureScalar(),
        integrator2.getInstantaneousPressureScalar(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_nested_sequence_close(
        "determinstic final average pressure tensor",
        integrator1.getAveragePressureTensor(),
        integrator2.getAveragePressureTensor(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_close(
        "determinstic final average pressure scalar",
        integrator1.getAveragePressureScalar(),
        integrator2.getAveragePressureScalar(),
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
    check_setters_and_getters()
    check_validation()
    check_short_propagation()
    check_deterministic_trajectories()

    print(
        "\033[32m"
        + "PASS: CudaLangevinPistonIntegrator Python API tests completed."
        + "\033[0m"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
