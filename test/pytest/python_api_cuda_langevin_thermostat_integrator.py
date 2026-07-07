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
    assert_sequence_close,
    assert_nested_sequence_close,
    assert_finite_temperature,
    expect_exception,
)

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
RANDOM_SEED: int = 314159
TEMPERATURE: float = 300.0
TIME_STEP: float = 0.001
REFERENCE_TEMPERATURE: float = 315.0
THERMOSTAT_FRICTION: float = 1.0
THERMOSTAT_SEED: int = 987654321
NUM_STEPS: int = 5
TOLERANCE: float = 1.0e-10
DETERMINISTIC_TOLERANCE: float = 0.0


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
    ctx.setRandomSeed(RANDOM_SEED)
    ctx.assignVelocitiesAtTemperature(TEMPERATURE)

    return ctx


def create_integrator() -> apo.CudaLangevinThermostatIntegrator:
    integrator = apo.CudaLangevinThermostatIntegrator(TIME_STEP)

    integrator.setTimeStep(TIME_STEP)
    integrator.setReferenceTemperature(REFERENCE_TEMPERATURE)
    integrator.setThermostatFriction(THERMOSTAT_FRICTION)
    integrator.setThermostatRngSeed(THERMOSTAT_SEED)

    return integrator


def check_setters_and_getters() -> None:
    print("Checking CudaLangevinThermostatIntegrator setters/getters...")

    integrator = create_integrator()

    assert_close(
        "CudaLangevinThermostatIntegrator.getReferenceTemperature",
        integrator.getReferenceTemperature(),
        REFERENCE_TEMPERATURE,
        TOLERANCE,
    )
    assert_close(
        "CudaLangevinThermostatIntegrator.getThermostatFriction",
        integrator.getThermostatFriction(),
        THERMOSTAT_FRICTION,
        TOLERANCE,
    )
    assert_close(
        "CudaLangevinThermostatIntegrator.getThermostatRngSeed",
        integrator.getThermostatRngSeed(),
        THERMOSTAT_SEED,
        TOLERANCE,
    )

    integrator.resetAverageTemperature()
    assert_close(
        "CudaLangevinThermostatIntegrator.getAverageTemperature after reset",
        integrator.getAverageTemperature(),
        0.0,
        TOLERANCE,
    )

    integrator.close()

    return


def check_validation() -> None:
    print("Checking CudaLangevinThermostatIntegrator validation...")

    integrator = create_integrator()

    expect_exception(
        "CudaLangevinThermostatIntegrator.setThermostatRngSeed rejects negative seed",
        ValueError,
        lambda: integrator.setThermostatRngSeed(-1),
    )
    expect_exception(
        "CudaLangevinThermostatIntegrator.setThermostatRngSeed rejects too-large seed",
        ValueError,
        lambda: integrator.setThermostatRngSeed(2**64),
    )
    expect_exception(
        "CudaLangevinThermostatIntegrator.setCharmmContext rejects non-CharmmContext",
        TypeError,
        lambda: integrator.setCharmmContext(object()),
    )
    expect_exception(
        "CudaLangevinThermostatIntegrator.propagate rejects negative step count",
        ValueError,
        lambda: integrator.propagate(-1),
    )
    expect_exception(
        "CudaLangevinThermostatIntegrator.propagate rejects missing CharmmContext",
        apo.ApoCharmmError,
        lambda: integrator.propagate(1),
    )
    expect_exception(
        "CudaLangevinThermostatIntegrator.initializeFromRestartFile rejects missing CharmmContext",
        apo.ApoCharmmError,
        lambda: integrator.initializeFromRestartFile("missing.rst"),
    )

    ctx = create_context()
    integrator.setCharmmContext(ctx)

    expect_exception(
        "CudaLangevinThermostatIntegrator.setCharmmContext rejects second context",
        apo.ApoCharmmError,
        lambda: integrator.setCharmmContext(ctx),
    )
    expect_exception(
        "CudaLangevinThermostatIntegrator.initializeFromRestartFile rejects missing restart file",
        apo.ApoCharmmError,
        lambda: integrator.initializeFromRestartFile("missing.rst"),
    )

    integrator.close()

    return


def check_short_propagation() -> None:
    print("Checking CudaLangevinThermostatIntegrator short propagation...")

    ctx = create_context()
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
        "CudaLangevinThermostatIntegrator.getInstantaneousTemperature",
        integrator.getInstantaneousTemperature(),
    )
    assert_finite_temperature(
        "CudaLangevinThermostatIntegrator.getAverageTemperature",
        integrator.getAverageTemperature(),
    )

    if final_coordinates_charges[0][:3] == initial_coordinates_charges[0][:3]:
        raise AssertionError("short propagation did not move atom 0 coordinates")

    integrator.close()

    return


def check_deterministic_trajectories() -> None:
    print("Checking deterministic Langevin thermostat trajectories...")

    ctx1 = create_context()
    ctx2 = create_context()

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
        "deterministic initial thermostat friction",
        integrator1.getThermostatFriction(),
        integrator2.getThermostatFriction(),
        DETERMINISTIC_TOLERANCE,
    )
    assert_equal(
        "deterministic initial thermostat seed",
        integrator1.getThermostatRngSeed(),
        integrator2.getThermostatRngSeed(),
    )

    initial_coordinates_charges1 = ctx1.getCoordinatesCharges()
    initial_coordinates_charges2 = ctx2.getCoordinatesCharges()

    integrator1.propagate(NUM_STEPS)
    integrator2.propagate(NUM_STEPS)

    final_coordinates_charges1 = ctx1.getCoordinatesCharges()
    final_coordinates_charges2 = ctx2.getCoordinatesCharges()

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
        + "PASS: CudaLangevinThermostatIntegrator Python API tests completed."
        + "\033[0m"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
