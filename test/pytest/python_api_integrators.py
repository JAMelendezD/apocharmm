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
    assert_close,
    assert_finite_temperature,
    assert_finite_scalar,
    assert_finite_sequence,
    assert_sequence_close,
    require_file,
)

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
RANDOM_SEED: int = 314159
TEMPERATURE: float = 300.0
TIME_STEP: float = 0.002
THERMOSTAT_FRICTION: float = 1.0
NOSE_HOOVER_PISTON_MASS: float = 500.0
LANGEVIN_PISTON_MASS: list[float] = [5.0]
LANGEVIN_PISTON_FRICTION: float = 5.0
REFERENCE_PRESSURE_TENSOR: list[list[float]] = [
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
]
NUM_STEPS: int = 10
TOLERANCE: float = 1.0e-10


def create_context(
    repo_root: Path,
) -> apo.CharmmContext:
    prm_path: str = require_file(repo_root / "test/data/toppar_water_ions.str")
    psf_path: str = require_file(repo_root / "test/data/waterbox.psf")
    crd_path: str = require_file(repo_root / "test/data/waterbox.crd")

    prm = apo.CharmmParameters(prm_path)
    psf = apo.CharmmPsf(psf_path)
    crd = apo.CharmmCrd(crd_path)

    fm = apo.ForceManager(psf, prm)
    fm.setBoxDimensions(BOX_DIMENSIONS)

    ctx = apo.CharmmContext(fm)
    ctx.setCoordinates(crd)
    ctx.setRandomSeedForVelocities(RANDOM_SEED)
    ctx.assignVelocitiesAtTemperature(TEMPERATURE)

    return ctx


def check_langevin_thermostat(repo_root: Path) -> None:
    print("Checking CudaLangevinThermostatIntegrator API...")
    ctx = create_context(repo_root)

    integrator = apo.CudaLangevinThermostatIntegrator(TIME_STEP)
    integrator.setTimeStep(TIME_STEP)
    integrator.setReferenceTemperature(TEMPERATURE)
    integrator.setThermostatFriction(THERMOSTAT_FRICTION)
    integrator.setThermostatRngSeed(RANDOM_SEED)

    assert_close(
        "CudaLangevinThermostatIntegrator.getReferenceTemperature",
        integrator.getReferenceTemperature(),
        TEMPERATURE,
        TOLERANCE,
    )
    assert_close(
        "CudaLangevinThermostatIntegrator.getThermostatFriction",
        integrator.getThermostatFriction(),
        THERMOSTAT_FRICTION,
        TOLERANCE,
    )
    if integrator.getThermostatRngSeed() != RANDOM_SEED:
        raise AssertionError(
            "CudaLangevinThermostatIntegrator.getThermostatRngSeed mismatch"
        )

    integrator.setCharmmContext(ctx)
    integrator.resetAverageTemperature()
    integrator.propagate(NUM_STEPS)

    assert_finite_temperature(
        "CharmmContext.computeTemperature", ctx.computeTemperature()
    )
    assert_finite_temperature(
        "CudaLangevinThermostatIntegrator.getAverageTemperature",
        integrator.getAverageTemperature(),
    )
    assert_finite_temperature(
        "CudaLangevinThermostatIntegrator.getInstantaneousTemperature",
        integrator.getInstantaneousTemperature(),
    )

    return


def check_nose_hoover(repo_root: Path) -> None:
    print("Checking CudaNoseHooverIntegrator API...")
    ctx = create_context(repo_root)

    integrator = apo.CudaNoseHooverIntegrator(TIME_STEP)
    integrator.setTimeStep(TIME_STEP)
    integrator.setReferenceTemperature(TEMPERATURE)
    integrator.setNoseHooverPistonMass(NOSE_HOOVER_PISTON_MASS)
    integrator.useOldTemperature(False)

    assert_close(
        "CudaNoseHooverIntegrator.getReferenceTemperature",
        integrator.getReferenceTemperature(),
        TEMPERATURE,
        TOLERANCE,
    )
    assert_close(
        "CudaNoseHooverIntegrator.getNoseHooverPistonMass",
        integrator.getNoseHooverPistonMass(),
        NOSE_HOOVER_PISTON_MASS,
        TOLERANCE,
    )

    integrator.setCharmmContext(ctx)
    integrator.resetAverageTemperature()
    integrator.propagate(NUM_STEPS)

    assert_finite_temperature(
        "CharmmContext.computeTemperature", ctx.computeTemperature()
    )
    assert_finite_temperature(
        "CudaNoseHooverIntegrator.getAverageTemperature",
        integrator.getAverageTemperature(),
    )
    assert_finite_temperature(
        "CudaNoseHooverIntegrator.getInstantaneousTemperature",
        integrator.getInstantaneousTemperature(),
    )

    return


def check_langevin_piston(repo_root: Path) -> None:
    print("Checking CudaLangevinPistonIntegrator API...")
    ctx = create_context(repo_root)

    integrator = apo.CudaLangevinPistonIntegrator(TIME_STEP)
    integrator.setTimeStep(TIME_STEP)
    integrator.useNoseHooverThermostat(True)
    integrator.setReferenceTemperature(TEMPERATURE)
    integrator.setNoseHooverPistonMass(NOSE_HOOVER_PISTON_MASS)
    integrator.useOldTemperature(False)
    integrator.setReferencePressure(REFERENCE_PRESSURE_TENSOR)
    integrator.setConstantSurfaceTension(False)
    integrator.setCrystalType(apo.CrystalType.CUBIC)
    integrator.setLangevinPistonMass(LANGEVIN_PISTON_MASS)
    integrator.setLangevinPistonFrictionSeed(RANDOM_SEED)
    integrator.setLangevinPistonFriction(LANGEVIN_PISTON_FRICTION)

    assert_close(
        "CudaLangevinPistonIntegrator.getReferenceTemperature",
        integrator.getReferenceTemperature(),
        TEMPERATURE,
        TOLERANCE,
    )
    assert_close(
        "CudaLangevinPistonIntegrator.getNoseHooverPistonMass",
        integrator.getNoseHooverPistonMass(),
        NOSE_HOOVER_PISTON_MASS,
        TOLERANCE,
    )
    for i in range(3):
        assert_sequence_close(
            "CudaLangevinPistonIntegrator.getReferencePressure",
            integrator.getReferencePressureTensor()[i],
            REFERENCE_PRESSURE_TENSOR[i],
            TOLERANCE,
        )
    if integrator.getCrystalType() != apo.CrystalType.CUBIC:
        raise AssertionError("CudaLangevinPistonIntegrator.getCrystalType mismatch")
    assert_sequence_close(
        "CudaLangevinPistonIntegrator.getLangevinPistonMass",
        integrator.getLangevinPistonMass(),
        LANGEVIN_PISTON_MASS,
        TOLERANCE,
    )
    if integrator.getLangevinPistonFrictionSeed() != RANDOM_SEED:
        raise AssertionError(
            "CudaLangevinPistonIntegrator.getLangevinPistonFrictionSeed mismatch"
        )

    integrator.setCharmmContext(ctx)
    integrator.resetAverages()
    integrator.propagate(NUM_STEPS)

    assert_finite_temperature(
        "CharmmContext.computeTemperature", ctx.computeTemperature()
    )
    assert_finite_temperature(
        "CudaLangevinPistonIntegrator.getAverageTemperature",
        integrator.getAverageTemperature(),
    )
    assert_finite_temperature(
        "CudaLangevinPistonIntegrator.getInstantaneousTemperature",
        integrator.getInstantaneousTemperature(),
    )
    for i in range(3):
        assert_finite_sequence(
            "CudaLangevinPistonIntegrator.getInstantaneousPressureTensor",
            integrator.getInstantaneousPressureTensor()[i],
        )
    assert_finite_scalar(
        "CudaLangevinPistonIntegrator.getInstantaneousPressureScalar",
        integrator.getInstantaneousPressureScalar(),
    )
    for i in range(3):
        assert_finite_sequence(
            "CudaLangevinPistonIntegrator.getAveragePressureTensor",
            integrator.getAveragePressureTensor()[i],
        )
    assert_finite_scalar(
        "CudaLangevinPistonIntegrator.getAveragePressureScalar",
        integrator.getAveragePressureScalar(),
    )

    return


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")

    check_langevin_thermostat(repo_root)
    check_nose_hoover(repo_root)
    check_langevin_piston(repo_root)

    print(
        "PASS: Integrator API setters, getters, base CudaIntegrator methods, context binding, and short propagation all completed."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
