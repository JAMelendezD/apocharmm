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

from python_api_test_helpers import assert_finite_temperature, require_file

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
RANDOM_SEED: int = 314159
TEMPERATURE: float = 300.0
TIME_STEP: float = 0.002
THERMOSTAT_FRICTION: float = 1.0
NUM_STEPS_WITH_RESTRAINT: int = 10000
NUM_STEPS_AFTER_UNSUBSCRIBE: int = 1000


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")

    parameter_path: str = require_file(repo_root / "test/data/toppar_water_ions.str")
    psf_path: str = require_file(repo_root / "test/data/nacl_pair.psf")
    crd_path: str = require_file(repo_root / "test/data/nacl_pair.cor")

    print("Loading input files...")
    prm = apo.CharmmParameters(parameter_path)
    psf = apo.CharmmPsf(psf_path)
    crd = apo.CharmmCrd(crd_path)

    print("Creating ForceManager and CharmmContext...")
    fm = apo.ForceManager(psf, prm)
    fm.setBoxDimensions(BOX_DIMENSIONS)

    ctx = apo.CharmmContext(fm)
    ctx.setCoordinates(crd)
    ctx.setRandomSeedForVelocities(RANDOM_SEED)
    ctx.assignVelocitiesAtTemperature(TEMPERATURE)

    print("Creating Langevin thermostat integrator...")
    integrator = apo.CudaLangevinThermostatIntegrator(TIME_STEP)
    integrator.setReferenceTemperature(TEMPERATURE)
    integrator.setThermostatFriction(THERMOSTAT_FRICTION)
    integrator.setThermostatRngSeed(RANDOM_SEED)
    integrator.setCharmmContext(ctx)

    initial_temperature: float = ctx.computeTemperature()
    assert_finite_temperature("initial", initial_temperature)
    print(f"Initial temperature: {initial_temperature:.6f} K")

    print("Creating atom selection")
    selector = apo.AtomSelector(psf)
    selection = selector.select("all")

    print("Creating and configuring harmonic restraint...")
    harm = apo.HarmonicRestraintForce(psf.getNumAtoms())
    harm.setSelection(selection)
    harm.setReferenceCoordinates(crd.getCoordinates())
    harm.setMasses(psf.getMasses())
    harm.setForceConstant(1.0)

    print("Subscribing harmonic restraint through ForceManager...")
    fm.subscribe(harm)

    print(f"Propagating {NUM_STEPS_WITH_RESTRAINT} steps with restraint...")
    integrator.propagate(NUM_STEPS_WITH_RESTRAINT)

    restrained_temperature: float = ctx.computeTemperature()
    assert_finite_temperature("restrained", restrained_temperature)
    print(f"Temperature after restrained propagation: {restrained_temperature:.6f} K")

    print("Updating harmonic restraint force constant...")
    harm.setForceConstant(0.5)

    print("Propagating one more step after updating restraint...")
    integrator.propagate(1)

    updated_temperature: float = ctx.computeTemperature()
    assert_finite_temperature("updated-restraint", updated_temperature)
    print(f"Temperature after restraint update: {updated_temperature:.6f} K")

    print("Unsubscribing harmonic restraint through ForceManager...")
    fm.unsubscribe(harm)

    print(f"Propagating {NUM_STEPS_AFTER_UNSUBSCRIBE} steps after unsubscribe...")
    integrator.propagate(NUM_STEPS_AFTER_UNSUBSCRIBE)

    final_temperature: float = ctx.computeTemperature()
    assert_finite_temperature("final", final_temperature)
    print(f"Final temperature: {final_temperature:.6f} K")

    print(
        "PASS: harmonic restraint API construction, setters, subscribe, update, unsubscribe, and propagation all completed."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
