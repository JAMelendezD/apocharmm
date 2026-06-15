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
    assert_file_created,
    assert_finite_temperature,
    require_file,
)

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
RANDOM_SEED: int = 314159
TEMPERATURE: float = 300.0
TIME_STEP: float = 0.002
THERMOSTAT_FRICTION: float = 1.0
INITIAL_REPORT_FREQUENCY: int = 1
REPORT_FREQUENCY: int = 2
NUM_STEPS: int = 4


def remove_if_exists(path: Path) -> None:
    if path.exists():
        path.unlink()
    return


def create_context(repo_root: Path, assign_velocities: bool) -> apo.CharmmContext:
    parameter_path: str = require_file(repo_root / "test/data/toppar_water_ions.str")
    psf_path: str = require_file(repo_root / "test/data/nacl_pair.psf")
    crd_path: str = require_file(repo_root / "test/data/nacl_pair.cor")

    prm = apo.CharmmParameters(parameter_path)
    psf = apo.CharmmPsf(psf_path)
    crd = apo.CharmmCrd(crd_path)

    fm = apo.ForceManager(psf, prm)
    fm.setBoxDimensions(BOX_DIMENSIONS)

    ctx = apo.CharmmContext(fm)
    ctx.setCoordinates(crd)

    if assign_velocities:
        ctx.setRandomSeedForVelocities(RANDOM_SEED)
        ctx.assignVelocitiesAtTemperature(TEMPERATURE)

    return ctx


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")
    output_dir: Path = repo_root / "test/pytest"
    output_dir.mkdir(parents=True, exist_ok=True)

    dcd_path: Path = output_dir / "tmp_python_api_subscriber.dcd"
    rst_path: Path = output_dir / "tmp_python_api_subscriber.rst"
    for path in (dcd_path, rst_path):
        remove_if_exists(path)

    print("Creating context and Langevin thermostat integrator...")
    ctx = create_context(repo_root, assign_velocities=True)
    integrator = apo.CudaLangevinThermostatIntegrator(TIME_STEP)
    integrator.setReferenceTemperature(TEMPERATURE)
    integrator.setThermostatFriction(THERMOSTAT_FRICTION)
    integrator.setThermostatRngSeed(RANDOM_SEED)
    integrator.setCharmmContext(ctx)

    print("Creating DCD and restart subscribers...")
    dcd = apo.DcdSubscriber(str(dcd_path), INITIAL_REPORT_FREQUENCY)
    rst = apo.RestartSubscriber(str(rst_path), INITIAL_REPORT_FREQUENCY)

    dcd.setReportFrequency(REPORT_FREQUENCY)
    rst.setReportFrequency(REPORT_FREQUENCY)

    if dcd.getReportFrequency() != REPORT_FREQUENCY:
        raise AssertionError("DcdSubscriber.getReportFrequency mismatch")

    if rst.getReportFrequency() != REPORT_FREQUENCY:
        raise AssertionError("RestartSubscriber.getReportFrequency mismatch")

    print("Subscribing outputs through CudaIntegrator.subscribe...")
    integrator.subscribe(dcd)
    integrator.subscribe(rst)

    print(f"Propagating {NUM_STEPS} steps to write subscriber outputs...")
    integrator.propagate(NUM_STEPS)
    assert_finite_temperature("post-subscriber", ctx.computeTemperature())

    print("Closing output objects before checking files...")
    integrator.close()
    dcd.close()
    rst.close()

    assert_file_created(dcd_path)
    assert_file_created(rst_path)

    print("Checking CudaIntegrator.initializeFromRestartFile...")
    rst_ctx = create_context(repo_root, assign_velocities=False)
    rst_integrator = apo.CudaLangevinThermostatIntegrator(TIME_STEP)
    rst_integrator.setReferenceTemperature(TEMPERATURE)
    rst_integrator.setThermostatFriction(THERMOSTAT_FRICTION)
    rst_integrator.setCharmmContext(rst_ctx)
    rst_integrator.initializeFromRestartFile(str(rst_path))
    rst_integrator.propagate(1)

    assert_finite_temperature("restart", rst_ctx.computeTemperature())

    print("Cleaning up subscriber output files...")
    rst_integrator.close()
    for path in (dcd_path, rst_path):
        remove_if_exists(path)

    print(
        "PASS: DcdSubscriber, RestartSubscriber, Subscriber report-frequency API, CudaIntegrator.subscribe, and restart initialization all completed."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
