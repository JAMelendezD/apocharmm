# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from __future__ import annotations

import math
from pathlib import Path
import sys

import apocharmm as apo

from python_api_test_helpers import expect_exception, require_file

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
TIME_STEP: float = 0.002


def check_primary_object_validation(
    parameter_path: str, psf: apo.CharmmPsf, prm: apo.CharmmParameters
) -> None:
    print("Checking primary object validation...")
    expect_exception(
        "ForceManager rejects a non-CharmmPsf psf argument",
        TypeError,
        lambda: apo.ForceManager(prm, prm),
    )
    expect_exception(
        "ForceManager rejects a non-CharmmParameters parameter argument",
        TypeError,
        lambda: apo.ForceManager(psf, psf),
    )

    return


def check_manager_context_validation(
    fm: apo.ForceManager, ctx: apo.CharmmContext, psf: apo.CharmmPsf
) -> None:
    print("Checking ForceManager and CharmmContext validation...")
    expect_exception(
        "ForceManager.setBoxDimensions rejects wrong length",
        apo.ApoCharmmError,
        lambda: fm.setBoxDimensions([1.0, 2.0]),
    )
    expect_exception(
        "ForceManager.setPeriodicBoundaryCondition rejects invalid integer",
        ValueError,
        lambda: fm.setPeriodicBoundaryCondition(999),
    )
    expect_exception(
        "CharmmContext rejects non-ForceManager construction",
        TypeError,
        lambda: apo.CharmmContext(psf),
    )
    expect_exception(
        "CharmmContext.setCoordinates rejects non-CharmmCrd",
        TypeError,
        lambda: ctx.setCoordinates(psf),
    )
    expect_exception(
        "CharmmContext.setPeriodicBoundaryCondition rejects invalid integer",
        ValueError,
        lambda: ctx.setPeriodicBoundaryCondition(999),
    )
    expect_exception(
        "CharmmContext.setRandomSeedForVelocities rejects negative seed",
        ValueError,
        lambda: ctx.setRandomSeedForVelocities(-1),
    )

    return


def check_selection_validation(
    selector: apo.AtomSelector, selection: apo.AtomSelection
) -> None:
    print("Checking atom-selection validation...")
    expect_exception(
        "AtomSelector rejects a non-string selection",
        TypeError,
        lambda: selector.select(1),
    )
    expect_exception(
        "AtomSelector rejects an empty selection string",
        ValueError,
        lambda: selector.select(""),
    )
    expect_exception(
        "AtomSelection.contains rejects bool atom index",
        TypeError,
        lambda: selection.contains(True),
    )
    expect_exception(
        "AtomSelection.contains reports C API out-of-range errors",
        apo.ApoCharmmError,
        lambda: selection.contains(selection.getNumAtoms()),
    )

    return


def check_harmonic_restraint_validation(
    fm: apo.ForceManager,
    psf: apo.CharmmPsf,
    crd: apo.CharmmCrd,
    selection: apo.AtomSelection,
) -> None:
    print("Checking HarmonicRestraintForce validation...")
    expect_exception(
        "HarmonicRestraintForce rejects zero atoms",
        ValueError,
        lambda: apo.HarmonicRestraintForce(0),
    )

    harm = apo.HarmonicRestraintForce(psf.getNumAtoms())
    expect_exception(
        "HarmonicRestraintForce.setSelection rejects non-AtomSelection",
        TypeError,
        lambda: harm.setSelection(psf),
    )
    expect_exception(
        "HarmonicRestraintForce.setForceConstant rejects non-finite value",
        apo.ApoCharmmError,
        lambda: harm.setForceConstant(math.inf),
    )
    expect_exception(
        "HarmonicRestraintForce.setForceConstant rejects negative value",
        apo.ApoCharmmError,
        lambda: harm.setForceConstant(-1.0),
    )
    expect_exception(
        "HarmonicRestraintForce.setReferenceCoordiantes rejects wrong coordinate length",
        ValueError,
        lambda: harm.setReferenceCoordinates([[0.0, 0.0]]),
    )
    expect_exception(
        "HarmonicRestraintForce.setBoxDimensions rejects wrong length",
        apo.ApoCharmmError,
        lambda: harm.setBoxDimensions([1.0, 2.0]),
    )
    expect_exception(
        "HarmonicRestraintForce.setBoxDimensions rejects non-finite entries",
        apo.ApoCharmmError,
        lambda: harm.setBoxDimensions([1.0, math.inf, 3.0]),
    )
    expect_exception(
        "HarmonicRestraintForce.setBoxDimensions rejects negative entries",
        apo.ApoCharmmError,
        lambda: harm.setBoxDimensions([1.0, -2.0, 3.0]),
    )
    expect_exception(
        "ForceManager.subscribe rejects non-subscribable objects",
        TypeError,
        lambda: fm.subscribe(object()),
    )
    expect_exception(
        "ForceManager.subscribe rejects non-string force tags",
        TypeError,
        lambda: fm.subscribe(harm, 1),
    )
    expect_exception(
        "ForceManager.subscribe rejects empty force tags",
        ValueError,
        lambda: fm.subscribe(harm, ""),
    )

    harm.setSelection(selection)
    harm.setReferenceCoordinates(crd.getCoordinates())
    harm.setMasses(psf.getMasses())
    harm.setForceConstant(1.0)
    fm.subscribe(harm)
    expect_exception(
        "ForceManager.subscribe rejects duplicate force objects",
        ValueError,
        lambda: fm.subscribe(harm),
    )
    fm.unsubscribe(harm)

    return


def check_integrator_validation(ctx: apo.CharmmContext, psf: apo.CharmmPsf) -> None:
    print("Checking integrator validation...")
    lt = apo.CudaLangevinThermostatIntegrator(TIME_STEP)
    expect_exception(
        "CudaIntegrator.setCharmmContext rejects non-CharmmContext",
        TypeError,
        lambda: lt.setCharmmContext(psf),
    )
    expect_exception(
        "CudaIntegrator.subscribe rejects non-Subscriber",
        TypeError,
        lambda: lt.subscribe(psf),
    )
    expect_exception(
        "CudaIntegrator.propagate rejects negative step count",
        ValueError,
        lambda: lt.propagate(-1),
    )

    expect_exception(
        "CudaLangevinThermostatIntegrator.setThermostatRngSeed rejects negative seed",
        ValueError,
        lambda: lt.setThermostatRngSeed(-1),
    )

    lp = apo.CudaLangevinPistonIntegrator(TIME_STEP)
    expect_exception(
        "CudaLangevinPistonIntegrator.setReferencePressure rejects wrong length",
        ValueError,
        lambda: lp.setReferencePressure([[1.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]),
    )
    expect_exception(
        "CudaLangevinPistonIntegrator.setReferencePressure rejects non-finite entries",
        apo.ApoCharmmError,
        lambda: lp.setReferencePressure(
            [[1.0, 0.0, 0.0], [0.0, math.inf, 0.0], [0.0, 0.0, 1.0]]
        ),
    )
    expect_exception(
        "CudaLangevinPistonIntegrator.setCrysalType rejects invalid integer",
        ValueError,
        lambda: lp.setCrystalType(999),
    )
    expect_exception(
        "CudaLangevinPistonIntegrator.setLangevinPistonMass rejects empty list",
        apo.ApoCharmmError,
        lambda: lp.setLangevinPistonMass([]),
    )
    expect_exception(
        "CudaLangevinPistonIntegrator.setLangevinPistonMass rejects too many entries",
        apo.ApoCharmmError,
        lambda: lp.setLangevinPistonMass([1.0, 2.0, 3.0, 4.0]),
    )
    expect_exception(
        "CudaLangevinPistonIntegrator.setLangevinPistonMass rejects non-finite entries",
        apo.ApoCharmmError,
        lambda: lp.setLangevinPistonMass([math.inf]),
    )
    expect_exception(
        "CudaLangevinPistonIntegrator.setLangevinPistonMass rejects negative entries",
        apo.ApoCharmmError,
        lambda: lp.setLangevinPistonMass([-1.0]),
    )
    expect_exception(
        "CudaLangevinPistonIntegrator.setLangevinPistonFrictionSeed rejects negative seed",
        ValueError,
        lambda: lp.setLangevinPistonFrictionSeed(-1),
    )

    return


def check_subscriber_validation(output_dir: Path) -> None:
    print("Checking subscriber validation...")
    dcd_path: Path = output_dir / "tmp_python_api_validation.dcd"
    rst_path: Path = output_dir / "tmp_python_api_validation.rst"

    expect_exception(
        "DcdSubscriber rejects bool report frequency",
        ValueError,
        lambda: apo.DcdSubscriber(str(dcd_path), True),
    )
    expect_exception(
        "Restartsubscriber rejects zero report frequency",
        ValueError,
        lambda: apo.RestartSubscriber(str(rst_path), 0),
    )

    dcd = apo.DcdSubscriber(str(dcd_path), 1)
    expect_exception(
        "Subscriber.setReportFrequency rejects bool",
        ValueError,
        lambda: dcd.setReportFrequency(True),
    )
    expect_exception(
        "Subscriber.setReportFrequency rejects zero",
        ValueError,
        lambda: dcd.setReportFrequency(0),
    )
    dcd.close()

    for path in (dcd_path, rst_path):
        if path.exists():
            path.unlink()

    return


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")
    output_dir: Path = repo_root / "test/pytest"
    output_dir.mkdir(parents=True, exist_ok=True)

    parameter_path: str = require_file(repo_root / "test/data/toppar_water_ions.str")
    psf_path: str = require_file(repo_root / "test/data/nacl_pair.psf")
    crd_path: str = require_file(repo_root / "test/data/nacl_pair.cor")

    print("Loading valid input objects for validation tests...")
    prm = apo.CharmmParameters(parameter_path)
    psf = apo.CharmmPsf(psf_path)
    crd = apo.CharmmCrd(crd_path)

    fm = apo.ForceManager(psf, prm)
    fm.setBoxDimensions(BOX_DIMENSIONS)
    ctx = apo.CharmmContext(fm)
    ctx.setCoordinates(crd)

    selector = apo.AtomSelector(psf)
    selection = selector.select("all")

    check_primary_object_validation(parameter_path, psf, prm)
    check_manager_context_validation(fm, ctx, psf)
    check_selection_validation(selector, selection)
    check_harmonic_restraint_validation(fm, psf, crd, selection)
    check_integrator_validation(ctx, psf)
    check_subscriber_validation(output_dir)

    print("PASS: Python API validation and error-checking paths all completed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
