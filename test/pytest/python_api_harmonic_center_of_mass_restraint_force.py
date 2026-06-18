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
    require_file,
)

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
SODIUM_REFERENCE: list[float] = [0.0, 0.0, 0.0]
CHLORIDE_REFERENCE: list[float] = [2.82, 2.82, 2.82]
CENTER_REFERENCE: list[float] = [1.41, 1.41, 1.41]


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")

    prm_path: str = require_file(repo_root / "test/data/toppar_water_ions.str")
    psf_path: str = require_file(repo_root / "test/data/nacl_pair.psf")
    crd_path: str = require_file(repo_root / "test/data/nacl_pair.cor")

    print("Loading input files...")
    prm = apo.CharmmParameters(prm_path)
    psf = apo.CharmmPsf(psf_path)
    crd = apo.CharmmCrd(crd_path)

    print("Creating atom selections...")
    selector = apo.AtomSelector(psf)
    all_selection = selector.select("all")
    sodium_selection = selector.select("type SOD")
    chloride_selection = selector.select("type CLA")

    print("Checking HMCM construction and simple setters...")
    hmcm = apo.HarmonicCenterOfMassRestraintForce(psf.getNumAtoms())
    assert_equal("HMCM default force tag", hmcm.default_force_tag, "hmcm")

    hmcm.setSelection(all_selection)
    hmcm.setForceConstant(2.5)
    hmcm.setReferencePosition(CENTER_REFERENCE)
    hmcm.setReferencePosition([1.0, 0.0, 0.0], [1, 0, 0])
    hmcm.setReferenceDistance(1.0)
    hmcm.setMasses(psf.getMasses())
    hmcm.setMassWeighting(True)
    hmcm.setMassWeighting(False)

    print("Checking HMCM input validation...")
    expect_exception(
        "HMCM rejects non-int atom count",
        TypeError,
        lambda: apo.HarmonicCenterOfMassRestraintForce(2.0),
    )
    expect_exception(
        "HMCM rejects zero atoms",
        ValueError,
        lambda: apo.HarmonicCenterOfMassRestraintForce(0),
    )
    expect_exception(
        "HMCM rejects non-selection object", TypeError, lambda: hmcm.setSelection("all")
    )
    expect_exception(
        "HMCM rejects negative force constant",
        apo.ApoCharmmError,
        lambda: hmcm.setForceConstant(-1.0),
    )
    expect_exception(
        "HMCM rejects reference position with wrong length",
        apo.ApoCharmmError,
        lambda: hmcm.setReferencePosition([0.0, 0.0]),
    )
    expect_exception(
        "HMCM rejects reference mask with wrong length",
        apo.ApoCharmmError,
        lambda: hmcm.setReferencePosition([0.0, 0.0, 0.0], [1, 1]),
    )
    expect_exception(
        "HMCM rejects inactive reference mask",
        apo.ApoCharmmError,
        lambda: hmcm.setReferencePosition([0.0, 0.0, 0.0], [0, 0, 0]),
    )
    expect_exception(
        "HMCM rejects negative reference distance",
        apo.ApoCharmmError,
        lambda: hmcm.setReferenceDistance(-1.0),
    )
    expect_exception(
        "HMCM rejects mass list with wrong length",
        apo.ApoCharmmError,
        lambda: hmcm.setMasses([1.0]),
    )

    print("Checking ForceManager subscription paths...")
    fm = apo.ForceManager(psf, prm)
    fm.setBoxDimensions(BOX_DIMENSIONS)

    fm.subscribe(hmcm)
    expect_exception(
        "ForceManager rejects duplicate HMCM subscription",
        ValueError,
        lambda: fm.subscribe(hmcm),
    )
    fm.unsubscribe(hmcm)

    sodium_hmcm = apo.HarmonicCenterOfMassRestraintForce(psf.getNumAtoms())
    sodium_hmcm.setSelection(sodium_selection)
    sodium_hmcm.setForceConstant(1.0)
    sodium_hmcm.setReferencePosition(SODIUM_REFERENCE)
    sodium_hmcm.setMasses(psf.getMasses())

    chloride_hmcm = apo.HarmonicCenterOfMassRestraintForce(psf.getNumAtoms())
    chloride_hmcm.setSelection(chloride_selection)
    chloride_hmcm.setForceConstant(1.0)
    chloride_hmcm.setReferencePosition(CHLORIDE_REFERENCE)
    chloride_hmcm.setMasses(psf.getMasses())

    fm.subscribe(sodium_hmcm, "sodium_hmcm")
    fm.subscribe(chloride_hmcm, "chloride_hmcm")

    print("Checking HMCM with a one-step dynamics smoke test...")
    ctx = apo.CharmmContext(fm)
    ctx.setCoordinates(crd)
    ctx.useHolonomicConstraints(True)
    ctx.setRandomSeedForVelocities(314159)
    ctx.assignVelocitiesAtTemperature(10.0)

    integrator = apo.CudaLangevinThermostatIntegrator(0.001)
    integrator.setCharmmContext(ctx)
    integrator.setReferenceTemperature(10.0)
    integrator.setThermostatFriction(1.0)
    integrator.setThermostatRngSeed(314159)
    integrator.propagate(1)

    fm.unsubscribe(sodium_hmcm)
    fm.unsubscribe(chloride_hmcm)

    integrator.propagate(1)

    print("PASS: HarmonicCenterOfMassRestraintForce Python API test completed.")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
