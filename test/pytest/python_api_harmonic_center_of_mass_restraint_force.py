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
    assert_finite_temperature,
    expect_exception,
)

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
RANDOM_SEED: int = 314159
TEMPERATURE: float = 300.0
TIME_STEP: float = 0.001


def create_system() -> tuple[
    apo.CharmmParameters,
    apo.CharmmPsf,
    apo.CharmmCrd,
    apo.ForceManager,
    apo.CharmmContext,
]:
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

    fm = ctx.getForceManager()

    return prm, psf, crd, fm, ctx


def create_configured_restraint(
    psf: apo.CharmmPsf,
) -> apo.HarmonicCenterOfMassRestraintForce:
    selector = apo.AtomSelector(psf)
    selection = selector.select("all")

    restraint = apo.HarmonicCenterOfMassRestraintForce(psf.getNumAtoms())
    restraint.setSelection(selection)
    restraint.setReferencePosition([1.41, 1.41, 1.41])
    restraint.setReferencePosition([1.41, 1.41, 1.41], [1, 1, 1])
    restraint.setReferenceDistance(0.0)
    restraint.setMasses(psf.getMasses())
    restraint.setMassWeighting(False)
    restraint.setForceConstant(1.0)

    return restraint


def check_construction_and_setters() -> None:
    print("Checking HarmonicCenterOfMassRestraintForce construction and setters...")

    _, psf, _, _, _ = create_system()

    restraint = create_configured_restraint(psf)

    assert_equal(
        "HarmonicCenterOfMassRestraintForce.default_force_tag",
        restraint.default_force_tag,
        "hmcm",
    )

    restraint.close()

    return


def check_subscription_and_short_propagation() -> None:
    print("Checking HarmonicCenterOfMassRestraintForce ForceManager subscription...")

    _, psf, crd, fm, ctx = create_system()

    restraint = create_configured_restraint(psf)

    fm.subscribe(restraint)

    expect_exception(
        "ForceManager rejects duplicate HarmonicCenterOfMassRestraintForce subscription",
        ValueError,
        lambda: fm.subscribe(restraint),
    )

    integrator = apo.CudaLangevinThermostatIntegrator(TIME_STEP)
    integrator.setReferenceTemperature(TEMPERATURE)
    integrator.setThermostatFriction(0.0)
    integrator.setThermostatRngSeed(RANDOM_SEED)
    integrator.setCharmmContext(ctx)

    integrator.propagate(10)

    assert_finite_temperature(
        "post harmonic center-of-mass restraint propagation", ctx.computeTemperature()
    )

    fm.unsubscribe(restraint)

    integrator.close()
    restraint.close()
    fm.close()
    ctx.close()

    return


def check_validation() -> None:
    print("Checking HarmonicCenterOfMassRestraintForce validation...")

    _, psf, _, fm, _ = create_system()
    restraint = apo.HarmonicCenterOfMassRestraintForce(psf.getNumAtoms())
    selector = apo.AtomSelector(psf)
    selection = selector.select("all")

    expect_exception(
        "HarmonicCenterOfMassRestraintForce rejects non-int num_atoms",
        TypeError,
        lambda: apo.HarmonicCenterOfMassRestraintForce(2.0),  # type: ignore[arg-type]
    )
    expect_exception(
        "HarmonicCenterOfMassRestraintForce rejects zero num_atoms",
        ValueError,
        lambda: apo.HarmonicCenterOfMassRestraintForce(0),
    )
    expect_exception(
        "HarmonicCenterOfMassRestraintForce rejects negative num_atoms",
        ValueError,
        lambda: apo.HarmonicCenterOfMassRestraintForce(-1),
    )
    expect_exception(
        "setSelection rejects non-AtomSelection",
        TypeError,
        lambda: restraint.setSelection(object()),  # type: ignore[arg-type]
    )
    expect_exception(
        "setForceConstant rejects negative force constant",
        apo.ApoCharmmError,
        lambda: restraint.setForceConstant(-1.0),
    )
    expect_exception(
        "setReferencePosition rejects wrong coordinate length",
        apo.ApoCharmmError,
        lambda: restraint.setReferencePosition([0.0, 0.0]),
    )
    expect_exception(
        "setReferencePosition rejects wrong reference mask length",
        apo.ApoCharmmError,
        lambda: restraint.setReferencePosition([0.0, 0.0, 0.0], [1, 1]),
    )
    expect_exception(
        "setReferencePosition rejects inactive reference mask",
        apo.ApoCharmmError,
        lambda: restraint.setReferencePosition([0.0, 0.0, 0.0], [0, 0, 0]),
    )
    expect_exception(
        "setReferencePosition rejects invalid reference mask value",
        apo.ApoCharmmError,
        lambda: restraint.setReferencePosition([0.0, 0.0, 0.0], [1, 2, 1]),
    )
    expect_exception(
        "setReferenceDistance rejects negative reference distance",
        apo.ApoCharmmError,
        lambda: restraint.setReferenceDistance(-1.0),
    )
    expect_exception(
        "setMasses rejects wrong length",
        apo.ApoCharmmError,
        lambda: restraint.setMasses([1.0]),
    )
    expect_exception(
        "subscribe rejects non-ForceManager",
        TypeError,
        lambda: restraint._subscribe_to_force_manager(object()),  # type: ignore[arg-type]
    )
    expect_exception(
        "subscribe rejects empty force tag",
        ValueError,
        lambda: restraint._subscribe_to_force_manager(fm, ""),
    )

    restraint.setSelection(selection)
    restraint.setReferencePosition([1.41, 1.41, 1.41])
    restraint.setReferenceDistance(0.0)
    restraint.setForceConstant(1.0)
    restraint.setMassWeighting(False)

    fm.subscribe(restraint, "custom-hmcm")
    fm.unsubscribe(restraint)

    restraint.close()
    fm.close()

    return


def check_close_invalidates_handle() -> None:
    print("Checking HarmonicCenterOfMassRestraintForce close invalidates handle...")

    restraint = apo.HarmonicCenterOfMassRestraintForce(2)
    restraint.close()

    expect_exception(
        "closed HarmonicCenterOfMassRestraintForce rejects setter",
        RuntimeError,
        lambda: restraint.setForceConstant(1.0),
    )

    return


def main(argc: int, argv: list[str]) -> int:
    check_construction_and_setters()
    check_subscription_and_short_propagation()
    check_validation()
    check_close_invalidates_handle()

    print(
        "\033[32m"
        + "PASS: HarmonicCenterOfMassRestraintForce Python API tests completed."
        + "\033[0m"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
