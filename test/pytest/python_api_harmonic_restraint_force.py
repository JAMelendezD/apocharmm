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
    ctx.setRandomSeedForVelocities(RANDOM_SEED)
    ctx.assignVelocitiesAtTemperature(TEMPERATURE)

    fm = ctx.getForceManager()

    return prm, psf, crd, fm, ctx


def create_configured_restraint(
    psf: apo.CharmmPsf, crd: apo.CharmmCrd
) -> apo.HarmonicRestraintForce:
    selector = apo.AtomSelector(psf)
    selection = selector.select("all")

    restraint = apo.HarmonicRestraintForce(psf.getNumAtoms())
    restraint.setSelection(selection)
    restraint.setReferenceCoordinates(crd.getCoordinates())
    restraint.setMasses(psf.getMasses())
    restraint.setForceConstant(10.0)
    restraint.setForceConstants([10.0 for _ in range(psf.getNumAtoms())])
    restraint.setBoxDimensions(BOX_DIMENSIONS)

    return restraint


def check_construction_and_setters() -> None:
    print("Checking HarmonicRestraintForce construction and setters...")

    _, psf, crd, _, _ = create_system()

    restraint = create_configured_restraint(psf, crd)

    assert_equal(
        "HarmonicRestraintForce.default_force_tag", restraint.default_force_tag, "harm"
    )

    restraint.close()

    return


def check_subscription_and_short_propagation() -> None:
    print("Checking HarmonicRestraintForce ForceManager subscription...")

    _, psf, crd, fm, ctx = create_system()

    restraint = create_configured_restraint(psf, crd)

    fm.subscribe(restraint)

    expect_exception(
        "ForceManager rejects duplicate HarmonicRestraintForce subscription",
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
        "post harmonic restraint propagation", ctx.computeTemperature()
    )

    fm.unsubscribe(restraint)

    integrator.close()
    restraint.close()
    fm.close()
    ctx.close()

    return


def check_validation() -> None:
    print("Checking HarmonicRestraintForce validation...")

    _, psf, crd, fm, _ = create_system()
    restraint = apo.HarmonicRestraintForce(psf.getNumAtoms())
    selector = apo.AtomSelector(psf)
    selection = selector.select("all")

    expect_exception(
        "HarmonicRestraintForce rejects non-int num_atoms",
        TypeError,
        lambda: apo.HarmonicRestraintForce(2.0),  # type: ignore[arg-type]
    )
    expect_exception(
        "HarmonicRestraintForce rejects zero num_atoms",
        ValueError,
        lambda: apo.HarmonicRestraintForce(0),
    )
    expect_exception(
        "HarmonicRestraintForce rejects negative num_atoms",
        ValueError,
        lambda: apo.HarmonicRestraintForce(-1),
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
        "setForceConstants rejects wrong length",
        apo.ApoCharmmError,
        lambda: restraint.setForceConstants([1.0]),
    )
    expect_exception(
        "setReferenceCoordinates rejects wrong coordinate length",
        ValueError,
        lambda: restraint.setReferenceCoordinates([[0.0, 0.0]]),
    )
    expect_exception(
        "setReferenceCoordinates rejects wrong atom count",
        apo.ApoCharmmError,
        lambda: restraint.setReferenceCoordinates([[0.0, 0.0, 0.0]]),
    )
    expect_exception(
        "setMasses rejects wrong length",
        apo.ApoCharmmError,
        lambda: restraint.setMasses([1.0]),
    )
    expect_exception(
        "setBoxDimensions rejects wrong length",
        apo.ApoCharmmError,
        lambda: restraint.setBoxDimensions([50.0, 50.0]),
    )
    expect_exception(
        "setBoxDimensions rejects negative dimension",
        apo.ApoCharmmError,
        lambda: restraint.setBoxDimensions([50.0, -1.0, 50.0]),
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
    restraint.setReferenceCoordinates(crd.getCoordinates())
    restraint.setMasses(psf.getMasses())
    restraint.setForceConstant(1.0)

    fm.subscribe(restraint, "custom-harmonic-restraint")
    fm.unsubscribe(restraint)

    restraint.close()
    fm.close()

    return


def check_close_invalidates_handle() -> None:
    print("Checking HarmonicRestraintForce close invalidates handle...")

    restraint = apo.HarmonicRestraintForce(2)
    restraint.close()

    expect_exception(
        "closed HarmonicRestraintForce rejects setter",
        RuntimeError,
        lambda: restraint.setForceConstant(1.0),
    )


def main(argc: int, argv: list[str]) -> int:
    check_construction_and_setters()
    check_subscription_and_short_propagation()
    check_validation()
    check_close_invalidates_handle()

    print(
        "\033[32m"
        + "PASS: HarmonicRestraintForce Python API tests completed."
        + "\033[0m"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
