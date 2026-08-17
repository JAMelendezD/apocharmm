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
import sys

import apocharmm as apo

from python_api_test_helpers import (
    get_data_path,
    require_file,
    assert_equal,
    assert_finite_temperature,
    expect_exception,
    expect_invalid_argument,
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

    expect_invalid_argument(
        "ForceManager rejects duplicate HarmonicRestraintForce subscription",
        lambda: fm.subscribe(restraint),
        "Force is already subscribed to this ForceManager",
        expected_context="ForceManager.subscribe(HarmonicRestraintForce)",
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

    expect_invalid_argument(
        "ForceManager rejects missing HarmonicRestraintForce unsubscription",
        lambda: fm.unsubscribe(restraint),
        "Force is not subscribed to this ForceManager",
        expected_context="ForceManager.unsubscribe(HarmonicRestraintForce)",
    )

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
    expect_invalid_argument(
        "HarmonicRestraintForce rejects zero num_atoms",
        lambda: apo.HarmonicRestraintForce(0),
        "Atom count must be positive; observed 0",
        expected_context="HarmonicRestraintForce construction",
    )
    expect_invalid_argument(
        "HarmonicRestraintForce rejects negative num_atoms",
        lambda: apo.HarmonicRestraintForce(-1),
        "Atom count must be positive; observed -1",
    )
    expect_exception(
        "HarmonicRestraintForce rejects out-of-range num_atoms",
        ValueError,
        lambda: apo.HarmonicRestraintForce(2**31),
    )
    expect_exception(
        "setSelection rejects non-AtomSelection",
        TypeError,
        lambda: restraint.setSelection(object()),  # type: ignore[arg-type]
    )
    force_constant_error = expect_invalid_argument(
        "setForceConstant rejects negative force constant",
        lambda: restraint.setForceConstant(-1.0),
        "Force constant must be non-negative",
        expected_context="HarmonicRestraintForce.setForceConstant(force_constant)",
    )
    assert_equal(
        "setForceConstant rendered native function occurrence count",
        force_constant_error.message.count(
            "apo_harmonic_restraint_force_set_force_constant"
        ),
        1,
    )
    expect_invalid_argument(
        "setForceConstant rejects non-finite force constant",
        lambda: restraint.setForceConstant(math.inf),
        "Force constant must be finite",
    )
    restraint.setForceConstant(0.0)
    expect_invalid_argument(
        "setForceConstants rejects wrong length",
        lambda: restraint.setForceConstants([1.0]),
        "Force-constant array size mismatch",
    )
    expect_invalid_argument(
        "setForceConstants rejects non-finite value",
        lambda: restraint.setForceConstants([1.0, math.inf]),
        "Force constant at index 1 must be finite",
    )
    expect_invalid_argument(
        "setForceConstants rejects negative value",
        lambda: restraint.setForceConstants([1.0, -1.0]),
        "Force constant at index 1 must be non-negative",
    )
    expect_exception(
        "setReferenceCoordinates rejects wrong coordinate length",
        ValueError,
        lambda: restraint.setReferenceCoordinates([[0.0, 0.0]]),
    )
    expect_invalid_argument(
        "setReferenceCoordinates rejects wrong atom count",
        lambda: restraint.setReferenceCoordinates([[0.0, 0.0, 0.0]]),
        "Reference-coordinate array size mismatch; expected 2, observed 1",
    )
    expect_invalid_argument(
        "setReferenceCoordinates rejects non-finite coordinate",
        lambda: restraint.setReferenceCoordinates(
            [[0.0, 0.0, 0.0], [0.0, math.inf, 0.0]]
        ),
        "Reference coordinate at atom index 1, Y component must be finite",
    )
    expect_invalid_argument(
        "setMasses rejects wrong length",
        lambda: restraint.setMasses([1.0]),
        "Mass array size mismatch; expected 2, observed 1",
    )
    expect_invalid_argument(
        "setMasses rejects non-finite mass",
        lambda: restraint.setMasses([1.0, math.inf]),
        "Mass at index 1 must be finite",
    )
    expect_invalid_argument(
        "setMasses rejects negative mass",
        lambda: restraint.setMasses([1.0, -1.0]),
        "Mass at index 1 must be non-negative",
    )
    restraint.setMasses([0.0, 1.0])
    expect_invalid_argument(
        "setBoxDimensions rejects wrong length",
        lambda: restraint.setBoxDimensions([50.0, 50.0]),
        "Box-dimension array size mismatch; expected 3, observed 2",
    )
    expect_invalid_argument(
        "setBoxDimensions rejects non-finite dimension",
        lambda: restraint.setBoxDimensions([50.0, math.inf, 50.0]),
        "Box dimension at index 1 must be finite",
    )
    expect_invalid_argument(
        "setBoxDimensions rejects zero dimension",
        lambda: restraint.setBoxDimensions([50.0, 0.0, 50.0]),
        "Box dimension at index 1 must be positive",
    )
    expect_invalid_argument(
        "setBoxDimensions rejects negative dimension",
        lambda: restraint.setBoxDimensions([50.0, -1.0, 50.0]),
        "Box dimension at index 1 must be positive",
    )
    restraint.setBoxDimensions(BOX_DIMENSIONS)
    expect_exception(
        "subscribe rejects non-ForceManager",
        TypeError,
        lambda: restraint._subscribe_to_force_manager(object()),  # type: ignore[arg-type]
    )
    expect_invalid_argument(
        "subscribe rejects empty force tag",
        lambda: fm.subscribe(restraint, ""),
        "Force tag must not be empty",
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
