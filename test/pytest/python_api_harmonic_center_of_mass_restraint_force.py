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

import apo_test_helpers as apo_test

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
    prm_path: str = apo_test.require_file(
        apo_test.get_toppar_dir() / "toppar_water_ions.str"
    )
    psf_path: str = apo_test.require_file(apo_test.get_data_dir() / "nacl_pair.psf")
    crd_path: str = apo_test.require_file(apo_test.get_data_dir() / "nacl_pair.cor")

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

    apo_test.assert_equal(
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

    apo_test.expect_invalid_argument(
        "ForceManager rejects duplicate HarmonicCenterOfMassRestraintForce subscription",
        lambda: fm.subscribe(restraint),
        "Force is already subscribed to this ForceManager",
        expected_context="ForceManager.subscribe(HarmonicCenterOfMassRestraintForce)",
    )

    integrator = apo.CudaLangevinThermostatIntegrator(TIME_STEP)
    integrator.setReferenceTemperature(TEMPERATURE)
    integrator.setThermostatFriction(0.0)
    integrator.setThermostatRngSeed(RANDOM_SEED)
    integrator.setCharmmContext(ctx)

    integrator.propagate(10)

    apo_test.assert_finite_temperature(
        "post harmonic center-of-mass restraint propagation", ctx.computeTemperature()
    )

    fm.unsubscribe(restraint)

    apo_test.expect_invalid_argument(
        "ForceManager rejects missing HarmonicCenterOfMassRestraintForce unsubscription",
        lambda: fm.unsubscribe(restraint),
        "Force is not subscribed to this ForceManager",
        expected_context=(
            "ForceManager.unsubscribe(HarmonicCenterOfMassRestraintForce)"
        ),
    )

    integrator.close()
    restraint.close()
    fm.close()
    ctx.close()

    return


def check_validation() -> None:
    print("Checking HarmonicCenterOfMassRestraintForce validation...")

    _, psf, _, fm, _ = create_system()
    num_atoms: int = psf.getNumAtoms()
    restraint = apo.HarmonicCenterOfMassRestraintForce(num_atoms)
    selector = apo.AtomSelector(psf)
    selection = selector.select("all")

    apo_test.expect_exception(
        "HarmonicCenterOfMassRestraintForce rejects non-int num_atoms",
        TypeError,
        lambda: apo.HarmonicCenterOfMassRestraintForce(2.0),  # type: ignore[arg-type]
    )
    apo_test.expect_invalid_argument(
        "HarmonicCenterOfMassRestraintForce rejects zero num_atoms",
        lambda: apo.HarmonicCenterOfMassRestraintForce(0),
        "Atom count must be positive; observed 0",
        expected_context="HarmonicCenterOfMassRestraintForce construction",
    )
    apo_test.expect_invalid_argument(
        "HarmonicCenterOfMassRestraintForce rejects negative num_atoms",
        lambda: apo.HarmonicCenterOfMassRestraintForce(-1),
        "Atom count must be positive; observed -1",
        expected_context="HarmonicCenterOfMassRestraintForce construction",
    )
    apo_test.expect_exception(
        "HarmonicCenterOfMassRestraintForce rejects out-of-range num_atoms",
        ValueError,
        lambda: apo.HarmonicCenterOfMassRestraintForce(2**31),
    )
    apo_test.expect_exception(
        "setSelection rejects non-AtomSelection",
        TypeError,
        lambda: restraint.setSelection(object()),  # type: ignore[arg-type]
    )

    force_constant_error = apo_test.expect_invalid_argument(
        "setForceConstant rejects negative force constant",
        lambda: restraint.setForceConstant(-1.0),
        "Force constant must be non-negative",
        expected_context=(
            "HarmonicCenterOfMassRestraintForce.setForceConstant(force_constant)"
        ),
    )
    apo_test.assert_equal(
        "setForceConstant rendered native function occurrence count",
        force_constant_error.message.count(
            "apo_harmonic_center_of_mass_restraint_force_set_force_constant"
        ),
        1,
    )
    apo_test.expect_invalid_argument(
        "setForceConstant rejects non-finite force constant",
        lambda: restraint.setForceConstant(math.inf),
        "Force constant must be finite",
        expected_context=(
            "HarmonicCenterOfMassRestraintForce.setForceConstant(force_constant)"
        ),
    )
    restraint.setForceConstant(0.0)

    apo_test.expect_invalid_argument(
        "setReferencePosition rejects wrong coordinate length",
        lambda: restraint.setReferencePosition([0.0, 0.0]),
        "Reference-position array size mismatch; expected 3, observed 2",
        expected_context=(
            "HarmonicCenterOfMassRestraintForce.setReferencePosition(reference_position)"
        ),
    )
    apo_test.expect_invalid_argument(
        "setReferencePosition rejects wrong reference mask length",
        lambda: restraint.setReferencePosition([0.0, 0.0, 0.0], [1, 1]),
        "Reference-mask array size mismatch; expected 3, observed 2",
    )
    apo_test.expect_invalid_argument(
        "setReferencePosition rejects inactive reference mask",
        lambda: restraint.setReferencePosition([0.0, 0.0, 0.0], [0, 0, 0]),
        "Reference mask must activate at least one coordinate",
    )
    apo_test.expect_invalid_argument(
        "setReferencePosition rejects invalid reference mask value",
        lambda: restraint.setReferencePosition([0.0, 0.0, 0.0], [1, 2, 1]),
        "Reference mask at index 1 must be 0 or 1; observed 2",
    )
    apo_test.expect_invalid_argument(
        "setReferencePosition rejects non-finite coordinate",
        lambda: restraint.setReferencePosition([0.0, math.inf, 0.0]),
        "Reference position at index 1 must be finite",
    )

    apo_test.expect_invalid_argument(
        "setReferenceDistance rejects negative reference distance",
        lambda: restraint.setReferenceDistance(-1.0),
        "Reference distance must be non-negative",
        expected_context=(
            "HarmonicCenterOfMassRestraintForce.setReferenceDistance(reference_distance)"
        ),
    )
    apo_test.expect_invalid_argument(
        "setReferenceDistance rejects non-finite reference distance",
        lambda: restraint.setReferenceDistance(math.inf),
        "Reference distance must be finite",
    )

    apo_test.expect_invalid_argument(
        "setMasses rejects wrong length",
        lambda: restraint.setMasses([1.0]),
        f"Mass array size mismatch; expected {num_atoms}, observed 1",
        expected_context="HarmonicCenterOfMassRestraintForce.setMasses(masses)",
    )

    nonfinite_masses: list[float] = [1.0 for _ in range(num_atoms)]
    nonfinite_masses[1] = math.inf
    apo_test.expect_invalid_argument(
        "setMasses rejects non-finite mass",
        lambda: restraint.setMasses(nonfinite_masses),
        "Mass at index 1 must be finite",
    )

    negative_masses: list[float] = [1.0 for _ in range(num_atoms)]
    negative_masses[1] = -1.0
    apo_test.expect_invalid_argument(
        "setMasses rejects negative mass",
        lambda: restraint.setMasses(negative_masses),
        "Mass at index 1 must be non-negative",
    )

    zero_masses: list[float] = [0.0 for _ in range(num_atoms)]
    apo_test.expect_invalid_argument(
        "setMasses rejects zero selected total weight",
        lambda: restraint.setMasses(zero_masses),
        "Selected atoms must have positive total weight",
    )

    restraint.setMassWeighting(False)
    apo_test.expect_invalid_argument(
        "setMassWeighting rejects zero selected total weight",
        lambda: restraint.setMassWeighting(True),
        "Selected atoms must have positive total weight",
        expected_context=("HarmonicCenterOfMassRestraintForce.setMassWeighting(flag)"),
    )
    restraint.setMasses(psf.getMasses())

    apo_test.expect_exception(
        "subscribe rejects non-ForceManager",
        TypeError,
        lambda: restraint._subscribe_to_force_manager(object()),  # type: ignore[arg-type]
    )
    apo_test.expect_invalid_argument(
        "subscribe rejects empty force tag",
        lambda: fm.subscribe(restraint, ""),
        "Force tag must not be empty",
        expected_context=("ForceManager.subscribe(HarmonicCenterOfMassRestraintForce)"),
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

    apo_test.expect_exception(
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
