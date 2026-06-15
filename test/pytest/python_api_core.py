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
    assert_equal,
    assert_sequence_close,
    require_file,
)

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]
EXPECTED_COORDINATES: list[list[float]] = [[0.0, 0.0, 0.0], [2.82, 2.82, 2.82]]
EXPECTED_CHARGES: list[float] = [1.0, -1.0]
EXPECTED_MASSES: list[float] = [22.9898, 35.45]
TOLERANCE: float = 1.0e-10


def assert_coordinates_close(
    label: str,
    observed: list[list[float]],
    expected: list[list[float]],
    tolerance: float,
) -> None:
    assert_equal(f"{label} atom count", len(observed), len(expected))

    for i in range(len(observed)):
        assert_sequence_close(f"{label}[{i}]", observed[i], expected[i], tolerance)

    return


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")

    parameter_path: str = require_file(repo_root / "test/data/toppar_water_ions.str")
    psf_path: str = require_file(repo_root / "test/data/nacl_pair.psf")
    crd_path: str = require_file(repo_root / "test/data/nacl_pair.cor")

    print("Loading input files...")
    prm = apo.CharmmParameters(parameter_path)
    prm_from_list = apo.CharmmParameters([parameter_path])
    psf = apo.CharmmPsf(psf_path)
    crd = apo.CharmmCrd(crd_path)

    print("Checking CharmmPsf scalar getters...")
    assert_equal("CharmmPsf.getNumAtoms", psf.getNumAtoms(), 2)
    assert_equal("CharmmPsf.getNumBonds", psf.getNumBonds(), 0)
    assert_equal("CharmmPsf.getNumAngles", psf.getNumAngles(), 0)
    assert_equal("CharmmPsf.getNumDihedrals", psf.getNumDihedrals(), 0)
    assert_equal("CharmmPsf.getNumImpropers", psf.getNumImpropers(), 0)
    assert_equal("CharmmPsf.getNumCrossTerms", psf.getNumCrossTerms(), 0)
    assert_close(
        "CharmmPsf.getNetCharge", psf.getNetCharge(), sum(EXPECTED_CHARGES), TOLERANCE
    )
    assert_close(
        "CharmmPsf.getTotalMass", psf.getTotalMass(), sum(EXPECTED_MASSES), TOLERANCE
    )
    assert_equal(
        "CharmmPsf.getFileName basename", Path(psf.getFileName()).name, "nacl_pair.psf"
    )

    print("Checking CharmmPsf array getters...")
    assert_equal(
        "CharmmPsf.getSegmentIdentifiers", psf.getSegmentIdentifiers(), ["NACL", "NACL"]
    )
    assert_equal("CharmmPsf.getResidueIdentifiers", psf.getResidueIdentifiers(), [1, 2])
    assert_equal("CharmmPsf.getResidueNames", psf.getResidueNames(), ["SOD", "CLA"])
    assert_equal("CharmmPsf.getAtomNames", psf.getAtomNames(), ["SOD", "CLA"])
    assert_equal("CharmmPsf.getAtomTypes", psf.getAtomTypes(), ["SOD", "CLA"])
    assert_sequence_close(
        "CharmmPsf.getCharges", psf.getCharges(), EXPECTED_CHARGES, TOLERANCE
    )
    assert_sequence_close(
        "CharmmPsf.getMasses", psf.getMasses(), EXPECTED_MASSES, TOLERANCE
    )

    print("Checking CharmmCrd getters...")
    assert_equal("CharmmCrd.getNumAtoms", crd.getNumAtoms(), psf.getNumAtoms())
    assert_coordinates_close(
        "CharmmCrd.getCoordinates",
        crd.getCoordinates(),
        EXPECTED_COORDINATES,
        TOLERANCE,
    )

    print("Checking ForceManager and CharmmContext getters/setters...")
    fm = apo.ForceManager(psf, prm)
    fm_from_list = apo.ForceManager(psf, prm_from_list)
    fm.setBoxDimensions(BOX_DIMENSIONS)
    fm.setPeriodicBoundaryCondition(apo.PeriodicBoundaryCondition.P1)
    assert_equal(
        "ForceManager.getPeriodicBoundaryCondition",
        fm.getPeriodicBoundaryCondition(),
        apo.PeriodicBoundaryCondition.P1,
    )

    ctx = apo.CharmmContext(fm)
    ctx.setCoordinates(crd)
    ctx.setPeriodicBoundaryCondition(apo.PeriodicBoundaryCondition.P1)
    ctx.useHolonomicConstraints(False)
    assert_equal("CharmmContext.getNumAtoms", ctx.getNumAtoms(), psf.getNumAtoms())
    assert_sequence_close(
        "CharmmContext.getBoxDimensions",
        ctx.getBoxDimensions(),
        BOX_DIMENSIONS,
        TOLERANCE,
    )
    assert_equal(
        "CharmmContext.getPeriodicBoundaryCondition",
        ctx.getPeriodicBoundaryCondition(),
        apo.PeriodicBoundaryCondition.P1,
    )

    coordinates_charges: list[list[float]] = ctx.getCoordinatesCharges()
    assert_equal(
        "CharmmContext.getCoordinatesCharges atom count",
        len(coordinates_charges),
        psf.getNumAtoms(),
    )
    for i in range(psf.getNumAtoms()):
        expected_coordinate_charge: list[float] = [
            EXPECTED_COORDINATES[i][0],
            EXPECTED_COORDINATES[i][1],
            EXPECTED_COORDINATES[i][2],
            EXPECTED_CHARGES[i],
        ]
        assert_sequence_close(
            f"CharmmContext.getCoordinatesCharges[{i}]",
            coordinates_charges[i],
            expected_coordinate_charge,
            TOLERANCE,
        )

    print("Checking AtomSelector and AtomSelection...")
    selector = apo.AtomSelector(psf)

    all_selection = selector.select("all")
    assert_equal(
        "All selection num atoms", all_selection.getNumAtoms(), psf.getNumAtoms()
    )
    assert_equal("All selection num selected", all_selection.getNumSelected(), 2)
    assert_equal("All selection atom indices", all_selection.getAtomIndices(), [0, 1])
    assert_equal("All selection contains atom 0", all_selection.contains(0), True)
    assert_equal("All selection contains atom 1", all_selection.contains(1), True)

    sodium_selection = selector.select("type SOD")
    assert_equal(
        "Sodium selection num atoms", sodium_selection.getNumAtoms(), psf.getNumAtoms()
    )
    assert_equal("Sodium selection num selected", sodium_selection.getNumSelected(), 1)
    assert_equal(
        "Sodium selection atom indices", sodium_selection.getAtomIndices(), [0]
    )
    assert_equal("Sodium selection contains atom 0", sodium_selection.contains(0), True)
    assert_equal(
        "Sodium selection contains atom 1", sodium_selection.contains(1), False
    )

    chloride_selection = selector.select("type CLA")
    assert_equal(
        "Chloride selection num atoms",
        chloride_selection.getNumAtoms(),
        psf.getNumAtoms(),
    )
    assert_equal(
        "Chloride selection num selected", chloride_selection.getNumSelected(), 1
    )
    assert_equal(
        "Chloride selection atom indices", chloride_selection.getAtomIndices(), [1]
    )
    assert_equal(
        "Chloride selection contains atom 0", chloride_selection.contains(0), False
    )
    assert_equal(
        "Chloride selection contains atom 1", chloride_selection.contains(1), True
    )

    print(
        "PASS: Core Python API objects, topology/coordinate getters, managers, context, enums, and selections all completed."
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
