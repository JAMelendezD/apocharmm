# BEGINLICENSE
#
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

import apo_test_helpers as apo_test

TOLERANCE: float = 1.0e-10
LINEAR_CHAIN_PSF: str = """PSF

       1 !NTITLE
 REMARKS generated CharmmPsf Python API unit test
       4 !NATOM
       1 SEG1     1 RES1 C1   CT1   -0.300000  12.0110           0
       2 SEG1     1 RES1 H1   HA1    0.100000   1.0080           0
       3 SEG1     2 RES2 C2   CT2    0.100000  12.0110           0
       4 SEG1     2 RES2 H2   HA2    0.100000   1.0080           0
       3 !NBOND: bonds
       1       2       2       3       3       4
       2 !NTHETA: angles
       1       2       3       2       3       4
       1 !NPHI: dihedrals
       1       2       3       4
       0 !NIMPHI: impropers
       0 !NDON: donors
       0 !NACC: acceptors
       0 !NCRTERM: cross-terms
"""
MALFORMED_PSF: str = """PSF

       1 !NTITLE
 REMARKS generated malformed CharmmPsf Python API unit test
       1 !NATOM
       1 SEG1     1 RES1 C1   CT1    0.000000  12.0110           0
"""
EXPECTED_SEGMENT_IDENTIFIERS: list[str] = ["SEG1", "SEG1", "SEG1", "SEG1"]
EXPECTED_RESIDUE_IDENTIFIERS: list[int] = [1, 1, 2, 2]
EXPECTED_RESIDUE_NAMES: list[str] = ["RES1", "RES1", "RES2", "RES2"]
EXPECTED_ATOM_NAMES: list[str] = ["C1", "H1", "C2", "H2"]
EXPECTED_ATOM_TYPES: list[str] = ["CT1", "HA1", "CT2", "HA2"]
EXPECTED_CHARGES: list[float] = [-0.3, 0.1, 0.1, 0.1]
EXPECTED_MASSES: list[float] = [12.011, 1.008, 12.011, 1.008]


def check_generated_linear_chain_psf(psf_path: Path) -> None:
    print("Checking generated CharmmPsf file...")

    apo_test.write_text_file(psf_path, LINEAR_CHAIN_PSF)
    psf = apo.CharmmPsf(str(psf_path))

    apo_test.assert_equal("CharmmPsf.getNumAtoms", psf.getNumAtoms(), 4)
    apo_test.assert_equal("CharmmPsf.getNumBonds", psf.getNumBonds(), 3)
    apo_test.assert_equal("CharmmPsf.getNumAngles", psf.getNumAngles(), 2)
    apo_test.assert_equal("CharmmPsf.getNumDihedrals", psf.getNumDihedrals(), 1)
    apo_test.assert_equal("CharmmPsf.getNumImpropers", psf.getNumImpropers(), 0)
    apo_test.assert_equal("CharmmPsf.getNumCrossTerms", psf.getNumCrossTerms(), 0)

    apo_test.assert_equal(
        "CharmmPsf.getSegmentIdentifiers",
        psf.getSegmentIdentifiers(),
        EXPECTED_SEGMENT_IDENTIFIERS,
    )
    apo_test.assert_equal(
        "CharmmPsf.getResidueIdentifiers",
        psf.getResidueIdentifiers(),
        EXPECTED_RESIDUE_IDENTIFIERS,
    )
    apo_test.assert_equal(
        "CharmmPsf.getResidueNames",
        psf.getResidueNames(),
        EXPECTED_RESIDUE_NAMES,
    )
    apo_test.assert_equal(
        "CharmmPsf.getAtomNames",
        psf.getAtomNames(),
        EXPECTED_ATOM_NAMES,
    )
    apo_test.assert_equal(
        "CharmmPsf.getAtomTypes",
        psf.getAtomTypes(),
        EXPECTED_ATOM_TYPES,
    )

    apo_test.assert_sequence_close(
        "CharmmPsf.getCharges",
        psf.getCharges(),
        EXPECTED_CHARGES,
        TOLERANCE,
    )
    apo_test.assert_sequence_close(
        "CharmmPsf.getMasses",
        psf.getMasses(),
        EXPECTED_MASSES,
        TOLERANCE,
    )
    apo_test.assert_close(
        "CharmmPsf.getNetCharge",
        psf.getNetCharge(),
        sum(EXPECTED_CHARGES),
        TOLERANCE,
    )
    apo_test.assert_close(
        "CharmmPsf.getTotalMass",
        psf.getTotalMass(),
        sum(EXPECTED_MASSES),
        TOLERANCE,
    )
    apo_test.assert_equal(
        "CharmmPsf.getFileName basename",
        Path(psf.getFileName()).name,
        psf_path.name,
    )

    return


def check_repository_nacl_pair_psf(repo_root: Path) -> None:
    print("Checking repository nacl_pair.psf...")

    psf_path: str = apo_test.require_file(repo_root / "test/data/nacl_pair.psf")
    psf = apo.CharmmPsf(psf_path)

    apo_test.assert_equal("nacl_pair getNumAtoms", psf.getNumAtoms(), 2)
    apo_test.assert_equal("nacl_pair getNumBonds", psf.getNumBonds(), 0)
    apo_test.assert_equal("nacl_pair getNumAngles", psf.getNumAngles(), 0)
    apo_test.assert_equal("nacl_pair getNumDihedrals", psf.getNumDihedrals(), 0)
    apo_test.assert_equal("nacl_pair getNumImpropers", psf.getNumImpropers(), 0)
    apo_test.assert_equal("nacl_pair getNumCrossTerms", psf.getNumCrossTerms(), 0)

    apo_test.assert_equal(
        "nacl_pair getSegmentIdentifiers",
        psf.getSegmentIdentifiers(),
        ["NACL", "NACL"],
    )
    apo_test.assert_equal(
        "nacl_pair getResidueIdentifiers",
        psf.getResidueIdentifiers(),
        [1, 2],
    )
    apo_test.assert_equal(
        "nacl_pair getResidueNames", psf.getResidueNames(), ["SOD", "CLA"]
    )
    apo_test.assert_equal("nacl_pair getAtomNames", psf.getAtomNames(), ["SOD", "CLA"])
    apo_test.assert_equal("nacl_pair getAtomTypes", psf.getAtomTypes(), ["SOD", "CLA"])

    expected_charges: list[float] = [1.0, -1.0]
    expected_masses: list[float] = [22.9898, 35.45]

    apo_test.assert_sequence_close(
        "nacl_pair getCharges",
        psf.getCharges(),
        expected_charges,
        TOLERANCE,
    )
    apo_test.assert_sequence_close(
        "nacl_pair getMasses",
        psf.getMasses(),
        expected_masses,
        TOLERANCE,
    )
    apo_test.assert_close("nacl_pair getNetCharge", psf.getNetCharge(), 0.0, TOLERANCE)
    apo_test.assert_close(
        "nacl_pair getTotalMass",
        psf.getTotalMass(),
        sum(expected_masses),
        TOLERANCE,
    )
    apo_test.assert_equal(
        "nacl_pair getFileName basename",
        Path(psf.getFileName()).name,
        "nacl_pair.psf",
    )

    return


def check_malformed_psf_file(psf_path: Path) -> None:
    print("Checking malformed CharmmPsf file error path...")

    apo_test.write_text_file(psf_path, MALFORMED_PSF)

    apo_test.expect_apo_error(
        "CharmmPsf rejects a malformed PSF file",
        lambda: apo.CharmmPsf(str(psf_path)),
        apo.APO_STATUS_RUNTIME_ERROR,
        f'Could not find BOND section in PSF "{psf_path}"',
        "CharmmPsf construction",
    )

    return


def check_missing_psf_file(psf_path: Path) -> None:
    print("Checking missing CharmmPsf file error path...")

    apo_test.remove_if_exists(psf_path)

    apo_test.expect_apo_error(
        "CharmmPsf rejects a missing PSF file",
        lambda: apo.CharmmPsf(str(psf_path)),
        apo.APO_STATUS_RUNTIME_ERROR,
        f'Failed to open file "{psf_path}"',
        "CharmmPsf construction",
    )

    apo_test.expect_apo_error(
        "CharmmPsf rejects an empty PSF path",
        lambda: apo.CharmmPsf(""),
        apo.APO_STATUS_INVALID_ARGUMENT,
        "apo_charmm_psf_create: PSF path is NULL or empty",
        "CharmmPsf construction",
    )

    return


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")
    output_dir: Path = repo_root / "test/pytest"
    output_dir.mkdir(parents=True, exist_ok=True)

    generated_psf_path: Path = output_dir / "tmp_python_api_charmm_psf_linear_chain.psf"
    malformed_psf_path: Path = output_dir / "tmp_python_api_charmm_psf_malformed.psf"
    missing_psf_path: Path = output_dir / "tmp_python_api_charmm_psf_missing.psf"

    generated_files: tuple[Path, ...] = (
        generated_psf_path,
        malformed_psf_path,
        missing_psf_path,
    )

    for path in generated_files:
        apo_test.remove_if_exists(path)

    try:
        check_generated_linear_chain_psf(generated_psf_path)
        check_repository_nacl_pair_psf(repo_root)
        check_malformed_psf_file(malformed_psf_path)
        check_missing_psf_file(missing_psf_path)
    finally:
        print("Cleaning up CharmmPsf Python API test files...")
        for path in generated_files:
            apo_test.remove_if_exists(path)

    print("\033[32m" + "PASS: CharmmPsf Python API tests completed." + "\033[0m")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
