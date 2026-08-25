# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from __future__ import annotations

import ctypes
from pathlib import Path
import sys

import apocharmm as apo

import apo_test_helpers as apo_test

NUM_ATOMS: int = 10
TEST_PSF_TEXT: str = """PSF

       1 !NTITLE
 REMARKS generated AtomSelector Python API test PSF
      10 !NATOM
       1 SEG1     1 ALA  N    NH1   -0.300000  14.0070           0
       2 SEG1     1 ALA  CA   CT1    0.100000  12.0110           0
       3 SEG1     1 ALA  CB   CT2    0.000000  12.0110           0
       4 SEG1     1 ALA  HB1  HA     0.100000   1.0080           0
       5 SEG2     2 GLY  N    NH1   -0.200000  14.0070           0
       6 SEG2     2 GLY  CA   CT1    0.100000  12.0110           0
       7 SEG2     2 GLY  HA1  HA     0.100000   1.0080           0
       8 WAT      3 TIP3 OH2  OT    -0.834000  15.9994           0
       9 WAT      3 TIP3 H1   HT     0.417000   1.0080           0
      10 WAT      3 TIP3 H2   HT     0.417000   1.0080           0
       7 !NBOND: bonds
       1       2       2       3       3       4       5       6
       6       7       8       9       8      10
       4 !NTHETA: angles
       1       2       3       2       3       4       5       6       7
       9       8      10
       1 !NPHI: dihedrals
       1       2       3       4
       0 !NIMPHI: impropers
       0 !NDON: donors
       0 !NACC: acceptors
       0 !NCRTERM: cross-terms
"""


def assert_selection(
    label: str,
    selector: apo.AtomSelector,
    selection_string: str,
    expected_indices: list[int],
) -> apo.AtomSelection:
    selection = selector.select(selection_string)

    apo_test.assert_equal(f"{label} num atoms", selection.getNumAtoms(), NUM_ATOMS)
    apo_test.assert_equal(
        f"{label} num selected", selection.getNumSelected(), len(expected_indices)
    )
    apo_test.assert_equal(
        f"{label} atom indices", selection.getAtomIndices(), expected_indices
    )

    expected_index_set: set[int] = set(expected_indices)
    for i in range(NUM_ATOMS):
        apo_test.assert_equal(
            f"{label} contains {i}", selection.contains(i), i in expected_index_set
        )

    return selection


def check_basic_selections(selector: apo.AtomSelector) -> None:
    print("Checking basic atom selections...")

    assert_selection("none", selector, "none", [])
    assert_selection("all", selector, "all", list(range(NUM_ATOMS)))

    assert_selection("bynum range", selector, "bynum 2:4", [1, 2, 3])
    assert_selection("reversed bynum range", selector, "bynu 4:2", [1, 2, 3])
    assert_selection("clamped bynum range", selector, "bynum -5:2", [0, 1])
    assert_selection("out-of-range bynum range", selector, "bynum 20:25", [])

    assert_selection("type CA", selector, "type CA", [1, 5])
    assert_selection("case-insensitive type", selector, "type ca", [1, 5])
    assert_selection("chemical CT1", selector, "chemical CT1", [1, 5])
    assert_selection("chemical CT%", selector, "chemical CT%", [1, 2, 5])

    assert_selection("segid SEG1", selector, "segid SEG1", [0, 1, 2, 3])
    assert_selection("segi SEG2", selector, "segi SEG2", [4, 5, 6])
    assert_selection("resid 1", selector, "resid 1", [0, 1, 2, 3])
    assert_selection("resi 1:2", selector, "resi 1:2", [0, 1, 2, 3, 4, 5, 6])
    assert_selection("resname GLY", selector, "resname GLY", [4, 5, 6])
    assert_selection("resn TIP3", selector, "resn TIP3", [7, 8, 9])
    assert_selection("atom SEG1 1 CA", selector, "atom SEG1 1 CA", [1])

    return


def check_wildcard_selections(selector: apo.AtomSelector) -> None:
    print("Checking wildcard atom selections...")

    assert_selection("type H*", selector, "type H*", [3, 6, 8, 9])
    assert_selection("type H+", selector, "type H+", [8, 9])
    assert_selection("type H#", selector, "type H#", [8, 9])
    assert_selection("type H%", selector, "type H%", [8, 9])
    assert_selection("type C%", selector, "type C%", [1, 2, 5])
    assert_selection("type C*", selector, "type C*", [1, 2, 5])

    return


def check_logical_and_prefix_selections(selector: apo.AtomSelector) -> None:
    print("Checking logical and prefix atom selections...")

    assert_selection("and", selector, "type CA .and. resn ALA", [1])
    assert_selection(
        "or and precedence",
        selector,
        "type CA .or. type N .and. resn ALA",
        [0, 1, 5],
    )
    assert_selection(
        "parentheses change precedence",
        selector,
        "(type CA .or. type N) .and. resn ALA",
        [0, 1],
    )
    assert_selection("not", selector, ".not. type H*", [0, 1, 2, 4, 5, 7])

    assert_selection("by residue", selector, ".byres. type HB1", [0, 1, 2, 3])
    assert_selection(
        "by residue with parentheses",
        selector,
        ".byres. (type CA .and. resn GLY)",
        [4, 5, 6],
    )
    assert_selection("by group", selector, ".bygroup. type HA1", [4, 5, 6])
    assert_selection("water by group", selector, ".bygroup. type H1", [7, 8, 9])

    assert_selection("bonded atom", selector, ".bonded. (atom SEG1 1 CA)", [0, 2])
    assert_selection("bonded bynum", selector, ".bonded. bynum 8", [8, 9])

    return


def check_validation(psf: apo.CharmmPsf, selector: apo.AtomSelector) -> None:
    print("Checking AtomSelector and AtomSelection validation...")

    apo_test.expect_exception(
        "AtomSelector rejects non-handle construction",
        TypeError,
        lambda: apo.AtomSelector(object()),  # type: ignore[arg-type]
    )
    apo_test.expect_exception(
        "AtomSelector rejects NULL-handle construction",
        TypeError,
        lambda: apo.AtomSelector(ctypes.c_void_p()),
    )
    apo_test.expect_exception(
        "AtomSelector.select rejects non-string selection",
        TypeError,
        lambda: selector.select(1),
    )
    invalid_selections: list[tuple[str, str]] = [
        ("", "selection_string is NULL or empty"),
        (".and. type CA", "Expected an atom selection at position 0"),
        (
            "type CA resn ALA",
            "Expected .AND., .OR., ')', or end of selection at position 8",
        ),
        ("(type CA .or. type N", "Found '(' without matching ')'"),
        (
            "type CA .or.",
            "Selection ended while expecting an atom selection at position 12",
        ),
        ("type", "Expected selection value after type at position 4"),
        ("all)", "Found ')' without matching '(' at position 3"),
        ("bynu A:C", "BYNU range requires integer atom numbers"),
        (".", "Unterminated dotted selection operator at position 0"),
    ]

    for selection_string, expected_diagnostic in invalid_selections:
        apo_test.expect_invalid_argument(
            f"AtomSelector reports invalid selection {selection_string!r}",
            lambda selection_string=selection_string: selector.select(selection_string),
            expected_diagnostic,
        )

    selection_error: apo.ApoCharmmError = apo_test.expect_invalid_argument(
        "AtomSelector reports unknown dotted operator",
        lambda: selector.select(".around. type CA"),
        'Unknown dotted atom selection operator ".around."',
    )
    apo_test.assert_equal(
        "AtomSelector.select Python operation context",
        selection_error.context,
        "AtomSelector.select(selection_string)",
    )
    apo_test.assert_equal(
        "AtomSelector.select rendered context occurrence count",
        selection_error.message.count(selection_error.context),
        1,
    )
    apo_test.assert_equal(
        "AtomSelector.select rendered native function occurrence count",
        selection_error.message.count("apo_atom_selector_select"),
        1,
    )
    apo_test.assert_equal(
        "AtomSelector.select rendered failed text count",
        selection_error.message.count("failed"),
        0,
    )

    selection: apo.AtomSelection = assert_selection(
        "validation selection", selector, "type CA", [1, 5]
    )

    apo_test.expect_exception(
        "AtomSelection.contains rejects bool atom index",
        TypeError,
        lambda: selection.contains(True),
    )
    apo_test.expect_exception(
        "AtomSelection.contains rejects non-int atom index",
        TypeError,
        lambda: selection.contains(1.0),
    )
    apo_test.expect_invalid_argument(
        "AtomSelection.contains rejects negative atom index",
        lambda: selection.contains(-1),
        f"Atom index is out of range; expected [0, {NUM_ATOMS}), observed -1",
    )
    apo_test.expect_exception(
        "AtomSelection.contains rejects atom index below int range",
        ValueError,
        lambda: selection.contains(-(2**31) - 1),
    )
    apo_test.expect_exception(
        "AtomSelection.contains rejects atom index above int range",
        ValueError,
        lambda: selection.contains(2**31),
    )

    contains_error = apo_test.expect_invalid_argument(
        "AtomSelection.contains reports out-of-range atom index",
        lambda: selection.contains(psf.getNumAtoms()),
        f"Atom index is out of range; expected [0, {NUM_ATOMS}), observed {NUM_ATOMS}",
    )
    apo_test.assert_equal(
        "AtomSelection.contains Python operation context",
        contains_error.context,
        "AtomSelection.contains(atom_index)",
    )
    apo_test.assert_equal(
        "AtomSelection.contains rendered context occurrence count",
        contains_error.message.count(contains_error.context),
        1,
    )
    apo_test.assert_equal(
        "AtomSelection.contains rendered native function occurrence count",
        contains_error.message.count("apo_atom_selection_contains"),
        1,
    )
    apo_test.assert_equal(
        "AtomSelection.contains rendered failed text count",
        contains_error.message.count("failed"),
        0,
    )

    return


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")
    output_dir: Path = repo_root / "test/pytest"
    output_dir.mkdir(parents=True, exist_ok=True)

    psf_path: Path = output_dir / "tmp_python_api_atom_selector.psf"
    apo_test.remove_if_exists(psf_path)

    try:
        psf_path.write_text(TEST_PSF_TEXT, encoding="utf-8")

        psf = apo.CharmmPsf(str(psf_path))
        selector = apo.AtomSelector(psf)

        check_basic_selections(selector)
        check_wildcard_selections(selector)
        check_logical_and_prefix_selections(selector)
        check_validation(psf, selector)
    finally:
        print("Cleaning up AtomSelector and AtomSelection Python API test files...")
        apo_test.remove_if_exists(psf_path)

    print(
        "\033[32m"
        + "PASS: AtomSelector and AtomSelection Python API tests completed."
        + "\033[0m"
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
