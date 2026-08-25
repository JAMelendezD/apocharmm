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

STANDARD_COORDINATES: list[list[float]] = [
    [1.25, -2.5, 3.75],
    [0.0, 4.5, -6.25],
    [-7.125, 8.0, 9.5],
]

EXTENDED_COORDINATES: list[list[float]] = [
    [-10.123456789, 20.25, -30.5],
    [1000.0, -2000.5, 3000.75],
    [0.125, 0.25, 0.5],
]


def format_standard_crd_line(
    atom_number: int,
    residue_number: int,
    residue_name: str,
    atom_name: str,
    x: float,
    y: float,
    z: float,
) -> str:
    return (
        f"{atom_number:5d}"
        f"{residue_number:5d} "
        f"{residue_name:<4} "
        f"{atom_name:<4}"
        f"{x:10.5f}"
        f"{y:10.5f}"
        f"{z:10.5f}"
        f" SEG1"
        f"{residue_number:5d}"
        f"{0.0:10.5f}\n"
    )


def format_extended_crd_line(
    atom_number: int,
    residue_number: int,
    residue_name: str,
    atom_name: str,
    x: float,
    y: float,
    z: float,
) -> str:
    return (
        f"{atom_number:10d}"
        f"{residue_number:10d} "
        f"{residue_name:<8} "
        f"{atom_name:<8}  "
        f"{x:20.10f}"
        f"{y:20.10f}"
        f"{z:20.10f}"
        f" SEG1"
        f"{residue_number:10d}"
        f"{0.0:20.10f}\n"
    )


def write_standard_crd(path: Path, coordinates: list[list[float]]) -> None:
    contents: list[str] = [
        "* generated standard CharmmCrd Python API test\n",
        "*\n",
        f"{len(coordinates):10d}\n",
    ]

    for i, xyz in enumerate(coordinates):
        residue_number: int = 1 if i < 2 else 2
        residue_name: str = "RS1" if residue_number == 1 else "RS2"
        atom_name: str = f"A{i + 1}"
        contents.append(
            format_standard_crd_line(
                i + 1,
                residue_number,
                residue_name,
                atom_name,
                xyz[0],
                xyz[1],
                xyz[2],
            )
        )

    path.write_text("".join(contents), encoding="utf-8")

    return


def write_extended_crd(path: Path, coordinates: list[list[float]]) -> None:
    contents: list[str] = [
        "* generated extended CharmmCrd Python API test\n",
        "*\n",
        f"{len(coordinates):10d} EXT\n",
    ]

    for i, xyz in enumerate(coordinates):
        residue_number: int = 1 if i < 2 else 2
        residue_name: str = "RES1" if residue_number == 1 else "RES2"
        atom_name: str = f"ATOM{i + 1}"
        contents.append(
            format_extended_crd_line(
                i + 1,
                residue_number,
                residue_name,
                atom_name,
                xyz[0],
                xyz[1],
                xyz[2],
            )
        )

    path.write_text("".join(contents), encoding="utf-8")

    return


def write_malformed_crd(path: Path) -> None:
    path.write_text(
        "* generated malformed CharmmCrd Python API test\n"
        "*\n"
        f"{1:10d}\n"
        "short\n",
        encoding="utf-8",
    )
    return


def check_standard_crd(path: Path) -> None:
    print("Checking standard CharmmCrd parsing...")

    write_standard_crd(path, STANDARD_COORDINATES)

    crd = apo.CharmmCrd(str(path))

    apo_test.assert_equal("standard CharmmCrd.getNumAtoms", crd.getNumAtoms(), 3)
    apo_test.assert_nested_sequence_close(
        "standard CharmmCrd.getCoordinates",
        crd.getCoordinates(),
        STANDARD_COORDINATES,
        TOLERANCE,
    )

    return


def check_extended_crd(path: Path) -> None:
    print("Checking extended CharmmCrd parsing...")

    write_extended_crd(path, EXTENDED_COORDINATES)

    crd = apo.CharmmCrd(str(path))

    apo_test.assert_equal("extended CharmmCrd.getNumAtoms", crd.getNumAtoms(), 3)
    apo_test.assert_nested_sequence_close(
        "extended CharmmCrd.getCoordinates",
        crd.getCoordinates(),
        EXTENDED_COORDINATES,
        TOLERANCE,
    )

    return


def check_repository_nacl_pair(repo_root: Path) -> None:
    print("Checking repository nacl_pair.cor...")

    crd_path: str = apo_test.require_file(repo_root / "test/data/nacl_pair.cor")
    crd = apo.CharmmCrd(crd_path)

    expected_coordinates: list[list[float]] = [
        [0.0, 0.0, 0.0],
        [2.82, 2.82, 2.82],
    ]

    apo_test.assert_equal("nacl_pair CharmmCrd.getNumAtoms", crd.getNumAtoms(), 2)
    apo_test.assert_nested_sequence_close(
        "nacl_pair CharmmCrd.getCoordinates",
        crd.getCoordinates(),
        expected_coordinates,
        TOLERANCE,
    )

    return


def check_missing_crd(path: Path) -> None:
    print("Checking missing CharmmCrd file error path...")

    apo_test.remove_if_exists(path)

    apo_test.expect_apo_error(
        "CharmmCrd rejects a missing coordinate file",
        lambda: apo.CharmmCrd(str(path)),
        apo.APO_STATUS_RUNTIME_ERROR,
        f'Failed to open file "{path}"',
        "CharmmCrd construction",
    )

    return


def check_malformed_crd(path: Path) -> None:
    print("Checking malformed CharmmCrd file error path...")

    write_malformed_crd(path)

    apo_test.expect_apo_error(
        "CharmmCrd rejects a malformed coordinate file",
        lambda: apo.CharmmCrd(str(path)),
        apo.APO_STATUS_RUNTIME_ERROR,
        (
            "Coordinate record 1 is truncated in CHARMM coordinate file "
            f'"{path}" at line 4'
        ),
        "CharmmCrd construction",
    )

    return


def check_invalid_atom_count_crd(path: Path) -> None:
    print("Checking invalid CharmmCrd atom-count error path...")

    path.write_text(
        "* generated invalid-count CharmmCrd Python API test\n" "*\n" "BAD\n",
        encoding="utf-8",
    )

    apo_test.expect_apo_error(
        "CharmmCrd rejects an invalid atom count",
        lambda: apo.CharmmCrd(str(path)),
        apo.APO_STATUS_RUNTIME_ERROR,
        (
            'Invalid atom count value "BAD" in CHARMM coordinate file '
            f'"{path}" at line 3'
        ),
        "CharmmCrd construction",
    )

    return


def check_empty_path() -> None:
    print("Checking empty CharmmCrd path error path...")

    apo_test.expect_apo_error(
        "CharmmCrd rejects an empty coordinate-file path",
        lambda: apo.CharmmCrd(""),
        apo.APO_STATUS_INVALID_ARGUMENT,
        "apo_charmm_crd_create: path is NULL or empty",
        "CharmmCrd construction",
    )

    return


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")
    output_dir: Path = repo_root / "test/pytest"
    output_dir.mkdir(parents=True, exist_ok=True)

    standard_crd_path: Path = output_dir / "tmp_python_api_charmm_crd_standard.cor"
    extended_crd_path: Path = output_dir / "tmp_python_api_charmm_crd_extended.cor"
    malformed_crd_path: Path = output_dir / "tmp_python_api_charmm_crd_malformed.cor"
    invalid_count_crd_path: Path = (
        output_dir / "tmp_python_api_charmm_crd_invalid_count.cor"
    )
    missing_crd_path: Path = output_dir / "tmp_python_api_charmm_crd_missing.cor"

    generated_files: tuple[Path, ...] = (
        standard_crd_path,
        extended_crd_path,
        malformed_crd_path,
        invalid_count_crd_path,
        missing_crd_path,
    )

    for path in generated_files:
        apo_test.remove_if_exists(path)

    try:
        check_standard_crd(standard_crd_path)
        check_extended_crd(extended_crd_path)
        check_repository_nacl_pair(repo_root)
        check_missing_crd(missing_crd_path)
        check_malformed_crd(malformed_crd_path)
        check_invalid_atom_count_crd(invalid_count_crd_path)
        check_empty_path()
    finally:
        print("Cleaning up CharmmCrd Python API test files...")
        for path in generated_files:
            apo_test.remove_if_exists(path)

    print("\033[32m" + "PASS: CharmmCrd Python API tests completed." + "\033[0m")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
