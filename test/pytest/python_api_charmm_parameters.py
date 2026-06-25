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
    write_text_file,
    require_file,
    remove_if_exists,
    expect_exception,
)

PARAMETER_TEXT: str = """* generated CharmmParameters Python API test file
*

NONBONDED
SOD      0.0      -0.0469      1.41075
CLA      0.0      -0.1500      2.27000

NBFIX
SOD      CLA      -0.083875    3.731

END
"""

SUPPLEMENTAL_PARAMETER_TEXT: str = """* generated supplemental CharmmParameters Python API test file
*

BONDS

END
"""

MALFORMED_PARAMETER_TEXT: str = """* generated malformed CharmmParameters Python API test file
*

NONBONDED
SOD      0.0      -0.0469      1.41075      99.0

END
"""

BOX_DIMENSIONS: list[float] = [50.0, 50.0, 50.0]


def check_single_file_construction(parameter_path: Path, repo_root: Path) -> None:
    print("Checking CharmmParameters single-file construction...")

    write_text_file(parameter_path, PARAMETER_TEXT)

    parameters = apo.CharmmParameters(str(parameter_path))

    # The Python wrapper does not expose CharmmParameters getters yet. Use the
    # object in ForceManager construction to verify that the handle is usable.
    psf_path: str = require_file(repo_root / "test/data/nacl_pair.psf")
    psf = apo.CharmmPsf(psf_path)

    fm = apo.ForceManager(psf, parameters)
    fm.setBoxDimensions(BOX_DIMENSIONS)

    return


def check_pathlike_construction(parameter_path: Path) -> None:
    print("Checking CharmmParameters Path-like construction...")

    write_text_file(parameter_path, PARAMETER_TEXT)

    apo.CharmmParameters(parameter_path)

    return


def check_file_list_construction(
    parameter_path: Path, supplemental_parameter_path: Path, repo_root: Path
) -> None:
    print("Checking CharmmParameters list/tuple construction...")

    write_text_file(parameter_path, PARAMETER_TEXT)
    write_text_file(supplemental_parameter_path, SUPPLEMENTAL_PARAMETER_TEXT)

    parameters_from_list = apo.CharmmParameters(
        [str(parameter_path), str(supplemental_parameter_path)]
    )
    parameters_from_tuple = apo.CharmmParameters(
        (parameter_path, supplemental_parameter_path)
    )

    psf_path: str = require_file(repo_root / "test/data/nacl_pair.psf")
    psf = apo.CharmmPsf(psf_path)

    fm_from_list = apo.ForceManager(psf, parameters_from_list)
    fm_from_tuple = apo.ForceManager(psf, parameters_from_tuple)

    fm_from_list.setBoxDimensions(BOX_DIMENSIONS)
    fm_from_tuple.setBoxDimensions(BOX_DIMENSIONS)

    return


def check_validation(
    missing_parameter_path: Path, malformed_parameter_path: Path
) -> None:
    print("Checking CharmmParameters validation...")

    remove_if_exists(missing_parameter_path)
    expect_exception(
        "CharmmParameters rejects a missing parameter file",
        apo.ApoCharmmError,
        lambda: apo.CharmmParameters(str(missing_parameter_path)),
    )

    write_text_file(malformed_parameter_path, MALFORMED_PARAMETER_TEXT)
    expect_exception(
        "CharmmParameters rejects a malformed parameter file",
        apo.ApoCharmmError,
        lambda: apo.CharmmParameters(str(malformed_parameter_path)),
    )

    expect_exception(
        "CharmmParameters rejects an empty parameter-file list",
        apo.ApoCharmmError,
        lambda: apo.CharmmParameters([]),
    )

    expect_exception(
        "CharmmParameters rejects a non-path argument",
        TypeError,
        lambda: apo.CharmmParameters(object()),
    )

    return


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")
    output_dir: Path = repo_root / "test/pytest"
    output_dir.mkdir(parents=True, exist_ok=True)

    parameter_path: Path = output_dir / "tmp_python_api_charmm_parameters.prm"
    pathlike_parameter_path: Path = (
        output_dir / "tmp_python_api_charmm_parameters_pathlike.prm"
    )
    supplemental_parameter_path: Path = (
        output_dir / "tmp_python_api_charmm_parameters_supplement.prm"
    )
    missing_parameter_path: Path = (
        output_dir / "tmp_python_api_charmm_parameters_missing.prm"
    )
    malformed_parameter_path: Path = (
        output_dir / "tmp_python_api_charmm_parameters_malformed.prm"
    )

    generated_files: tuple[Path, ...] = (
        parameter_path,
        pathlike_parameter_path,
        supplemental_parameter_path,
        missing_parameter_path,
        malformed_parameter_path,
    )

    for path in generated_files:
        remove_if_exists(path)

    try:
        check_single_file_construction(parameter_path, repo_root)
        check_pathlike_construction(pathlike_parameter_path)
        check_file_list_construction(
            parameter_path, supplemental_parameter_path, repo_root
        )
        check_validation(missing_parameter_path, malformed_parameter_path)
    finally:
        print("Cleaning up CharmmParameters Python API test files...")
        for path in generated_files:
            remove_if_exists(path)

    print("\033[32m" + "PASS: CharmmParameters Python API tests completed." + "\033[0m")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
