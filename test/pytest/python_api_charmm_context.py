# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

import sys

import apocharmm as apo

from python_api_test_helpers import get_data_path, assert_equal, assert_sequence_close

BOX_DIMENSIONS: tuple[float] = (40.0, 41.0, 42.0)


def main(argc: int, argv: list[str]) -> int:
    data_path: Path = get_data_path()

    prm = apo.CharmmParameters(data_path / "toppar_water_ions.str")
    psf = apo.CharmmPsf(data_path / "nacl_pair.psf")
    crd = apo.CharmmCrd(data_path / "nacl_pair.cor")

    ctx = apo.CharmmContext(psf, prm)
    ctx.setBoxDimensions(BOX_DIMENSIONS)

    assert_sequence_close(
        "context-first box dimensions", ctx.getBoxDimensions(), BOX_DIMENSIONS, 0.0
    )

    ctx.setCoordinates(crd)

    assert_equal("context-first num atoms", ctx.getNumAtoms(), psf.getNumAtoms())

    xyzq = ctx.getCoordinatesCharges()
    xyzm = ctx.getVelocityMass()

    assert_equal("coordinate row count", len(xyzq), psf.getNumAtoms())
    assert_equal("velocity row count", len(xyzm), psf.getNumAtoms())

    print("\033[32m" + "PASS: CharmmContext Python API tests completed." + "\033[0m")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
