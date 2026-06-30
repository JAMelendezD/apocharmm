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

from python_api_test_helpers import (
    get_data_path,
    assert_equal,
    assert_close,
    assert_sequence_close,
)

BOX_DIMENSIONS: tuple[float] = (40.0, 41.0, 42.0)
TOLERANCE: float = 1.0e-7


def main(argc: int, argv: list[str]) -> int:
    data_path: Path = get_data_path()

    prm = apo.CharmmParameters(data_path / "toppar_water_ions.str")
    psf = apo.CharmmPsf(data_path / "nacl_pair.psf")
    crd = apo.CharmmCrd(data_path / "nacl_pair.cor")

    ctx = apo.CharmmContext(psf, prm)
    ctx.setBoxDimensions(BOX_DIMENSIONS)

    ctx.setKappa(0.45)
    ctx.setCutoff(9.0)
    ctx.setCtonnb(7.5)
    ctx.setCtofnb(8.5)
    ctx.setFFTGrid([32, 34, 36])
    ctx.setPmeSplineOrder(6)
    ctx.setVdwType(apo.VdwType.DBEXP)

    assert_close("context kappa", ctx.getKappa(), 0.45, TOLERANCE)
    assert_close("context cutoff", ctx.getCutoff(), 9.0, TOLERANCE)
    assert_close("context ctonnb", ctx.getCtonnb(), 7.5, TOLERANCE)
    assert_close("context ctofnb", ctx.getCtofnb(), 8.5, TOLERANCE)
    assert_equal("context FFT grid", ctx.getFFTGrid(), (32, 34, 36))
    assert_equal("context PME spline order", ctx.getPmeSplineOrder(), 6)
    assert_equal("context VdW type", ctx.getVdwType(), apo.VdwType.DBEXP)

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
