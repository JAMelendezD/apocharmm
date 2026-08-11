# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

import apocharmm as apo
import sys


def main(argc: int, argv: list[str]) -> int:
    # Input variables
    box_dims = [50.0, 50.0, 50.0]
    temperature = 300.0

    # Load CHARMM force field
    prm = apo.CharmmParameters("test/data/toppar_water_ions.str")

    # Load PSF and coordinates
    psf = apo.CharmmPsf("test/data/waterbox.psf")
    crd = apo.CharmmCrd("test/data/waterbox.crd")

    ctx = apo.CharmmContext(psf, prm)
    ctx.setBoxDimensions(box_dims)
    ctx.setCoordinates(crd)

    # Calculate and print potential energy
    ctx.calculatePotentialEnergy(True)

    return 0


if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
