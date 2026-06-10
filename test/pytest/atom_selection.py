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


def main(argc, argv):
    # Load PSF
    psf = apo.CharmmPsf("test/data/1lvz.psf")

    selector = apo.AtomSelector(psf)

    all_atoms = selector.select("all")
    assert all_atoms.getNumAtoms() == psf.getNumAtoms()
    assert all_atoms.getNumSelected() == psf.getNumAtoms()
    assert all_atoms.getAtomIndices() == list(range(psf.getNumAtoms()))

    ca_atoms = selector.select("type CA")
    ca_indices = ca_atoms.getAtomIndices()
    for atom_index in ca_indices:
        assert ca_atoms.contains(atom_index)

    return 0


if __name__ == "__main__":
    main(len(sys.argv), sys.argv)
