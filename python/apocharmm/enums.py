# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from enum import IntEnum


class PeriodicBoundaryCondition(IntEnum):
    NONE = 0
    P1 = 1
    P21 = 2


class CrystalType(IntEnum):
    NONE = 0
    CUBIC = 1
    TETRAGONAL = 2
    ORTHORHOMBIC = 3
