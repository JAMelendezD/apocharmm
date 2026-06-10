# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from collections.abc import Sequence
import ctypes

from ._base import _ApoObject
from ._lib import lib
from .enums import PeriodicBoundaryCondition
from .error import check_status

from .charmm_parameters import CharmmParameters
from .charmm_psf import CharmmPsf

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    lib().apo_force_manager_create.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_create.restype = ctypes.c_int

    lib().apo_force_manager_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_force_manager_destroy.restype = None

    lib().apo_force_manager_set_box_dimensions.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
    ]
    lib().apo_force_manager_set_box_dimensions.restype = ctypes.c_int

    lib().apo_force_manager_set_periodic_boundary_condition.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
    ]
    lib().apo_force_manager_set_periodic_boundary_condition.restype = ctypes.c_int

    lib().apo_force_manager_get_periodic_boundary_condition.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_periodic_boundary_condition.restype = ctypes.c_int

    _prototypes_initialized = True

    return


class ForceManager(_ApoObject):
    _destroy_function_name = "apo_force_manager_destroy"

    def __init__(self, psf: CharmmPsf, parameters: CharmmParameters) -> None:
        _initialize_prototypes()
        super().__init__()

        if not isinstance(psf, CharmmPsf):
            raise TypeError("ForceManager expects a CharmmPsf")

        if not isinstance(parameters, CharmmParameters):
            raise TypeError("ForceManager expects a CharmmParameters")

        handle = ctypes.c_void_p()

        status = lib().apo_force_manager_create(
            ctypes.byref(handle), psf.handle, parameters.handle
        )

        check_status(status, "ForceManager construction failed")

        if handle.value is None:
            raise RuntimeError(
                "apo_force_manager_create returned success but produced a NULL handle"
            )

        self._handle = handle

        self._psf: CharmmPsf = psf
        self._parameters: CharmmParameters = parameters

        return

    def setBoxDimensions(self, box_dimensions: Sequence[float]) -> None:
        _initialize_prototypes()

        if len(box_dimensions) != 3:
            raise ValueError(
                "ForceManager.setBoxDimensions expects exactly three box_dimensions"
            )

        x: float = float(box_dimensions[0])
        y: float = float(box_dimensions[1])
        z: float = float(box_dimensions[2])

        c_x: ctypes.c_double = ctypes.c_double(x)
        c_y: ctypes.c_double = ctypes.c_double(y)
        c_z: ctypes.c_double = ctypes.c_double(z)

        status = lib().apo_force_manager_set_box_dimensions(self.handle, c_x, c_y, c_z)

        check_status(status, "ForceManager.setBoxDimensions() failed")

        return

    def setPeriodicBoundaryCondition(
        self, pbc: PeriodicBoundaryCondition | int
    ) -> None:
        _initialize_prototypes()

        if isinstance(pbc, bool):
            raise TypeError("pbc must be a PeriodicBoundaryCondition")

        try:
            pbc_value: PeriodicBoundaryCondition = PeriodicBoundaryCondition(pbc)
        except ValueError as exc:
            raise ValueError(f"invalid pbc: {pbc!r}") from exc

        c_pbc: ctypes.c_int = ctypes.c_int(int(pbc_value))

        status = lib().apo_force_manager_set_periodic_boundary_condition(
            self.handle, c_pbc
        )

        check_status(
            status,
            "ForceManager.setPeriodicBoundaryCondition(pbc) failed",
        )

        return

    def getPeriodicBoundaryCondition(self) -> PeriodicBoundaryCondition:
        _initialize_prototypes()

        pbc = ctypes.c_int()

        status = lib().apo_force_manager_get_periodic_boundary_condition(
            ctypes.byref(pbc), self.handle
        )

        check_status(
            status,
            "ForceManager.getPeriodicBoundaryCondition() failed",
        )

        return PeriodicBoundaryCondition(pbc.value)
