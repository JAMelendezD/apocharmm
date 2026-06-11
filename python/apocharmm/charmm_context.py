# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

import ctypes

from ._base import _ApoObject
from ._lib import lib
from .enums import PeriodicBoundaryCondition
from .error import check_status

from .charmm_crd import CharmmCrd
from .force_manager import ForceManager

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    lib().apo_charmm_context_create.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_create.restype = ctypes.c_int

    lib().apo_charmm_context_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_charmm_context_destroy.restype = None

    lib().apo_charmm_context_set_coordinates.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_set_coordinates.restype = ctypes.c_int

    lib().apo_charmm_context_set_periodic_boundary_condition.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
    ]
    lib().apo_charmm_context_set_periodic_boundary_condition.restype = ctypes.c_int

    lib().apo_charmm_context_set_random_seed_for_velocities.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint64,
    ]
    lib().apo_charmm_context_set_random_seed_for_velocities.restype = ctypes.c_int

    lib().apo_charmm_context_use_holonomic_constraints.argtypes = [
        ctypes.c_void_p,
        ctypes.c_bool,
    ]
    lib().apo_charmm_context_use_holonomic_constraints.restype = ctypes.c_int

    lib().apo_charmm_context_get_num_atoms.argtypes = [
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_num_atoms.restype = ctypes.c_int

    lib().apo_charmm_context_get_coordinates_charges.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_coordinates_charges.restype = ctypes.c_int

    lib().apo_charmm_context_get_box_dimensions.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_box_dimensions.restype = ctypes.c_int

    lib().apo_charmm_context_get_periodic_boundary_condition.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_periodic_boundary_condition.restype = ctypes.c_int

    lib().apo_charmm_context_assign_velocities_at_temperature.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib().apo_charmm_context_assign_velocities_at_temperature.restype = ctypes.c_int

    lib().apo_charmm_context_compute_temperature.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_compute_temperature.restype = ctypes.c_int

    _prototypes_initialized = True

    return


class CharmmContext(_ApoObject):
    _destroy_function_name = "apo_charmm_context_destroy"

    def __init__(self, force_manager: ForceManager) -> None:
        _initialize_prototypes()
        super().__init__()

        if not isinstance(force_manager, ForceManager):
            raise TypeError("CharmmContext expects a ForceManager")

        handle: ctypes.c_void_p = ctypes.c_void_p()

        status = lib().apo_charmm_context_create(
            ctypes.byref(handle), force_manager.handle
        )

        check_status(status, "CharmmContext construction failed")

        if handle.value is None:
            raise RuntimeError(
                "apo_charmm_context_create returned success but produced a NULL handle"
            )

        self._handle = handle
        self._force_manager: ForceManager = force_manager

        return

    def setCoordinates(self, crd: CharmmCrd) -> None:
        _initialize_prototypes()

        if not isinstance(crd, CharmmCrd):
            raise TypeError("CharmmContext.setCoordinates expects a CharmmCrd")

        status = lib().apo_charmm_context_set_coordinates(self.handle, crd.handle)

        check_status(status, "CharmmContext.setCoordinates(crd) failed")

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

        status = lib().apo_charmm_context_set_periodic_boundary_condition(
            self.handle, c_pbc
        )

        check_status(status, "CharmmContext.setPeriodicBoundaryCondition(pbc) failed")

        return

    def setRandomSeedForVelocities(self, seed: int) -> None:
        _initialize_prototypes()

        if isinstance(seed, bool) or seed < 0 or seed > 2**64 - 1:
            raise ValueError("seed must fit in uint64_t")

        c_seed: ctypes.c_uint64 = ctypes.c_uint64(seed)

        status = lib().apo_charmm_context_set_random_seed_for_velocities(
            self.handle, c_seed
        )

        check_status(status, "CharmmContext.setRandomSeedForVelocities(seed) failed")

        return

    def useHolonomicConstraints(self, use_holonomic_constraints: bool) -> None:
        _initialize_prototypes()

        c_use_holonomic_constraints: ctypes.c_bool = ctypes.c_bool(
            use_holonomic_constraints
        )

        status = lib().apo_charmm_context_use_holonomic_constraints(
            self.handle, c_use_holonomic_constraints
        )

        check_status(
            status,
            "CharmmContext.useHolonomicConstraints(useHolonomicConstraints) failed",
        )

        return

    def getNumAtoms(self) -> int:
        _initialize_prototypes()

        num_atoms = ctypes.c_size_t()

        status = lib().apo_charmm_context_get_num_atoms(
            ctypes.byref(num_atoms), self.handle
        )

        check_status(status, "CharmmContext.getNumAtoms() failed")

        return int(num_atoms.value)

    def getCoordinatesCharges(self) -> list[list[float]]:
        _initialize_prototypes()

        num_atoms: int = self.getNumAtoms()
        buffer_len: int = 4 * num_atoms
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(buffer_len)

        buffer_type = ctypes.c_double * buffer_len
        buffer = buffer_type()

        status = lib().apo_charmm_context_get_coordinates_charges(
            buffer, c_buffer_len, self.handle
        )

        check_status(status, "CharmmContext.getCoordinatesCharges() failed")

        coordinates_charges: list[list[float]] = []
        for i in range(num_atoms):
            coordinates_charges.append(
                [
                    float(buffer[i * 4 + 0]),
                    float(buffer[i * 4 + 1]),
                    float(buffer[i * 4 + 2]),
                    float(buffer[i * 4 + 3]),
                ]
            )

        return coordinates_charges

    def getBoxDimensions(self) -> tuple[float, float, float]:
        _initialize_prototypes()

        c_x: ctypes.c_double = ctypes.c_double()
        c_y: ctypes.c_double = ctypes.c_double()
        c_z: ctypes.c_double = ctypes.c_double()

        status = lib().apo_charmm_context_get_box_dimensions(
            ctypes.byref(c_x), ctypes.byref(c_y), ctypes.byref(c_z), self.handle
        )

        check_status(status, "CharmmContext.getBoxDimensions() failed")

        return (float(c_x.value), float(c_y.value), float(c_z.value))

    def getPeriodicBoundaryCondition(self) -> PeriodicBoundaryCondition:
        _initialize_prototypes()

        pbc = ctypes.c_int()

        status = lib().apo_charmm_context_get_periodic_boundary_condition(
            ctypes.byref(pbc), self.handle
        )

        check_status(status, "CharmmContext.getPeriodicBoundaryCondition() failed")

        return PeriodicBoundaryCondition(pbc.value)

    def assignVelocitiesAtTemperature(self, temperature: float) -> None:
        _initialize_prototypes()

        c_temperature: ctypes.c_double = ctypes.c_double(temperature)

        status = lib().apo_charmm_context_assign_velocities_at_temperature(
            self.handle, c_temperature
        )

        check_status(
            status, "CharmmContext.assignVelocitiesAtTemperature(temperature) failed"
        )

        return

    def computeTemperature(self) -> float:
        _initialize_prototypes()

        temperature = ctypes.c_double()

        status = lib().apo_charmm_context_compute_temperature(
            ctypes.byref(temperature), self.handle
        )

        check_status(status, "CharmmContext.computeTemperature() failed")

        return float(temperature.value)
