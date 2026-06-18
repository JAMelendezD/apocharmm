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
from .atom_selection import AtomSelection
from .error import check_status
from .force_manager import ForceManager

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    lib().apo_harmonic_center_of_mass_restraint_force_create.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_int,
    ]
    lib().apo_harmonic_center_of_mass_restraint_force_create.restype = ctypes.c_int

    lib().apo_harmonic_center_of_mass_restraint_force_destroy.argtypes = [
        ctypes.c_void_p
    ]
    lib().apo_harmonic_center_of_mass_restraint_force_destroy.restype = None

    lib().apo_harmonic_center_of_mass_restraint_force_set_selection.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib().apo_harmonic_center_of_mass_restraint_force_set_selection.restype = (
        ctypes.c_int
    )

    lib().apo_harmonic_center_of_mass_restraint_force_set_force_constant.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib().apo_harmonic_center_of_mass_restraint_force_set_force_constant.restype = (
        ctypes.c_int
    )

    lib().apo_harmonic_center_of_mass_restraint_force_set_reference_position.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
    ]
    lib().apo_harmonic_center_of_mass_restraint_force_set_reference_position.restype = (
        ctypes.c_int
    )

    lib().apo_harmonic_center_of_mass_restraint_force_set_reference_distance.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib().apo_harmonic_center_of_mass_restraint_force_set_reference_distance.restype = (
        ctypes.c_int
    )

    lib().apo_harmonic_center_of_mass_restraint_force_set_masses.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib().apo_harmonic_center_of_mass_restraint_force_set_masses.restype = ctypes.c_int

    lib().apo_harmonic_center_of_mass_restraint_force_set_mass_weighting.argtypes = [
        ctypes.c_void_p,
        ctypes.c_bool,
    ]
    lib().apo_harmonic_center_of_mass_restraint_force_set_mass_weighting.restype = (
        ctypes.c_int
    )

    lib().apo_force_manager_subscribe_harmonic_center_of_mass_restraint_force.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_char_p,
    ]
    lib().apo_force_manager_subscribe_harmonic_center_of_mass_restraint_force.restype = (
        ctypes.c_int
    )

    lib().apo_force_manager_unsubscribe_harmonic_center_of_mass_restraint_force.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_unsubscribe_harmonic_center_of_mass_restraint_force.restype = (
        ctypes.c_int
    )

    _prototypes_initialized = True

    return


class HarmonicCenterOfMassRestraintForce(_ApoObject):
    _destroy_function_name = "apo_harmonic_center_of_mass_restraint_force_destroy"

    def __init__(self, num_atoms: int) -> None:
        super().__init__()

        if not isinstance(num_atoms, int):
            raise TypeError("num_atoms must be an int")

        if num_atoms <= 0 or num_atoms > 2**31 - 1:
            raise ValueError("num_atoms must fit in positive int")

        _initialize_prototypes()

        handle: ctypes.c_void_p = ctypes.c_void_p()
        c_num_atoms: ctypes.c_int = ctypes.c_int(num_atoms)

        status = lib().apo_harmonic_center_of_mass_restraint_force_create(
            ctypes.byref(handle), c_num_atoms
        )

        check_status(status, "HarmonicCenterOfMassRestraintForce construction failed")

        if handle.value is None:
            raise RuntimeError(
                "apo_harmonic_center_of_mass_restraint_force_create returned success but produced a NULL handle"
            )

        self._handle = handle
        self._default_force_tag: str = "hmcm"

        return

    @property
    def default_force_tag(self) -> str:
        return self._default_force_tag

    def setSelection(self, selection: AtomSelection) -> None:
        _initialize_prototypes()

        if not isinstance(selection, AtomSelection):
            raise TypeError(
                "HarmonicCenterOfMassRestraintForce.setSelection expects an AtomSelection"
            )

        status = lib().apo_harmonic_center_of_mass_restraint_force_set_selection(
            self.handle, selection.handle
        )

        check_status(
            status, "HarmonicCenterOfMassRestraintForce.setSelection(selection) failed"
        )

        return

    def setForceConstant(self, force_constant: float) -> None:
        _initialize_prototypes()

        c_force_constant: ctypes.c_double = ctypes.c_double(float(force_constant))

        status = lib().apo_harmonic_center_of_mass_restraint_force_set_force_constant(
            self.handle, c_force_constant
        )

        check_status(
            status, "HarmonicCenterOfMassRestraintForceConstant(force_constant) failed"
        )

        return

    def setReferencePosition(
        self,
        reference_position: Sequence[float],
        reference_mask: Sequence[int] | None = None,
    ) -> None:
        _initialize_prototypes()

        reference_position_values: list[float] = [
            float(value) for value in reference_position
        ]

        if reference_mask is None:
            reference_mask_values: list[int] = [1, 1, 1]
        else:
            reference_mask_values: list[int] = [int(value) for value in reference_mask]

        c_buffer_type1 = ctypes.c_double * len(reference_position_values)
        c_buffer1 = c_buffer_type1(*reference_position_values)
        c_buffer_len1: ctypes.c_size_t = ctypes.c_size_t(len(reference_position_values))

        c_buffer_type2 = ctypes.c_int * len(reference_mask_values)
        c_buffer2 = c_buffer_type2(*reference_mask_values)
        c_buffer_len2: ctypes.c_size_t = ctypes.c_size_t(len(reference_mask_values))

        status = (
            lib().apo_harmonic_center_of_mass_restraint_force_set_reference_position(
                self.handle, c_buffer1, c_buffer_len1, c_buffer2, c_buffer_len2
            )
        )

        check_status(
            status,
            "HarmonicCenterOfMassRestraintForce.setReferencePosition(reference_position) failed",
        )

        return

    def setReferenceDistance(self, reference_distance: float) -> None:
        _initialize_prototypes()

        c_reference_distance: ctypes.c_double = ctypes.c_double(
            float(reference_distance)
        )

        status = (
            lib().apo_harmonic_center_of_mass_restraint_force_set_reference_distance(
                self.handle, c_reference_distance
            )
        )

        check_status(
            status,
            "HarmonicCenterOfMassRestraintForce.setReferenceDistance(reference_distance) failed",
        )

        return

    def setMasses(self, masses: Sequence[float]) -> None:
        _initialize_prototypes()

        mass_values: list[float] = [float(value) for value in masses]

        c_buffer_type = ctypes.c_double * len(mass_values)
        c_buffer = c_buffer_type(*mass_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(mass_values))

        status = lib().apo_harmonic_center_of_mass_restraint_force_set_masses(
            self.handle, c_buffer, c_buffer_len
        )

        check_status(
            status, "HarmonicCenterOfMassRestraintForce.setMasses(masses) failed"
        )

        return

    def setMassWeighting(self, flag: bool) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        status = lib().apo_harmonic_center_of_mass_restraint_force_set_mass_weighting(
            self.handle, c_flag
        )

        check_status(
            status, "HarmonicCenterOfMassRestraintForce.setMassWeighting(flag) failed"
        )

        return

    def _subscribe_to_force_manager(
        self, force_manager: ForceManager, force_tag: str | None = None
    ) -> None:
        _initialize_prototypes()

        if not isinstance(force_manager, ForceManager):
            raise TypeError("force_manager must be a ForceManager")

        if force_tag is None:
            force_tag_value: str = self._default_force_tag
        else:
            if not isinstance(force_tag, str):
                raise TypeError("force_tag must be a str")

            force_tag_value = force_tag

        if force_tag_value == "":
            raise ValueError("force_tag must not be empty")

        encoded_force_tag: bytes = force_tag_value.encode("utf-8")
        c_force_tag: ctypes.c_char_p = ctypes.c_char_p(encoded_force_tag)

        status = (
            lib().apo_force_manager_subscribe_harmonic_center_of_mass_restraint_force(
                force_manager.handle, self.handle, c_force_tag
            )
        )

        check_status(
            status, "ForceManager.subscribe(HarmonicCenterOfMassRestraintForce) failed"
        )

        return

    def _unsubscribe_from_force_manager(self, force_manager: ForceManager) -> None:
        _initialize_prototypes()

        status = (
            lib().apo_force_manager_unsubscribe_harmonic_center_of_mass_restraint_force(
                force_manager.handle, self.handle
            )
        )

        check_status(
            status,
            "ForceManager.unsubscribe(HarmonicCenterOfMassRestraintForce) failed",
        )

        return
