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
import math

from ._base import _ApoObject
from ._lib import encode_path, lib
from .atom_selection import AtomSelection
from .error import check_status
from .force_manager import ForceManager

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    lib().apo_harmonic_restraint_force_create.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_int,
    ]
    lib().apo_harmonic_restraint_force_create.restype = ctypes.c_int

    lib().apo_harmonic_restraint_force_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_harmonic_restraint_force_destroy.restype = None

    lib().apo_harmonic_restraint_force_set_selection.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib().apo_harmonic_restraint_force_set_selection.restype = ctypes.c_int

    lib().apo_harmonic_restraint_force_set_force_constant.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib().apo_harmonic_restraint_force_set_force_constant.restype = ctypes.c_int

    lib().apo_harmonic_restraint_force_set_force_constants.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib().apo_harmonic_restraint_force_set_force_constants.restype = ctypes.c_int

    lib().apo_harmonic_restraint_force_set_reference_coordinates.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib().apo_harmonic_restraint_force_set_reference_coordinates.restype = ctypes.c_int

    lib().apo_harmonic_restraint_force_set_masses.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib().apo_harmonic_restraint_force_set_masses.restype = ctypes.c_int

    lib().apo_harmonic_restraint_force_set_box_dimensions.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib().apo_harmonic_restraint_force_set_box_dimensions.restype = ctypes.c_int

    lib().apo_force_manager_subscribe_harmonic_restraint_force.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
        ctypes.c_char_p,
    ]
    lib().apo_force_manager_subscribe_harmonic_restraint_force.restype = ctypes.c_int

    lib().apo_force_manager_unsubscribe_harmonic_restraint_force.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_unsubscribe_harmonic_restraint_force.restype = ctypes.c_int

    _prototypes_initialized = True

    return


class HarmonicRestraintForce(_ApoObject):
    _destroy_function_name = "apo_harmonic_restraint_force_destroy"

    def __init__(self, num_atoms: int) -> None:
        super().__init__()

        if not isinstance(num_atoms, int):
            raise TypeError("num_atoms must be an int")

        if num_atoms <= 0 or num_atoms > 2**31 - 1:
            raise ValueError("num_atoms must fit in positive int")

        _initialize_prototypes()

        handle: ctypes.c_void_p = ctypes.c_void_p()
        c_num_atoms: ctypes.c_int = ctypes.c_int(num_atoms)

        status = lib().apo_harmonic_restraint_force_create(
            ctypes.byref(handle), c_num_atoms
        )

        check_status(status, "HarmonicRestraintForce construction failed")

        if handle.value is None:
            raise RuntimeError(
                "apo_harmonic_restraint_force_create returned success but produced a NULL handle"
            )

        self._handle = handle
        self._default_force_tag: str = "harm"

        return

    @property
    def default_force_tag(self) -> str:
        return self._default_force_tag

    def setSelection(self, selection: AtomSelection) -> None:
        _initialize_prototypes()

        if not isinstance(selection, AtomSelection):
            raise TypeError(
                "HarmonicRestraintForce.setSelection expects an AtomSelection"
            )

        status = lib().apo_harmonic_restraint_force_set_selection(
            self.handle, selection.handle
        )

        check_status(status, "HarmonicRestraintForce.setSelection(selection) failed")

        return

    def setForceConstant(self, force_constant: float) -> None:
        _initialize_prototypes()

        c_force_constant: ctypes.c_double = ctypes.c_double(float(force_constant))

        status = lib().apo_harmonic_restraint_force_set_force_constant(
            self.handle, c_force_constant
        )

        check_status(
            status, "HarmonicRestraintForce.setForceConstant(force_constant) failed"
        )

        return

    def setForceConstants(self, force_constants: Sequence[float]) -> None:
        _initialize_prototypes()

        force_constant_values: list[float] = [float(value) for value in force_constants]

        c_buffer_type = ctypes.c_double * len(force_constant_values)
        c_buffer = c_buffer_type(*force_constant_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(force_constant_values))

        status = lib().apo_harmonic_restraint_force_set_force_constants(
            self.handle, c_buffer, c_buffer_len
        )

        check_status(
            status, "HarmonicRestraintForce.setForceConstants(force_constants) failed"
        )

        return

    def setReferenceCoordinates(
        self, reference_coordinates: Sequence[Sequence[float]]
    ) -> None:
        _initialize_prototypes()

        flattened_coordinates: list[float] = []

        for i, coordinate in enumerate(reference_coordinates):
            coordinate_values: list[float] = [float(value) for value in coordinate]

            if len(coordinate_values) != 3:
                raise ValueError(
                    f"reference_coordinates[{i}] must contain exactly 3 elements"
                )

            flattened_coordinates.extend(coordinate_values)

        c_buffer_type = ctypes.c_double * len(flattened_coordinates)
        c_buffer = c_buffer_type(*flattened_coordinates)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(flattened_coordinates))

        status = lib().apo_harmonic_restraint_force_set_reference_coordinates(
            self.handle, c_buffer, c_buffer_len
        )

        check_status(
            status,
            "HarmonicRestraintForce.setReferenceCoordinates(reference_coordinates) failed",
        )

        return

    def setMasses(self, masses: Sequence[float]) -> None:
        _initialize_prototypes()

        mass_values: list[float] = [float(value) for value in masses]

        c_buffer_type = ctypes.c_double * len(mass_values)
        c_buffer = c_buffer_type(*mass_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(mass_values))

        status = lib().apo_harmonic_restraint_force_set_masses(
            self.handle, c_buffer, c_buffer_len
        )

        check_status(status, "HarmonicRestraintForce.setMasses(masses) failed")

        return

    def setBoxDimensions(self, box_dimensions: Sequence[float]) -> None:
        _initialize_prototypes()

        box_values: list[float] = [float(value) for value in box_dimensions]

        c_buffer_type = ctypes.c_double * len(box_values)
        c_buffer = c_buffer_type(*box_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(box_values))

        status = lib().apo_harmonic_restraint_force_set_box_dimensions(
            self.handle, c_buffer, c_buffer_len
        )

        check_status(
            status, "HarmonicRestraintForce.setBoxDimensions(box_dimensions) failed"
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

        status = lib().apo_force_manager_subscribe_harmonic_restraint_force(
            force_manager.handle, self.handle, c_force_tag
        )

        check_status(status, "ForceManager.subscribe(HarmonicRestraintForce) failed")

        return

    def _unsubscribe_from_force_manager(self, force_manager: ForceManager) -> None:
        _initialize_prototypes()

        status = lib().apo_force_manager_unsubscribe_harmonic_restraint_force(
            force_manager.handle, self.handle
        )

        check_status(status, "ForceManager.unsubscribe(HarmonicRestraintForce) failed")

        return
