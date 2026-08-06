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
from .error import check_status, configure_status_function
from .force_manager import ForceManager

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_harmonic_restraint_force_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_int],
        "HarmonicRestraintForce construction",
    )

    lib().apo_harmonic_restraint_force_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_harmonic_restraint_force_destroy.restype = None

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_selection,
        [ctypes.c_void_p, ctypes.c_void_p],
        "HarmonicRestraintForce.setSelection(selection)",
    )

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_force_constant,
        [ctypes.c_void_p, ctypes.c_double],
        "HarmonicRestraintForce.setForceConstant(force_constant)",
    )

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_force_constants,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "HarmonicRestraintForce.setForceConstants(force_constants)",
    )

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_reference_coordinates,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "HarmonicRestraintForce.setReferenceCoordinates(reference_coordinates)",
    )

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_masses,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "HarmonicRestraintForce.setMasses(masses)",
    )

    configure_status_function(
        lib().apo_harmonic_restraint_force_set_box_dimensions,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "HarmonicRestraintForce.setBoxDimensions(box_dimensions)",
    )

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

        if num_atoms < -(2**31) or num_atoms > 2**31 - 1:
            raise ValueError("num_atoms must fit in int")

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

        lib().apo_harmonic_restraint_force_set_selection(self.handle, selection.handle)

        return

    def setForceConstant(self, force_constant: float) -> None:
        _initialize_prototypes()

        c_force_constant: ctypes.c_double = ctypes.c_double(float(force_constant))

        lib().apo_harmonic_restraint_force_set_force_constant(
            self.handle, c_force_constant
        )

        return

    def setForceConstants(self, force_constants: Sequence[float]) -> None:
        _initialize_prototypes()

        force_constant_values: list[float] = [float(value) for value in force_constants]

        c_buffer_type = ctypes.c_double * len(force_constant_values)
        c_buffer = c_buffer_type(*force_constant_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(force_constant_values))

        lib().apo_harmonic_restraint_force_set_force_constants(
            self.handle, c_buffer, c_buffer_len
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

        lib().apo_harmonic_restraint_force_set_reference_coordinates(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setMasses(self, masses: Sequence[float]) -> None:
        _initialize_prototypes()

        mass_values: list[float] = [float(value) for value in masses]

        c_buffer_type = ctypes.c_double * len(mass_values)
        c_buffer = c_buffer_type(*mass_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(mass_values))

        lib().apo_harmonic_restraint_force_set_masses(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setBoxDimensions(self, box_dimensions: Sequence[float]) -> None:
        _initialize_prototypes()

        box_values: list[float] = [float(value) for value in box_dimensions]

        c_buffer_type = ctypes.c_double * len(box_values)
        c_buffer = c_buffer_type(*box_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(box_values))

        lib().apo_harmonic_restraint_force_set_box_dimensions(
            self.handle, c_buffer, c_buffer_len
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
