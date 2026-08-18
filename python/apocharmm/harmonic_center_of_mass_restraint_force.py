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
from .error import configure_status_function
from .force_manager import ForceManager

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_int],
        "HarmonicCenterOfMassRestraintForce construction",
    )

    lib().apo_harmonic_center_of_mass_restraint_force_destroy.argtypes = [
        ctypes.c_void_p
    ]
    lib().apo_harmonic_center_of_mass_restraint_force_destroy.restype = None

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_selection,
        [ctypes.c_void_p, ctypes.c_void_p],
        "HarmonicCenterOfMassRestraintForce.setSelection(selection)",
    )

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_force_constant,
        [ctypes.c_void_p, ctypes.c_double],
        "HarmonicCenterOfMassRestraintForce.setForceConstant(force_constant)",
    )

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_reference_position,
        [
            ctypes.c_void_p,
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t,
            ctypes.POINTER(ctypes.c_int),
            ctypes.c_size_t,
        ],
        "HarmonicCenterOfMassRestraintForce.setReferencePosition(reference_position)",
    )

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_reference_distance,
        [ctypes.c_void_p, ctypes.c_double],
        "HarmonicCenterOfMassRestraintForce.setReferenceDistance(reference_distance)",
    )

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_masses,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "HarmonicCenterOfMassRestraintForce.setMasses(masses)",
    )

    configure_status_function(
        lib().apo_harmonic_center_of_mass_restraint_force_set_mass_weighting,
        [ctypes.c_void_p, ctypes.c_bool],
        "HarmonicCenterOfMassRestraintForce.setMassWeighting(flag)",
    )

    configure_status_function(
        lib().apo_force_manager_subscribe_harmonic_center_of_mass_restraint_force,
        [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_char_p],
        "ForceManager.subscribe(HarmonicCenterOfMassRestraintForce)",
    )

    configure_status_function(
        lib().apo_force_manager_unsubscribe_harmonic_center_of_mass_restraint_force,
        [ctypes.c_void_p, ctypes.c_void_p],
        "ForceManager.unsubscribe(HarmonicCenterOfMassRestraintForce)",
    )

    _prototypes_initialized = True

    return


class HarmonicCenterOfMassRestraintForce(_ApoObject):
    _destroy_function_name = "apo_harmonic_center_of_mass_restraint_force_destroy"

    def __init__(self, num_atoms: int) -> None:
        super().__init__()

        if not isinstance(num_atoms, int):
            raise TypeError("num_atoms must be an int")

        if num_atoms < -(2**31) or num_atoms > 2**31 - 1:
            raise ValueError("num_atoms must fit in int")

        _initialize_prototypes()

        handle: ctypes.c_void_p = ctypes.c_void_p()
        c_num_atoms: ctypes.c_int = ctypes.c_int(num_atoms)

        lib().apo_harmonic_center_of_mass_restraint_force_create(
            ctypes.byref(handle), c_num_atoms
        )

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

        lib().apo_harmonic_center_of_mass_restraint_force_set_selection(
            self.handle, selection.handle
        )

        return

    def setForceConstant(self, force_constant: float) -> None:
        _initialize_prototypes()

        c_force_constant: ctypes.c_double = ctypes.c_double(float(force_constant))

        lib().apo_harmonic_center_of_mass_restraint_force_set_force_constant(
            self.handle, c_force_constant
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

        lib().apo_harmonic_center_of_mass_restraint_force_set_reference_position(
            self.handle, c_buffer1, c_buffer_len1, c_buffer2, c_buffer_len2
        )

        return

    def setReferenceDistance(self, reference_distance: float) -> None:
        _initialize_prototypes()

        c_reference_distance: ctypes.c_double = ctypes.c_double(
            float(reference_distance)
        )

        lib().apo_harmonic_center_of_mass_restraint_force_set_reference_distance(
            self.handle, c_reference_distance
        )

        return

    def setMasses(self, masses: Sequence[float]) -> None:
        _initialize_prototypes()

        mass_values: list[float] = [float(value) for value in masses]

        c_buffer_type = ctypes.c_double * len(mass_values)
        c_buffer = c_buffer_type(*mass_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(mass_values))

        lib().apo_harmonic_center_of_mass_restraint_force_set_masses(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setMassWeighting(self, flag: bool) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        lib().apo_harmonic_center_of_mass_restraint_force_set_mass_weighting(
            self.handle, c_flag
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

        lib().apo_force_manager_subscribe_harmonic_center_of_mass_restraint_force(
            force_manager.handle, self.handle, c_force_tag
        )

        return

    def _unsubscribe_from_force_manager(self, force_manager: ForceManager) -> None:
        _initialize_prototypes()

        lib().apo_force_manager_unsubscribe_harmonic_center_of_mass_restraint_force(
            force_manager.handle, self.handle
        )

        return
