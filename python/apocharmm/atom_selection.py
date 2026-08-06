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
from .error import configure_status_function

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    lib().apo_atom_selection_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_atom_selection_destroy.restype = None

    configure_status_function(
        lib().apo_atom_selection_get_num_atoms,
        [ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p],
        "AtomSelection.getNumAtoms()",
    )

    configure_status_function(
        lib().apo_atom_selection_get_num_selected,
        [ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p],
        "AtomSelection.getNumSelected()",
    )

    configure_status_function(
        lib().apo_atom_selection_get_atom_indices,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.c_void_p],
        "AtomSelection.getAtomIndices()",
    )

    configure_status_function(
        lib().apo_atom_selection_contains,
        [ctypes.POINTER(ctypes.c_bool), ctypes.c_void_p, ctypes.c_int],
        "AtomSelection.contains(atom_index)",
    )

    _prototypes_initialized = True

    return


class AtomSelection(_ApoObject):
    _destroy_function_name = "apo_atom_selection_destroy"

    def __init__(self, handle: ctypes.c_void_p) -> None:
        _initialize_prototypes()
        super().__init__()

        if not isinstance(handle, ctypes.c_void_p):
            raise TypeError("AtomSelection expects a ctypes.c_void_p handle")

        if handle.value is None:
            raise RuntimeError("AtomSelection construction failed: NULL handle")

        self._handle = handle

        return

    def getNumAtoms(self) -> int:
        _initialize_prototypes()

        c_num_atoms: ctypes.c_size_t = ctypes.c_size_t()

        lib().apo_atom_selection_get_num_atoms(ctypes.byref(c_num_atoms), self.handle)

        return int(c_num_atoms.value)

    def getNumSelected(self) -> int:
        _initialize_prototypes()

        c_num_selected: ctypes.c_size_t = ctypes.c_size_t()

        lib().apo_atom_selection_get_num_selected(
            ctypes.byref(c_num_selected), self.handle
        )

        return int(c_num_selected.value)

    def getAtomIndices(self) -> list[int]:
        _initialize_prototypes()

        num_selected: int = self.getNumSelected()

        c_buffer_type = ctypes.c_int * num_selected
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_selected)

        lib().apo_atom_selection_get_atom_indices(c_buffer, c_buffer_len, self.handle)

        atom_indices: list[int] = []
        for i in range(num_selected):
            atom_indices.append(int(c_buffer[i]))

        return atom_indices

    def contains(self, atom_index: int) -> bool:
        _initialize_prototypes()

        if not isinstance(atom_index, int) or isinstance(atom_index, bool):
            raise TypeError("atom_index must be an int")

        if atom_index < -(2**31) or atom_index > 2**31 - 1:
            raise ValueError("atom_index must fit in int")

        c_is_selected: ctypes.c_bool = ctypes.c_bool()
        c_atom_index: ctypes.c_int = ctypes.c_int(atom_index)

        lib().apo_atom_selection_contains(
            ctypes.byref(c_is_selected), self.handle, c_atom_index
        )

        return bool(c_is_selected.value)
