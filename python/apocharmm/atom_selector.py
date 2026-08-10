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
from .atom_selection import AtomSelection
from .charmm_psf import CharmmPsf
from .error import configure_status_function

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_atom_selector_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p],
        "AtomSelector construction",
    )

    lib().apo_atom_selector_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_atom_selector_destroy.restype = None

    configure_status_function(
        lib().apo_atom_selector_select,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p, ctypes.c_char_p],
        "AtomSelector.select(selection_string)",
    )

    _prototypes_initialized = True

    return


class AtomSelector(_ApoObject):
    _destroy_function_name = "apo_atom_selector_destroy"

    def __init__(self, psf: CharmmPsf) -> None:
        _initialize_prototypes()
        super().__init__()

        if not isinstance(psf, CharmmPsf):
            raise TypeError("AtomSelector expects a CharmmPsf")

        handle: ctypes.c_void_p = ctypes.c_void_p()

        lib().apo_atom_selector_create(ctypes.byref(handle), psf.handle)

        if handle.value is None:
            raise RuntimeError(
                "apo_atom_selector_create returned success but produced a NULL handle"
            )

        self._handle = handle
        self._psf: CharmmPsf = psf

        return

    def select(self, selection_string: str) -> AtomSelection:
        _initialize_prototypes()

        if not isinstance(selection_string, str):
            raise TypeError("selection_string must be a str")

        handle: ctypes.c_void_p = ctypes.c_void_p()

        encoded_selection: bytes = selection_string.encode("utf-8")
        c_selection_string: ctypes.c_char_p = ctypes.c_char_p(encoded_selection)

        lib().apo_atom_selector_select(
            ctypes.byref(handle), self.handle, c_selection_string
        )

        if handle.value is None:
            raise RuntimeError(
                "apo_atom_selector_select returned success but produced a NULL handle"
            )

        return AtomSelection(handle)
