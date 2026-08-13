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
from ._lib import encode_path, lib
from ._types import FilePath
from .error import configure_status_function

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_charmm_crd_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p],
        "CharmmCrd construction",
    )

    lib().apo_charmm_crd_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_charmm_crd_destroy.restype = None

    configure_status_function(
        lib().apo_charmm_crd_get_num_atoms,
        [ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p],
        "CharmmCrd.getNumAtoms()",
    )

    configure_status_function(
        lib().apo_charmm_crd_get_coordinates,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmCrd.getCoordinates()",
    )

    _prototypes_initialized = True

    return


class CharmmCrd(_ApoObject):
    _destroy_function_name = "apo_charmm_crd_destroy"

    def __init__(self, path: FilePath) -> None:
        _initialize_prototypes()
        super().__init__()

        handle: ctypes.c_void_p = ctypes.c_void_p()

        encoded_path: bytes = encode_path(path)
        c_path: ctypes.c_char_p = ctypes.c_char_p(encoded_path)

        lib().apo_charmm_crd_create(ctypes.byref(handle), c_path)

        if handle.value is None:
            raise RuntimeError(
                "apo_charmm_crd_create returned success but produced a NULL handle"
            )

        self._handle = handle

        return

    def getNumAtoms(self) -> int:
        _initialize_prototypes()

        c_num_atoms = ctypes.c_size_t()

        lib().apo_charmm_crd_get_num_atoms(ctypes.byref(c_num_atoms), self.handle)

        return int(c_num_atoms.value)

    def getCoordinates(self) -> list[list[float]]:
        _initialize_prototypes()

        num_atoms: int = self.getNumAtoms()

        c_buffer_type = ctypes.c_double * (num_atoms * 3)
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_atoms * 3)

        lib().apo_charmm_crd_get_coordinates(c_buffer, c_buffer_len, self.handle)

        coordinates: list[list[float]] = []
        for i in range(num_atoms):
            coordinates.append(
                [
                    float(c_buffer[i * 3 + 0]),
                    float(c_buffer[i * 3 + 1]),
                    float(c_buffer[i * 3 + 2]),
                ]
            )

        return coordinates
