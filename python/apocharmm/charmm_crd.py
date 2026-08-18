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
    """
    @anchor python_charmm_crd
    @brief Owns coordinates parsed from one CHARMM CRD or COR file.

    Construct the wrapper from a `str`, `bytes`, `os.PathLike[str]`, or
    `os.PathLike[bytes]` path. The path is converted with `os.fsencode`, copied
    through the C ABI, and not retained. Native parsing completes synchronously
    before construction returns.

    `getCoordinates()` returns a newly allocated nested Python list in
    atom-major `[x, y, z]` order. Returned values never alias native storage.
    The wrapper owns one C handle. `close()`, `destroy()`, and context-manager
    exit release it idempotently. Accessing a getter after closure raises
    `RuntimeError`.

    The wrapper provides no internal synchronization. Do not overlap `close()`
    with another method call from a different thread.

    @see coordinates
    """

    _destroy_function_name = "apo_charmm_crd_destroy"

    def __init__(self, path: FilePath) -> None:
        """
        @brief Constructs an owning wrapper from one CHARMM coordinate file.

        Standard and extended CHARMM fixed-width coordinate records are
        supported. File-name extensions are not validated. The encoded path and
        file contents are not retained after native construction finishes.

        @param[in] path Path supplied as `str`, `bytes`,
        `os.PathLike[str]`, or `os.PathLike[bytes]`. Its encoded representation
        must be nonempty.

        @throws TypeError If `path` cannot be converted by `os.fsencode`.
        @throws ApoCharmmError If the encoded path is empty, the file cannot be
        read, the atom count or a coordinate record is malformed, a coordinate
        is non-finite, or native host allocation fails.
        @throws RuntimeError If `APOCHARMM_LIBRARY_PATH` is unset or empty, or
        native construction reports success but produces a NULL handle.
        @throws OSError If the configured apoCHARMM C ABI shared library cannot
        be loaded.

        @post On success, this wrapper owns a live coordinate handle.
        @warning An encoded path containing an embedded null byte is truncated
        at that byte by the C-string boundary instead of being rejected.
        """
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
        """
        @brief Returns the parsed atom count.

        @return Non-negative dimensionless Python `int`.
        @throws RuntimeError If this wrapper is closed.
        @throws ApoCharmmError If the native handle is invalid or reports a
        negative atom count.
        """
        _initialize_prototypes()

        c_num_atoms = ctypes.c_size_t()

        lib().apo_charmm_crd_get_num_atoms(ctypes.byref(c_num_atoms), self.handle)

        return int(c_num_atoms.value)

    def getCoordinates(self) -> list[list[float]]:
        """
        @brief Returns a copied nested list of Cartesian coordinates.

        The result contains one three-element list per atom in `[x, y, z]`
        order. Values are Python `float` objects copied from native
        double-precision host storage and use angstroms. An empty coordinate set
        returns an empty list.

        @return Newly allocated `list[list[float]]` with exact shape `(N, 3)`.
        @throws RuntimeError If this wrapper is closed.
        @throws ApoCharmmError If the native count or coordinate getter rejects
        the handle or output buffer.
        @throws OverflowError If `3 * N` cannot be represented as a ctypes
        array length.
        @throws MemoryError If the temporary C buffer or returned Python lists
        cannot be allocated.

        @post Mutating the returned lists does not modify native coordinates.
        """
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
