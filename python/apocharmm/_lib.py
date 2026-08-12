# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

import ctypes
import os

from ._types import FilePath

_library: ctypes.CDLL | None = None


def lib() -> ctypes.CDLL:
    """
    @brief Returns the cached apoCHARMM C ABI shared-library wrapper.

    The first call loads the path named by `APOCHARMM_LIBRARY_PATH` and
    configures `apo_last_error()` to return a raw address. Keeping that result
    as `ctypes.c_void_p` lets `check_status()` distinguish NULL from an empty C
    string and copy the diagnostic before its native lifetime ends. Later calls
    return the same `ctypes.CDLL` object.

    @return The process-local cached `ctypes.CDLL` wrapper. The wrapper is owned
    by this module and must not be closed by callers.
    @throws RuntimeError If `APOCHARMM_LIBRARY_PATH` is unset or empty.
    @throws OSError If the configured shared library cannot be loaded.

    @post `apo_last_error` accepts no arguments and returns a raw pointer value.
    """
    global _library

    if _library is not None:
        return _library

    library_path = os.environ.get("APOCHARMM_LIBRARY_PATH")
    if library_path is None or library_path == "":
        raise RuntimeError("APOCHARMM_LIBRARY_PATH must point to libapocharmm_c.so")

    library = ctypes.CDLL(library_path)

    library.apo_last_error.argtypes = []
    library.apo_last_error.restype = ctypes.c_void_p

    _library = library

    return library


def encode_path(path: FilePath) -> bytes:
    return os.fsencode(path)
