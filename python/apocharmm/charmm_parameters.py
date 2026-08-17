# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

"""
@brief Provides the owning Python wrapper for CHARMM parameter sets.

`CharmmParameters` parses one or more CHARMM `.prm` or `.str` files through
apoCHARMM's C ABI and supplies the resulting native object to managers and
contexts.

@see charmm_parameters
"""

import ctypes

from ._base import _ApoObject
from ._lib import encode_path, lib
from ._types import FilePaths
from .error import configure_status_function

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_charmm_parameters_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p],
        "CharmmParameters construction",
    )

    configure_status_function(
        lib().apo_charmm_parameters_create_from_files,
        [
            ctypes.POINTER(ctypes.c_void_p),
            ctypes.POINTER(ctypes.c_char_p),
            ctypes.c_size_t,
        ],
        "CharmmParameters construction",
    )

    lib().apo_charmm_parameters_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_charmm_parameters_destroy.restype = None

    _prototypes_initialized = True

    return


class CharmmParameters(_ApoObject):
    """
    @anchor python_charmm_parameters
    @brief Owns a parsed native CHARMM parameter set.

    Construct the wrapper from one path or from an ordered `list` or `tuple` of
    paths. Each path may be `str`, `bytes`, `os.PathLike[str]`, or
    `os.PathLike[bytes]`. Arbitrary sequence types are not treated as multiple
    files. Path values are converted with `os.fsencode`, copied through the C
    ABI, and are not retained by the wrapper. File-name extensions are not
    validated.

    The wrapper owns one C handle. `close()` and context-manager exit release
    that handle and are idempotent. Accessing `handle`, or passing the closed
    wrapper to an API that accesses it, raises `RuntimeError`. A native
    ForceManager created from this wrapper retains independent shared ownership
    of the parameter set.

    The Python layer currently exposes construction and lifetime only; native
    parameter maps and incremental file loading are not exposed. Native status
    failures are raised as `ApoCharmmError`. The wrapper provides no internal
    synchronization for concurrent host access or closure.

    @see charmm_parameters
    """

    _destroy_function_name = "apo_charmm_parameters_destroy"

    def __init__(self, paths: FilePaths) -> None:
        """
        @brief Constructs an owning wrapper from one or more parameter files.

        A `list` or `tuple` is interpreted as an ordered file list. Every other
        accepted path object is interpreted as one file. All path values are
        encoded and copied before native parsing; this object does not retain
        the original Python path objects. File-name extensions are not
        validated.

        @param[in] paths One path, or a non-empty `list` or `tuple` of paths.
        Each path must be `str`, `bytes`, `os.PathLike[str]`, or
        `os.PathLike[bytes]` and must encode to a non-empty native path.

        @throws TypeError If `paths`, or an element of a `list` or `tuple`,
        cannot be converted by `os.fsencode`.
        @throws ApoCharmmError If the native C ABI rejects an empty path or file
        list, cannot open or parse a file, or reports another native
        construction failure.
        @throws RuntimeError If `APOCHARMM_LIBRARY_PATH` is unset or empty, or
        if the native function reports success but returns a NULL handle.
        @throws OSError If the configured apoCHARMM C ABI shared library cannot
        be loaded.

        @post On success, this wrapper owns a live native parameter handle.
        @warning An encoded path containing an embedded null byte is currently
        truncated at that byte by the C-string boundary instead of being
        rejected; see @ref charmm_parameters.
        """
        _initialize_prototypes()
        super().__init__()

        handle: ctypes.c_void_p = ctypes.c_void_p()

        if isinstance(paths, (list, tuple)):
            encoded_paths: list[bytes] = [encode_path(path) for path in paths]
            num_paths: int = len(encoded_paths)
            c_num_paths: ctypes.c_size_t = ctypes.c_size_t(num_paths)

            path_array_type = ctypes.c_char_p * len(encoded_paths)
            path_array = path_array_type(*encoded_paths)

            lib().apo_charmm_parameters_create_from_files(
                ctypes.byref(handle), path_array, c_num_paths
            )

            function_name = "apo_charmm_parameters_create_from_files"
        else:
            encoded_path: bytes = encode_path(paths)
            c_path: ctypes.c_char_p = ctypes.c_char_p(encoded_path)

            lib().apo_charmm_parameters_create(ctypes.byref(handle), c_path)

            function_name = "apo_charmm_parameters_create"

        if handle.value is None:
            raise RuntimeError(
                "{} returned success but produced a NULL handle".format(function_name)
            )

        self._handle = handle

        return
