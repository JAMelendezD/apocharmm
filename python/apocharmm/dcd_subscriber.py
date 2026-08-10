# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

import ctypes

from ._lib import encode_path, lib
from ._types import FilePath
from .error import configure_status_function
from .subscriber import Subscriber

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_dcd_subscriber_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p],
        "DcdSubscriber construction",
    )

    configure_status_function(
        lib().apo_dcd_subscriber_create_with_report_frequency,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_char_p, ctypes.c_int],
        "DcdSubscriber construction",
    )

    lib().apo_dcd_subscriber_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_dcd_subscriber_destroy.restype = None

    configure_status_function(
        lib().apo_dcd_subscriber_as_subscriber,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p],
        "DcdSubscriber base-subscriber conversion",
    )

    _prototypes_initialized = True

    return


class DcdSubscriber(Subscriber):
    _destroy_function_name = "apo_dcd_subscriber_destroy"

    def __init__(self, path: FilePath, report_frequency: int | None = None) -> None:
        _initialize_prototypes()
        super().__init__()

        handle: ctypes.c_void_p = ctypes.c_void_p()

        encoded_path: bytes = encode_path(path)
        c_path: ctypes.c_char_p = ctypes.c_char_p(encoded_path)

        if report_frequency is None:
            function_name = "apo_dcd_subscriber_create"
            lib().apo_dcd_subscriber_create(ctypes.byref(handle), c_path)
        else:
            if report_frequency < -(2**31) or report_frequency > 2**31 - 1:
                raise ValueError("report_frequency must fit in int")

            c_report_frequency: ctypes.c_int = ctypes.c_int(report_frequency)
            function_name = "apo_dcd_subscriber_create_with_report_frequency"

            lib().apo_dcd_subscriber_create_with_report_frequency(
                ctypes.byref(handle), c_path, c_report_frequency
            )

        if handle.value is None:
            raise RuntimeError(
                "{} returned success but produced a NULL handle".format(function_name)
            )

        self._handle = handle

        subscriber_handle: ctypes.c_void_p = ctypes.c_void_p()

        lib().apo_dcd_subscriber_as_subscriber(
            ctypes.byref(subscriber_handle), self.handle
        )

        if subscriber_handle.value is None:
            raise RuntimeError(
                "apo_dcd_subscriber_as_subscriber returned success but produced a NULL handle"
            )

        self._subscriber_handle = subscriber_handle
        self._path: FilePath = path

        return

    def close(self) -> None:
        super().close()
        self._subscriber_handle = ctypes.c_void_p()
        return
