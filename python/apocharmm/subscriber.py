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

    configure_status_function(
        lib().apo_subscriber_set_report_frequency,
        [ctypes.c_void_p, ctypes.c_int],
        "Subscriber.setReportFrequency(report_frequency)",
    )

    configure_status_function(
        lib().apo_subscriber_get_report_frequency,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "Subscriber.getReportFrequency()",
    )

    _prototypes_initialized = True

    return


class Subscriber(_ApoObject):
    def __init__(self) -> None:
        super().__init__()
        self._subscriber_handle: ctypes.c_void_p = ctypes.c_void_p()
        return

    @property
    def subscriber_handle(self) -> ctypes.c_void_p:
        if self._subscriber_handle is None or self._subscriber_handle.value is None:
            raise RuntimeError("apoCHARMM subscriber object has been destroyed")
        return self._subscriber_handle

    def setReportFrequency(self, report_frequency: int) -> None:
        _initialize_prototypes()

        if report_frequency < -(2**31) or report_frequency > 2**31 - 1:
            raise ValueError("report_frequency must fit in int")

        c_report_frequency: ctypes.c_int = ctypes.c_int(report_frequency)

        lib().apo_subscriber_set_report_frequency(
            self.subscriber_handle, c_report_frequency
        )

        return

    def getReportFrequency(self) -> int:
        _initialize_prototypes()

        c_report_frequency = ctypes.c_int()

        lib().apo_subscriber_get_report_frequency(
            ctypes.byref(c_report_frequency), self.subscriber_handle
        )

        return int(c_report_frequency.value)
