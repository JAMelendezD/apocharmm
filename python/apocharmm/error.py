# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

import ctypes

APO_STATUS_OK: int = 0
APO_STATUS_INVALID_ARGUMENT: int = 1
APO_STATUS_RUNTIME_ERROR: int = 2
APO_STATUS_CUDA_ERROR: int = 3
APO_STATUS_NOT_INITIALIZED: int = 4
APO_STATUS_NOT_IMPLEMENTED: int = 5


_STATUS_NAMES: dict[int, str] = {
    APO_STATUS_OK: "APO_STATUS_OK",
    APO_STATUS_INVALID_ARGUMENT: "APO_STATUS_INVALID_ARGUMENT",
    APO_STATUS_RUNTIME_ERROR: "APO_STATUS_RUNTIME_ERROR",
    APO_STATUS_CUDA_ERROR: "APO_STATUS_CUDA_ERROR",
    APO_STATUS_NOT_INITIALIZED: "APO_STATUS_NOT_INITIALIZED",
    APO_STATUS_NOT_IMPLEMENTED: "APO_STATUS_NOT_IMPLEMENTED",
}
_UNKNOWN_STATUS_NAME: str = "APO_STATUS_UNKNOWN"
_NATIVE_DIAGNOSTIC_FALLBACK: str = "Unknown apoCHARMM C API error"


class ApoCharmmError(RuntimeError):
    def __init__(self, status: int, context: str, native_diagnostic: str) -> None:
        self._status: int = int(status)
        self._status_name: str = _STATUS_NAMES.get(self._status, _UNKNOWN_STATUS_NAME)
        self._context: str = context
        self._native_diagnostic: str = native_diagnostic
        self._message: str = "{} [{} ({})]: {}".format(
            self._context, self._status_name, self._status, self._native_diagnostic
        )
        super().__init__(self._message)
        return

    @property
    def status(self) -> int:
        return self._status

    @property
    def status_name(self) -> str:
        return self._status_name

    @property
    def context(self) -> str:
        return self._context

    @property
    def native_diagnostic(self) -> str:
        return self._native_diagnostic

    @property
    def message(self) -> str:
        return self._message


def check_status(status: int, context: str) -> None:
    if status == APO_STATUS_OK:
        return

    from ._lib import lib

    message_address: int | None = lib().apo_last_error()
    if message_address is None:
        native_diagnostic: str = _NATIVE_DIAGNOSTIC_FALLBACK
    else:
        message_bytes: bytes = ctypes.string_at(message_address)
        native_diagnostic = message_bytes.decode("utf-8", errors="replace")
        if native_diagnostic == "":
            native_diagnostic = _NATIVE_DIAGNOSTIC_FALLBACK

    raise ApoCharmmError(status, context, native_diagnostic)
