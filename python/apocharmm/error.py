# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from collections.abc import Callable, Sequence
import ctypes
from typing import Protocol

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


class _StatusFunction(Protocol):
    argtypes: Sequence[object] | None
    restype: object
    errcheck: Callable[[int, "_StatusFunction", tuple[object, ...]], int] | None


class ApoCharmmError(RuntimeError):
    def __init__(self, status: int, context: str, native_diagnostic: str) -> None:
        self._status: int = int(status)
        self._status_name: str = _STATUS_NAMES.get(self._status, _UNKNOWN_STATUS_NAME)
        self._context: str = context
        self._native_diagnostic: str = native_diagnostic

        rendered_native_diagnostic: str = native_diagnostic.replace(
            "\r\n", "\n"
        ).replace("\r", "\n")
        header: str = "{} [{} ({})]".format(
            self._context, self._status_name, self._status
        )

        separator: str = ": "
        if "\n" in rendered_native_diagnostic:
            separator = "" if rendered_native_diagnostic.startswith("\n") else "\n"

        self._message: str = "{}{}{}".format(
            header, separator, rendered_native_diagnostic
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


def configure_status_function(
    function: _StatusFunction, argtypes: Sequence[object], context: str
) -> None:
    def raise_on_error(
        result: int, _function: _StatusFunction, _arguments: tuple[object, ...]
    ) -> int:
        check_status(result, context)
        return result

    function.argtypes = list(argtypes)
    function.restype = ctypes.c_int
    function.errcheck = raise_on_error

    return
