# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from __future__ import annotations

from collections.abc import Callable
import ctypes
from pathlib import Path
import sys
from typing import cast

import apocharmm as apo
import apocharmm._lib as apo_lib
from apocharmm.error import check_status, configure_status_function

import apo_test_helpers as apo_test

NATIVE_DIAGNOSTIC_FALLBACK: str = "Unknown apoCHARMM C API error"


class _FakeErrorLibrary:
    def __init__(self, message_address: int | None) -> None:
        self._message_address: int | None = message_address
        return

    def apo_last_error(self) -> int | None:
        return self._message_address


def capture_check_status_error(
    status: int, context: str, native_message: bytes | None
) -> apo.ApoCharmmError:
    message_buffer: ctypes.Array[ctypes.c_char] | None = None
    message_address: int | None = None

    if native_message is not None:
        message_buffer = ctypes.create_string_buffer(native_message)
        message_address = ctypes.addressof(message_buffer)

    fake_library = _FakeErrorLibrary(message_address)
    original_lib: Callable[[], ctypes.CDLL] = apo_lib.lib

    def fake_lib() -> ctypes.CDLL:
        return cast(ctypes.CDLL, fake_library)

    apo_lib.lib = fake_lib
    try:
        error = apo_test.expect_exception(
            f"check_status rejects status {status}",
            apo.ApoCharmmError,
            lambda: check_status(status, context),
        )
    finally:
        apo_lib.lib = original_lib

    return error


def check_exception_model() -> None:
    print("Checking ApoCharmError fields and rendered message...")

    context: str = "Python operation failed"
    native_diagnostic: str = "apo_native_operation: native diagnostic"
    expected_message: str = (
        "Python operation failed [APO_STATUS_INVALID_ARGUMENT (1)]: apo_native_operation: native diagnostic"
    )

    error = apo.ApoCharmmError(
        apo.APO_STATUS_INVALID_ARGUMENT, context, native_diagnostic
    )

    apo_test.assert_equal(
        "ApoCharmmError RuntimeError subclass",
        issubclass(apo.ApoCharmmError, RuntimeError),
        True,
    )
    apo_test.assert_equal(
        "ApoCharmmError RuntimeError instance", isinstance(error, RuntimeError), True
    )
    apo_test.assert_equal(
        "Known numeric status", error.status, apo.APO_STATUS_INVALID_ARGUMENT
    )
    apo_test.assert_equal(
        "Known status name", error.status_name, "APO_STATUS_INVALID_ARGUMENT"
    )
    apo_test.assert_equal("Python operation context", error.context, context)
    apo_test.assert_equal(
        "Native diagnostic", error.native_diagnostic, native_diagnostic
    )
    apo_test.assert_equal("Stored rendered message", error.message, expected_message)
    apo_test.assert_equal("Exception string", str(error), expected_message)
    apo_test.assert_equal(
        "Rendered context occurence count", error.message.count(context), 1
    )

    unknown_status: int = 987654321
    unknown_error = apo.ApoCharmmError(
        unknown_status, "Unknown-status operation", "unknown-status diagnostic"
    )

    apo_test.assert_equal(
        "Unknown numeric status", unknown_error.status, unknown_status
    )
    apo_test.assert_equal(
        "Unknown status name", unknown_error.status_name, "APO_STATUS_UNKNOWN"
    )
    apo_test.assert_equal(
        "Unknown rendered message",
        unknown_error.message,
        "Unknown-status operation [APO_STATUS_UNKNOWN (987654321)]: unknown-status diagnostic",
    )

    return


def check_native_message_handling() -> None:
    print("Checking native diagnostic fallback and decoding...")

    null_error = capture_check_status_error(
        apo.APO_STATUS_RUNTIME_ERROR, "Null native message operation", None
    )
    apo_test.assert_equal(
        "Null native message fallback",
        null_error.native_diagnostic,
        NATIVE_DIAGNOSTIC_FALLBACK,
    )
    apo_test.assert_equal(
        "Null native message rendered fallback",
        null_error.message,
        "Null native message operation [APO_STATUS_RUNTIME_ERROR (2)]: Unknown apoCHARMM C API error",
    )

    empty_error = capture_check_status_error(
        apo.APO_STATUS_RUNTIME_ERROR, "Empty native message operation", b""
    )
    apo_test.assert_equal(
        "Empty native message fallback",
        empty_error.native_diagnostic,
        NATIVE_DIAGNOSTIC_FALLBACK,
    )
    apo_test.assert_equal(
        "Empty native message rendered fallback",
        empty_error.message,
        "Empty native message operation [APO_STATUS_RUNTIME_ERROR (2)]: Unknown apoCHARMM C API error",
    )

    invalid_utf8_error = capture_check_status_error(
        apo.APO_STATUS_RUNTIME_ERROR,
        "Invalid UTF-8 operation",
        b"native \xff diagnostic",
    )
    apo_test.assert_equal(
        "Invalid UTF-8 replacement decoding",
        invalid_utf8_error.native_diagnostic,
        "native \ufffd diagnostic",
    )
    apo_test.assert_equal(
        "Invalid UTF-8 rendered message",
        invalid_utf8_error.message,
        "Invalid UTF-8 operation [APO_STATUS_RUNTIME_ERROR (2)]: native \ufffd diagnostic",
    )

    return


def check_configured_status_function() -> None:
    print("Checking configured apo_status function behavior...")

    success_function = ctypes.CFUNCTYPE(ctypes.c_int)(lambda: apo.APO_STATUS_OK)
    configure_status_function(success_function, [], "Configured success operation")

    apo_test.assert_equal("Configured status argtypes", success_function.argtypes, [])
    apo_test.assert_equal(
        "Configured status restype", success_function.restype, ctypes.c_int
    )
    apo_test.assert_equal(
        "Configured status errcheck installed",
        callable(success_function.errcheck),
        True,
    )
    apo_test.assert_equal(
        "Configured status success result", success_function(), apo.APO_STATUS_OK
    )

    unknown_status: int = 987654321
    context: str = "Configured status operation"
    message_buffer: ctypes.Array[ctypes.c_char] = ctypes.create_string_buffer(
        b"configured native diagnostic"
    )
    fake_library = _FakeErrorLibrary(ctypes.addressof(message_buffer))
    original_lib: Callable[[], ctypes.CDLL] = apo_lib.lib

    def fake_lib() -> ctypes.CDLL:
        return cast(ctypes.CDLL, fake_library)

    status_function = ctypes.CFUNCTYPE(ctypes.c_int)(lambda: unknown_status)
    configure_status_function(status_function, [], context)

    apo_lib.lib = fake_lib
    try:
        error = apo_test.expect_exception(
            "Configured status function preserves unknown status",
            apo.ApoCharmmError,
            lambda: status_function(),
        )
    finally:
        apo_lib.lib = original_lib

    apo_test.assert_equal(
        "Configured unknown numeric status", error.status, unknown_status
    )
    apo_test.assert_equal(
        "Configured unknown status name", error.status_name, "APO_STATUS_UNKNOWN"
    )
    apo_test.assert_equal("Configured operation context", error.context, context)
    apo_test.assert_equal(
        "Configured native diagnostic",
        error.native_diagnostic,
        "configured native diagnostic",
    )
    apo_test.assert_equal(
        "Configured rendered message",
        error.message,
        "Configured status operation [APO_STATUS_UNKNOWN (987654321)]: configured native diagnostic",
    )

    return


def check_python_and_native_failures(repo_root: Path) -> None:
    print("Checking Python validation and native wrapper failures...")

    psf_path: str = apo_test.require_file(repo_root / "test/data/nacl_pair.psf")
    psf = apo.CharmmPsf(psf_path)
    selector = apo.AtomSelector(psf)

    type_error = apo_test.expect_exception(
        "AtomSelector.select preserves Python TypeError",
        TypeError,
        lambda: selector.select(cast(str, 1)),
    )
    apo_test.assert_equal(
        "Python type-check exception type", type(type_error), TypeError
    )

    empty_selection_error: apo.ApoCharmmError = apo_test.expect_exception(
        "AtomSelector.select maps an empty selection to invalid argument",
        apo.ApoCharmmError,
        lambda: selector.select(""),
    )
    apo_test.assert_equal(
        "empty selection exception type",
        type(empty_selection_error),
        apo.ApoCharmmError,
    )
    apo_test.assert_equal(
        "empty selection numeric status",
        empty_selection_error.status,
        apo.APO_STATUS_INVALID_ARGUMENT,
    )
    apo_test.assert_equal(
        "empty selection status name",
        empty_selection_error.status_name,
        "APO_STATUS_INVALID_ARGUMENT",
    )
    apo_test.assert_equal(
        "empty selection context",
        empty_selection_error.context,
        "AtomSelector.select(selection_string)",
    )
    apo_test.assert_equal(
        "empty selection native diagnostic",
        empty_selection_error.native_diagnostic,
        "apo_atom_selector_select: selection_string is NULL or empty",
    )

    invalid_range_error: apo.ApoCharmmError = apo_test.expect_exception(
        "AtomSelector.select maps an invalid BYNU range to invalid argument",
        apo.ApoCharmmError,
        lambda: selector.select("bynu A:C"),
    )
    apo_test.assert_equal(
        "invalid BYNU range exception type",
        type(invalid_range_error),
        apo.ApoCharmmError,
    )
    apo_test.assert_equal(
        "invalid BYNU range numeric status",
        invalid_range_error.status,
        apo.APO_STATUS_INVALID_ARGUMENT,
    )
    apo_test.assert_equal(
        "invalid BYNU range status name",
        invalid_range_error.status_name,
        "APO_STATUS_INVALID_ARGUMENT",
    )
    apo_test.assert_equal(
        "invalid BYNU range context",
        invalid_range_error.context,
        "AtomSelector.select(selection_string)",
    )
    apo_test.assert_equal(
        "invalid BYNU range diagnostic contains raw message",
        "BYNU range requires integer atom numbers"
        in invalid_range_error.native_diagnostic,
        True,
    )
    apo_test.assert_equal(
        "invalid BYNU range rendered failed text count",
        invalid_range_error.message.count("failed"),
        0,
    )

    invalid_argument_error = apo_test.expect_exception(
        "CharmmParameters maps a native invalid argument",
        apo.ApoCharmmError,
        lambda: apo.CharmmParameters([]),
    )
    apo_test.assert_equal(
        "native invalid-argument exception type",
        type(invalid_argument_error),
        apo.ApoCharmmError,
    )
    apo_test.assert_equal(
        "native invalid-argument numeric status",
        invalid_argument_error.status,
        apo.APO_STATUS_INVALID_ARGUMENT,
    )
    apo_test.assert_equal(
        "native invalid-argument status name",
        invalid_argument_error.status_name,
        "APO_STATUS_INVALID_ARGUMENT",
    )
    apo_test.assert_equal(
        "native invalid-argument context",
        invalid_argument_error.context,
        "CharmmParameters construction",
    )
    apo_test.assert_equal(
        "end-to-end native diagnostic",
        invalid_argument_error.native_diagnostic,
        "apo_charmm_parameters_create_from_files: "
        "paths must contain at least one parameter file",
    )
    apo_test.assert_equal(
        "end-to-end rendered context occurrence count",
        invalid_argument_error.message.count(invalid_argument_error.context),
        1,
    )
    apo_test.assert_equal(
        "end-to-end rendered native function occurrence count",
        invalid_argument_error.message.count("apo_charmm_parameters_create_from_files"),
        1,
    )
    apo_test.assert_equal(
        "end-to-end rendered failed text count",
        invalid_argument_error.message.count("failed"),
        0,
    )

    return


def main(argc: int, argv: list[str]) -> int:
    repo_root: Path = Path(argv[1]) if argc > 1 else Path(".")

    check_exception_model()
    check_native_message_handling()
    check_configured_status_function()
    check_python_and_native_failures(repo_root)

    print("\033[32m" + "PASS: Python error API tests completed." + "\033[0m")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(len(sys.argv), sys.argv))
