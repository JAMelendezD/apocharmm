# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# ENDLICENSE

import ctypes
import unittest

from unittest.mock import patch

from apocharmm._base import _ApoObject


class _FailingLibrary:
    def apo_test_destroy(self, handle: ctypes.c_void_p) -> None:
        del handle
        raise RuntimeError("destroy failed")


class _TestObject(_ApoObject):
    _destroy_function_name = "apo_test_destroy"


class TestApoObjectCleanup(unittest.TestCase):
    def test_close_invalidates_handle_before_destroy_failure(self) -> None:
        obj = _TestObject()
        obj._handle = ctypes.c_void_p(1)

        with patch("apocharmm._base.lib", return_value=_FailingLibrary()):
            with self.assertRaisesRegex(RuntimeError, "destroy failed"):
                obj.close()

        self.assertIsNotNone(obj._handle)
        self.assertIsNone(obj._handle.value)

    def test_del_suppresses_destroy_failure(self) -> None:
        obj = _TestObject()
        obj._handle = ctypes.c_void_p(1)

        with patch("apocharmm._base.lib", return_value=_FailingLibrary()):
            self.assertIsNone(obj.__del__())

        self.assertIsNotNone(obj._handle)
        self.assertIsNone(obj._handle.value)


if __name__ == "__main__":
    unittest.main()
