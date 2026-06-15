# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

from collections.abc import Sequence
import ctypes
from typing import Protocol, cast

from ._base import _ApoObject
from ._lib import lib
from .enums import PeriodicBoundaryCondition
from .error import check_status

from .charmm_parameters import CharmmParameters
from .charmm_psf import CharmmPsf

_prototypes_initialized: bool = False


class _SubscribableForce(Protocol):
    def _subscribe_to_force_manager(
        self, force_manager: "ForceManager", force_tag: str | None = None
    ) -> None:
        return

    def _unsubscribe_from_force_manager(self, force_manager: "ForceManager") -> None:
        return


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    lib().apo_force_manager_create.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_create.restype = ctypes.c_int

    lib().apo_force_manager_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_force_manager_destroy.restype = None

    lib().apo_force_manager_set_box_dimensions.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib().apo_force_manager_set_box_dimensions.restype = ctypes.c_int

    lib().apo_force_manager_set_periodic_boundary_condition.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
    ]
    lib().apo_force_manager_set_periodic_boundary_condition.restype = ctypes.c_int

    lib().apo_force_manager_get_periodic_boundary_condition.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_periodic_boundary_condition.restype = ctypes.c_int

    _prototypes_initialized = True

    return


class ForceManager(_ApoObject):
    _destroy_function_name = "apo_force_manager_destroy"

    def __init__(self, psf: CharmmPsf, parameters: CharmmParameters) -> None:
        _initialize_prototypes()
        super().__init__()

        if not isinstance(psf, CharmmPsf):
            raise TypeError("ForceManager expects a CharmmPsf")

        if not isinstance(parameters, CharmmParameters):
            raise TypeError("ForceManager expects a CharmmParameters")

        handle = ctypes.c_void_p()

        status = lib().apo_force_manager_create(
            ctypes.byref(handle), psf.handle, parameters.handle
        )

        check_status(status, "ForceManager construction failed")

        if handle.value is None:
            raise RuntimeError(
                "apo_force_manager_create returned success but produced a NULL handle"
            )

        self._handle = handle

        self._psf: CharmmPsf = psf
        self._parameters: CharmmParameters = parameters
        self._subscribed_forces: list[_SubscribableForce] = []

        return

    def close(self) -> None:
        super().close()

        if hasattr(self, "_subscribed_forces"):
            self._subscribed_forces = []

        return

    def subscribe(
        self, force: _SubscribableForce, force_tag: str | None = None
    ) -> None:
        _initialize_prototypes()

        if force_tag is not None and not isinstance(force_tag, str):
            raise TypeError("force_tag must be a str")

        if force_tag == "":
            raise ValueError("force_tag must not be empty")

        subscribe_method = getattr(force, "_subscribe_to_force_manager", None)
        if not callable(subscribe_method):
            raise TypeError(
                "ForceManager.subscribe expects an object with _subscribe_to_force_manager(force_manager, force_tag)"
            )

        if any(existing_force is force for existing_force in self._subscribed_forces):
            raise ValueError("force is already subscribed to this ForceManager")

        subscribable_force: _SubscribableForce = cast(_SubscribableForce, force)
        subscribable_force._subscribe_to_force_manager(self, force_tag)

        self._subscribed_forces.append(subscribable_force)

        return

    def unsubscribe(self, force: _SubscribableForce) -> None:
        _initialize_prototypes()

        unsubscribe_method = getattr(force, "_unsubscribe_from_force_manager", None)
        if not callable(unsubscribe_method):
            raise TypeError(
                "ForceManager.unsubscribe expects an object with _unsubscribe_from_force_manager(force_manager)"
            )

        subscribable_force: _SubscribableForce = cast(_SubscribableForce, force)
        subscribable_force._unsubscribe_from_force_manager(self)

        for index, existing_force in enumerate(self._subscribed_forces):
            if existing_force is force:
                del self._subscribed_forces[index]
                break

        return

    def setBoxDimensions(self, box_dimensions: Sequence[float]) -> None:
        _initialize_prototypes()

        box_values: list[float] = [float(value) for value in box_dimensions]

        c_buffer_type = ctypes.c_double * len(box_values)
        c_buffer = c_buffer_type(*box_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(box_values))

        status = lib().apo_force_manager_set_box_dimensions(
            self.handle, c_buffer, c_buffer_len
        )

        check_status(status, "ForceManager.setBoxDimensions() failed")

        return

    def setPeriodicBoundaryCondition(
        self, pbc: PeriodicBoundaryCondition | int
    ) -> None:
        _initialize_prototypes()

        try:
            pbc_value: PeriodicBoundaryCondition = PeriodicBoundaryCondition(pbc)
        except ValueError as exc:
            raise ValueError(f"invalid pbc: {pbc!r}") from exc

        c_pbc: ctypes.c_int = ctypes.c_int(int(pbc_value))

        status = lib().apo_force_manager_set_periodic_boundary_condition(
            self.handle, c_pbc
        )

        check_status(
            status,
            "ForceManager.setPeriodicBoundaryCondition(pbc) failed",
        )

        return

    def getPeriodicBoundaryCondition(self) -> PeriodicBoundaryCondition:
        _initialize_prototypes()

        c_pbc = ctypes.c_int()

        status = lib().apo_force_manager_get_periodic_boundary_condition(
            ctypes.byref(c_pbc), self.handle
        )

        check_status(
            status,
            "ForceManager.getPeriodicBoundaryCondition() failed",
        )

        return PeriodicBoundaryCondition(c_pbc.value)
