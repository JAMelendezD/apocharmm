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
from .enums import PeriodicBoundaryCondition, VdwType
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

    lib().apo_force_manager_set_kappa.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib().apo_force_manager_set_kappa.restype = ctypes.c_int

    lib().apo_force_manager_set_cutoff.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib().apo_force_manager_set_cutoff.restype = ctypes.c_int

    lib().apo_force_manager_set_ctonnb.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib().apo_force_manager_set_ctonnb.restype = ctypes.c_int

    lib().apo_force_manager_set_ctofnb.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib().apo_force_manager_set_ctofnb.restype = ctypes.c_int

    lib().apo_force_manager_set_fft_grid.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
    ]
    lib().apo_force_manager_set_fft_grid.restype = ctypes.c_int

    lib().apo_force_manager_set_pme_spline_order.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
    ]
    lib().apo_force_manager_set_pme_spline_order.restype = ctypes.c_int

    lib().apo_force_manager_set_periodic_boundary_condition.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
    ]
    lib().apo_force_manager_set_periodic_boundary_condition.restype = ctypes.c_int

    lib().apo_force_manager_set_vdw_type.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib().apo_force_manager_set_vdw_type.restype = ctypes.c_int

    lib().apo_force_manager_set_print_energy_decomposition.argtypes = [
        ctypes.c_void_p,
        ctypes.c_bool,
    ]
    lib().apo_force_manager_set_print_energy_decomposition.restype = ctypes.c_int

    lib().apo_force_manager_get_num_atoms.argtypes = [
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_num_atoms.restype = ctypes.c_int

    lib().apo_force_manager_is_initialized.argtypes = [
        ctypes.POINTER(ctypes.c_bool),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_is_initialized.restype = ctypes.c_int

    lib().apo_force_manager_get_box_dimensions.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_box_dimensions.restype = ctypes.c_int

    lib().apo_force_manager_get_kappa.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_kappa.restype = ctypes.c_int

    lib().apo_force_manager_get_cutoff.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_cutoff.restype = ctypes.c_int

    lib().apo_force_manager_get_ctonnb.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_ctonnb.restype = ctypes.c_int

    lib().apo_force_manager_get_ctofnb.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_ctofnb.restype = ctypes.c_int

    lib().apo_force_manager_get_fft_grid.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_fft_grid.restype = ctypes.c_int

    lib().apo_force_manager_get_pme_spline_order.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_pme_spline_order.restype = ctypes.c_int

    lib().apo_force_manager_get_periodic_boundary_condition.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_periodic_boundary_condition.restype = ctypes.c_int

    lib().apo_force_manager_get_vdw_type.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_get_vdw_type.restype = ctypes.c_int

    lib().apo_force_manager_is_composite.argtypes = [
        ctypes.POINTER(ctypes.c_bool),
        ctypes.c_void_p,
    ]
    lib().apo_force_manager_is_composite.restype = ctypes.c_int

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

        self._psf: CharmmPsf | None = psf
        self._parameters: CharmmParameters | None = parameters
        self._subscribed_forces: list[_SubscribableForce] = []

        return

    @classmethod
    def _from_handle(cls, handle: ctypes.c_void_p) -> "ForceManager":
        _initialize_prototypes()

        if handle.value is None:
            raise RuntimeError("ForceManager._from_handle received a NULL handle")

        obj = cls.__new__(cls)
        _ApoObject.__init__(obj)

        obj._handle = handle
        obj._psf = None
        obj._parameters = None
        obj._subscribed_forces = []

        return obj

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

        index_to_remove: int | None = None
        for index, existing_force in enumerate(self._subscribed_forces):
            if existing_force is force:
                index_to_remove = index
                break

        if index_to_remove is None:
            raise ValueError("force is not subscribed to this ForceManager")

        subscribable_force: _SubscribableForce = cast(_SubscribableForce, force)
        subscribable_force._unsubscribe_from_force_manager(self)

        del self._subscribed_forces[index_to_remove]

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

        check_status(status, "ForceManager.setBoxDimensions(box_dimensions) failed")

        return

    def setKappa(self, kappa: float) -> None:
        _initialize_prototypes()

        c_kappa: ctypes.c_double = ctypes.c_double(kappa)

        status = lib().apo_force_manager_set_kappa(self.handle, c_kappa)

        check_status(status, "ForceManager.setKappa(kappa) failed")

        return

    def setCutoff(self, cutoff: float) -> None:
        _initialize_prototypes()

        c_cutoff: ctypes.c_double = ctypes.c_double(cutoff)

        status = lib().apo_force_manager_set_cutoff(self.handle, c_cutoff)

        check_status(status, "ForceManager.setCutoff(cutoff) failed")

        return

    def setCtonnb(self, ctonnb: float) -> None:
        _initialize_prototypes()

        c_ctonnb: ctypes.c_double = ctypes.c_double(ctonnb)

        status = lib().apo_force_manager_set_ctonnb(self.handle, c_ctonnb)

        check_status(status, "ForceManager.setCtonnb(ctonnb) failed")

        return

    def setCtofnb(self, ctofnb: float) -> None:
        _initialize_prototypes()

        c_ctofnb: ctypes.c_double = ctypes.c_double(ctofnb)

        status = lib().apo_force_manager_set_ctofnb(self.handle, c_ctofnb)

        check_status(status, "ForceManager.setCtofnb(ctofnb) failed")

        return

    def setFFTGrid(self, grid: Sequence[int]) -> None:
        _initialize_prototypes()

        grid_values: list[int] = [int(value) for value in grid]

        c_buffer_type = ctypes.c_int * len(grid_values)
        c_buffer = c_buffer_type(*grid_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(grid_values))

        status = lib().apo_force_manager_set_fft_grid(
            self.handle, c_buffer, c_buffer_len
        )

        check_status(status, "ForceManager.setFFTGrid(grid) failed")

        return

    def setPmeSplineOrder(self, order: int) -> None:
        _initialize_prototypes()

        c_order: ctypes.c_int = ctypes.c_int(order)

        status = lib().apo_force_manager_set_pme_spline_order(self.handle, c_order)

        check_status(status, "ForceManager.setPmeSplineOrder(order) failed")

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

    def setVdwType(self, vdw_type: VdwType | int) -> None:
        _initialize_prototypes()

        try:
            vdw_type_value: VdwType = VdwType(vdw_type)
        except ValueError as exc:
            raise ValueError(f"invalid vdw_type: {vdw_type!r}") from exc

        c_vdw_type: ctypes.c_int = ctypes.c_int(int(vdw_type_value))

        status = lib().apo_force_manager_set_vdw_type(self.handle, c_vdw_type)

        check_status(status, "ForceManager.setVdwType(vdw_type) failed")

        return

    def setPrintEnergyDecomposition(self, flag: bool = True) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        status = lib().apo_force_manager_set_print_energy_decomposition(
            self.handle, c_flag
        )

        check_status(status, "ForceManager.setPrintEnergyDecomposition(flag) failed")

        return

    def getNumAtoms(self) -> int:
        _initialize_prototypes()

        c_num_atoms: ctypes.c_size_t = ctypes.c_size_t()

        status = lib().apo_force_manager_get_num_atoms(
            ctypes.byref(c_num_atoms), self.handle
        )

        check_status(status, "ForceManager.getNumAtoms() failed")

        return int(c_num_atoms.value)

    def isInitialized(self) -> bool:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool()

        status = lib().apo_force_manager_is_initialized(
            ctypes.byref(c_flag), self.handle
        )

        check_status(status, "ForceManager.isInitialized() failed")

        return bool(c_flag.value)

    def getBoxDimensions(self) -> tuple[float, float, float]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_double * 3
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(3)

        status = lib().apo_force_manager_get_box_dimensions(
            c_buffer, c_buffer_len, self.handle
        )

        check_status(status, "ForceManager.getBoxDimensions() failed")

        return (float(c_buffer[0]), float(c_buffer[1]), float(c_buffer[2]))

    def getKappa(self) -> float:
        _initialize_prototypes()

        c_kappa: ctypes.c_double = ctypes.c_double()

        status = lib().apo_force_manager_get_kappa(ctypes.byref(c_kappa), self.handle)

        check_status(status, "ForceManager.getKappa() failed")

        return float(c_kappa.value)

    def getCutoff(self) -> float:
        _initialize_prototypes()

        c_cutoff: ctypes.c_double = ctypes.c_double()

        status = lib().apo_force_manager_get_cutoff(ctypes.byref(c_cutoff), self.handle)

        check_status(status, "ForceManager.getCutoff() failed")

        return float(c_cutoff.value)

    def getCtonnb(self) -> float:
        _initialize_prototypes()

        c_ctonnb: ctypes.c_double = ctypes.c_double()

        status = lib().apo_force_manager_get_ctonnb(ctypes.byref(c_ctonnb), self.handle)

        check_status(status, "ForceManager.getCtonnb() failed")

        return float(c_ctonnb.value)

    def getCtofnb(self) -> float:
        _initialize_prototypes()

        c_ctofnb: ctypes.c_double = ctypes.c_double()

        status = lib().apo_force_manager_get_ctofnb(ctypes.byref(c_ctofnb), self.handle)

        check_status(status, "ForceManager.getCtofnb() failed")

        return float(c_ctofnb.value)

    def getFFTGrid(self) -> tuple[int, int, int]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_int * 3
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(3)

        status = lib().apo_force_manager_get_fft_grid(
            c_buffer, c_buffer_len, self.handle
        )

        check_status(status, "ForceManager.getFFTGrid() failed")

        return (int(c_buffer[0]), int(c_buffer[1]), int(c_buffer[2]))

    def getPmeSplineOrder(self) -> int:
        _initialize_prototypes()

        c_order: ctypes.c_int = ctypes.c_int()

        status = lib().apo_force_manager_get_pme_spline_order(
            ctypes.byref(c_order), self.handle
        )

        check_status(status, "ForceManager.getPmeSplineOrder() failed")

        return float(c_order.value)

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

    def getVdwType(self) -> VdwType:
        _initialize_prototypes()

        c_vdw_type: ctypes.c_int = ctypes.c_int()

        status = lib().apo_force_manager_get_vdw_type(
            ctypes.byref(c_vdw_type), self.handle
        )

        check_status(status, "ForceManager.getVdwType() failed")

        return VdwType(c_vdw_type.value)

    def isComposite(self) -> bool:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool()

        status = lib().apo_force_manager_is_composite(ctypes.byref(c_flag), self.handle)

        check_status(status, "ForceManager.isComposite() failed")

        return bool(c_flag.value)
