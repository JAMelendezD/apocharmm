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
from .error import configure_status_function

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

    configure_status_function(
        lib().apo_force_manager_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p, ctypes.c_void_p],
        "ForceManager construction",
    )

    lib().apo_force_manager_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_force_manager_destroy.restype = None

    configure_status_function(
        lib().apo_force_manager_set_box_dimensions,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "ForceManager.setBoxDimensions(box_dimensions)",
    )

    configure_status_function(
        lib().apo_force_manager_set_kappa,
        [ctypes.c_void_p, ctypes.c_double],
        "ForceManager.setKappa(kappa)",
    )

    configure_status_function(
        lib().apo_force_manager_set_cutoff,
        [ctypes.c_void_p, ctypes.c_double],
        "ForceManager.setCutoff(cutoff)",
    )

    configure_status_function(
        lib().apo_force_manager_set_ctonnb,
        [ctypes.c_void_p, ctypes.c_double],
        "ForceManager.setCtonnb(ctonnb)",
    )

    configure_status_function(
        lib().apo_force_manager_set_ctofnb,
        [ctypes.c_void_p, ctypes.c_double],
        "ForceManager.setCtofnb",
    )

    configure_status_function(
        lib().apo_force_manager_set_fft_grid,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_size_t],
        "ForceManager.setFFTGrid(grid)",
    )

    configure_status_function(
        lib().apo_force_manager_set_pme_spline_order,
        [ctypes.c_void_p, ctypes.c_int],
        "ForceManager.setPmeSplineOrder(order)",
    )

    configure_status_function(
        lib().apo_force_manager_set_periodic_boundary_condition,
        [ctypes.c_void_p, ctypes.c_int],
        "ForceManager.setPeriodicBoundaryCondition(pbc)",
    )

    configure_status_function(
        lib().apo_force_manager_set_vdw_type,
        [ctypes.c_void_p, ctypes.c_int],
        "ForceManager.setVdwType(vdw_type)",
    )

    configure_status_function(
        lib().apo_force_manager_set_print_energy_decomposition,
        [ctypes.c_void_p, ctypes.c_bool],
        "ForceManager.setPrintEnergyDecomposition(flag)",
    )

    configure_status_function(
        lib().apo_force_manager_get_num_atoms,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "ForceManager.getNumAtoms()",
    )

    configure_status_function(
        lib().apo_force_manager_is_initialized,
        [ctypes.POINTER(ctypes.c_bool), ctypes.c_void_p],
        "ForceManager.isInitialized()",
    )

    configure_status_function(
        lib().apo_force_manager_get_box_dimensions,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_void_p],
        "ForceManager.getBoxDimensions()",
    )

    configure_status_function(
        lib().apo_force_manager_get_kappa,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "ForceManager.getKappa()",
    )

    configure_status_function(
        lib().apo_force_manager_get_cutoff,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "ForceManager.getCutoff",
    )

    configure_status_function(
        lib().apo_force_manager_get_ctonnb,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "ForceManager.getCtonnb",
    )

    configure_status_function(
        lib().apo_force_manager_get_ctofnb,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "ForceManager.getCtofnb",
    )

    configure_status_function(
        lib().apo_force_manager_get_fft_grid,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.c_void_p],
        "ForceManager.getFFTGrid()",
    )

    configure_status_function(
        lib().apo_force_manager_get_pme_spline_order,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "ForceManager.getPmeSplineOrder()",
    )

    configure_status_function(
        lib().apo_force_manager_get_periodic_boundary_condition,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "ForceManager.getPeriodicBoundaryCondition()",
    )

    configure_status_function(
        lib().apo_force_manager_get_vdw_type,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "ForceManager.getVdwType()",
    )

    configure_status_function(
        lib().apo_force_manager_is_composite,
        [ctypes.POINTER(ctypes.c_bool), ctypes.c_void_p],
        "ForceManager.isComposite()",
    )

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

        lib().apo_force_manager_create(
            ctypes.byref(handle), psf.handle, parameters.handle
        )

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

        is_native_force: bool = isinstance(force, _ApoObject)

        if force_tag == "" and not is_native_force:
            raise ValueError("force_tag must not be empty")

        subscribe_method = getattr(force, "_subscribe_to_force_manager", None)
        if not callable(subscribe_method):
            raise TypeError(
                "ForceManager.subscribe expects an object with _subscribe_to_force_manager(force_manager, force_tag)"
            )

        if not is_native_force and any(
            existing_force is force for existing_force in self._subscribed_forces
        ):
            raise ValueError("force is already subscribed to this ForceManager")

        subscribable_force: _SubscribableForce = cast(_SubscribableForce, force)
        subscribable_force._subscribe_to_force_manager(self, force_tag)

        self._subscribed_forces.append(subscribable_force)

        return

    def unsubscribe(self, force: _SubscribableForce) -> None:
        _initialize_prototypes()

        is_native_force: bool = isinstance(force, _ApoObject)

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

        if index_to_remove is None and not is_native_force:
            raise ValueError("force is not subscribed to this ForceManager")

        subscribable_force: _SubscribableForce = cast(_SubscribableForce, force)
        subscribable_force._unsubscribe_from_force_manager(self)

        if index_to_remove is not None:
            del self._subscribed_forces[index_to_remove]

        return

    def setBoxDimensions(self, box_dimensions: Sequence[float]) -> None:
        _initialize_prototypes()

        box_values: list[float] = [float(value) for value in box_dimensions]

        c_buffer_type = ctypes.c_double * len(box_values)
        c_buffer = c_buffer_type(*box_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(box_values))

        lib().apo_force_manager_set_box_dimensions(self.handle, c_buffer, c_buffer_len)

        return

    def setKappa(self, kappa: float) -> None:
        _initialize_prototypes()

        c_kappa: ctypes.c_double = ctypes.c_double(kappa)

        lib().apo_force_manager_set_kappa(self.handle, c_kappa)

        return

    def setCutoff(self, cutoff: float) -> None:
        _initialize_prototypes()

        c_cutoff: ctypes.c_double = ctypes.c_double(cutoff)

        lib().apo_force_manager_set_cutoff(self.handle, c_cutoff)

        return

    def setCtonnb(self, ctonnb: float) -> None:
        _initialize_prototypes()

        c_ctonnb: ctypes.c_double = ctypes.c_double(ctonnb)

        lib().apo_force_manager_set_ctonnb(self.handle, c_ctonnb)

        return

    def setCtofnb(self, ctofnb: float) -> None:
        _initialize_prototypes()

        c_ctofnb: ctypes.c_double = ctypes.c_double(ctofnb)

        lib().apo_force_manager_set_ctofnb(self.handle, c_ctofnb)

        return

    def setFFTGrid(self, grid: Sequence[int]) -> None:
        _initialize_prototypes()

        grid_values: list[int] = [int(value) for value in grid]

        c_buffer_type = ctypes.c_int * len(grid_values)
        c_buffer = c_buffer_type(*grid_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(grid_values))

        lib().apo_force_manager_set_fft_grid(self.handle, c_buffer, c_buffer_len)

        return

    def setPmeSplineOrder(self, order: int) -> None:
        _initialize_prototypes()

        c_order: ctypes.c_int = ctypes.c_int(order)

        lib().apo_force_manager_set_pme_spline_order(self.handle, c_order)

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

        lib().apo_force_manager_set_periodic_boundary_condition(self.handle, c_pbc)

        return

    def setVdwType(self, vdw_type: VdwType | int) -> None:
        _initialize_prototypes()

        try:
            vdw_type_value: VdwType = VdwType(vdw_type)
        except ValueError as exc:
            raise ValueError(f"invalid vdw_type: {vdw_type!r}") from exc

        c_vdw_type: ctypes.c_int = ctypes.c_int(int(vdw_type_value))

        lib().apo_force_manager_set_vdw_type(self.handle, c_vdw_type)

        return

    def setPrintEnergyDecomposition(self, flag: bool = True) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        lib().apo_force_manager_set_print_energy_decomposition(self.handle, c_flag)

        return

    def getNumAtoms(self) -> int:
        _initialize_prototypes()

        c_num_atoms: ctypes.c_int = ctypes.c_int()

        lib().apo_force_manager_get_num_atoms(ctypes.byref(c_num_atoms), self.handle)

        return int(c_num_atoms.value)

    def isInitialized(self) -> bool:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool()

        lib().apo_force_manager_is_initialized(ctypes.byref(c_flag), self.handle)

        return bool(c_flag.value)

    def getBoxDimensions(self) -> tuple[float, float, float]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_double * 3
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(3)

        lib().apo_force_manager_get_box_dimensions(c_buffer, c_buffer_len, self.handle)

        return (float(c_buffer[0]), float(c_buffer[1]), float(c_buffer[2]))

    def getKappa(self) -> float:
        _initialize_prototypes()

        c_kappa: ctypes.c_double = ctypes.c_double()

        lib().apo_force_manager_get_kappa(ctypes.byref(c_kappa), self.handle)

        return float(c_kappa.value)

    def getCutoff(self) -> float:
        _initialize_prototypes()

        c_cutoff: ctypes.c_double = ctypes.c_double()

        lib().apo_force_manager_get_cutoff(ctypes.byref(c_cutoff), self.handle)

        return float(c_cutoff.value)

    def getCtonnb(self) -> float:
        _initialize_prototypes()

        c_ctonnb: ctypes.c_double = ctypes.c_double()

        lib().apo_force_manager_get_ctonnb(ctypes.byref(c_ctonnb), self.handle)

        return float(c_ctonnb.value)

    def getCtofnb(self) -> float:
        _initialize_prototypes()

        c_ctofnb: ctypes.c_double = ctypes.c_double()

        lib().apo_force_manager_get_ctofnb(ctypes.byref(c_ctofnb), self.handle)

        return float(c_ctofnb.value)

    def getFFTGrid(self) -> tuple[int, int, int]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_int * 3
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(3)

        lib().apo_force_manager_get_fft_grid(c_buffer, c_buffer_len, self.handle)

        return (int(c_buffer[0]), int(c_buffer[1]), int(c_buffer[2]))

    def getPmeSplineOrder(self) -> int:
        _initialize_prototypes()

        c_order: ctypes.c_int = ctypes.c_int()

        lib().apo_force_manager_get_pme_spline_order(ctypes.byref(c_order), self.handle)

        return float(c_order.value)

    def getPeriodicBoundaryCondition(self) -> PeriodicBoundaryCondition:
        _initialize_prototypes()

        c_pbc = ctypes.c_int()

        lib().apo_force_manager_get_periodic_boundary_condition(
            ctypes.byref(c_pbc), self.handle
        )

        return PeriodicBoundaryCondition(c_pbc.value)

    def getVdwType(self) -> VdwType:
        _initialize_prototypes()

        c_vdw_type: ctypes.c_int = ctypes.c_int()

        lib().apo_force_manager_get_vdw_type(ctypes.byref(c_vdw_type), self.handle)

        return VdwType(c_vdw_type.value)

    def isComposite(self) -> bool:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool()

        lib().apo_force_manager_is_composite(ctypes.byref(c_flag), self.handle)

        return bool(c_flag.value)
