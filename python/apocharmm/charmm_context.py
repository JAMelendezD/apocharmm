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

from ._base import _ApoObject
from ._lib import lib
from .enums import PeriodicBoundaryCondition, VdwType
from .error import check_status

from .charmm_crd import CharmmCrd
from .charmm_parameters import CharmmParameters
from .charmm_psf import CharmmPsf
from .force_manager import ForceManager

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    lib().apo_charmm_context_create.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_create.restype = ctypes.c_int

    lib().apo_charmm_context_create_from_psf_parameters.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_create_from_psf_parameters.restype = ctypes.c_int

    lib().apo_charmm_context_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_charmm_context_destroy.restype = None

    lib().apo_charmm_context_set_box_dimensions.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
    ]
    lib().apo_charmm_context_set_box_dimensions.restype = ctypes.c_int

    lib().apo_charmm_context_set_kappa.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib().apo_charmm_context_set_kappa.restype = ctypes.c_int

    lib().apo_charmm_context_set_cutoff.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib().apo_charmm_context_set_cutoff.restype = ctypes.c_int

    lib().apo_charmm_context_set_ctonnb.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib().apo_charmm_context_set_ctonnb.restype = ctypes.c_int

    lib().apo_charmm_context_set_ctofnb.argtypes = [ctypes.c_void_p, ctypes.c_double]
    lib().apo_charmm_context_set_ctofnb.restype = ctypes.c_int

    lib().apo_charmm_context_set_fft_grid.argtypes = [
        ctypes.c_void_p,
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
    ]
    lib().apo_charmm_context_set_fft_grid.restype = ctypes.c_int

    lib().apo_charmm_context_set_pme_spline_order.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
    ]
    lib().apo_charmm_context_set_pme_spline_order.restype = ctypes.c_int

    lib().apo_charmm_context_set_periodic_boundary_condition.argtypes = [
        ctypes.c_void_p,
        ctypes.c_int,
    ]
    lib().apo_charmm_context_set_periodic_boundary_condition.restype = ctypes.c_int

    lib().apo_charmm_context_set_vdw_type.argtypes = [ctypes.c_void_p, ctypes.c_int]
    lib().apo_charmm_context_set_vdw_type.restype = ctypes.c_int

    lib().apo_charmm_context_set_coordinates.argtypes = [
        ctypes.c_void_p,
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_set_coordinates.restype = ctypes.c_int

    lib().apo_charmm_context_set_random_seed.argtypes = [
        ctypes.c_void_p,
        ctypes.c_uint64,
    ]
    lib().apo_charmm_context_set_random_seed.restype = ctypes.c_int

    lib().apo_charmm_context_use_holonomic_constraints.argtypes = [
        ctypes.c_void_p,
        ctypes.c_bool,
    ]
    lib().apo_charmm_context_use_holonomic_constraints.restype = ctypes.c_int

    lib().apo_charmm_context_get_num_atoms.argtypes = [
        ctypes.POINTER(ctypes.c_size_t),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_num_atoms.restype = ctypes.c_int

    lib().apo_charmm_context_get_coordinates_charges.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_coordinates_charges.restype = ctypes.c_int

    lib().apo_charmm_context_get_velocity_mass.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_velocity_mass.restype = ctypes.c_int

    lib().apo_charmm_context_get_box_dimensions.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_box_dimensions.restype = ctypes.c_int

    lib().apo_charmm_context_get_kappa.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_kappa.restype = ctypes.c_int

    lib().apo_charmm_context_get_cutoff.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_cutoff.restype = ctypes.c_int

    lib().apo_charmm_context_get_ctonnb.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_ctonnb.restype = ctypes.c_int

    lib().apo_charmm_context_get_ctofnb.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_ctofnb.restype = ctypes.c_int

    lib().apo_charmm_context_get_fft_grid.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_size_t,
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_fft_grid.restype = ctypes.c_int

    lib().apo_charmm_context_get_pme_spline_order.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_pme_spline_order.restype = ctypes.c_int

    lib().apo_charmm_context_get_periodic_boundary_condition.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_periodic_boundary_condition.restype = ctypes.c_int

    lib().apo_charmm_context_get_vdw_type.argtypes = [
        ctypes.POINTER(ctypes.c_int),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_vdw_type.restype = ctypes.c_int

    lib().apo_charmm_context_get_force_manager.argtypes = [
        ctypes.POINTER(ctypes.c_void_p),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_get_force_manager.restype = ctypes.c_int

    lib().apo_charmm_context_assign_velocities_at_temperature.argtypes = [
        ctypes.c_void_p,
        ctypes.c_double,
    ]
    lib().apo_charmm_context_assign_velocities_at_temperature.restype = ctypes.c_int

    lib().apo_charmm_context_compute_temperature.argtypes = [
        ctypes.POINTER(ctypes.c_double),
        ctypes.c_void_p,
    ]
    lib().apo_charmm_context_compute_temperature.restype = ctypes.c_int

    _prototypes_initialized = True

    return


class CharmmContext(_ApoObject):
    _destroy_function_name = "apo_charmm_context_destroy"

    def __init__(
        self,
        force_manager_or_psf: ForceManager | CharmmPsf,
        parameters: CharmmParameters | None = None,
    ) -> None:
        _initialize_prototypes()
        super().__init__()

        handle: ctypes.c_void_p = ctypes.c_void_p()

        if isinstance(force_manager_or_psf, ForceManager) and parameters is None:
            status = lib().apo_charmm_context_create(
                ctypes.byref(handle), force_manager_or_psf.handle
            )
            error_context = "CharmmContext construction from ForceManager failed"
            self._force_manager: ForceManager | None = force_manager_or_psf
            self._psf: CharmmPsf | None = None
            self._parameters: CharmmParameters | None = None

        elif isinstance(force_manager_or_psf, CharmmPsf) and isinstance(
            parameters, CharmmParameters
        ):
            status = lib().apo_charmm_context_create_from_psf_parameters(
                ctypes.byref(handle), force_manager_or_psf.handle, parameters.handle
            )
            error_context = (
                "CharmmContext construction from CharmmPsf/CharmmParameters failed"
            )
            self._force_manager = None
            self._psf = force_manager_or_psf
            self._parameters = parameters

        else:
            raise TypeError(
                "CharmmContext expects either ForceManager or (CharmmPsf, CharmmParameters)"
            )

        check_status(status, error_context)

        if handle.value is None:
            raise RuntimeError(
                "CharmmContext construction returned success but produced a NULL handle"
            )

        self._handle = handle

        return

    def setBoxDimensions(self, box_dimensions: Sequence[float]) -> None:
        _initialize_prototypes()

        box_values: list[float] = [float(value) for value in box_dimensions]

        c_buffer_type = ctypes.c_double * len(box_values)
        c_buffer = c_buffer_type(*box_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(box_values))

        status = lib().apo_charmm_context_set_box_dimensions(
            self.handle, c_buffer, c_buffer_len
        )

        check_status(status, "CharmmContext.setBoxDimensions(box_dimensions) failed")

        return

    def setKappa(self, kappa: float) -> None:
        _initialize_prototypes()

        c_kappa: ctypes.c_double = ctypes.c_double(kappa)

        status = lib().apo_charmm_context_set_kappa(self.handle, c_kappa)

        check_status(status, "CharmmContext.setKappa(kappa) failed")

        return

    def setCutoff(self, cutoff: float) -> None:
        _initialize_prototypes()

        c_cutoff: ctypes.c_double = ctypes.c_double(cutoff)

        status = lib().apo_charmm_context_set_cutoff(self.handle, c_cutoff)

        check_status(status, "CharmmContext.setCutoff(cutoff) failed")

        return

    def setCtonnb(self, ctonnb: float) -> None:
        _initialize_prototypes()

        c_ctonnb: ctypes.c_double = ctypes.c_double(ctonnb)

        status = lib().apo_charmm_context_set_ctonnb(self.handle, c_ctonnb)

        check_status(status, "CharmmContext.setCtonnb(ctonnb) failed")

        return

    def setCtofnb(self, ctofnb: float) -> None:
        _initialize_prototypes()

        c_ctofnb: ctypes.c_double = ctypes.c_double(ctofnb)

        status = lib().apo_charmm_context_set_ctofnb(self.handle, c_ctofnb)

        check_status(status, "CharmmContext.setCtofnb(ctofnb) failed")

        return

    def setFFTGrid(self, grid: Sequence[int]) -> None:
        _initialize_prototypes()

        grid_values: list[int] = [int(value) for value in grid]

        c_buffer_type = ctypes.c_int * len(grid_values)
        c_buffer = c_buffer_type(*grid_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(grid_values))

        status = lib().apo_charmm_context_set_fft_grid(
            self.handle, c_buffer, c_buffer_len
        )

        check_status(status, "CharmmContext.setFFTGrid(grid) failed")

        return

    def setPmeSplineOrder(self, order: int) -> None:
        _initialize_prototypes()

        c_order: ctypes.c_int = ctypes.c_int(order)

        status = lib().apo_charmm_context_set_pme_spline_order(self.handle, c_order)

        check_status(status, "CharmmContext.setPmeSplineOrder(order) failed")

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

        status = lib().apo_charmm_context_set_periodic_boundary_condition(
            self.handle, c_pbc
        )

        check_status(status, "CharmmContext.setPeriodicBoundaryCondition(pbc) failed")

        return

    def setVdwType(self, vdw_type: VdwType | int) -> None:
        _initialize_prototypes()

        try:
            vdw_type_value: VdwType = VdwType(vdw_type)
        except ValueError as exc:
            raise ValueError(f"invalid vdw_type: {vdw_type!r}") from exc

        c_vdw_type: ctypes.c_int = ctypes.c_int(int(vdw_type_value))

        status = lib().apo_charmm_context_set_vdw_type(self.handle, c_vdw_type)

        check_status(status, "CharmmContext.setVdwType(vdw_type) failed")

        return

    def setCoordinates(self, crd: CharmmCrd) -> None:
        _initialize_prototypes()

        if not isinstance(crd, CharmmCrd):
            raise TypeError("CharmmContext.setCoordinates expects a CharmmCrd")

        status = lib().apo_charmm_context_set_coordinates(self.handle, crd.handle)

        check_status(status, "CharmmContext.setCoordinates(crd) failed")

        return

    def setRandomSeed(self, seed: int) -> None:
        _initialize_prototypes()

        if seed < 0 or seed > 2**64 - 1:
            raise ValueError("seed must fit in uint64_t")

        c_seed: ctypes.c_uint64 = ctypes.c_uint64(seed)

        status = lib().apo_charmm_context_set_random_seed(self.handle, c_seed)

        check_status(status, "CharmmContext.setRandomSeed(seed) failed")

        return

    def useHolonomicConstraints(self, flag: bool) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        status = lib().apo_charmm_context_use_holonomic_constraints(self.handle, c_flag)

        check_status(status, "CharmmContext.useHolonomicConstraints(flag) failed")

        return

    def getNumAtoms(self) -> int:
        _initialize_prototypes()

        num_atoms: ctypes.c_size_t = ctypes.c_size_t()

        status = lib().apo_charmm_context_get_num_atoms(
            ctypes.byref(num_atoms), self.handle
        )

        check_status(status, "CharmmContext.getNumAtoms() failed")

        return int(num_atoms.value)

    def getCoordinatesCharges(self) -> list[list[float]]:
        _initialize_prototypes()

        num_atoms: int = self.getNumAtoms()

        c_buffer_type = ctypes.c_double * (num_atoms * 4)
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_atoms * 4)

        status = lib().apo_charmm_context_get_coordinates_charges(
            c_buffer, c_buffer_len, self.handle
        )

        check_status(status, "CharmmContext.getCoordinatesCharges() failed")

        xyzq: list[list[float]] = []
        for i in range(num_atoms):
            xyzq.append(
                [
                    float(c_buffer[i * 4 + 0]),
                    float(c_buffer[i * 4 + 1]),
                    float(c_buffer[i * 4 + 2]),
                    float(c_buffer[i * 4 + 3]),
                ]
            )

        return xyzq

    def getVelocityMass(self) -> list[list[float]]:
        _initialize_prototypes()

        num_atoms: int = self.getNumAtoms()

        c_buffer_type = ctypes.c_double * (num_atoms * 4)
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_atoms * 4)

        status = lib().apo_charmm_context_get_velocity_mass(
            c_buffer, c_buffer_len, self.handle
        )

        check_status(status, "CharmmContext.getVelocityMass() failed")

        xyzm: list[list[float]] = []
        for i in range(num_atoms):
            xyzm.append(
                [
                    float(c_buffer[i * 4 + 0]),
                    float(c_buffer[i * 4 + 1]),
                    float(c_buffer[i * 4 + 2]),
                    float(c_buffer[i * 4 + 3]),
                ]
            )

        return xyzm

    def getBoxDimensions(self) -> tuple[float, float, float]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_double * 3
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(3)

        status = lib().apo_charmm_context_get_box_dimensions(
            c_buffer, c_buffer_len, self.handle
        )

        check_status(status, "CharmmContext.getBoxDimensions() failed")

        return (float(c_buffer[0]), float(c_buffer[1]), float(c_buffer[2]))

    def getKappa(self) -> float:
        _initialize_prototypes()

        c_kappa: ctypes.c_double = ctypes.c_double()

        status = lib().apo_charmm_context_get_kappa(ctypes.byref(c_kappa), self.handle)

        check_status(status, "CharmmContext.getKappa() failed")

        return float(c_kappa.value)

    def getCutoff(self) -> float:
        _initialize_prototypes()

        c_cutoff: ctypes.c_double = ctypes.c_double()

        status = lib().apo_charmm_context_get_cutoff(
            ctypes.byref(c_cutoff), self.handle
        )

        check_status(status, "CharmmContext.getCutoff() failed")

        return float(c_cutoff.value)

    def getCtonnb(self) -> float:
        _initialize_prototypes()

        c_ctonnb: ctypes.c_double = ctypes.c_double()

        status = lib().apo_charmm_context_get_ctonnb(
            ctypes.byref(c_ctonnb), self.handle
        )

        check_status(status, "CharmmContext.getCtonnb() failed")

        return float(c_ctonnb.value)

    def getCtofnb(self) -> float:
        _initialize_prototypes()

        c_ctofnb: ctypes.c_double = ctypes.c_double()

        status = lib().apo_charmm_context_get_ctofnb(
            ctypes.byref(c_ctofnb), self.handle
        )

        check_status(status, "CharmmContext.getCtofnb() failed")

        return float(c_ctofnb.value)

    def getFFTGrid(self) -> tuple[int, int, int]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_int * 3
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(3)

        status = lib().apo_charmm_context_get_fft_grid(
            c_buffer, c_buffer_len, self.handle
        )

        check_status(status, "CharmmContext.getFFTGrid() failed")

        return (int(c_buffer[0]), int(c_buffer[1]), int(c_buffer[2]))

    def getPmeSplineOrder(self) -> int:
        _initialize_prototypes()

        c_order: ctypes.c_int = ctypes.c_int()

        status = lib().apo_charmm_context_get_pme_spline_order(
            ctypes.byref(c_order), self.handle
        )

        check_status(status, "CharmmContext.getPmeSplineOrder() failed")

        return int(c_order.value)

    def getPeriodicBoundaryCondition(self) -> PeriodicBoundaryCondition:
        _initialize_prototypes()

        c_pbc = ctypes.c_int()

        status = lib().apo_charmm_context_get_periodic_boundary_condition(
            ctypes.byref(c_pbc), self.handle
        )

        check_status(status, "CharmmContext.getPeriodicBoundaryCondition() failed")

        return PeriodicBoundaryCondition(c_pbc.value)

    def getVdwType(self) -> VdwType:
        _initialize_prototypes()

        c_vdw_type: ctypes.c_int = ctypes.c_int()

        status = lib().apo_charmm_context_get_vdw_type(
            ctypes.byref(c_vdw_type), self.handle
        )

        check_status(status, "CharmmContext.getVdwType() failed")

        return VdwType(c_vdw_type.value)

    def getForceManager(self) -> ForceManager:
        _initialize_prototypes()

        if self._force_manager is not None:
            return self._force_manager

        fm_handle: ctypes.c_void_p = ctypes.c_void_p()

        status = lib().apo_charmm_context_get_force_manager(
            ctypes.byref(fm_handle), self.handle
        )

        check_status(status, "CharmmContext.getForceManager() failed")

        self._force_manager = ForceManager._from_handle(fm_handle)

        return self._force_manager

    def assignVelocitiesAtTemperature(self, temperature: float) -> None:
        _initialize_prototypes()

        c_temperature: ctypes.c_double = ctypes.c_double(temperature)

        status = lib().apo_charmm_context_assign_velocities_at_temperature(
            self.handle, c_temperature
        )

        check_status(
            status, "CharmmContext.assignVelocitiesAtTemperature(temperature) failed"
        )

        return

    def computeTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        status = lib().apo_charmm_context_compute_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        check_status(status, "CharmmContext.computeTemperature() failed")

        return float(c_temperature.value)
