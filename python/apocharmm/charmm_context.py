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
from ._lib import encode_path, lib
from ._types import FilePath
from .enums import PeriodicBoundaryCondition, VdwType
from .error import configure_status_function

from .charmm_crd import CharmmCrd
from .charmm_parameters import CharmmParameters
from .charmm_psf import CharmmPsf
from .force_manager import ForceManager


def _flatten_rows(
    values: Sequence[Sequence[float]], width: int, argument_name: str
) -> list[float]:
    flattened: list[float] = []

    for index, row in enumerate(values):
        row_values: list[float] = [float(value) for value in row]
        if len(row_values) != width:
            raise ValueError(
                f"{argument_name}[{index}] must contain exactly {width} values; observed {len(row_values)}"
            )
        flattened.extend(row_values)

    return flattened


_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_charmm_context_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p],
        "CharmmContext construction from ForceManager",
    )

    configure_status_function(
        lib().apo_charmm_context_create_from_psf_parameters,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p, ctypes.c_void_p],
        "CharmmContext construction from CharmmPsf/CharmmParameters",
    )

    lib().apo_charmm_context_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_charmm_context_destroy.restype = None

    configure_status_function(
        lib().apo_charmm_context_set_prm,
        [ctypes.c_void_p, ctypes.c_void_p],
        "CharmmContext.setPrm(parameters)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_psf,
        [ctypes.c_void_p, ctypes.c_void_p],
        "CharmmContext.setPsf(psf)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_force_manager,
        [ctypes.c_void_p, ctypes.c_void_p],
        "CharmmContext.setForceManager(force_manager)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_coordinates_charges,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "CharmmContext.setCoordinatesCharges(coordinates_charges)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_coordinates_from_array,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "CharmmContext.setCoordinates(coordinates)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_coordinates,
        [ctypes.c_void_p, ctypes.c_void_p],
        "CharmmContext.setCoordinates(crd)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_charges,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "CharmmContext.setCharges(charges)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_velocities_inverse_masses,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "CharmmContext.setVelocitiesInverseMasses(velocities_inverse_masses)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_velocities,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "CharmmContext.setVelocities(velocities)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_velocities_from_charmm_velocity_file,
        [ctypes.c_void_p, ctypes.c_char_p],
        "CharmmContext.setVelocitiesFromCHARMMVelocityFile(file_name)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_masses,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "CharmmContext.setMasses(masses)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_temperature,
        [ctypes.c_void_p, ctypes.c_double],
        "CharmmContext.setTemperature(temperature)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_periodic_boundary_condition,
        [ctypes.c_void_p, ctypes.c_int],
        "CharmmContext.setPeriodicBoundaryCondition(pbc)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_box_dimensions,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "CharmmContext.setBoxDimensions(box_dimensions)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_random_seed,
        [ctypes.c_void_p, ctypes.c_uint64],
        "CharmmContext.setRandomSeed(seed)",
    )

    configure_status_function(
        lib().apo_charmm_context_use_holonomic_constraints,
        [ctypes.c_void_p, ctypes.c_bool],
        "CharmmContext.useHolonomicConstraints(flag)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_kappa,
        [ctypes.c_void_p, ctypes.c_double],
        "CharmmContext.setKappa(kappa)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_cutoff,
        [ctypes.c_void_p, ctypes.c_double],
        "CharmmContext.setCutoff(cutoff)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_ctonnb,
        [ctypes.c_void_p, ctypes.c_double],
        "CharmmContext.setCtonnb(ctonnb)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_ctofnb,
        [ctypes.c_void_p, ctypes.c_double],
        "CharmmContext.setCtofnb(ctofnb)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_fft_grid,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_int), ctypes.c_size_t],
        "CharmmContext.setFFTGrid(grid)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_pme_spline_order,
        [ctypes.c_void_p, ctypes.c_int],
        "CharmmContext.setPmeSplineOrder(order)",
    )

    configure_status_function(
        lib().apo_charmm_context_set_vdw_type,
        [ctypes.c_void_p, ctypes.c_int],
        "CharmmContext.setVdwType(vdw_type)",
    )

    configure_status_function(
        lib().apo_charmm_context_get_num_atoms,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "CharmmContext.getNumAtoms()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_num_degrees_of_freedom,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "CharmmContext.getNumDegreesOfFreedom()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_coordinates_charges,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmContext.getCoordinatesCharges()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_velocity_mass,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmContext.getVelocityMass()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_periodic_boundary_condition,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "CharmmContext.getPeriodicBoundaryCondition()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_box_dimensions,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmContext.getBoxDimensions()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_random_seed,
        [ctypes.POINTER(ctypes.c_uint64), ctypes.c_void_p],
        "CharmmContext.getRandomSeed()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_volume,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CharmmContext.getVolume()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_kappa,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CharmmContext.getKappa()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_cutoff,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CharmmContext.getCutoff()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_ctonnb,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CharmmContext.getCtonnb()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_ctofnb,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CharmmContext.getCtofnb()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_fft_grid,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_size_t, ctypes.c_void_p],
        "CharmmContext.getFFTGrid()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_pme_spline_order,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "CharmmContext.getPmeSplineOrder()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_vdw_type,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "CharmmContext.getVdwType()",
    )

    configure_status_function(
        lib().apo_charmm_context_get_force_manager,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p],
        "CharmmContext.getForceManager()",
    )

    configure_status_function(
        lib().apo_charmm_context_assign_velocities_at_temperature,
        [ctypes.c_void_p, ctypes.c_double],
        "CharmmContext.assignVelocitiesAtTemperature(temperature)",
    )

    configure_status_function(
        lib().apo_charmm_context_compute_temperature,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CharmmContext.computeTemperature()",
    )

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
            lib().apo_charmm_context_create(
                ctypes.byref(handle), force_manager_or_psf.handle
            )
            self._force_manager: ForceManager | None = force_manager_or_psf
            self._psf: CharmmPsf | None = None
            self._parameters: CharmmParameters | None = None

        elif isinstance(force_manager_or_psf, CharmmPsf) and isinstance(
            parameters, CharmmParameters
        ):
            lib().apo_charmm_context_create_from_psf_parameters(
                ctypes.byref(handle), force_manager_or_psf.handle, parameters.handle
            )
            self._force_manager = None
            self._psf = force_manager_or_psf
            self._parameters = parameters

        else:
            raise TypeError(
                "CharmmContext expects either ForceManager or (CharmmPsf, CharmmParameters)"
            )

        if handle.value is None:
            raise RuntimeError(
                "CharmmContext construction returned success but produced a NULL handle"
            )

        self._handle = handle

        return

    def setPrm(self, parameters: CharmmParameters) -> None:
        _initialize_prototypes()

        if not isinstance(parameters, CharmmParameters):
            raise TypeError("CharmmContext.setPrm expects CharmmParameters")

        lib().apo_charmm_context_set_prm(self.handle, parameters.handle)
        self._parameters = parameters

        return

    def setPsf(self, psf: CharmmPsf) -> None:
        _initialize_prototypes()

        if not isinstance(psf, CharmmPsf):
            raise TypeError("CharmmContext.setPsf expects CharmmPsf")

        lib().apo_charmm_context_set_psf(self.handle, psf.handle)
        self._psf = psf

        return

    def setForceManager(self, force_manager: ForceManager) -> None:
        _initialize_prototypes()

        if not isinstance(force_manager, ForceManager):
            raise TypeError("CharmmContext.setForceManager expects ForceManager")

        lib().apo_charmm_context_set_force_manager(self.handle, force_manager.handle)
        self._force_manager = force_manager
        self._psf = None
        self._parameters = None

        return

    def setCoordinatesCharges(
        self, coordinates_charges: Sequence[Sequence[float]]
    ) -> None:
        _initialize_prototypes()

        values: list[float] = _flatten_rows(
            coordinates_charges, 4, "coordinates_charges"
        )
        c_buffer_type = ctypes.c_double * len(values)
        c_buffer = c_buffer_type(*values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(values))

        lib().apo_charmm_context_set_coordinates_charges(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setCoordinates(
        self, coordinates: CharmmCrd | Sequence[Sequence[float]]
    ) -> None:
        _initialize_prototypes()

        if isinstance(coordinates, CharmmCrd):
            lib().apo_charmm_context_set_coordinates(self.handle, coordinates.handle)
            return

        values: list[float] = _flatten_rows(coordinates, 3, "coordinates")
        c_buffer_type = ctypes.c_double * len(values)
        c_buffer = c_buffer_type(*values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(values))

        lib().apo_charmm_context_set_coordinates_from_array(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setCharges(self, charges: Sequence[float]) -> None:
        _initialize_prototypes()

        values: list[float] = [float(value) for value in charges]
        c_buffer_type = ctypes.c_double * len(values)
        c_buffer = c_buffer_type(*values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(values))

        lib().apo_charmm_context_set_charges(self.handle, c_buffer, c_buffer_len)

        return

    def setVelocitiesInverseMasses(
        self, velocities_inverse_masses: Sequence[Sequence[float]]
    ) -> None:
        _initialize_prototypes()

        values: list[float] = _flatten_rows(
            velocities_inverse_masses, 4, "velocities_inverse_masses"
        )
        c_buffer_type = ctypes.c_double * len(values)
        c_buffer = c_buffer_type(*values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(values))

        lib().apo_charmm_context_set_velocities_inverse_masses(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setVelocities(self, velocities: Sequence[Sequence[float]]) -> None:
        _initialize_prototypes()

        values: list[float] = _flatten_rows(velocities, 3, "velocities")
        c_buffer_type = ctypes.c_double * len(values)
        c_buffer = c_buffer_type(*values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(values))

        lib().apo_charmm_context_set_velocities(self.handle, c_buffer, c_buffer_len)

        return

    def setVelocitiesFromCHARMMVelocityFile(self, path: FilePath) -> None:
        _initialize_prototypes()

        encoded_path: bytes = encode_path(path)
        c_path: ctypes.c_char_p = ctypes.c_char_p(encoded_path)

        lib().apo_charmm_context_set_velocities_from_charmm_velocity_file(
            self.handle, c_path
        )

        return

    def setMasses(self, masses: Sequence[float]) -> None:
        _initialize_prototypes()

        values: list[float] = [float(value) for value in masses]
        c_buffer_type = ctypes.c_double * len(values)
        c_buffer = c_buffer_type(*values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(values))

        lib().apo_charmm_context_set_masses(self.handle, c_buffer, c_buffer_len)

        return

    def setTemperature(self, temperature: float) -> None:
        _initialize_prototypes()

        c_temperature: ctypes.c_double = ctypes.c_double(temperature)

        lib().apo_charmm_context_set_temperature(self.handle, c_temperature)

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

        lib().apo_charmm_context_set_periodic_boundary_condition(self.handle, c_pbc)

        return

    def setBoxDimensions(self, box_dimensions: Sequence[float]) -> None:
        _initialize_prototypes()

        values: list[float] = [float(value) for value in box_dimensions]

        c_buffer_type = ctypes.c_double * len(values)
        c_buffer = c_buffer_type(*values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(values))

        lib().apo_charmm_context_set_box_dimensions(self.handle, c_buffer, c_buffer_len)

        return

    def setRandomSeed(self, seed: int) -> None:
        _initialize_prototypes()

        if seed < 0 or seed > 2**64 - 1:
            raise ValueError("seed must fit in uint64_t")

        c_seed: ctypes.c_uint64 = ctypes.c_uint64(seed)

        lib().apo_charmm_context_set_random_seed(self.handle, c_seed)

        return

    def useHolonomicConstraints(self, flag: bool) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        lib().apo_charmm_context_use_holonomic_constraints(self.handle, c_flag)

        return

    def setKappa(self, kappa: float) -> None:
        _initialize_prototypes()

        c_kappa: ctypes.c_double = ctypes.c_double(kappa)

        lib().apo_charmm_context_set_kappa(self.handle, c_kappa)

        return

    def setCutoff(self, cutoff: float) -> None:
        _initialize_prototypes()

        c_cutoff: ctypes.c_double = ctypes.c_double(cutoff)

        lib().apo_charmm_context_set_cutoff(self.handle, c_cutoff)

        return

    def setCtonnb(self, ctonnb: float) -> None:
        _initialize_prototypes()

        c_ctonnb: ctypes.c_double = ctypes.c_double(ctonnb)

        lib().apo_charmm_context_set_ctonnb(self.handle, c_ctonnb)

        return

    def setCtofnb(self, ctofnb: float) -> None:
        _initialize_prototypes()

        c_ctofnb: ctypes.c_double = ctypes.c_double(ctofnb)

        lib().apo_charmm_context_set_ctofnb(self.handle, c_ctofnb)

        return

    def setFFTGrid(self, grid: Sequence[int]) -> None:
        _initialize_prototypes()

        values: list[int] = [int(value) for value in grid]
        c_buffer_type = ctypes.c_int * len(values)
        c_buffer = c_buffer_type(*values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(values))

        lib().apo_charmm_context_set_fft_grid(self.handle, c_buffer, c_buffer_len)

        return

    def setPmeSplineOrder(self, order: int) -> None:
        _initialize_prototypes()

        c_order: ctypes.c_int = ctypes.c_int(order)

        lib().apo_charmm_context_set_pme_spline_order(self.handle, c_order)

        return

    def setVdwType(self, vdw_type: VdwType | int) -> None:
        _initialize_prototypes()

        try:
            vdw_type_value: VdwType = VdwType(vdw_type)
        except ValueError as exc:
            raise ValueError(f"invalid vdw_type: {vdw_type!r}") from exc

        c_vdw_type: ctypes.c_int = ctypes.c_int(int(vdw_type_value))

        lib().apo_charmm_context_set_vdw_type(self.handle, c_vdw_type)

        return

    def getNumAtoms(self) -> int:
        _initialize_prototypes()

        num_atoms: ctypes.c_int = ctypes.c_int()

        lib().apo_charmm_context_get_num_atoms(ctypes.byref(num_atoms), self.handle)

        return int(num_atoms.value)

    def getNumDegreesOfFreedom(self) -> int:
        _initialize_prototypes()

        ndegf: ctypes.c_int = ctypes.c_int()

        lib().apo_charmm_context_get_num_degrees_of_freedom(
            ctypes.byref(ndegf), self.handle
        )

        return int(ndegf.value)

    def getCoordinatesCharges(self) -> list[list[float]]:
        _initialize_prototypes()

        num_atoms: int = self.getNumAtoms()

        c_buffer_type = ctypes.c_double * (num_atoms * 4)
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(num_atoms * 4)

        lib().apo_charmm_context_get_coordinates_charges(
            c_buffer, c_buffer_len, self.handle
        )

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

        lib().apo_charmm_context_get_velocity_mass(c_buffer, c_buffer_len, self.handle)

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

    def getPeriodicBoundaryCondition(self) -> PeriodicBoundaryCondition:
        _initialize_prototypes()

        c_pbc = ctypes.c_int()

        lib().apo_charmm_context_get_periodic_boundary_condition(
            ctypes.byref(c_pbc), self.handle
        )

        return PeriodicBoundaryCondition(c_pbc.value)

    def getBoxDimensions(self) -> tuple[float, float, float]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_double * 3
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(3)

        lib().apo_charmm_context_get_box_dimensions(c_buffer, c_buffer_len, self.handle)

        return (float(c_buffer[0]), float(c_buffer[1]), float(c_buffer[2]))

    def getRandomSeed(self) -> int:
        _initialize_prototypes()

        c_seed: ctypes.c_uint64 = ctypes.c_uint64()

        lib().apo_charmm_context_get_random_seed(ctypes.byref(c_seed), self.handle)

        return int(c_seed.value)

    def getVolume(self) -> float:
        _initialize_prototypes()

        c_volume: ctypes.c_double = ctypes.c_double()

        lib().apo_charmm_context_get_volume(ctypes.byref(c_volume), self.handle)

        return float(c_volume.value)

    def getKappa(self) -> float:
        _initialize_prototypes()

        c_kappa: ctypes.c_double = ctypes.c_double()

        lib().apo_charmm_context_get_kappa(ctypes.byref(c_kappa), self.handle)

        return float(c_kappa.value)

    def getCutoff(self) -> float:
        _initialize_prototypes()

        c_cutoff: ctypes.c_double = ctypes.c_double()

        lib().apo_charmm_context_get_cutoff(ctypes.byref(c_cutoff), self.handle)

        return float(c_cutoff.value)

    def getCtonnb(self) -> float:
        _initialize_prototypes()

        c_ctonnb: ctypes.c_double = ctypes.c_double()

        lib().apo_charmm_context_get_ctonnb(ctypes.byref(c_ctonnb), self.handle)

        return float(c_ctonnb.value)

    def getCtofnb(self) -> float:
        _initialize_prototypes()

        c_ctofnb: ctypes.c_double = ctypes.c_double()

        lib().apo_charmm_context_get_ctofnb(ctypes.byref(c_ctofnb), self.handle)

        return float(c_ctofnb.value)

    def getFFTGrid(self) -> tuple[int, int, int]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_int * 3
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(3)

        lib().apo_charmm_context_get_fft_grid(c_buffer, c_buffer_len, self.handle)

        return (int(c_buffer[0]), int(c_buffer[1]), int(c_buffer[2]))

    def getPmeSplineOrder(self) -> int:
        _initialize_prototypes()

        c_order: ctypes.c_int = ctypes.c_int()

        lib().apo_charmm_context_get_pme_spline_order(
            ctypes.byref(c_order), self.handle
        )

        return int(c_order.value)

    def getVdwType(self) -> VdwType:
        _initialize_prototypes()

        c_vdw_type: ctypes.c_int = ctypes.c_int()

        lib().apo_charmm_context_get_vdw_type(ctypes.byref(c_vdw_type), self.handle)

        return VdwType(c_vdw_type.value)

    def getForceManager(self) -> ForceManager:
        _initialize_prototypes()

        if self._force_manager is not None:
            return self._force_manager

        fm_handle: ctypes.c_void_p = ctypes.c_void_p()

        lib().apo_charmm_context_get_force_manager(ctypes.byref(fm_handle), self.handle)

        self._force_manager = ForceManager._from_handle(fm_handle)

        return self._force_manager

    def assignVelocitiesAtTemperature(self, temperature: float) -> None:
        _initialize_prototypes()

        c_temperature: ctypes.c_double = ctypes.c_double(temperature)

        lib().apo_charmm_context_assign_velocities_at_temperature(
            self.handle, c_temperature
        )

        return

    def computeTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        lib().apo_charmm_context_compute_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        return float(c_temperature.value)
