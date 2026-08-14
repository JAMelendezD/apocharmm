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

from ._lib import lib
from .enums import CrystalType
from .error import configure_status_function
from .cuda_integrator import CudaIntegrator

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_double],
        "CudaLangevinPistonIntegrator construction",
    )

    lib().apo_cuda_langevin_piston_integrator_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_cuda_langevin_piston_integrator_destroy.restype = None

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_use_nose_hoover_thermostat,
        [ctypes.c_void_p, ctypes.c_bool],
        "CudaLangevinPistonIntegrator.useNoseHooverThermostat(flag)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_set_reference_temperature,
        [ctypes.c_void_p, ctypes.c_double],
        "CudaLangevinPistonIntegrator.setReferenceTemperature(temperature)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_set_nose_hoover_piston_mass,
        [ctypes.c_void_p, ctypes.c_double],
        "CudaLangevinPistonIntegrator.setNoseHooverPistonMass(mass)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_use_old_temperature,
        [ctypes.c_void_p, ctypes.c_bool],
        "CudaLangevinPistonIntegrator.useOldTemperature(flag)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_set_reference_pressure,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "CudaLangevinPistonIntegrator.setReferencePressure(pressure_tensor)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_set_constant_surface_tension,
        [ctypes.c_void_p, ctypes.c_bool],
        "CudaLangevinPistonIntegrator.setConstantSurfaceTension(flag)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_set_crystal_type,
        [ctypes.c_void_p, ctypes.c_int],
        "CudaLangevinPistonIntegrator.setCrystalType(crystal_type)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_mass,
        [ctypes.c_void_p, ctypes.POINTER(ctypes.c_double), ctypes.c_size_t],
        "CudaLangevinPistonIntegrator.setLangevinPistonMass(mass)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_friction_seed,
        [ctypes.c_void_p, ctypes.c_uint64],
        "CudaLangevinPistonIntegrator.setLangevinPistonFrictionSeed(seed)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_friction,
        [ctypes.c_void_p, ctypes.c_double],
        "CudaLangevinPistonIntegrator.setLangevinPistonFriction(friction)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_reset_averages,
        [ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.resetAverages()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_reference_temperature,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getReferenceTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_nose_hoover_piston_mass,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getNoseHooverPistonMass()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_average_temperature,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getAverageTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_reference_pressure_tensor,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getReferencePressureTensor()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_crystal_type,
        [ctypes.POINTER(ctypes.c_int), ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getCrystalType()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_mass,
        [
            ctypes.POINTER(ctypes.c_size_t),
            ctypes.POINTER(ctypes.c_double),
            ctypes.c_size_t,
            ctypes.c_void_p,
        ],
        "CudaLangevinPistonIntegrator.getLangevinPistonMass()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_tensor,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getInstantaneousPressureTensor()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_scalar,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getInstantaneousPressureScalar()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_average_pressure_tensor,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_size_t, ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getAveragePressureTensor()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_average_pressure_scalar,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getAveragePressureScalar()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_friction_seed,
        [ctypes.POINTER(ctypes.c_uint64), ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getLangevinPistonFrictionSeed()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_get_instantaneous_temperature,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaLangevinPistonIntegrator.getInstantaneousTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_piston_integrator_as_cuda_integrator,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p],
        "CudaLangevinPistonIntegrator base-integrator conversion",
    )

    _prototypes_initialized = True

    return


class CudaLangevinPistonIntegrator(CudaIntegrator):
    _destroy_function_name = "apo_cuda_langevin_piston_integrator_destroy"

    def __init__(self, time_step: float) -> None:
        _initialize_prototypes()
        super().__init__()

        handle: ctypes.c_void_p = ctypes.c_void_p()
        c_time_step: ctypes.c_double = ctypes.c_double(time_step)

        lib().apo_cuda_langevin_piston_integrator_create(
            ctypes.byref(handle), c_time_step
        )

        if handle.value is None:
            raise RuntimeError(
                "apo_cuda_langevin_piston_integrator_create returned success but produced a NULL handle"
            )

        self._handle = handle

        integrator_handle: ctypes.c_void_p = ctypes.c_void_p()

        lib().apo_cuda_langevin_piston_integrator_as_cuda_integrator(
            ctypes.byref(integrator_handle), self.handle
        )

        if integrator_handle.value is None:
            raise RuntimeError(
                "apo_cuda_langevin_piston_integrator_as_cuda_integrator returned success but produced a NULL handle"
            )

        self._integrator_handle = integrator_handle

        return

    def useNoseHooverThermostat(self, flag: bool = True) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        lib().apo_cuda_langevin_piston_integrator_use_nose_hoover_thermostat(
            self.handle, c_flag
        )

        return

    def setReferenceTemperature(self, temperature: float) -> None:
        _initialize_prototypes()

        c_temperature: ctypes.c_double = ctypes.c_double(temperature)

        lib().apo_cuda_langevin_piston_integrator_set_reference_temperature(
            self.handle, c_temperature
        )

        return

    def setNoseHooverPistonMass(self, mass: float) -> None:
        _initialize_prototypes()

        c_mass: ctypes.c_double = ctypes.c_double(mass)

        lib().apo_cuda_langevin_piston_integrator_set_nose_hoover_piston_mass(
            self.handle, c_mass
        )

        return

    def useOldTemperature(self, flag: bool = True) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        lib().apo_cuda_langevin_piston_integrator_use_old_temperature(
            self.handle, c_flag
        )

        return

    def setReferencePressure(self, pressure_tensor: Sequence[Sequence[float]]) -> None:
        _initialize_prototypes()

        flattened_pressure_tensor: list[float] = []

        for pressure in pressure_tensor:
            flattened_pressure_tensor.extend(float(value) for value in pressure)

        c_buffer_type = ctypes.c_double * len(flattened_pressure_tensor)
        c_buffer = c_buffer_type(*flattened_pressure_tensor)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(flattened_pressure_tensor))

        lib().apo_cuda_langevin_piston_integrator_set_reference_pressure(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setConstantSurfaceTension(self, flag: bool = True) -> None:
        _initialize_prototypes()

        c_flag: ctypes.c_bool = ctypes.c_bool(flag)

        lib().apo_cuda_langevin_piston_integrator_set_constant_surface_tension(
            self.handle, c_flag
        )

        return

    def setCrystalType(self, crystal_type: CrystalType | int) -> None:
        _initialize_prototypes()

        c_crystal_type: ctypes.c_int = ctypes.c_int(int(crystal_type))

        lib().apo_cuda_langevin_piston_integrator_set_crystal_type(
            self.handle, c_crystal_type
        )

        return

    def setLangevinPistonMass(self, masses: Sequence[float]) -> None:
        _initialize_prototypes()

        mass_values: list[float] = [float(value) for value in masses]

        c_buffer_type = ctypes.c_double * len(mass_values)
        c_buffer = c_buffer_type(*mass_values)
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(len(mass_values))

        lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_mass(
            self.handle, c_buffer, c_buffer_len
        )

        return

    def setLangevinPistonFrictionSeed(self, seed: int) -> None:
        _initialize_prototypes()

        if seed < 0 or seed > 2**64 - 1:
            raise ValueError("seed must fit in uint64_t")

        c_seed: ctypes.c_uint64 = ctypes.c_uint64(seed)

        lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_friction_seed(
            self.handle, c_seed
        )

        return

    def setLangevinPistonFriction(self, friction: float) -> None:
        _initialize_prototypes()

        c_friction: ctypes.c_double = ctypes.c_double(friction)

        lib().apo_cuda_langevin_piston_integrator_set_langevin_piston_friction(
            self.handle, c_friction
        )

        return

    def resetAverages(self) -> None:
        _initialize_prototypes()

        lib().apo_cuda_langevin_piston_integrator_reset_averages(self.handle)

        return

    def getReferenceTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        lib().apo_cuda_langevin_piston_integrator_get_reference_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        return float(c_temperature.value)

    def getNoseHooverPistonMass(self) -> float:
        _initialize_prototypes()

        c_mass = ctypes.c_double()

        lib().apo_cuda_langevin_piston_integrator_get_nose_hoover_piston_mass(
            ctypes.byref(c_mass), self.handle
        )

        return float(c_mass.value)

    def getAverageTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        lib().apo_cuda_langevin_piston_integrator_get_average_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        return float(c_temperature.value)

    def getReferencePressureTensor(self) -> list[list[float]]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_double * 9
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(9)

        lib().apo_cuda_langevin_piston_integrator_get_reference_pressure_tensor(
            c_buffer, c_buffer_len, self.handle
        )

        pressure_tensor: list[list[float]] = []
        for i in range(3):
            pressure_tensor.append(
                [
                    float(c_buffer[i * 3 + 0]),
                    float(c_buffer[i * 3 + 1]),
                    float(c_buffer[i * 3 + 2]),
                ]
            )

        return pressure_tensor

    def getCrystalType(self) -> CrystalType:
        _initialize_prototypes()

        c_crystal_type = ctypes.c_int()

        lib().apo_cuda_langevin_piston_integrator_get_crystal_type(
            ctypes.byref(c_crystal_type), self.handle
        )

        return CrystalType(c_crystal_type.value)

    def getLangevinPistonMass(self) -> list[float]:
        _initialize_prototypes()

        c_num_mass: ctypes.c_size_t = ctypes.c_size_t()

        c_buffer_type = ctypes.c_double * 3
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(3)

        lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_mass(
            ctypes.byref(c_num_mass), c_buffer, c_buffer_len, self.handle
        )

        num_mass: int = int(c_num_mass.value)

        masses: list[float] = []
        for i in range(num_mass):
            masses.append(float(c_buffer[i]))

        return masses

    def getInstantaneousPressureTensor(self) -> list[list[float]]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_double * 9
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(9)

        lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_tensor(
            c_buffer, c_buffer_len, self.handle
        )

        pressure_tensor: list[list[float]] = []
        for i in range(3):
            pressure_tensor.append(
                [
                    float(c_buffer[i * 3 + 0]),
                    float(c_buffer[i * 3 + 1]),
                    float(c_buffer[i * 3 + 2]),
                ]
            )

        return pressure_tensor

    def getInstantaneousPressureScalar(self) -> float:
        _initialize_prototypes()

        c_pressure_scalar = ctypes.c_double()

        lib().apo_cuda_langevin_piston_integrator_get_instantaneous_pressure_scalar(
            ctypes.byref(c_pressure_scalar), self.handle
        )

        return float(c_pressure_scalar.value)

    def getAveragePressureTensor(self) -> list[list[float]]:
        _initialize_prototypes()

        c_buffer_type = ctypes.c_double * 9
        c_buffer = c_buffer_type()
        c_buffer_len: ctypes.c_size_t = ctypes.c_size_t(9)

        lib().apo_cuda_langevin_piston_integrator_get_average_pressure_tensor(
            c_buffer, c_buffer_len, self.handle
        )

        pressure_tensor: list[list[float]] = []
        for i in range(3):
            pressure_tensor.append(
                [
                    float(c_buffer[i * 3 + 0]),
                    float(c_buffer[i * 3 + 1]),
                    float(c_buffer[i * 3 + 2]),
                ]
            )

        return pressure_tensor

    def getAveragePressureScalar(self) -> float:
        _initialize_prototypes()

        c_pressure_scalar = ctypes.c_double()

        lib().apo_cuda_langevin_piston_integrator_get_average_pressure_scalar(
            ctypes.byref(c_pressure_scalar), self.handle
        )

        return float(c_pressure_scalar.value)

    def getLangevinPistonFrictionSeed(self) -> int:
        _initialize_prototypes()

        c_seed = ctypes.c_uint64()

        lib().apo_cuda_langevin_piston_integrator_get_langevin_piston_friction_seed(
            ctypes.byref(c_seed), self.handle
        )

        return int(c_seed.value)

    def getInstantaneousTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        lib().apo_cuda_langevin_piston_integrator_get_instantaneous_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        return float(c_temperature.value)
