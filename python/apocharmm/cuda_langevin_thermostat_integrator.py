# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

import ctypes

from ._lib import lib
from .error import configure_status_function
from .cuda_integrator import CudaIntegrator

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_create,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_double],
        "CudaLangevinThermostatIntegrator construction",
    )

    lib().apo_cuda_langevin_thermostat_integrator_destroy.argtypes = [ctypes.c_void_p]
    lib().apo_cuda_langevin_thermostat_integrator_destroy.restype = None

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_set_reference_temperature,
        [ctypes.c_void_p, ctypes.c_double],
        "CudaLangevinThermostatIntegrator.setReferenceTemperature(temperature)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_set_thermostat_friction,
        [ctypes.c_void_p, ctypes.c_double],
        "CudaLangevinThermostatIntegrator.setThermostatFriction(friction)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_set_thermostat_rng_seed,
        [ctypes.c_void_p, ctypes.c_uint64],
        "CudaLangevinThermostatIntegrator.setThermostatRngSeed(seed)",
    )

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_reset_average_temperature,
        [ctypes.c_void_p],
        "CudaLangevinThermostatIntegrator.resetAverageTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_get_reference_temperature,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaLangevinThermostatIntegrator.getReferenceTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_get_thermostat_friction,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaLangevinThermostatIntegrator.getThermostatFriction()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_get_thermostat_rng_seed,
        [ctypes.POINTER(ctypes.c_uint64), ctypes.c_void_p],
        "CudaLangevinThermostatIntegrator.getThermostatRngSeed()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_get_average_temperature,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaLangevinThermostatIntegrator.getAverageTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_get_instantaneous_temperature,
        [ctypes.POINTER(ctypes.c_double), ctypes.c_void_p],
        "CudaLangevinThermostatIntegrator.getInstantaneousTemperature()",
    )

    configure_status_function(
        lib().apo_cuda_langevin_thermostat_integrator_as_cuda_integrator,
        [ctypes.POINTER(ctypes.c_void_p), ctypes.c_void_p],
        "CudaLangevinThermostatIntegrator base-integrator conversion",
    )

    _prototypes_initialized = True

    return


class CudaLangevinThermostatIntegrator(CudaIntegrator):
    _destroy_function_name = "apo_cuda_langevin_thermostat_integrator_destroy"

    def __init__(self, time_step: float) -> None:
        _initialize_prototypes()
        super().__init__()

        handle: ctypes.c_void_p = ctypes.c_void_p()
        c_time_step: ctypes.c_double = ctypes.c_double(time_step)

        lib().apo_cuda_langevin_thermostat_integrator_create(
            ctypes.byref(handle), c_time_step
        )

        if handle.value is None:
            raise RuntimeError(
                "apo_cuda_langevin_thermostat_integrator_create returned success but produced a NULL handle"
            )

        self._handle = handle

        integrator_handle: ctypes.c_void_p = ctypes.c_void_p()

        lib().apo_cuda_langevin_thermostat_integrator_as_cuda_integrator(
            ctypes.byref(integrator_handle), self.handle
        )

        if integrator_handle.value is None:
            raise RuntimeError(
                "apo_cuda_langevin_thermostat_integrator_as_cuda_integrator returned success but produced a NULL handle"
            )

        self._integrator_handle = integrator_handle

        return

    def setReferenceTemperature(self, temperature: float) -> None:
        _initialize_prototypes()

        c_temperature: ctypes.c_double = ctypes.c_double(temperature)

        lib().apo_cuda_langevin_thermostat_integrator_set_reference_temperature(
            self.handle, c_temperature
        )

        return

    def setThermostatFriction(self, friction: float) -> None:
        _initialize_prototypes()

        c_friction: ctypes.c_double = ctypes.c_double(friction)

        lib().apo_cuda_langevin_thermostat_integrator_set_thermostat_friction(
            self.handle, c_friction
        )

        return

    def setThermostatRngSeed(self, seed: int) -> None:
        _initialize_prototypes()

        if seed < 0 or seed > 2**64 - 1:
            raise ValueError("seed must fit in uint64_t")

        c_seed: ctypes.c_uint64 = ctypes.c_uint64(seed)

        lib().apo_cuda_langevin_thermostat_integrator_set_thermostat_rng_seed(
            self.handle, c_seed
        )

        return

    def resetAverageTemperature(self) -> None:
        _initialize_prototypes()

        lib().apo_cuda_langevin_thermostat_integrator_reset_average_temperature(
            self.handle
        )

        return

    def getReferenceTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        lib().apo_cuda_langevin_thermostat_integrator_get_reference_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        return float(c_temperature.value)

    def getThermostatFriction(self) -> float:
        _initialize_prototypes()

        c_friction = ctypes.c_double()

        lib().apo_cuda_langevin_thermostat_integrator_get_thermostat_friction(
            ctypes.byref(c_friction), self.handle
        )

        return float(c_friction.value)

    def getThermostatRngSeed(self) -> int:
        _initialize_prototypes()

        c_seed = ctypes.c_uint64()

        lib().apo_cuda_langevin_thermostat_integrator_get_thermostat_rng_seed(
            ctypes.byref(c_seed), self.handle
        )

        return int(c_seed.value)

    def getAverageTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        lib().apo_cuda_langevin_thermostat_integrator_get_average_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        return float(c_temperature.value)

    def getInstantaneousTemperature(self) -> float:
        _initialize_prototypes()

        c_temperature = ctypes.c_double()

        lib().apo_cuda_langevin_thermostat_integrator_get_instantaneous_temperature(
            ctypes.byref(c_temperature), self.handle
        )

        return float(c_temperature.value)
