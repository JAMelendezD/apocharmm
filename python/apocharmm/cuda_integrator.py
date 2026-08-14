# BEGINLICENSE
# This file is part of apoCHARMM, which is distributed under the BSD 3-clause
# license, as described in the LICENSE file in the top level directory of this
# project.
#
# Author: James E. Gonzales II
#
# ENDLICENSE

import ctypes
import signal
import threading

from ._base import _ApoObject
from ._lib import encode_path, lib
from ._types import FilePath
from .error import configure_status_function

from .charmm_context import CharmmContext
from .subscriber import Subscriber

_prototypes_initialized: bool = False


def _initialize_prototypes() -> None:
    global _prototypes_initialized

    if _prototypes_initialized:
        return

    configure_status_function(
        lib().apo_cuda_integrator_set_time_step,
        [ctypes.c_void_p, ctypes.c_double],
        "CudaIntegrator.setTimeStep(time_step)",
    )

    configure_status_function(
        lib().apo_cuda_integrator_set_charmm_context,
        [ctypes.c_void_p, ctypes.c_void_p],
        "CudaIntegrator.setCharmmContext(context)",
    )

    configure_status_function(
        lib().apo_cuda_integrator_subscribe,
        [ctypes.c_void_p, ctypes.c_void_p],
        "CudaIntegrator.subscribe(subscriber)",
    )

    configure_status_function(
        lib().apo_cuda_integrator_unsubscribe,
        [ctypes.c_void_p, ctypes.c_void_p],
        "CudaIntegrator.unsubscribe(subscriber)",
    )

    configure_status_function(
        lib().apo_cuda_integrator_propagate,
        [ctypes.c_void_p, ctypes.c_int],
        "CudaIntegrator.propagate(num_steps)",
    )

    configure_status_function(
        lib().apo_cuda_integrator_initialize_from_restart_file,
        [ctypes.c_void_p, ctypes.c_char_p],
        "CudaIntegrator.initializeFromRestartFile(path)",
    )

    _prototypes_initialized = True

    return


class CudaIntegrator(_ApoObject):
    def __init__(self) -> None:
        super().__init__()

        self._integrator_handle: ctypes.c_void_p = ctypes.c_void_p()
        self._context: CharmmContext | None = None
        self._subscribers: list[Subscriber] = []

        return

    @property
    def integrator_handle(self) -> ctypes.c_void_p:
        if self._integrator_handle is None or self._integrator_handle.value is None:
            raise RuntimeError("apoCHARMM integrator object has been destroyed")

        return self._integrator_handle

    def close(self) -> None:
        super().close()

        self._integrator_handle = ctypes.c_void_p()
        self._context = None
        self._subscribers = []

        return

    def setTimeStep(self, time_step: float) -> None:
        _initialize_prototypes()

        c_time_step: ctypes.c_double = ctypes.c_double(time_step)

        lib().apo_cuda_integrator_set_time_step(self.integrator_handle, c_time_step)

        return

    def setCharmmContext(self, context: CharmmContext) -> None:
        _initialize_prototypes()

        if not isinstance(context, CharmmContext):
            raise TypeError("CudaIntegrator.setCharmmContext expects a CharmmContext")

        lib().apo_cuda_integrator_set_charmm_context(
            self.integrator_handle, context.handle
        )

        self._context = context

        return

    def subscribe(self, subscriber: Subscriber) -> None:
        _initialize_prototypes()

        if not isinstance(subscriber, Subscriber):
            raise TypeError("CudaIntegrator.subscribe expects a Subscriber")

        lib().apo_cuda_integrator_subscribe(
            self.integrator_handle, subscriber.subscriber_handle
        )

        self._subscribers.append(subscriber)

        return

    def unsubscribe(self, subscriber: Subscriber) -> None:
        _initialize_prototypes()

        if not isinstance(subscriber, Subscriber):
            raise TypeError("CudaIntegrator.unsubscribe expects a Subscriber")

        lib().apo_cuda_integrator_unsubscribe(
            self.integrator_handle, subscriber.subscriber_handle
        )

        self._subscribers.remove(subscriber)

        return

    def propagate(self, num_steps: int) -> None:
        _initialize_prototypes()

        if num_steps < -(2**31) or num_steps > 2**31 - 1:
            raise ValueError("num_steps must fit in int")

        c_num_steps: ctypes.c_int = ctypes.c_int(num_steps)

        restore_signal_handler: bool = (
            threading.current_thread() is threading.main_thread()
        )
        previous_sigint_handler = None

        if restore_signal_handler:
            previous_sigint_handler = signal.getsignal(signal.SIGINT)
            signal.signal(signal.SIGINT, signal.SIG_DFL)

        try:
            lib().apo_cuda_integrator_propagate(self.integrator_handle, c_num_steps)
        finally:
            if restore_signal_handler:
                signal.signal(signal.SIGINT, previous_sigint_handler)

        return

    def initializeFromRestartFile(self, path: FilePath) -> None:
        _initialize_prototypes()

        encoded_path: bytes = encode_path(path)
        c_path: ctypes.c_char_p = ctypes.c_char_p(encoded_path)

        lib().apo_cuda_integrator_initialize_from_restart_file(
            self.integrator_handle, c_path
        )

        return
