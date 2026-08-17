"""Periodic linked-cell k-nearest search.

Arrays are DLPack. Pass any ``__dlpack__()`` object (numpy, torch,
jax, cupy, metatomic). The result is a DLPack int32 tensor; use
``numpy.from_dlpack`` or ``torch.from_dlpack``.
"""

from linkcell._lib import __version__, gpu_available, knearest

__all__ = ["__version__", "gpu_available", "knearest"]
