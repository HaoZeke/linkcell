"""Periodic linked-cell k-nearest search.

Arrays are DLPack. Pass any ``__dlpack__()`` object (numpy, torch,
jax, cupy, metatomic) for ``xyz`` and ``cell``. The result is
``(indices, dist2)``; use ``numpy.from_dlpack`` or
``torch.from_dlpack`` on each.
"""

from linkcell._lib import __version__, gpu_available, knearest

__all__ = ["__version__", "gpu_available", "knearest"]
