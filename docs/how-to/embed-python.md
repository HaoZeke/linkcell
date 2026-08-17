# Embed linkcell from Python

The Python module wraps the C waist with DLPack (dlpk). It does not
take `numpy.ndarray` as the ABI. Any object with `__dlpack__()` is
valid input. The return value implements `__dlpack__()`; consume it
with `numpy.from_dlpack` or `torch.from_dlpack`.

```python
import numpy as np
import linkcell

xyz = np.ascontiguousarray([[0.2, 0.0, 0.0], [9.4, 0.0, 0.0]], dtype=np.float64)
cell = np.ascontiguousarray([10.0, 10.0, 10.0], dtype=np.float64)
nn = np.from_dlpack(linkcell.knearest(xyz, cell, k=1))
assert int(nn[0, 0]) == 1
```

`xyz` is float64 `(n, 3)`. `cell` is float64 `(3,)` (ortho lengths),
`(3, 3)` lattice rows (vesin / `lc_cell` order), or `(4, 3)` rows plus
origin. Optional `mask` is length `n`. The GIL is detached for the
rayon CPU walk.

CUDA `xyz` (`kDLCUDA`, including a `torch.Tensor` on GPU) is passed
by device pointer into `lc_gpu_knearest`. `cell` stays host. The
return capsule is on the same CUDA device; `torch.from_dlpack`
takes it without a host bounce. `linkcell.gpu_available()` is true
when the driver and nvrtc load (gpulite, no CUDA SDK at build
time).

Wheels:

- limited ABI: one `abi3` wheel per platform, CPython 3.12+
- free-threaded set: `cp313t` and `cp314t` (no stable ABI on `t`
  until PyO3 `abi3t` / CPython 3.15)
