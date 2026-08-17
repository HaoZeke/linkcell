# linkcell

<p align="center">
  <img src="assets/branding/linkcell-logo-light.svg" width="360" alt="linkcell">
</p>

Periodic **linked-cell k-nearest** neighbour search for molecular simulations.

vesin builds cutoff pair lists. nanoflann builds Euclidean KD-trees without a
minimum-image convention. This crate is the piece those two leave open: the
linked-cell walk of Allen and Tildesley (*Computer Simulation of Liquids*),
a k-heap per source, shells expanded until the k-th neighbour cannot sit
outside the visited cube. The optional gpulite path runs that walk on a
CUDA device (`linkcell::gpu::Workspace`); pair lists stay on the device.

It is a LODE library. The Rust crate is the implementation. The C ABI
(`lc_*`) is the hourglass waist, the same shape as
[readcon-core](https://github.com/lode-org/readcon-core). C++ is a RAII
header over that ABI.

## Install

Rust:

```
cargo add linkcell
```

C and C++ consumers take the `staticlib` (`--features capi`, on by default)
plus `include/linkcell.h` or `include/linkcell.hpp`. Meson, CMake, and
pkg-config all install that archive and those headers.

Python takes the same search through DLPack (dlpk). Any `__dlpack__()`
object (numpy, torch, jax, cupy) is a valid `xyz` / `cell`. A CUDA
`xyz` (`torch.Tensor` on GPU) stays on device; `torch.from_dlpack`
consumes the result. `cell` is host.
Wheels: one CPython 3.12 limited-ABI (`abi3`) artifact per platform, and
one free-threaded set (`cp313t` / `cp314t`).

```
pip install linkcell
```

```python
import numpy as np
import linkcell

xyz = np.array([[0.2, 0.0, 0.0], [9.4, 0.0, 0.0]], dtype=np.float64)
cell = np.array([10.0, 10.0, 10.0], dtype=np.float64)
nn = np.from_dlpack(linkcell.knearest(xyz, cell, 1))
```

### Meson

```
meson setup build
meson compile -C build
meson install -C build
```

As a wrap, Meson exposes `linkcell_dep`:

```
[wrap-git]
url = https://github.com/d-SEAMS/linkcell.git
revision = v0.3.0
depth = 1

[provide]
linkcell = linkcell_dep
```

```meson
linkcell_dep = dependency('linkcell', fallback: ['linkcell', 'linkcell_dep'])
```

### CMake

```
cmake -B build -DCMAKE_INSTALL_PREFIX=$PREFIX
cmake --build build
cmake --install build
```

```cmake
find_package(linkcell 0.2 REQUIRED)
target_link_libraries(app PRIVATE linkcell::linkcell)
```

In the same build tree the target is `linkcell::linkcell`.

### pkg-config

```
pkg-config --cflags --libs linkcell
```

Both Meson and CMake write `linkcell.pc` (Libs includes the Rust
sysroot: pthread, dl, m on Linux).

## Rust

```rust
use linkcell::{knearest, Cell};

let sim = Cell::ortho(10.0, 10.0, 10.0)?;
let sheared = Cell::from_vectors(
    [10.0, 0.0, 0.0],
    [5.0, 8.66, 0.0],
    [0.0, 0.0, 10.0],
    [0.0, 0.0, 0.0],
)?;
let xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
let rows = knearest(&xyz, &sim, 1, None, None)?;
assert_eq!(rows[0].indices, vec![1]);
```

`mask[i] == false` removes a point as both a source and a candidate.
`cell_hint` is the target cell edge; `None` uses 3.0 in the box units.

## C

```c
#include "linkcell.h"

double xyz[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
lc_cell box = lc_cell_ortho(10.0, 10.0, 10.0);
int out[2];
if (lc_knearest(xyz, 2, &box, 1, NULL, 0.0, out) != 0) {
  return 1;
}
```

`n` and `k` are `size_t`. `out` has length `n * k`. Unused slots are `-1`.
Neighbours of source `i` are `out[i*k + 0 ..]`, nearest first.

## C++

```cpp
#include "linkcell.hpp"

const linkcell::Cell box = linkcell::Cell::ortho(10.0, 10.0, 10.0);
const double xyz[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
const linkcell::Neighbours nn = linkcell::knearest(xyz, 2, box, 1);
```

`nn` owns the packed `n * k` buffer. Unused slots are `-1`.
`nn.neighbour(i, j)` is the j-th neighbour of i. Failure throws
`linkcell::Error`.

The optional device walk reads occupancy from the environment:
`LINKCELL_TPP` (threads per particle) and `LINKCELL_BLOCK` (CUDA
block size). Unset, the library picks a pair that maximises particles
per block under the device thread and 48 KiB shared-memory limits.
d-SEAMS writes the same keys from `SEAMS_CONFIG` / `seams --tpp`.

## Docs

| Quadrant | Page |
| --- | --- |
| Tutorial | [Two points](docs/tutorials/two-points.md) |
| How-to | [Embed from C](docs/how-to/embed-c.md), [Embed from C++](docs/how-to/embed-cpp.md), [pkg-config](docs/how-to/pkg-config.md) |
| Reference | [C ABI](docs/reference/c-api.md), [Algorithm](docs/reference/algorithm.md) |
| Explanation | [MIC and cells](docs/explanation/mic-and-cells.md) |

Rust API: [docs.rs/linkcell](https://docs.rs/linkcell). Map: [docs/index.md](docs/index.md).

## Design

- The cell is a general parallelepiped: three lattice vectors plus an
  origin. Orthorhombic boxes are `Cell::ortho` / `lc_cell_ortho`.
  Binning is in fractional space, so a sheared dump is not treated as
  orthogonal.
- Points fold into the primary cell once. Each source then walks
  Chebyshev shells of neighbour cells. Pair distances are a Cartesian
  subtract plus that cell's lattice translation (`dist2_shifted` and
  `lattice_shift`), the vesin / LAMMPS ghost construction. Orthorhombic
  boxes use three independent wraps and skip the two Hinv matvecs.
- One lattice shift per unique cell is wrong unless every wrap of that
  cell is visited. The walk visits integer cell offsets, so each wrap
  of a bin is a separate visit.
- The search does not take a cutoff. A cell-size hint only sets the bin
  width. Shells grow until the k-heap is exact.
- vesin remains the right library for a *cutoff* pair list.

## License

MIT
