# linkcell

Periodic **linked-cell k-nearest** neighbour search for molecular simulations.

vesin builds cutoff pair lists. nanoflann builds Euclidean KD-trees without a
minimum-image convention. This crate is the piece those two leave open: the
linked-cell walk of Allen and Tildesley (*Computer Simulation of Liquids*),
a k-heap per source, shells expanded until the k-th neighbour cannot sit
outside the visited cube.

It is a LODE library. The Rust crate is the implementation. The C ABI
(`lc_*`) is the hourglass waist, the same shape as [readcon-core](https://github.com/lode-org/readcon-core).
C++ is a RAII header over that ABI.

## Install

```
cargo add linkcell
```

C / C++ consumers build the `staticlib` (`--features capi`, on by default)
and include `include/linkcell.h` or `include/linkcell.hpp`.

## Rust

```rust
use linkcell::{knearest, OrthoBox};

let sim = OrthoBox::new(10.0, 10.0, 10.0)?;
let xyz = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]];
let rows = knearest(&xyz, &sim, 1, None, None)?;
assert_eq!(rows[0].indices, vec![1]);
```

`mask[i] == false` removes a point as both a source and a candidate.
`cell_hint` is the target cell edge; `None` uses 3.0 in the box units.

## C

```c
double xyz[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
lc_box box = {.lx = 10, .ly = 10, .lz = 10};
int out[2];
lc_knearest(xyz, 2, &box, 1, NULL, 0.0, out);
```

`out` has length `n * k`. Unused slots are `-1`. Neighbours of source `i`
are `out[i*k + 0 ..]`, nearest first.

## Design

- Orthorhombic periodic boxes only. Triclinic cells are a later crate
  revision, not a silent approximation.
- Distances are minimum-image, the same folding a LAMMPS dump box uses.
- The search does not take a cutoff. A cell-size hint only sets the bin
  width. Shells grow until the k-heap is exact.
- vesin remains the right library for a *cutoff* pair list.

## License

MIT
