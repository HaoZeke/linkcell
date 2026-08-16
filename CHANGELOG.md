# Changelog

## Unreleased

## 0.2.1

- `Error::MaskLen` when `mask` is `Some` and `mask.len() != n`.
- `Error::Overflow` when `n * k` does not fit a slice. `BufferSize`
  remains a caller `out` whose length is not `n * k`.
- `lc_last_error` is written only by `lc_knearest`. `lc_version` does
  not clear the slot. An interior NUL in a message no longer drops the
  string to `NULL`.
- C++ `Neighbours` owns the packed `n * k` buffer. `knearest_into`
  takes `out_len`. There is no 5-argument `knearest(..., int *out)`.
- CI checks `include/linkcell.h` against cbindgen 0.29.4.
- Branding: sheared linked cells, k=4 neighbours, periodic wrap
  (`assets/branding/`).

## 0.2.0

ABI break. Existing 0.1.x C and C++ callers must rebuild.

- `lc_knearest` takes `size_t n` and `size_t k`. The 0.1.x `int` counts
  overflowed on large frames.
- `lc_last_error` is thread-local. Concurrent searches no longer share
  one process-wide slot. `lc_version` stays process-static.
- C++ `linkcell::knearest` writes a packed `n * k` index buffer (unused
  slots `-1`) and returns a `Neighbours` view. It no longer returns
  `std::vector<std::vector<int>>`. Failure throws `linkcell::Error`.
- Rust `Error::Empty` is an empty point list only. A wrong-length
  `knearest_into` buffer is `Error::BufferSize`. An overflowing
  linked-cell mesh is `Error::TooManyCells`.
- The per-source k-set keeps one image of each neighbour (the
  nearest). The same particle visited through two wraps does not
  occupy two slots.
- `knearest_brute` on a sheared cell takes the 27-image minimum.
  The single parallelepiped wrap is not the Wigner-Seitz cell of a
  60-degree hex prism.

## 0.1.2

Fold once, then Cartesian pair distances plus a lattice shift
(the vesin / LAMMPS ghost trick). Sources run in parallel. The C
ABI writes indices into the caller buffer.

## 0.1.1

Orthorhombic boxes use the three-wrap minimum image, not two
3x3 matvecs. The C ABI reads packed xyz in place. The k-heap for
k <= 16 stays on the stack.

## 0.1.0

Periodic linked-cell k-nearest neighbour search. Rust crate, C ABI
(`lc_*`), C++ header. The cell is a general parallelepiped;
orthorhombic boxes are `Cell::ortho` / `lc_cell_ortho`.

Installable from Meson (`linkcell_dep`, `pkg.generate`), CMake
(`find_package(linkcell)`, `linkcell::linkcell`), and pkg-config
(`linkcell.pc`).
