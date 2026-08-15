# Changelog

## 0.1.0

Periodic linked-cell k-nearest neighbour search. Rust crate, C ABI
(`lc_*`), C++ header. The cell is a general parallelepiped;
orthorhombic boxes are `Cell::ortho` / `lc_cell_ortho`.

Installable from Meson (`linkcell_dep`, `pkg.generate`), CMake
(`find_package(linkcell)`, `linkcell::linkcell`), and pkg-config
(`linkcell.pc`).
