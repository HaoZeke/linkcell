#!/usr/bin/env bash
# Build, test, and install linkcell with Meson and CMake, then check
# pkg-config and find_package against the CMake prefix.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

MESON_BUILD="${LINKCELL_MESON_BUILD:-$ROOT/build-meson}"
CMAKE_BUILD="${LINKCELL_CMAKE_BUILD:-$ROOT/build-cmake}"
MESON_DEST="${LINKCELL_MESON_DEST:-$ROOT/destdir-meson}"
CMAKE_PREFIX="${LINKCELL_CMAKE_PREFIX:-$ROOT/prefix-cmake}"

echo "=== cargo test ==="
cargo test --locked --manifest-path "$ROOT/Cargo.toml"

echo "=== meson ==="
meson setup "$MESON_BUILD" --wipe
meson compile -C "$MESON_BUILD"
meson test -C "$MESON_BUILD" --print-errorlogs
rm -rf "$MESON_DEST"
meson install -C "$MESON_BUILD" --destdir "$MESON_DEST"

prefix_meson=""
for cand in \
  "$MESON_DEST/usr/local" \
  "$MESON_DEST/usr" \
  "$MESON_DEST"; do
  if [[ -f "$cand/include/linkcell.h" ]]; then
    prefix_meson="$cand"
    break
  fi
done
if [[ -z "$prefix_meson" ]]; then
  echo "meson install did not write include/linkcell.h under $MESON_DEST"
  find "$MESON_DEST" -type f | head
  exit 1
fi
test -f "$prefix_meson/include/linkcell.hpp"

pc_meson=""
for cand in \
  "$prefix_meson/lib/pkgconfig/linkcell.pc" \
  "$prefix_meson/lib/x86_64-linux-gnu/pkgconfig/linkcell.pc" \
  "$prefix_meson/share/pkgconfig/linkcell.pc"; do
  if [[ -f "$cand" ]]; then
    pc_meson="$(dirname "$cand")"
    break
  fi
done
if [[ -z "$pc_meson" ]]; then
  echo "meson install did not write linkcell.pc"
  find "$MESON_DEST" -name 'linkcell.pc'
  exit 1
fi
export PKG_CONFIG_PATH="$pc_meson${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
pkg-config --exists --print-errors linkcell
echo "meson pkg-config: $(pkg-config --modversion linkcell)"
echo "meson cflags: $(pkg-config --cflags linkcell)"
echo "meson libs: $(pkg-config --libs linkcell)"
unset PKG_CONFIG_PATH

echo "=== cmake ==="
cmake -S "$ROOT" -B "$CMAKE_BUILD" -DCMAKE_INSTALL_PREFIX="$CMAKE_PREFIX"
cmake --build "$CMAKE_BUILD"
ctest --test-dir "$CMAKE_BUILD" --output-on-failure
rm -rf "$CMAKE_PREFIX"
cmake --install "$CMAKE_BUILD"
test -f "$CMAKE_PREFIX/include/linkcell.h"
test -f "$CMAKE_PREFIX/include/linkcell.hpp"
test -f "$CMAKE_PREFIX/lib/liblinkcell.a" \
  -o -f "$CMAKE_PREFIX/lib64/liblinkcell.a"
test -f "$CMAKE_PREFIX/lib/pkgconfig/linkcell.pc" \
  -o -f "$CMAKE_PREFIX/lib64/pkgconfig/linkcell.pc"
test -f "$CMAKE_PREFIX/lib/cmake/linkcell/linkcellConfig.cmake" \
  -o -f "$CMAKE_PREFIX/lib64/cmake/linkcell/linkcellConfig.cmake"

pc_cmake=""
for cand in \
  "$CMAKE_PREFIX/lib/pkgconfig" \
  "$CMAKE_PREFIX/lib64/pkgconfig"; do
  if [[ -f "$cand/linkcell.pc" ]]; then
    pc_cmake="$cand"
    break
  fi
done
export PKG_CONFIG_PATH="$pc_cmake"
pkg-config --exists --print-errors linkcell
echo "cmake pkg-config: $(pkg-config --modversion linkcell)"
echo "cmake cflags: $(pkg-config --cflags linkcell)"
echo "cmake libs: $(pkg-config --libs linkcell)"

cmake -S "$ROOT/tests/cmake-consumer" -B "$CMAKE_BUILD/find-consumer" \
  -DCMAKE_PREFIX_PATH="$CMAKE_PREFIX" \
  -DLINKCELL_EXAMPLE="$ROOT/examples/two_points.cpp"
cmake --build "$CMAKE_BUILD/find-consumer"
"$CMAKE_BUILD/find-consumer/find_lc"

echo "OK"
