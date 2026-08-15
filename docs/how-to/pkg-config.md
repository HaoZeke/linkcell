# Use pkg-config

Meson and CMake both write `linkcell.pc` at install. The file names
the archive and the Rust sysroot the staticlib still needs.

## Prefix

After `meson install` or `cmake --install` into `$PREFIX`:

```
export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig:$PREFIX/lib64/pkgconfig
pkg-config --exists --print-errors linkcell
pkg-config --modversion linkcell
pkg-config --cflags --libs linkcell
```

On Debian multiarch the `.pc` may sit under
`$PREFIX/lib/x86_64-linux-gnu/pkgconfig`. Put that directory on
`PKG_CONFIG_PATH` instead.

`--cflags` is `-I$PREFIX/include`. `--libs` is `-L$PREFIX/lib
-llinkcell` plus the sysroot: `-ldl -lpthread -lm` on Linux,
`-lpthread -lm` on macOS.

## Compile

```
cc app.c $(pkg-config --cflags --libs linkcell) -o app
c++ -std=c++17 app.cpp $(pkg-config --cflags --libs linkcell) -o app
```

The search is a static archive. The order is the object, then
`--libs`. Do not drop the sysroot flags and add them by hand; the
`.pc` is the list.

## Check

`scripts/check-install.sh` in this repository installs both Meson and
CMake prefixes and runs `pkg-config --exists linkcell` against each.
