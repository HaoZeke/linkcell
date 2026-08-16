# Embed linkcell from C

Link the installed `staticlib` and include `linkcell.h`. The C ABI is
the hourglass waist; every other language wraps `lc_*`.

`lc_knearest` takes packed `double` xyz (`n` triples), `size_t n`,
`size_t k`, and a caller-owned `int` buffer of length `n * k`. Unused
slots are `-1`. Returns 0 on success. On failure, read `lc_last_error()`
on the same thread and do not free the pointer.

## Meson wrap

`subprojects/linkcell.wrap`:

```
[wrap-git]
url = https://github.com/d-SEAMS/linkcell.git
revision = v0.2.3
depth = 1

[provide]
linkcell = linkcell_dep
```

```meson
project('app', 'c', default_options: ['c_std=c11'])
linkcell_dep = dependency('linkcell', fallback: ['linkcell', 'linkcell_dep'])
executable('app', 'app.c', dependencies: linkcell_dep)
```

`linkcell_dep` already carries the Rust sysroot (pthread, dl, m on
Linux).

## CMake

After `cmake --install` into `$PREFIX`:

```cmake
cmake_minimum_required(VERSION 3.22)
project(app LANGUAGES C)
find_package(linkcell 0.2 REQUIRED)
add_executable(app app.c)
target_link_libraries(app PRIVATE linkcell::linkcell)
```

```
cmake -S . -B build -DCMAKE_PREFIX_PATH=$PREFIX
cmake --build build
```

In the linkcell build tree the same imported name is
`linkcell::linkcell`.

## Call

```c
#include "linkcell.h"

double xyz[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
lc_cell box = lc_cell_ortho(10.0, 10.0, 10.0);
int out[2];
if (lc_knearest(xyz, 2, &box, 1, NULL, 0.0, out) != 0) {
  fprintf(stderr, "%s\n", lc_last_error());
}
```

A sheared box is three lattice vectors plus an origin in `lc_cell`
(`ax..cz`, `ox..oz`), the same row layout vesin uses.

Full field list: [C ABI](../reference/c-api.md). A first program:
[two points](../tutorials/two-points.md).
