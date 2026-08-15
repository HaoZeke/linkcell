# Embed linkcell from C++

`include/linkcell.hpp` is a C++17 RAII header over `linkcell.h`. It
does not add a second library. Link the same `staticlib` a C consumer
links.

`linkcell::knearest` writes packed indices into a caller-owned
`int *` of length `n * k` and returns a `Neighbours` view over that
buffer. It does not return `std::vector<std::vector<int>>`. Unused
slots are `-1`. Failure throws `linkcell::Error` (`std::runtime_error`
with the C status and the thread-local `lc_last_error` text).

## Meson wrap

Same wrap as C. Ask Meson for C++:

```meson
project('app', 'cpp', default_options: ['cpp_std=c++17'])
linkcell_dep = dependency('linkcell', fallback: ['linkcell', 'linkcell_dep'])
executable('app', 'app.cpp', dependencies: linkcell_dep)
```

## CMake

```cmake
cmake_minimum_required(VERSION 3.22)
project(app LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
find_package(linkcell 0.2 REQUIRED)
add_executable(app app.cpp)
target_link_libraries(app PRIVATE linkcell::linkcell)
```

`tests/cmake-consumer/` in this repository is that pattern against an
installed prefix.

## Call

```cpp
#include "linkcell.hpp"

const linkcell::Cell box = linkcell::Cell::ortho(10.0, 10.0, 10.0);
const double xyz[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
int out[2];
const linkcell::Neighbours nn = linkcell::knearest(xyz, 2, box, 1, out);
```

A general parallelepiped is `linkcell::Cell::from_vectors(a, b, c,
origin)`. `Cell::raw()` is the `lc_cell` the ABI consumes.

`mask` is an optional `const int *` of length `n` (nonzero keeps the
point). `cell_hint <= 0` selects the default bin edge.

Signatures: [C ABI](../reference/c-api.md). First program:
[two points](../tutorials/two-points.md).
