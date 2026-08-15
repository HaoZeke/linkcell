# C ABI

Header: `include/linkcell.h`. Prefix `lc_`. Caller owns every
buffer. C++ (`include/linkcell.hpp`) is a header-only wrap of these
entry points.

## Cell

```c
typedef struct lc_cell {
  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
  double ox, oy, oz;
} lc_cell;

lc_cell lc_cell_ortho(double lx, double ly, double lz);
```

Lattice vectors are rows `a`, `b`, `c` (the same order vesin uses)
plus an origin. `lc_cell_ortho` builds a diagonal box at the origin.

## Search

```c
int lc_knearest(const double *xyz, size_t n, const lc_cell *simbox,
                size_t k, const int *mask, double cell_hint,
                int *out_nn);
```

| Argument | Contract |
| --- | --- |
| `xyz` | `n` packed `x y z` triples. Not null. |
| `n` | Point count. `size_t`. Zero is an error. |
| `simbox` | Periodic parallelepiped. Not null. |
| `k` | Neighbours per source. `size_t`. Zero is an error. |
| `mask` | `NULL` (keep every point) or `n` ints, nonzero to include. A zero drops the point as both source and candidate. |
| `cell_hint` | Target bin edge. `<= 0` selects the default (3.0 in box units). |
| `out_nn` | Caller-owned, length `n * k`. Unused slots are `-1`. Neighbours of source `i` are `out_nn[i*k + t]`, nearest first. |

Returns 0 on success, nonzero on failure.

## Errors and version

```c
const char *lc_last_error(void);
const char *lc_version(void);
```

Both return pointers the caller must not free. `lc_last_error` is
thread-local and `NULL` after a successful `lc_*` call on that thread.
The pointer is invalid after the next `lc_*` call on the same thread.
Distinct searches may run concurrently; each thread reads its own
slot. `lc_version` is process-static.

## C++ wrap

```cpp
namespace linkcell {
class Error : public std::runtime_error { /* status() is the C return */ };

struct Cell {
  std::array<double, 3> a, b, c, origin;
  static Cell ortho(double lx, double ly, double lz);
  static Cell from_vectors(std::array<double, 3> a, std::array<double, 3> b,
                           std::array<double, 3> c,
                           std::array<double, 3> origin = {});
  lc_cell raw() const;
};

class Neighbours {
  int neighbour(std::size_t i, std::size_t j) const; /* or -1 */
  std::size_t n() const;
  std::size_t k() const;
  const int *data() const;
};

Neighbours knearest(const double *xyz, std::size_t n, const Cell &cell,
                    std::size_t k, int *out,
                    const int *mask = nullptr, double cell_hint = 0.0);

const char *version();
}
```

`knearest` writes the same packed `n * k` buffer as `lc_knearest` and
returns a non-owning `Neighbours` view. It does not return
`std::vector<std::vector<int>>`. An overload without `out` allocates
`std::vector<int>`. Failure throws `linkcell::Error`. Requires C++17.

## Rust errors (same crate)

`linkcell::Error`:

| Variant | Meaning |
| --- | --- |
| `ZeroK` | `k == 0` |
| `BadBox` | A box length is not strictly positive, or H is singular |
| `Empty` | The point list is empty |
| `BufferSize` | `knearest_into` / `out` length is not `n * k` |

`Error::Empty` is not a wrong-length buffer.

Walk and stop rule: [algorithm](algorithm.md).
