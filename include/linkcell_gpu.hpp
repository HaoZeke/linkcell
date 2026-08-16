#ifndef LINKCELL_GPU_HPP
#define LINKCELL_GPU_HPP

#include "linkcell.hpp"

#include <cstddef>

/** @file linkcell_gpu.hpp
 *  @brief Device-resident periodic k-nearest search.
 *
 *  Same contract as `knearest_into`: packed `n * k` indices, nearest
 *  first, unused slots `-1`. The walk is the host algorithm (fold,
 *  bin, Chebyshev shells until the k-th neighbour cannot sit outside
 *  the visited cube). `xyz` and `out` are CUDA device pointers.
 *  `cell` is host. Orthorhombic cells only.
 *
 *  The workspace keeps the bin arrays across calls (grow on demand).
 *  Built only when the meson `with_gpulite` feature is on.
 */

namespace linkcell {
namespace gpu {

#if defined(LINKCELL_HAS_GPULITE)

bool available();

class Workspace {
public:
  Workspace();
  ~Workspace();
  Workspace(Workspace &&) noexcept;
  Workspace &operator=(Workspace &&) noexcept;
  Workspace(const Workspace &) = delete;
  Workspace &operator=(const Workspace &) = delete;

  /// Device `xyz` is `n * 3` doubles. Device `out` is `n * k` ints.
  void knearest_into(const double *xyz, std::size_t n, const Cell &cell,
                     std::size_t k, int *out, std::size_t out_len,
                     const int *mask = nullptr, double cell_hint = 0.0);

  /// `nFrames` systems that share one orthorhombic cell. `xyz` is
  /// frame-major `nFrames * n * 3`, `out` is `nFrames * n * k`.
  void knearest_into_many(const double *xyz, std::size_t n,
                          std::size_t nFrames, const Cell &cell, std::size_t k,
                          int *out, std::size_t out_len,
                          const int *mask = nullptr, double cell_hint = 0.0);

private:
  struct Impl;
  Impl *impl_;
};

#else

inline bool available() { return false; }

#endif

} // namespace gpu
} // namespace linkcell

#endif
