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
 *  `cell` is a host `Cell`, or a device-packed lattice (3 / 9 / 12
 *  doubles) inverted on device. General parallelepiped. `k <= 16`.
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
  /// Device `out_d2` is `n * k` doubles or nullptr.
  void knearest_into(const double *xyz, std::size_t n, const Cell &cell,
                     std::size_t k, int *out, std::size_t out_len,
                     const int *mask = nullptr, double cell_hint = 0.0,
                     double *out_d2 = nullptr);

  /// `nFrames` systems that share one cell. `xyz` is
  /// frame-major `nFrames * n * 3`, `out` / `out_d2` are
  /// `nFrames * n * k`. Kernels enqueue on `queue()`; `wait` is
  /// the stream barrier.
  void knearest_into_many(const double *xyz, std::size_t n,
                          std::size_t nFrames, const Cell &cell, std::size_t k,
                          int *out, std::size_t out_len,
                          const int *mask = nullptr, double cell_hint = 0.0,
                          bool wait = true, const double *frameBox = nullptr,
                          double *out_d2 = nullptr);

  /// Device-packed cell: 3 ortho lengths, 9 lattice rows, or 12
  /// rows plus origin. Invert, widths, and H/Hinv stay on device.
  /// The host reads four launch ints (`nx`, `ny`, `nz`, `nC`).
  void knearest_into_many_dcell(const double *xyz, std::size_t n,
                                std::size_t nFrames, const double *cell,
                                int cell_n, std::size_t k, int *out,
                                std::size_t out_len, const int *mask = nullptr,
                                double cell_hint = 0.0, bool wait = true,
                                double *out_d2 = nullptr);

  /// Persistent CUDA stream for this workspace. All device work is
  /// enqueued here; call `wait()` before reading device results.
  void *queue();
  void wait();

private:
  struct Impl;
  Impl *impl_;
  void launchWalk(const double *xyz, std::size_t n, std::size_t nFrames,
                  std::size_t k, int *out, const int *mask, bool wait,
                  double *out_d2, int nx, int ny, int nz, int nC,
                  int maxReach);
};

#else

inline bool available() { return false; }

#endif

} // namespace gpu
} // namespace linkcell

#endif
