#ifndef LINKCELL_HPP
#define LINKCELL_HPP

#if defined(__cplusplus) && __cplusplus < 201703L
#error "linkcell.hpp requires C++17 or later"
#endif

#include "linkcell.h"

#include <array>
#include <stdexcept>
#include <string>
#include <vector>

namespace linkcell {

struct Cell {
  std::array<double, 3> a{};
  std::array<double, 3> b{};
  std::array<double, 3> c{};
  std::array<double, 3> origin{};

  static Cell ortho(double lx, double ly, double lz) {
    Cell cell;
    cell.a = {lx, 0.0, 0.0};
    cell.b = {0.0, ly, 0.0};
    cell.c = {0.0, 0.0, lz};
    return cell;
  }

  lc_cell raw() const {
    return lc_cell{a[0], a[1], a[2], b[0], b[1], b[2],
                   c[0], c[1], c[2], origin[0], origin[1], origin[2]};
  }
};

/// k-nearest neighbour indices, nearest first. Empty row if the point
/// was masked or isolated.
inline std::vector<std::vector<int>>
knearest(const std::vector<std::array<double, 3>> &xyz, const Cell &box, int k,
         const int *mask = nullptr, double cell_hint = 0.0) {
  const int n = static_cast<int>(xyz.size());
  if (n == 0) {
    return {};
  }
  std::vector<double> packed(static_cast<std::size_t>(n) * 3);
  for (int i = 0; i < n; i++) {
    packed[static_cast<std::size_t>(i) * 3 + 0] = xyz[static_cast<std::size_t>(i)][0];
    packed[static_cast<std::size_t>(i) * 3 + 1] = xyz[static_cast<std::size_t>(i)][1];
    packed[static_cast<std::size_t>(i) * 3 + 2] = xyz[static_cast<std::size_t>(i)][2];
  }
  std::vector<int> out(static_cast<std::size_t>(n) * static_cast<std::size_t>(k),
                       -1);
  const lc_cell raw = box.raw();
  if (lc_knearest(packed.data(), n, &raw, k, mask, cell_hint, out.data()) !=
      0) {
    const char *msg = lc_last_error();
    throw std::runtime_error(msg ? msg : "linkcell: knearest failed");
  }
  std::vector<std::vector<int>> rows(static_cast<std::size_t>(n));
  for (int i = 0; i < n; i++) {
    for (int t = 0; t < k; t++) {
      const int j = out[static_cast<std::size_t>(i) * static_cast<std::size_t>(k) +
                        static_cast<std::size_t>(t)];
      if (j >= 0) {
        rows[static_cast<std::size_t>(i)].push_back(j);
      }
    }
  }
  return rows;
}

inline const char *version() { return lc_version(); }

} // namespace linkcell

#endif
