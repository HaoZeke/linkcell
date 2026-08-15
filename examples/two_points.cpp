#include "linkcell.hpp"

#include <cstddef>
#include <exception>
#include <iostream>

// Neighbours of source i sit at nn.neighbour(i, j). Unused slots are -1.
static int run(const char *label, const double *xyz, std::size_t n,
               const linkcell::Cell &box, std::size_t k) {
  const linkcell::Neighbours nn = linkcell::knearest(xyz, n, box, k);
  std::cout << label << "\n";
  for (std::size_t i = 0; i < n; ++i) {
    std::cout << i << " ->";
    for (std::size_t t = 0; t < k; ++t) {
      const int j = nn.neighbour(i, t);
      if (t == 0 && j < 0) {
        std::cerr << "unexpected neighbour rows\n";
        return 1;
      }
      std::cout << " " << j;
    }
    std::cout << "\n";
  }
  return 0;
}

int main() {
  const double xyz[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  const std::size_t n = 2;
  const std::size_t k = 1;
  const linkcell::Cell ortho = linkcell::Cell::ortho(10.0, 10.0, 10.0);
  const linkcell::Cell sheared = linkcell::Cell::from_vectors(
      {10.0, 0.0, 0.0}, {5.0, 8.66, 0.0}, {0.0, 0.0, 10.0});

  try {
    if (run("ortho", xyz, n, ortho, k) != 0 ||
        run("sheared", xyz, n, sheared, k) != 0) {
      return 1;
    }
  } catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    return 1;
  }
  return 0;
}
