#include "linkcell.hpp"

#include <iostream>

int main() {
  const linkcell::Cell box = linkcell::Cell::ortho(10.0, 10.0, 10.0);
  const std::vector<std::array<double, 3>> xyz{{0.0, 0.0, 0.0},
                                               {1.0, 0.0, 0.0}};
  const auto rows = linkcell::knearest(xyz, box, 1);
  if (rows.size() != 2 || rows[0].empty() || rows[1].empty()) {
    std::cerr << "unexpected neighbour rows\n";
    return 1;
  }
  std::cout << "0 -> " << rows[0][0] << "\n1 -> " << rows[1][0] << "\n";
  return 0;
}
