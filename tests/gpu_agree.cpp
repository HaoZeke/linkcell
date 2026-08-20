#include "linkcell.hpp"
#include "linkcell_gpu.hpp"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

#ifdef LINKCELL_HAS_GPULITE
#include "gpulite_compat.hpp"
#endif

namespace {

int fail(const char *msg) {
  std::cerr << msg << '\n';
  return 1;
}

#ifdef LINKCELL_HAS_GPULITE
int runCase(const char *name, const std::vector<double> &xyz, std::size_t n,
            const linkcell::Cell &cell, std::size_t k) {
  std::vector<int> host(n * k, -2);
  linkcell::knearest_into(xyz.data(), n, cell, k, host.data(), host.size());

  auto &rt = gpulite::CUDART::instance();
  void *dxyz = nullptr;
  void *dout = nullptr;
  const std::size_t xyzBytes = n * 3 * sizeof(double);
  const std::size_t outBytes = n * k * sizeof(int);
  if (rt.cudaMalloc(&dxyz, xyzBytes) != cudaSuccess) {
    return fail("cudaMalloc xyz");
  }
  if (rt.cudaMalloc(&dout, outBytes) != cudaSuccess) {
    rt.cudaFree(dxyz);
    return fail("cudaMalloc out");
  }
  rt.cudaMemcpy(dxyz, xyz.data(), xyzBytes, cudaMemcpyHostToDevice);
  linkcell::gpu::Workspace ws;
  ws.knearest_into(static_cast<const double *>(dxyz), n, cell, k,
                   static_cast<int *>(dout), n * k);
  std::vector<int> dev(n * k, -3);
  rt.cudaMemcpy(dev.data(), dout, outBytes, cudaMemcpyDeviceToHost);
  rt.cudaFree(dxyz);
  rt.cudaFree(dout);

  for (std::size_t i = 0; i < n * k; ++i) {
    if (host[i] != dev[i]) {
      std::cerr << name << " mismatch at " << i << ": host " << host[i]
                << " device " << dev[i] << '\n';
      return 1;
    }
  }
  return 0;
}
#endif

} // namespace

int main() {
#ifndef LINKCELL_HAS_GPULITE
  return 0;
#else
  if (!linkcell::gpu::available()) {
    std::cerr << "skip: no CUDA/nvrtc\n";
    return 0;
  }
  const linkcell::Cell ortho = linkcell::Cell::ortho(10.0, 10.0, 10.0);
  const std::vector<double> face = {0.2, 1.0, 1.0, 9.7, 1.0, 1.0};
  if (runCase("ortho-face", face, 2, ortho, 1) != 0) {
    return 1;
  }

  const linkcell::Cell hex = linkcell::Cell::from_vectors(
      {10.0, 0.0, 0.0}, {5.0, 8.660254037844386, 0.0}, {0.0, 0.0, 10.0});
  const std::vector<double> hexFace = {0.2, 0.1, 1.0, 9.7, 0.1, 1.0};
  if (runCase("hex-prism-face", hexFace, 2, hex, 1) != 0) {
    return 1;
  }

  const linkcell::Cell tri = linkcell::Cell::from_vectors(
      {10.0, 0.0, 0.0}, {3.0, 9.0, 0.0}, {1.0, 2.0, 8.0});
  std::vector<double> pts;
  pts.reserve(48 * 3);
  std::uint64_t s = 7;
  auto frac = [&s]() {
    s = s * 6364136223846793005ULL + 1;
    return static_cast<double>(s >> 11) / static_cast<double>(1ULL << 53);
  };
  for (int i = 0; i < 48; ++i) {
    const double u = frac();
    const double v = frac();
    const double w = frac();
    pts.push_back(10.0 * u + 3.0 * v + 1.0 * w);
    pts.push_back(9.0 * v + 2.0 * w);
    pts.push_back(8.0 * w);
  }
  if (runCase("triclinic-48-k4", pts, 48, tri, 4) != 0) {
    return 1;
  }
  return 0;
#endif
}
