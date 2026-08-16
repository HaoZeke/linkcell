#include "linkcell_gpu.hpp"

#ifdef LINKCELL_HAS_GPULITE

#include <gpulite/gpulite.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

using gpulite::CUDART;
using gpulite::KernelFactory;
using gpulite::NVRTC;
using gpulite::dim3;

namespace linkcell {
namespace gpu {
namespace {

constexpr int kMaxK = 16;
constexpr int kMaxCells = 262144;

const char *kKernels = R"CUDA(
__device__ inline int remEuclid(int a, int n) {
  int r = a % n;
  return r < 0 ? r + n : r;
}
__device__ inline int divEuclid(int a, int n) {
  int q = a / n;
  int r = a % n;
  if (r < 0) --q;
  return q;
}

extern "C" __global__ void bin_atoms(const double* __restrict__ xyz,
    const int* __restrict__ mask, int n, double lx, double ly, double lz,
    double ox, double oy, double oz, int nx, int ny, int nz,
    int* __restrict__ cellOf, int* __restrict__ cellCount,
    double* __restrict__ folded) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  if (mask && mask[i] == 0) {
    cellOf[i] = -1;
    folded[i * 3 + 0] = 0;
    folded[i * 3 + 1] = 0;
    folded[i * 3 + 2] = 0;
    return;
  }
  double fx = (xyz[i * 3 + 0] - ox) / lx;
  double fy = (xyz[i * 3 + 1] - oy) / ly;
  double fz = (xyz[i * 3 + 2] - oz) / lz;
  fx -= floor(fx); fy -= floor(fy); fz -= floor(fz);
  folded[i * 3 + 0] = fx * lx + ox;
  folded[i * 3 + 1] = fy * ly + oy;
  folded[i * 3 + 2] = fz * lz + oz;
  int cx = (int)(fx * nx); if (cx < 0) cx = 0; if (cx >= nx) cx = nx - 1;
  int cy = (int)(fy * ny); if (cy < 0) cy = 0; if (cy >= ny) cy = ny - 1;
  int cz = (int)(fz * nz); if (cz < 0) cz = 0; if (cz >= nz) cz = nz - 1;
  const int cid = (cz * ny + cy) * nx + cx;
  cellOf[i] = cid;
  atomicAdd(cellCount + cid, 1);
}

extern "C" __global__ void prefix_cells(int* __restrict__ cellCount,
    int* __restrict__ cellOff, int nC) {
  extern __shared__ int shared[];
  const int tid = threadIdx.x;
  const int nthreads = blockDim.x;
  const int chunk = (nC + nthreads - 1) / nthreads;
  const int start = tid * chunk;
  const int end = start + chunk < nC ? start + chunk : nC;
  int local = 0;
  for (int c = start; c < end; ++c) {
    cellOff[c] = local;
    local += cellCount[c];
    cellCount[c] = 0;
  }
  shared[tid] = local;
  __syncthreads();
  if (tid == 0) {
    int acc = 0;
    for (int t = 0; t < nthreads; ++t) {
      const int v = shared[t];
      shared[t] = acc;
      acc += v;
    }
    cellOff[nC] = acc;
  }
  __syncthreads();
  const int off = shared[tid];
  for (int c = start; c < end; ++c) cellOff[c] += off;
}

extern "C" __global__ void scatter_atoms(const int* __restrict__ cellOf,
    int* __restrict__ cellCount, const int* __restrict__ cellOff,
    int* __restrict__ order, int n) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  const int cid = cellOf[i];
  if (cid < 0) return;
  const int slot = atomicAdd(cellCount + cid, 1);
  order[cellOff[cid] + slot] = i;
}

extern "C" __global__ void knearest_shells(const double* __restrict__ folded,
    const int* __restrict__ cellOf, const int* __restrict__ cellOff,
    const int* __restrict__ order, int n, int nx, int ny, int nz,
    double lx, double ly, double lz, double cellMin, int k, int maxReach,
    int* __restrict__ out) {
  const int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i >= n) return;
  for (int t = 0; t < k; ++t) out[i * k + t] = -1;
  const int cid = cellOf[i];
  if (cid < 0) return;
  const int nxy = nx * ny;
  const int iz = cid / nxy;
  const int iy = (cid % nxy) / nx;
  const int ix = cid % nx;
  const double ixp = folded[i * 3 + 0];
  const double iyp = folded[i * 3 + 1];
  const double izp = folded[i * 3 + 2];
  double bestR2[16];
  int bestJ[16];
  int found = 0;
  for (int reach = 1; reach <= maxReach; ++reach) {
    for (int dx = -reach; dx <= reach; ++dx) {
      for (int dy = -reach; dy <= reach; ++dy) {
        for (int dz = -reach; dz <= reach; ++dz) {
          if (reach > 1 && dx != -reach && dx != reach && dy != -reach &&
              dy != reach && dz != -reach && dz != reach) {
            continue;
          }
          const int jx = ix + dx;
          const int jy = iy + dy;
          const int jz = iz + dz;
          const int ncx = remEuclid(jx, nx);
          const int ncy = remEuclid(jy, ny);
          const int ncz = remEuclid(jz, nz);
          const double sx = (double)divEuclid(jx, nx) * lx;
          const double sy = (double)divEuclid(jy, ny) * ly;
          const double sz = (double)divEuclid(jz, nz) * lz;
          const int ncid = (ncz * ny + ncy) * nx + ncx;
          const int a0 = cellOff[ncid];
          const int a1 = cellOff[ncid + 1];
          for (int s = a0; s < a1; ++s) {
            const int j = order[s];
            if (j == i) continue;
            const double rx = folded[j * 3 + 0] - ixp + sx;
            const double ry = folded[j * 3 + 1] - iyp + sy;
            const double rz = folded[j * 3 + 2] - izp + sz;
            const double r2 = rx * rx + ry * ry + rz * rz;
            if (r2 <= 0.0) continue;
            int have = -1;
            for (int t = 0; t < found; ++t) {
              if (bestJ[t] == j) {
                have = t;
                break;
              }
            }
            if (have >= 0) {
              if (r2 < bestR2[have]) bestR2[have] = r2;
              continue;
            }
            if (found < k) {
              bestR2[found] = r2;
              bestJ[found] = j;
              ++found;
            } else {
              int w = 0;
              for (int t = 1; t < k; ++t) {
                if (bestR2[t] > bestR2[w]) w = t;
              }
              if (r2 < bestR2[w]) {
                bestR2[w] = r2;
                bestJ[w] = j;
              }
            }
          }
        }
      }
    }
    if (found >= k) {
      double worst = bestR2[0];
      for (int t = 1; t < found; ++t) {
        if (bestR2[t] > worst) worst = bestR2[t];
      }
      const double bound = (double)reach * cellMin;
      if (worst <= bound * bound) break;
    }
  }
  for (int a = 1; a < found; ++a) {
    const double key = bestR2[a];
    const int id = bestJ[a];
    int p = a;
    while (p > 0 && (bestR2[p - 1] > key ||
                     (bestR2[p - 1] == key && bestJ[p - 1] > id))) {
      bestR2[p] = bestR2[p - 1];
      bestJ[p] = bestJ[p - 1];
      --p;
    }
    bestR2[p] = key;
    bestJ[p] = id;
  }
  for (int t = 0; t < found; ++t) out[i * k + t] = bestJ[t];
}
)CUDA";

void checkCuda(cudaError_t st, const char *what) {
  if (st != cudaSuccess) {
    throw Error(std::string(what) + ": " +
                CUDART::instance().cudaGetErrorString(st));
  }
}

struct DevPtr {
  void *p = nullptr;
  ~DevPtr() {
    if (p != nullptr && CUDART::loaded()) {
      CUDART::instance().cudaFree(p);
    }
  }
  void reset() {
    if (p != nullptr && CUDART::loaded()) {
      CUDART::instance().cudaFree(p);
    }
    p = nullptr;
  }
};

void growOne(DevPtr &d, std::size_t bytes, const char *what) {
  d.reset();
  checkCuda(CUDART::instance().cudaMalloc(&d.p, bytes), what);
}

bool cellIsOrtho(const Cell &cell) {
  constexpr double tol = 1.0e-12;
  return std::fabs(cell.a[1]) < tol && std::fabs(cell.a[2]) < tol &&
         std::fabs(cell.b[0]) < tol && std::fabs(cell.b[2]) < tol &&
         std::fabs(cell.c[0]) < tol && std::fabs(cell.c[1]) < tol &&
         cell.a[0] > 0.0 && cell.b[1] > 0.0 && cell.c[2] > 0.0;
}

} // namespace

bool available() {
  return CUDART::loaded() && NVRTC::loaded();
}

struct Workspace::Impl {
  int capN = 0;
  int capC = 0;
  DevPtr dcellOf, dcellCount, dcellOff, dorder, dfolded;

  void ensure(int n, int nC) {
    if (n > capN) {
      capN = n;
      growOne(dcellOf, static_cast<std::size_t>(capN) * sizeof(int), "cellOf");
      growOne(dorder, static_cast<std::size_t>(capN) * sizeof(int), "order");
      growOne(dfolded, static_cast<std::size_t>(capN) * 3 * sizeof(double),
              "folded");
    }
    if (nC > capC) {
      capC = nC;
      growOne(dcellCount, static_cast<std::size_t>(capC) * sizeof(int),
              "cellCount");
      growOne(dcellOff, static_cast<std::size_t>(capC + 1) * sizeof(int),
              "cellOff");
    }
  }
};

Workspace::Workspace() : impl_(new Impl) {}
Workspace::~Workspace() { delete impl_; }
Workspace::Workspace(Workspace &&o) noexcept : impl_(o.impl_) {
  o.impl_ = nullptr;
}
Workspace &Workspace::operator=(Workspace &&o) noexcept {
  if (this != &o) {
    delete impl_;
    impl_ = o.impl_;
    o.impl_ = nullptr;
  }
  return *this;
}

void Workspace::knearest_into(const double *xyz, std::size_t n, const Cell &cell,
                              std::size_t k, int *out, std::size_t out_len,
                              const int *mask, double cell_hint) {
  if (!available()) {
    throw Error("CUDA driver or nvrtc not loaded");
  }
  if (n == 0) {
    throw Error("no points");
  }
  if (k == 0) {
    throw Error("k must be at least 1");
  }
  if (k > static_cast<std::size_t>(kMaxK)) {
    throw Error("device k-nearest supports k <= 16");
  }
  if (out_len != n * k) {
    throw Error("out buffer length must be n * k");
  }
  if (xyz == nullptr || out == nullptr) {
    throw Error("null pointer");
  }
  if (!cellIsOrtho(cell)) {
    throw Error("device k-nearest is orthorhombic only");
  }

  double edge = cell_hint;
  if (!(edge > 0.0)) {
    edge = 3.0;
  }
  const double lx = cell.a[0];
  const double ly = cell.b[1];
  const double lz = cell.c[2];
  edge = std::min(edge, std::min(lx, std::min(ly, lz)));
  int nx = static_cast<int>(std::floor(lx / edge));
  int ny = static_cast<int>(std::floor(ly / edge));
  int nz = static_cast<int>(std::floor(lz / edge));
  if (nx < 1) {
    nx = 1;
  }
  if (ny < 1) {
    ny = 1;
  }
  if (nz < 1) {
    nz = 1;
  }
  const int nC = nx * ny * nz;
  if (nC <= 0 || nC > kMaxCells) {
    throw Error("too many cells");
  }
  const double cellMin =
      std::min(lx / static_cast<double>(nx),
               std::min(ly / static_cast<double>(ny), lz / static_cast<double>(nz)));
  const int maxReach = std::max(nx, std::max(ny, nz)) / 2 + 1;
  const int nI = static_cast<int>(n);
  const int kI = static_cast<int>(k);

  impl_->ensure(nI, nC);
  auto &rt = CUDART::instance();
  checkCuda(rt.cudaMemset(impl_->dcellCount.p, 0,
                          static_cast<std::size_t>(nC) * sizeof(int)),
            "zero cells");

  auto &factory = KernelFactory::instance(0);
  const std::vector<std::string> opt{"-std=c++17"};
  auto *kBin = factory.create("bin_atoms", kKernels, "linkcell.cu", opt);
  auto *kPref = factory.create("prefix_cells", kKernels, "linkcell.cu", opt);
  auto *kScat = factory.create("scatter_atoms", kKernels, "linkcell.cu", opt);
  auto *kNear = factory.create("knearest_shells", kKernels, "linkcell.cu", opt);

  const int block = 128;
  const int grid = (nI + block - 1) / block;
  void *xyzP = const_cast<double *>(xyz);
  void *maskV = const_cast<int *>(mask);
  void *outP = out;
  int nCv = nC;
  int maxR = maxReach;
  double cmin = cellMin;
  double ox = cell.origin[0];
  double oy = cell.origin[1];
  double oz = cell.origin[2];
  {
    std::vector<void *> a = {&xyzP, &maskV, &nI, &lx, &ly, &lz, &ox, &oy, &oz,
                             &nx,   &ny,    &nz, &impl_->dcellOf.p,
                             &impl_->dcellCount.p, &impl_->dfolded.p};
    kBin->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
  }
  {
    std::vector<void *> a = {&impl_->dcellCount.p, &impl_->dcellOff.p, &nCv};
    kPref->launch(dim3(1), dim3(128), 128 * sizeof(int), nullptr, a, true);
  }
  {
    std::vector<void *> a = {&impl_->dcellOf.p, &impl_->dcellCount.p,
                             &impl_->dcellOff.p, &impl_->dorder.p, &nI};
    kScat->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
  }
  {
    std::vector<void *> a = {&impl_->dfolded.p, &impl_->dcellOf.p,
                             &impl_->dcellOff.p, &impl_->dorder.p,
                             &nI, &nx, &ny, &nz, &lx, &ly, &lz, &cmin, &kI,
                             &maxR, &outP};
    kNear->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
  }
}

} // namespace gpu
} // namespace linkcell

#endif
