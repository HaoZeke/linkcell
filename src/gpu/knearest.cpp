#include "linkcell_gpu.hpp"

#ifdef LINKCELL_HAS_GPULITE

#include <gpulite/gpulite.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
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
#define TPP 8
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

__device__ inline void heapPush(double* bestR2, int* bestJ, int* found, int k,
    double r2, int j) {
  int have = -1;
  for (int t = 0; t < *found; ++t) {
    if (bestJ[t] == j) {
      have = t;
      break;
    }
  }
  if (have >= 0) {
    if (r2 < bestR2[have]) bestR2[have] = r2;
    return;
  }
  if (*found < k) {
    bestR2[*found] = r2;
    bestJ[*found] = j;
    ++(*found);
    return;
  }
  int w = 0;
  for (int t = 1; t < k; ++t) {
    if (bestR2[t] > bestR2[w]) w = t;
  }
  if (r2 < bestR2[w]) {
    bestR2[w] = r2;
    bestJ[w] = j;
  }
}

extern "C" __global__ void bin_atoms(const double* __restrict__ xyz,
    const int* __restrict__ mask, int n, int nFrames, double lx, double ly,
    double lz, double ox, double oy, double oz, int nx, int ny, int nz,
    int nC, int* __restrict__ cellOf, int* __restrict__ cellCount,
    double* __restrict__ folded) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n * nFrames) return;
  const int f = tid / n;
  const int i = tid % n;
  const int id = f * n + i;
  if (mask && mask[id] == 0) {
    cellOf[id] = -1;
    folded[id * 3 + 0] = 0;
    folded[id * 3 + 1] = 0;
    folded[id * 3 + 2] = 0;
    return;
  }
  double fx = (xyz[id * 3 + 0] - ox) / lx;
  double fy = (xyz[id * 3 + 1] - oy) / ly;
  double fz = (xyz[id * 3 + 2] - oz) / lz;
  fx -= floor(fx); fy -= floor(fy); fz -= floor(fz);
  folded[id * 3 + 0] = fx * lx + ox;
  folded[id * 3 + 1] = fy * ly + oy;
  folded[id * 3 + 2] = fz * lz + oz;
  int cx = (int)(fx * nx); if (cx < 0) cx = 0; if (cx >= nx) cx = nx - 1;
  int cy = (int)(fy * ny); if (cy < 0) cy = 0; if (cy >= ny) cy = ny - 1;
  int cz = (int)(fz * nz); if (cz < 0) cz = 0; if (cz >= nz) cz = nz - 1;
  const int cid = (cz * ny + cy) * nx + cx;
  cellOf[id] = cid;
  atomicAdd(cellCount + f * nC + cid, 1);
}

// Tiled Blelloch exclusive scan, one frame per block (HOOMD / CUB DeviceScan).
extern "C" __global__ void prefix_cells(int* __restrict__ cellCount,
    int* __restrict__ cellOff, int nC, int nFrames) {
  extern __shared__ int sh[];
  const int f = blockIdx.x;
  if (f >= nFrames) return;
  const int tid = threadIdx.x;
  const int nthreads = blockDim.x;
  int* counts = cellCount + f * nC;
  int* offs = cellOff + f * (nC + 1);
  int carry = 0;
  for (int base = 0; base < nC; base += nthreads) {
    const int i = base + tid;
    const int val = (i < nC) ? counts[i] : 0;
    sh[tid] = val;
    __syncthreads();
    for (int d = 1; d < nthreads; d <<= 1) {
      int add = 0;
      if (tid >= d) add = sh[tid - d];
      __syncthreads();
      sh[tid] += add;
      __syncthreads();
    }
    // sh[tid] is inclusive; exclusive is sh[tid] - val
    if (i < nC) {
      offs[i] = carry + sh[tid] - val;
      counts[i] = 0;
    }
    __syncthreads();
    if (tid == nthreads - 1) {
      sh[0] = carry + sh[tid];
    }
    __syncthreads();
    carry = sh[0];
    __syncthreads();
  }
  if (tid == 0) offs[nC] = carry;
}

// HOOMD / vesin: store coordinates in cell order for coalesced stencil reads.
extern "C" __global__ void scatter_atoms(const double* __restrict__ folded,
    const int* __restrict__ cellOf, int* __restrict__ cellCount,
    const int* __restrict__ cellOff, int* __restrict__ order,
    double* __restrict__ sorted, int* __restrict__ home, int n, int nFrames,
    int nC) {
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid >= n * nFrames) return;
  const int f = tid / n;
  const int i = tid % n;
  const int id = f * n + i;
  const int cid = cellOf[id];
  if (cid < 0) {
    home[id] = -1;
    return;
  }
  const int slot = atomicAdd(cellCount + f * nC + cid, 1);
  const int dest = f * n + cellOff[f * (nC + 1) + cid] + slot;
  order[dest] = i;
  home[id] = dest;
  sorted[dest * 3 + 0] = folded[id * 3 + 0];
  sorted[dest * 3 + 1] = folded[id * 3 + 1];
  sorted[dest * 3 + 2] = folded[id * 3 + 2];
}

// HOOMD NeighborListGPUBinned: TPP threads split the cell occupants.
// Cabana VerletLayout2D: packed n*k output. Host stop: worst <= (reach*cell_min)^2.
extern "C" __global__ void knearest_shells(const double* __restrict__ sorted,
    const int* __restrict__ cellOf, const int* __restrict__ cellOff,
    const int* __restrict__ order, const int* __restrict__ home,
    const int* __restrict__ sdx, const int* __restrict__ sdy,
    const int* __restrict__ sdz, const int* __restrict__ reachOff,
    int n, int nFrames, int nC, int nx, int ny, int nz,
    double lx, double ly, double lz, double cellMin, int k, int maxReach,
    int* __restrict__ out) {
  const int lane = threadIdx.x % TPP;
  const int gid = (blockIdx.x * blockDim.x + threadIdx.x) / TPP;
  const int active = gid < n * nFrames;
  const int f = active ? gid / n : 0;
  const int i = active ? gid % n : 0;
  const int id = f * n + i;
  const int cid = active ? cellOf[id] : -1;
  if (active && lane == 0) {
    for (int t = 0; t < k; ++t) out[id * k + t] = -1;
  }
  const int nxy = nx * ny;
  const int iz = cid >= 0 ? cid / nxy : 0;
  const int iy = cid >= 0 ? (cid % nxy) / nx : 0;
  const int ix = cid >= 0 ? cid % nx : 0;
  double ixp = 0, iyp = 0, izp = 0;
  if (cid >= 0) {
    const int dest = home[id];
    ixp = sorted[dest * 3 + 0];
    iyp = sorted[dest * 3 + 1];
    izp = sorted[dest * 3 + 2];
  }
  double bestR2[16];
  int bestJ[16];
  int found = 0;
  extern __shared__ unsigned char raw[];
  const int ppb = blockDim.x / TPP;
  const int pslot = threadIdx.x / TPP;
  double* sh_d2 = (double*)raw;
  int* sh_j = (int*)(sh_d2 + ppb * TPP * 16);
  int* sh_n = (int*)(sh_j + ppb * TPP * 16);
  int* sh_stop = sh_n + ppb * TPP;
  for (int reach = 1; reach <= maxReach; ++reach) {
    const int s0 = reachOff[reach];
    const int s1 = reachOff[reach + 1];
    if (cid >= 0) {
      for (int st = s0; st < s1; ++st) {
        const int dx = sdx[st];
        const int dy = sdy[st];
        const int dz = sdz[st];
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
        const int a0 = f * n + cellOff[f * (nC + 1) + ncid];
        const int a1 = f * n + cellOff[f * (nC + 1) + ncid + 1];
        for (int s = a0 + lane; s < a1; s += TPP) {
          const int j = order[s];
          if (j == i) continue;
          const double rx = sorted[s * 3 + 0] - ixp + sx;
          const double ry = sorted[s * 3 + 1] - iyp + sy;
          const double rz = sorted[s * 3 + 2] - izp + sz;
          const double r2 = rx * rx + ry * ry + rz * rz;
          if (r2 <= 0.0) continue;
          heapPush(bestR2, bestJ, &found, k, r2, j);
        }
      }
    }
    const int base = (pslot * TPP + lane) * 16;
    sh_n[pslot * TPP + lane] = found;
    for (int t = 0; t < 16; ++t) {
      sh_d2[base + t] = t < found ? bestR2[t] : 0;
      sh_j[base + t] = t < found ? bestJ[t] : -1;
    }
    __syncwarp();
    if (lane == 0) {
      found = 0;
      for (int L = 0; L < TPP; ++L) {
        const int nL = sh_n[pslot * TPP + L];
        const int bL = (pslot * TPP + L) * 16;
        for (int t = 0; t < nL; ++t) {
          heapPush(bestR2, bestJ, &found, k, sh_d2[bL + t], sh_j[bL + t]);
        }
      }
      int stop = 0;
      if (found >= k) {
        double worst = bestR2[0];
        for (int t = 1; t < found; ++t) {
          if (bestR2[t] > worst) worst = bestR2[t];
        }
        const double bound = (double)reach * cellMin;
        if (worst <= bound * bound) stop = 1;
      }
      sh_n[pslot * TPP] = found;
      sh_stop[pslot] = stop;
      for (int t = 0; t < found; ++t) {
        sh_d2[pslot * TPP * 16 + t] = bestR2[t];
        sh_j[pslot * TPP * 16 + t] = bestJ[t];
      }
    }
    __syncwarp();
    found = sh_n[pslot * TPP];
    for (int t = 0; t < found; ++t) {
      bestR2[t] = sh_d2[pslot * TPP * 16 + t];
      bestJ[t] = sh_j[pslot * TPP * 16 + t];
    }
    if (sh_stop[pslot]) break;
  }
  if (lane != 0) return;
  for (int a = 1; a < found; ++a) {
    const double key = bestR2[a];
    const int idj = bestJ[a];
    int p = a;
    while (p > 0 && (bestR2[p - 1] > key ||
                     (bestR2[p - 1] == key && bestJ[p - 1] > idj))) {
      bestR2[p] = bestR2[p - 1];
      bestJ[p] = bestJ[p - 1];
      --p;
    }
    bestR2[p] = key;
    bestJ[p] = idj;
  }
  for (int t = 0; t < found; ++t) out[id * k + t] = bestJ[t];
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

void buildStencil(int maxReach, std::vector<int> &dx, std::vector<int> &dy,
                  std::vector<int> &dz, std::vector<int> &reachOff) {
  dx.clear();
  dy.clear();
  dz.clear();
  reachOff.assign(static_cast<std::size_t>(maxReach) + 2, 0);
  for (int r = 1; r <= maxReach; ++r) {
    reachOff[static_cast<std::size_t>(r)] = static_cast<int>(dx.size());
    for (int ix = -r; ix <= r; ++ix) {
      for (int iy = -r; iy <= r; ++iy) {
        for (int iz = -r; iz <= r; ++iz) {
          if (r > 1 && std::abs(ix) != r && std::abs(iy) != r &&
              std::abs(iz) != r) {
            continue;
          }
          dx.push_back(ix);
          dy.push_back(iy);
          dz.push_back(iz);
        }
      }
    }
  }
  reachOff[static_cast<std::size_t>(maxReach) + 1] = static_cast<int>(dx.size());
}

struct Workspace::Impl {
  int capN = 0;
  int capC = 0;
  int capSt = 0;
  int cachedReach = -1;
  DevPtr dcellOf, dcellCount, dcellOff, dorder, dfolded, dsorted, dhome;
  DevPtr dsdx, dsdy, dsdz, dreachOff;

  void ensure(int nTot, int nCtot, int nOff) {
    if (nTot > capN) {
      capN = nTot;
      growOne(dcellOf, static_cast<std::size_t>(capN) * sizeof(int), "cellOf");
      growOne(dorder, static_cast<std::size_t>(capN) * sizeof(int), "order");
      growOne(dhome, static_cast<std::size_t>(capN) * sizeof(int), "home");
      growOne(dfolded, static_cast<std::size_t>(capN) * 3 * sizeof(double),
              "folded");
      growOne(dsorted, static_cast<std::size_t>(capN) * 3 * sizeof(double),
              "sorted");
    }
    if (nCtot > capC) {
      capC = nCtot;
      growOne(dcellCount, static_cast<std::size_t>(capC) * sizeof(int),
              "cellCount");
      growOne(dcellOff, static_cast<std::size_t>(nOff) * sizeof(int),
              "cellOff");
    }
  }

  void ensureStencil(int maxReach) {
    if (cachedReach == maxReach && dsdx.p != nullptr) {
      return;
    }
    std::vector<int> dx, dy, dz, off;
    buildStencil(maxReach, dx, dy, dz, off);
    const int nSt = static_cast<int>(dx.size());
    if (nSt > capSt) {
      capSt = nSt;
      growOne(dsdx, static_cast<std::size_t>(capSt) * sizeof(int), "sdx");
      growOne(dsdy, static_cast<std::size_t>(capSt) * sizeof(int), "sdy");
      growOne(dsdz, static_cast<std::size_t>(capSt) * sizeof(int), "sdz");
    }
    growOne(dreachOff, off.size() * sizeof(int), "reachOff");
    auto &rt = CUDART::instance();
    checkCuda(rt.cudaMemcpy(dsdx.p, dx.data(), dx.size() * sizeof(int),
                            cudaMemcpyHostToDevice),
              "HtoD sdx");
    checkCuda(rt.cudaMemcpy(dsdy.p, dy.data(), dy.size() * sizeof(int),
                            cudaMemcpyHostToDevice),
              "HtoD sdy");
    checkCuda(rt.cudaMemcpy(dsdz.p, dz.data(), dz.size() * sizeof(int),
                            cudaMemcpyHostToDevice),
              "HtoD sdz");
    checkCuda(rt.cudaMemcpy(dreachOff.p, off.data(), off.size() * sizeof(int),
                            cudaMemcpyHostToDevice),
              "HtoD reachOff");
    cachedReach = maxReach;
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
  knearest_into_many(xyz, n, 1, cell, k, out, out_len, mask, cell_hint);
}

void Workspace::knearest_into_many(const double *xyz, std::size_t n,
                                   std::size_t nFrames, const Cell &cell,
                                   std::size_t k, int *out, std::size_t out_len,
                                   const int *mask, double cell_hint) {
  if (!available()) {
    throw Error("CUDA driver or nvrtc not loaded");
  }
  if (n == 0 || nFrames == 0) {
    throw Error("no points");
  }
  if (k == 0) {
    throw Error("k must be at least 1");
  }
  if (k > static_cast<std::size_t>(kMaxK)) {
    throw Error("device k-nearest supports k <= 16");
  }
  if (out_len != n * k * nFrames) {
    throw Error("out buffer length must be nFrames * n * k");
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
  double lx = cell.a[0];
  double ly = cell.b[1];
  double lz = cell.c[2];
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
  int nI = static_cast<int>(n);
  int nF = static_cast<int>(nFrames);
  int kI = static_cast<int>(k);
  const int nTot = nI * nF;
  const int nCtot = nC * nF;
  const int nOff = (nC + 1) * nF;

  impl_->ensure(nTot, nCtot, nOff);
  impl_->ensureStencil(maxReach);
  auto &rt = CUDART::instance();
  checkCuda(rt.cudaMemset(impl_->dcellCount.p, 0,
                          static_cast<std::size_t>(nCtot) * sizeof(int)),
            "zero cells");

  auto &factory = KernelFactory::instance(0);
  const std::vector<std::string> opt{"-std=c++17"};
  auto *kBin = factory.create("bin_atoms", kKernels, "linkcell.cu", opt);
  auto *kPref = factory.create("prefix_cells", kKernels, "linkcell.cu", opt);
  auto *kScat = factory.create("scatter_atoms", kKernels, "linkcell.cu", opt);
  auto *kNear = factory.create("knearest_shells", kKernels, "linkcell.cu", opt);

  const int block = 128;
  const int grid = (nTot + block - 1) / block;
  const int tppGrid = (nTot * 8 + block - 1) / block;
  const int ppb = block / 8;
  const std::size_t sh = static_cast<std::size_t>(ppb) * 8 * 16 * sizeof(double) +
                         static_cast<std::size_t>(ppb) * 8 * 16 * sizeof(int) +
                         static_cast<std::size_t>(ppb) * 8 * sizeof(int) +
                         static_cast<std::size_t>(ppb) * sizeof(int);
  void *xyzP = const_cast<double *>(xyz);
  void *maskV = const_cast<int *>(mask);
  void *outP = out;
  int nCv = nC;
  int maxR = maxReach;
  double cmin = cellMin;
  double ox = cell.origin[0];
  double oy = cell.origin[1];
  double oz = cell.origin[2];
  auto launchArgs = [](void **raw, std::size_t n) {
    return std::vector<void *>(raw, raw + n);
  };
  {
    void *raw[] = {&xyzP,
                   &maskV,
                   &nI,
                   &nF,
                   &lx,
                   &ly,
                   &lz,
                   &ox,
                   &oy,
                   &oz,
                   &nx,
                   &ny,
                   &nz,
                   &nCv,
                   &impl_->dcellOf.p,
                   &impl_->dcellCount.p,
                   &impl_->dfolded.p};
    auto a = launchArgs(raw, sizeof(raw) / sizeof(raw[0]));
    kBin->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
  }
  {
    void *raw[] = {&impl_->dcellCount.p, &impl_->dcellOff.p, &nCv, &nF};
    auto a = launchArgs(raw, sizeof(raw) / sizeof(raw[0]));
    kPref->launch(dim3(nF), dim3(128), 128 * sizeof(int), nullptr, a, true);
  }
  {
    void *raw[] = {&impl_->dfolded.p,    &impl_->dcellOf.p, &impl_->dcellCount.p,
                   &impl_->dcellOff.p,   &impl_->dorder.p,  &impl_->dsorted.p,
                   &impl_->dhome.p,      &nI,               &nF,
                   &nCv};
    auto a = launchArgs(raw, sizeof(raw) / sizeof(raw[0]));
    kScat->launch(dim3(grid), dim3(block), 0, nullptr, a, true);
  }
  {
    void *raw[] = {&impl_->dsorted.p,   &impl_->dcellOf.p, &impl_->dcellOff.p,
                   &impl_->dorder.p,    &impl_->dhome.p,   &impl_->dsdx.p,
                   &impl_->dsdy.p,      &impl_->dsdz.p,    &impl_->dreachOff.p,
                   &nI,                 &nF,               &nCv,
                   &nx,                 &ny,               &nz,
                   &lx,                 &ly,               &lz,
                   &cmin,               &kI,               &maxR,
                   &outP};
    auto a = launchArgs(raw, sizeof(raw) / sizeof(raw[0]));
    kNear->launch(dim3(tppGrid), dim3(block), sh, nullptr, a, true);
  }
}

} // namespace gpu
} // namespace linkcell

#endif
