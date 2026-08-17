#include "linkcell_gpu.h"
#include "linkcell_gpu.hpp"

#if defined(LINKCELL_HAS_GPULITE)
#include <gpulite/gpulite.hpp>
#endif

#include <new>
#include <string>

namespace {

thread_local std::string g_err;

void set_err(const char *msg) { g_err = msg ? msg : ""; }

} // namespace

extern "C" {

int lc_gpu_available(void) {
#if defined(LINKCELL_HAS_GPULITE)
  return linkcell::gpu::available() ? 1 : 0;
#else
  return 0;
#endif
}

const char *lc_gpu_last_error(void) {
  return g_err.empty() ? nullptr : g_err.c_str();
}

struct lc_gpu_workspace {
#if defined(LINKCELL_HAS_GPULITE)
  linkcell::gpu::Workspace ws;
#endif
};

lc_gpu_workspace *lc_gpu_workspace_new(void) {
#if defined(LINKCELL_HAS_GPULITE)
  if (!linkcell::gpu::available()) {
    set_err("CUDA driver or nvrtc not loaded");
    return nullptr;
  }
  try {
    return new lc_gpu_workspace();
  } catch (const std::exception &e) {
    set_err(e.what());
    return nullptr;
  }
#else
  set_err("built without gpulite");
  return nullptr;
#endif
}

void lc_gpu_workspace_free(lc_gpu_workspace *ws) { delete ws; }

#if defined(LINKCELL_HAS_GPULITE)
static linkcell::Cell cell_from_c(const lc_cell *simbox) {
  return linkcell::Cell(*simbox);
}
#endif

int lc_gpu_knearest(lc_gpu_workspace *ws, const double *xyz, size_t n,
                    const struct lc_cell *simbox, size_t k, const int *mask,
                    double cell_hint, int *out_nn) {
  return lc_gpu_knearest_many(ws, xyz, n, 1, simbox, k, mask, cell_hint,
                              out_nn, 1);
}

int lc_gpu_knearest_many(lc_gpu_workspace *ws, const double *xyz, size_t n,
                         size_t n_frames, const struct lc_cell *simbox,
                         size_t k, const int *mask, double cell_hint,
                         int *out_nn, int wait) {
#if defined(LINKCELL_HAS_GPULITE)
  if (ws == nullptr || xyz == nullptr || simbox == nullptr ||
      out_nn == nullptr) {
    set_err("null pointer");
    return 1;
  }
  try {
    auto cell = cell_from_c(simbox);
    const size_t out_len = n * k * n_frames;
    ws->ws.knearest_into_many(xyz, n, n_frames, cell, k, out_nn, out_len, mask,
                              cell_hint, wait != 0, nullptr);
    set_err(nullptr);
    return 0;
  } catch (const std::exception &e) {
    set_err(e.what());
    return 1;
  }
#else
  (void)ws;
  (void)xyz;
  (void)n;
  (void)n_frames;
  (void)simbox;
  (void)k;
  (void)mask;
  (void)cell_hint;
  (void)out_nn;
  (void)wait;
  set_err("built without gpulite");
  return 1;
#endif
}

void *lc_gpu_queue(lc_gpu_workspace *ws) {
#if defined(LINKCELL_HAS_GPULITE)
  return ws == nullptr ? nullptr : ws->ws.queue();
#else
  (void)ws;
  return nullptr;
#endif
}

void lc_gpu_wait(lc_gpu_workspace *ws) {
#if defined(LINKCELL_HAS_GPULITE)
  if (ws != nullptr) {
    ws->ws.wait();
  }
#else
  (void)ws;
#endif
}

int lc_gpu_alloc(void **ptr, size_t bytes) {
#if defined(LINKCELL_HAS_GPULITE)
  if (ptr == nullptr) {
    set_err("null pointer");
    return 1;
  }
  if (!linkcell::gpu::available()) {
    set_err("CUDA driver or nvrtc not loaded");
    return 1;
  }
  try {
    auto &rt = gpulite::CUDART::instance();
    auto st = rt.cudaMalloc(ptr, bytes);
    if (st != cudaSuccess) {
      set_err(rt.cudaGetErrorString(st));
      return 1;
    }
    set_err(nullptr);
    return 0;
  } catch (const std::exception &e) {
    set_err(e.what());
    return 1;
  }
#else
  (void)ptr;
  (void)bytes;
  set_err("built without gpulite");
  return 1;
#endif
}

void lc_gpu_free(void *ptr) {
#if defined(LINKCELL_HAS_GPULITE)
  if (ptr != nullptr && gpulite::CUDART::loaded()) {
    gpulite::CUDART::instance().cudaFree(ptr);
  }
#else
  (void)ptr;
#endif
}

int lc_gpu_fill_i32(void *ptr, int value, size_t n) {
#if defined(LINKCELL_HAS_GPULITE)
  if (ptr == nullptr && n != 0) {
    set_err("null pointer");
    return 1;
  }
  if (!linkcell::gpu::available()) {
    set_err("CUDA driver or nvrtc not loaded");
    return 1;
  }
  try {
    auto &rt = gpulite::CUDART::instance();
    if (value == 0 || value == -1) {
      const unsigned char byte = value == 0 ? 0 : 0xFF;
      auto st = rt.cudaMemset(ptr, byte, n * sizeof(int));
      if (st != cudaSuccess) {
        set_err(rt.cudaGetErrorString(st));
        return 1;
      }
      set_err(nullptr);
      return 0;
    }
    set_err("lc_gpu_fill_i32 only supports 0 and -1");
    return 1;
  } catch (const std::exception &e) {
    set_err(e.what());
    return 1;
  }
#else
  (void)ptr;
  (void)value;
  (void)n;
  set_err("built without gpulite");
  return 1;
#endif
}

} // extern "C"
