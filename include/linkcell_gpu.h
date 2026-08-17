#ifndef LINKCELL_GPU_H
#define LINKCELL_GPU_H

#pragma once

#include "linkcell.h"

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Device-resident k-nearest. Same packed n*k int layout as lc_knearest.
 * xyz and out_nn are CUDA device pointers. simbox is a host lc_cell.
 * lc_gpu_knearest_many_dcell takes a device-packed cell (3, 9, or 12
 * doubles) and inverts H on device.
 * Built when meson with_gpulite is on; otherwise lc_gpu_available is 0
 * and the other entry points fail.
 *
 * Distinct workspaces may run concurrently. Each workspace owns one
 * CUDA stream (lc_gpu_queue). Call lc_gpu_wait before reading device
 * results if wait was 0 on the search.
 */

typedef struct lc_gpu_workspace lc_gpu_workspace;

int lc_gpu_available(void);
const char *lc_gpu_last_error(void);

lc_gpu_workspace *lc_gpu_workspace_new(void);
void lc_gpu_workspace_free(lc_gpu_workspace *ws);

int lc_gpu_knearest(lc_gpu_workspace *ws, const double *xyz, size_t n,
                    const struct lc_cell *simbox, size_t k, const int *mask,
                    double cell_hint, int *out_nn, double *out_d2);

int lc_gpu_knearest_many(lc_gpu_workspace *ws, const double *xyz, size_t n,
                         size_t n_frames, const struct lc_cell *simbox,
                         size_t k, const int *mask, double cell_hint,
                         int *out_nn, double *out_d2, int wait);

int lc_gpu_knearest_many_dcell(lc_gpu_workspace *ws, const double *xyz,
                               size_t n, size_t n_frames, const double *cell,
                               int cell_n, size_t k, const int *mask,
                               double cell_hint, int *out_nn, double *out_d2,
                               int wait);

void *lc_gpu_queue(lc_gpu_workspace *ws);
void lc_gpu_wait(lc_gpu_workspace *ws);

int lc_gpu_alloc(void **ptr, size_t bytes);
void lc_gpu_free(void *ptr);
int lc_gpu_fill_i32(void *ptr, int value, size_t n);
int lc_gpu_memcpy(void *dst, const void *src, size_t bytes, int kind);

#ifdef __cplusplus
}
#endif

#endif
