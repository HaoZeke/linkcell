#ifndef LINKCELL_H
#define LINKCELL_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * linkcell C API.
 *
 * Memory ownership
 *   - lc_knearest writes into caller-owned out_nn (length n*k).
 *   - lc_last_error and lc_version return process-static pointers. Do not free.
 *
 * Thread safety
 *   - Distinct searches may run concurrently.
 *   - lc_last_error is a process-wide slot; do not read it from two
 *     threads after concurrent failures.
 */

typedef struct lc_box {
  double lx;
  double ly;
  double lz;
  double xlo;
  double ylo;
  double zlo;
} lc_box;

const char *lc_version(void);
const char *lc_last_error(void);

int lc_knearest(const double *xyz, int n, const lc_box *simbox, int k,
                const int *mask, double cell_hint, int *out_nn);

#ifdef __cplusplus
}
#endif

#endif
