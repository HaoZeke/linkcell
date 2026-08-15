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

typedef struct lc_cell {
  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
  double ox, oy, oz;
} lc_cell;

static inline lc_cell lc_cell_ortho(double lx, double ly, double lz) {
  lc_cell c;
  c.ax = lx;
  c.ay = 0.0;
  c.az = 0.0;
  c.bx = 0.0;
  c.by = ly;
  c.bz = 0.0;
  c.cx = 0.0;
  c.cy = 0.0;
  c.cz = lz;
  c.ox = 0.0;
  c.oy = 0.0;
  c.oz = 0.0;
  return c;
}

const char *lc_version(void);
const char *lc_last_error(void);

int lc_knearest(const double *xyz, int n, const lc_cell *simbox, int k,
                const int *mask, double cell_hint, int *out_nn);

#ifdef __cplusplus
}
#endif

#endif
