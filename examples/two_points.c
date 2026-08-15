#include "linkcell.h"

#include <stdio.h>
#include <stdlib.h>

/* Neighbours of source i sit at out[i * k + 0 .. k). Unused slots are -1. */
static int run(const char *label, const double *xyz, size_t n,
               const lc_cell *box, size_t k, int *out) {
  if (lc_knearest(xyz, n, box, k, NULL, 0.0, out) != 0) {
    fprintf(stderr, "%s\n", lc_last_error());
    return 1;
  }
  printf("%s\n", label);
  for (size_t i = 0; i < n; i++) {
    printf("%zu ->", i);
    for (size_t t = 0; t < k; t++) {
      printf(" %d", out[i * k + t]);
    }
    printf("\n");
  }
  return 0;
}

int main(void) {
  const size_t n = 2;
  const size_t k = 1;
  const double xyz[] = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  const lc_cell ortho = lc_cell_ortho(10.0, 10.0, 10.0);
  const lc_cell sheared = {
      .ax = 10.0,
      .ay = 0.0,
      .az = 0.0,
      .bx = 5.0,
      .by = 8.66,
      .bz = 0.0,
      .cx = 0.0,
      .cy = 0.0,
      .cz = 10.0,
      .ox = 0.0,
      .oy = 0.0,
      .oz = 0.0,
  };

  int *out = malloc(n * k * sizeof *out);
  if (out == NULL) {
    return 1;
  }
  if (run("ortho", xyz, n, &ortho, k, out) != 0 ||
      run("sheared", xyz, n, &sheared, k, out) != 0) {
    free(out);
    return 1;
  }
  free(out);
  return 0;
}
