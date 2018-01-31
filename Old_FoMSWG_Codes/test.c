#include "fishercalc.c"

int main(int argc, char **argv) {
  FILE *fp;
  COSMOPARAM p;
  double **F, **F_small;

  F = dmatrix(0, NPARAMTOT-1, 0, NPARAMTOT-1);
  F_small = dmatrix(0, NWPARAM-1, 0, NWPARAM-1);
  set_default_cosmology(&p);

  get_firas_fisher_matrix(F, 0x1);

#ifdef WL_VERSION
  get_wl_fisher_matrix(F, argv[1], &p, 0x18); // stat only
#endif
#ifdef BAO_VERSION
  get_bao_fisher_matrix2(F, argv[1], &p, 0x4);
#endif

  write_swg_fisher(F, &p,argv[2]);
  return(0);

}

