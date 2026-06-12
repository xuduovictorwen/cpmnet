#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Fortran calls */
extern void F77_NAME(cpm_predict)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cpm_predict_cdf)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cpm_predict_psr)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(cpm_predict_quantile)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ormll_beta_ccd_enet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(ormll_path_enet)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
  {"cpm_predict",         (DL_FUNC) &F77_NAME(cpm_predict),         12},
  {"cpm_predict_cdf",     (DL_FUNC) &F77_NAME(cpm_predict_cdf),     12},
  {"cpm_predict_psr",     (DL_FUNC) &F77_NAME(cpm_predict_psr),     10},
  {"cpm_predict_quantile", (DL_FUNC) &F77_NAME(cpm_predict_quantile), 13},
  {"ormll_beta_ccd_enet", (DL_FUNC) &F77_NAME(ormll_beta_ccd_enet), 15},
  {"ormll_path_enet",     (DL_FUNC) &F77_NAME(ormll_path_enet),     21},
  {NULL, NULL, 0}
};

void R_init_cpmnet(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
