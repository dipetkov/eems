
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP rEEMSplots_tiles2contours(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rEEMSplots_tiles2contours_standardize(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"rEEMSplots_tiles2contours",             (DL_FUNC) &rEEMSplots_tiles2contours,             5},
    {"rEEMSplots_tiles2contours_standardize", (DL_FUNC) &rEEMSplots_tiles2contours_standardize, 5},
    {NULL, NULL, 0}
};

void R_init_rEEMSplots(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
