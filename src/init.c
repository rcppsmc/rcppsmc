#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP blockpfGaussianOpt(SEXP, SEXP, SEXP);
extern SEXP pfLineartBS(SEXP, SEXP, SEXP, SEXP);
extern SEXP pfNonlinBS(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"blockpfGaussianOpt", (DL_FUNC) &blockpfGaussianOpt, 3},
    {"pfLineartBS",        (DL_FUNC) &pfLineartBS,        4},
    {"pfNonlinBS",         (DL_FUNC) &pfNonlinBS,         2},
    {NULL, NULL, 0}
};

void R_init_RcppSMC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
