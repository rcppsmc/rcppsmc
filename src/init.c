#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP blockpfGaussianOpt_impl(SEXP, SEXP, SEXP);
extern SEXP pfLineartBS_impl(SEXP, SEXP, SEXP, SEXP);
extern SEXP pfNonlinBS_impl(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"blockpfGaussianOpt_impl", (DL_FUNC) &blockpfGaussianOpt_impl, 3},
    {"pfLineartBS_impl",        (DL_FUNC) &pfLineartBS_impl,        4},
    {"pfNonlinBS_impl",         (DL_FUNC) &pfNonlinBS_impl,         2},
    {NULL, NULL, 0}
};

void R_init_RcppSMC(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
