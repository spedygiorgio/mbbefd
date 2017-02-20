#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP mbbefd_g2a(SEXP, SEXP);
extern SEXP mbbefd_rmbbefdC(SEXP, SEXP, SEXP);
extern SEXP mbbefd_rMBBEFDC(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"mbbefd_g2a",      (DL_FUNC) &mbbefd_g2a,      2},
    {"mbbefd_rmbbefdC", (DL_FUNC) &mbbefd_rmbbefdC, 3},
    {"mbbefd_rMBBEFDC", (DL_FUNC) &mbbefd_rMBBEFDC, 3},
    {NULL, NULL, 0}
};

void R_init_mbbefd(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
