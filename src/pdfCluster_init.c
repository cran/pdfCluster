#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void apply_valley_measure(void *, void *, void *, void *);
extern void c_kepdfN(void *, void *, void *, void *, void *, void *, void *, void *);
extern void c_kepdft7(void *, void *, void *, void *, void *, void *, void *, void *);
extern void valley_measure(void *, void *, void *);

/* .Call calls */
extern SEXP g_components(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"apply_valley_measure", (DL_FUNC) &apply_valley_measure, 4},
    {"c_kepdfN",             (DL_FUNC) &c_kepdfN,             8},
    {"c_kepdft7",            (DL_FUNC) &c_kepdft7,            8},
    {"valley_measure",       (DL_FUNC) &valley_measure,       3},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"g_components", (DL_FUNC) &g_components, 2},
    {NULL, NULL, 0}
};

void R_init_pdfCluster(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
