/* C implementation of hamming.count(B) and the band walk.
 *
 * Why long double: Hamming numbers get denser in log space as n grows (spacing is
 * about 2.45 / B^2), and from n ~ 1e13 it is smaller than the rounding error of
 * double arithmetic, which then misclassifies neighbouring numbers. Here every cell's
 * remainder is computed directly as B - k*ln5 - j*ln3 in 80-bit x87 long double
 * (error ~ 1e-13 in log space, about 400x below the spacing at n = 1e16), with no
 * accumulated drift from repeated subtraction.
 *
 * count() and band_terms() share cell_top(), so they classify every cell identically;
 * the result relies on that consistency. A 1e-11 guard (in log2 units) is kept from the
 * pure-Python version; it shifts the threshold uniformly, which is harmless.
 *
 * The pure-Python fallback in hamming.py still uses double and is only reliable up to
 * about n = 1e12. Build with `python setup.py build_ext --inplace`. Do not compile with
 * -ffast-math. Loops are parallelized with OpenMP (set OMP_NUM_THREADS to limit).
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#if LDBL_MANT_DIG < 64
#error "_hamming_count needs an 80-bit long double (x86-64); use the pure-Python fallback"
#endif

static const long double LN2 = 0.69314718055994530941723212145817656807550013436026L;
static const long double LN3 = 1.09861228866810969139524523692252570464749055782275L;
static const long double LN5 = 1.60943791243410037460075933322618763952560135426852L;
static const long double INV_LN2 = 1.44269504088896340735992468100189213742664595415299L;
static const long double GUARD = 1e-11L;

/* Largest i with i*ln2 <= rem (+ guard); rem must be >= 0. */
static inline long
cell_top(long double rem)
{
    return (long)(rem * INV_LN2 + GUARD);
}

/* Cells (j, k) with k*ln5 + j*ln3 <= B; calls body with `rem` = B - k*ln5 - j*ln3. */
#define PARALLEL_MIN_ROWS 2048

static PyObject *
hamming_count(PyObject *self, PyObject *arg)
{
    double Bd = PyFloat_AsDouble(arg);
    if (Bd == -1.0 && PyErr_Occurred())
        return NULL;
    if (Bd < 0)
        return PyLong_FromLong(0);

    const long double B = Bd;
    const long nk = (long)(B / LN5) + 2;
    unsigned long long total = 0;

    Py_BEGIN_ALLOW_THREADS
    #pragma omp parallel for schedule(dynamic, 16) reduction(+:total) if(nk > PARALLEL_MIN_ROWS)
    for (long k = 0; k < nk; k++) {
        const long double remk = B - (long double)k * LN5;
        if (remk < 0)
            continue;
        for (long j = 0;; j++) {
            const long double rem = remk - (long double)j * LN3;
            if (rem < 0)
                break;
            total += (unsigned long long)cell_top(rem) + 1;
        }
    }
    Py_END_ALLOW_THREADS

    return PyLong_FromUnsignedLongLong(total);
}

typedef struct { long j, k, i; } Term;

/* band_terms(lo, hi): the (j, k, i) exponents of every number 2^i 3^j 5^k that count(hi)
 * includes and count(lo) does not. */
static PyObject *
hamming_band_terms(PyObject *self, PyObject *args)
{
    double lod, hid;
    if (!PyArg_ParseTuple(args, "dd", &lod, &hid))
        return NULL;

    const long double lo = lod, hi = hid;
    Term *terms = NULL;
    size_t nterms = 0, cap = 0;
    int oom = 0;

    if (hi >= 0) {
        const long nk = (long)(hi / LN5) + 2;
        Py_BEGIN_ALLOW_THREADS
        #pragma omp parallel for schedule(dynamic, 16) if(nk > PARALLEL_MIN_ROWS)
        for (long k = 0; k < nk; k++) {
            const long double remk_hi = hi - (long double)k * LN5;
            const long double remk_lo = lo - (long double)k * LN5;
            if (remk_hi < 0)
                continue;
            for (long j = 0;; j++) {
                const long double rem_hi = remk_hi - (long double)j * LN3;
                if (rem_hi < 0)
                    break;
                const long double rem_lo = remk_lo - (long double)j * LN3;
                const long top_hi = cell_top(rem_hi);
                const long top_lo = rem_lo >= 0 ? cell_top(rem_lo) : -1;
                for (long i = top_lo + 1; i <= top_hi; i++) {
                    #pragma omp critical(hamming_terms)
                    {
                        if (nterms == cap) {
                            size_t ncap = cap ? cap * 2 : 1024;
                            Term *grown = realloc(terms, ncap * sizeof(Term));
                            if (grown) { terms = grown; cap = ncap; } else { oom = 1; }
                        }
                        if (!oom) {
                            terms[nterms].j = j;
                            terms[nterms].k = k;
                            terms[nterms].i = i;
                            nterms++;
                        }
                    }
                }
            }
        }
        Py_END_ALLOW_THREADS
    }

    if (oom) {
        free(terms);
        return PyErr_NoMemory();
    }
    PyObject *out = PyList_New(0);
    for (size_t t = 0; out != NULL && t < nterms; t++) {
        PyObject *tup = Py_BuildValue("(lll)", terms[t].j, terms[t].k, terms[t].i);
        if (tup == NULL || PyList_Append(out, tup) < 0) {
            Py_XDECREF(tup);
            Py_CLEAR(out);
            break;
        }
        Py_DECREF(tup);
    }
    free(terms);
    return out;
}

static PyMethodDef methods[] = {
    {"count", hamming_count, METH_O,
     "count(B): number of Hamming numbers whose natural log is <= B."},
    {"band_terms", hamming_band_terms, METH_VARARGS,
     "band_terms(lo, hi): (j, k, i) exponents of Hamming numbers with lo < ln <= hi."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT, "_hamming_count", NULL, -1, methods
};

PyMODINIT_FUNC
PyInit__hamming_count(void)
{
    (void)LN2;
    return PyModule_Create(&module);
}
