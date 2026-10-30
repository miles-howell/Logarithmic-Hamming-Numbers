/* Exact rank check: exact_count(I, J, K) = number of (i, j, k) >= 0 with
 * 2^i 3^j 5^k <= 2^I 3^J 5^K, i.e. the 1-indexed rank of that Hamming number.
 *
 * Independent of the float logic in hamming.py: for every (j, k) cell it decides
 *     i <= t,  t = I + (J-j)*log2(3) + (K-k)*log2(5)
 * with 80-bit long double. Any cell whose t is within MARGIN of an integer is redone
 * in 128-bit quad precision, and if even that cannot decide it the call raises
 * ValueError rather than guess. (By unique factorization t is never an exact integer
 * except on the cell (J, K) itself, which is counted directly.)
 *
 * Needs x86-64 (80-bit long double) and libquadmath. Parallelized with OpenMP over k.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <float.h>
#include <math.h>
#include <quadmath.h>

#if LDBL_MANT_DIG < 64
#error "_hamming_exact needs an 80-bit long double (x86-64)"
#endif

#define MARGIN  1e-9L
#define QMARGIN 1e-20Q

static __float128 Q3, Q5;   /* log2(3), log2(5) */
static long double R3, R5;

static PyObject *
exact_count(PyObject *self, PyObject *args)
{
    long I, J, K;
    if (!PyArg_ParseTuple(args, "lll", &I, &J, &K))
        return NULL;
    if (I < 0 || J < 0 || K < 0) {
        PyErr_SetString(PyExc_ValueError, "exponents must be non-negative");
        return NULL;
    }

    const double ln3 = log(3.0), ln5 = log(5.0);
    const double L = I * log(2.0) + J * ln3 + K * ln5;
    const long kmax = (long)(L / ln5) + 2;
    unsigned long long total = 0;
    int ambiguous = 0;

    Py_BEGIN_ALLOW_THREADS
    #pragma omp parallel for schedule(dynamic, 16) reduction(+:total)
    for (long k = 0; k <= kmax; k++) {
        const double rem = L - k * ln5;
        const long jmax = rem > 0 ? (long)(rem / ln3) + 2 : 1;
        const long double base = (long double)I + (long double)(K - k) * R5;
        for (long j = 0; j <= jmax; j++) {
            long double t = base + (long double)(J - j) * R3;
            long double fl = floorl(t);
            long double frac = t - fl;
            if (frac < MARGIN || frac > 1.0L - MARGIN) {
                if (j == J && k == K) {          /* the target itself: i <= I exactly */
                    total += (unsigned long long)I + 1;
                    continue;
                }
                __float128 tq = (__float128)I + (__float128)(J - j) * Q3
                                              + (__float128)(K - k) * Q5;
                __float128 fq = floorq(tq);
                __float128 fr = tq - fq;
                if (fr < QMARGIN || fr > 1.0Q - QMARGIN) {
                    #pragma omp atomic write
                    ambiguous = 1;
                    continue;
                }
                fl = (long double)fq;
            }
            if (fl >= 0)
                total += (unsigned long long)fl + 1;
        }
    }
    Py_END_ALLOW_THREADS

    if (ambiguous) {
        PyErr_SetString(PyExc_ValueError,
                        "a cell was too close to an integer to decide even in quad precision");
        return NULL;
    }
    return PyLong_FromUnsignedLongLong(total);
}

static PyMethodDef methods[] = {
    {"exact_count", exact_count, METH_VARARGS,
     "exact_count(I, J, K): exact rank of 2**I * 3**J * 5**K among Hamming numbers."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT, "_hamming_exact", NULL, -1, methods
};

PyMODINIT_FUNC
PyInit__hamming_exact(void)
{
    Q3 = log2q(3.0Q);
    Q5 = log2q(5.0Q);
    R3 = (long double)Q3;
    R5 = (long double)Q5;
    return PyModule_Create(&module);
}
