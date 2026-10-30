/* C implementation of hamming.count(B).
 *
 * Mirrors the pure-Python loop operation for operation (same double
 * arithmetic, same repeated subtraction, same 1e-11 guard) so both
 * backends return identical counts. Build with `python setup.py build_ext --inplace`.
 * Do not compile with -ffast-math or FMA contraction: that could change rounding.
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <math.h>
#include <stdint.h>

static PyObject *
hamming_count(PyObject *self, PyObject *arg)
{
    double B = PyFloat_AsDouble(arg);
    if (B == -1.0 && PyErr_Occurred())
        return NULL;
    if (B < 0)
        return PyLong_FromLong(0);

    const double a = log(2.0), b = log(3.0), c = log(5.0);
    const double inv_a = 1.0 / a;
    uint64_t total = 0;

    Py_BEGIN_ALLOW_THREADS
    double rem_k = B;
    while (rem_k >= 0) {
        double rem_j = rem_k;
        while (rem_j >= 0) {
            total += (uint64_t)(rem_j * inv_a + 1e-11) + 1;
            rem_j -= b;
        }
        rem_k -= c;
    }
    Py_END_ALLOW_THREADS

    return PyLong_FromUnsignedLongLong(total);
}

/* band_terms(lo, hi): the (j, k, i) exponents of every number 2^i 3^j 5^k with
 * lo < ln(value) <= hi. Mirrors the pure-Python extraction loop cell for cell. */
static PyObject *
hamming_band_terms(PyObject *self, PyObject *args)
{
    double lo, hi;
    if (!PyArg_ParseTuple(args, "dd", &lo, &hi))
        return NULL;

    const double b = log(3.0), c = log(5.0);
    const double inv_a = 1.0 / log(2.0);
    PyObject *out = PyList_New(0);
    if (out == NULL)
        return NULL;

    long k = 0;
    double rem_k_hi = hi, rem_k_lo = lo;
    while (rem_k_hi >= 0) {
        long j = 0;
        double rem_j_hi = rem_k_hi, rem_j_lo = rem_k_lo;
        while (rem_j_hi >= 0) {
            long top_hi = (long)(rem_j_hi * inv_a + 1e-11);
            long top_lo = rem_j_lo >= 0 ? (long)(rem_j_lo * inv_a + 1e-11) : -1;
            for (long i = top_lo + 1; i <= top_hi; i++) {
                PyObject *t = Py_BuildValue("(lll)", j, k, i);
                if (t == NULL || PyList_Append(out, t) < 0) {
                    Py_XDECREF(t);
                    Py_DECREF(out);
                    return NULL;
                }
                Py_DECREF(t);
            }
            j++;
            rem_j_hi -= b;
            rem_j_lo -= b;
        }
        k++;
        rem_k_hi -= c;
        rem_k_lo -= c;
    }
    return out;
}

static PyMethodDef methods[] = {
    {"band_terms", hamming_band_terms, METH_VARARGS,
     "band_terms(lo, hi): (j, k, i) exponents of Hamming numbers with lo < ln <= hi."},
    {"count", hamming_count, METH_O,
     "count(B): number of Hamming numbers whose natural log is <= B."},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
    PyModuleDef_HEAD_INIT, "_hamming_count", NULL, -1, methods
};

PyMODINIT_FUNC
PyInit__hamming_count(void)
{
    return PyModule_Create(&module);
}
