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

static PyMethodDef methods[] = {
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
