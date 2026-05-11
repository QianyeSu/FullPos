#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_20_API_VERSION
#include <numpy/arrayobject.h>

#include <limits.h>

extern void fullpos_pressure_ppq_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *levels,
    double *output,
    int *ierr);

extern void fullpos_pressure_ppuv_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *u_values,
    const double *v_values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *levels,
    double *u_output,
    double *v_output,
    int *ierr);

extern void fullpos_pressure_ppt_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *levels,
    double *output,
    int *ierr);

extern void fullpos_column_pressure_ppq_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *target_pressures,
    double *output,
    int *ierr);

extern void fullpos_column_pressure_ppuv_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *u_values,
    const double *v_values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *target_pressures,
    double *u_output,
    double *v_output,
    int *ierr);

extern void fullpos_column_pressure_ppt_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *target_pressures,
    double *output,
    int *ierr);

extern void fullpos_apache_column_pressure_tuvq_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *temperature,
    const double *u_values,
    const double *v_values,
    const double *q_values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *target_pressures,
    const double *target_surface_pressure,
    const int *enable_lescale,
    double *t_output,
    double *u_output,
    double *v_output,
    double *q_output,
    int *ierr);

extern void fullpos_hybrid_pressure_ppq_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *target_ak,
    const double *target_bk,
    double *output,
    int *ierr);

extern void fullpos_hybrid_pressure_ppuv_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *u_values,
    const double *v_values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *target_ak,
    const double *target_bk,
    double *u_output,
    double *v_output,
    int *ierr);

extern void fullpos_hybrid_pressure_ppt_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *target_ak,
    const double *target_bk,
    double *output,
    int *ierr);

extern void fullpos_theta_pressures_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *temperature,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *theta_levels,
    double *output,
    int *ierr);

extern void fullpos_eta_pressures_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *eta_levels,
    double *output,
    int *ierr);

extern void fullpos_temperature_pressures_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *temperature,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *temperature_levels,
    double *output,
    int *ierr);

extern void fullpos_height_above_orography_pressures_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *temperature,
    const double *specific_humidity,
    const int *has_specific_humidity,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *orography_geopotential,
    const double *height_levels,
    double *output,
    int *ierr);

extern void fullpos_height_above_sea_pressures_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *temperature,
    const double *specific_humidity,
    const int *has_specific_humidity,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *orography_geopotential,
    const double *height_levels,
    double *output,
    int *ierr);

extern void fullpos_potential_vorticity_pressures_c(
    const int *ncol,
    const int *nlev,
    const int *nout,
    const double *pv_values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *coriolis,
    const double *pv_levels,
    double *output,
    int *ierr);

extern void fullpos_diagnose_potential_vorticity_c(
    const int *ncol,
    const int *nlev,
    const double *u_values,
    const double *v_values,
    const double *temperature,
    const double *relative_vorticity,
    const double *temperature_meridional_gradient,
    const double *temperature_zonal_gradient,
    const double *surface_pressure_meridional_gradient,
    const double *surface_pressure_zonal_gradient,
    const double *kappa_values,
    const double *ak,
    const double *bk,
    const double *ps,
    const double *coriolis,
    double *potential_vorticity,
    double *potential_temperature,
    int *ierr);

extern void fullpos_gprcp_kappa_c(
    const int *ncol,
    const int *nlev,
    const double *q_values,
    double *kappa_values,
    int *ierr);

static PyArrayObject *as_c_array(PyObject *obj, int typenum, int min_ndim, int max_ndim)
{
    PyArrayObject *arr = (PyArrayObject *)PyArray_FROM_OTF(obj, typenum, NPY_ARRAY_IN_ARRAY);
    if (!arr) {
        return NULL;
    }
    int ndim = PyArray_NDIM(arr);
    if (ndim < min_ndim || ndim > max_ndim) {
        PyErr_Format(PyExc_ValueError, "expected array with %d..%d dimensions, got %d", min_ndim, max_ndim, ndim);
        Py_DECREF(arr);
        return NULL;
    }
    return arr;
}

static PyArrayObject *as_fortran_array(PyObject *obj, int typenum, int min_ndim, int max_ndim)
{
    PyArrayObject *arr = (PyArrayObject *)PyArray_FROM_OTF(obj, typenum, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_ALIGNED);
    if (!arr) {
        return NULL;
    }
    int ndim = PyArray_NDIM(arr);
    if (ndim < min_ndim || ndim > max_ndim) {
        PyErr_Format(PyExc_ValueError, "expected array with %d..%d dimensions, got %d", min_ndim, max_ndim, ndim);
        Py_DECREF(arr);
        return NULL;
    }
    return arr;
}

static int check_common(PyArrayObject *values, PyArrayObject *ak, PyArrayObject *bk, PyArrayObject *ps, PyArrayObject *levels)
{
    if (PyArray_NDIM(values) != 2) {
        PyErr_SetString(PyExc_ValueError, "values must be a 2D float64 array with shape (ncol, nlev)");
        return -1;
    }
    if (PyArray_NDIM(ak) != 1 || PyArray_NDIM(bk) != 1 || PyArray_NDIM(ps) != 1 || PyArray_NDIM(levels) != 1) {
        PyErr_SetString(PyExc_ValueError, "ak, bk, ps, and levels must be 1D float64 arrays");
        return -1;
    }
    if (PyArray_DIM(values, 0) > INT_MAX || PyArray_DIM(values, 1) > INT_MAX || PyArray_DIM(levels, 0) > INT_MAX) {
        PyErr_SetString(PyExc_ValueError, "input dimensions exceed the native int interface");
        return -1;
    }
    if (PyArray_DIM(ak, 0) != PyArray_DIM(values, 1) + 1 || PyArray_DIM(bk, 0) != PyArray_DIM(values, 1) + 1) {
        PyErr_SetString(PyExc_ValueError, "ak and bk must have nlev + 1 half-level coefficients");
        return -1;
    }
    if (PyArray_DIM(ps, 0) != PyArray_DIM(values, 0)) {
        PyErr_SetString(PyExc_ValueError, "ps length must match values.shape[0]");
        return -1;
    }
    return 0;
}

static int check_column_pressures(PyArrayObject *values, PyArrayObject *target_pressures)
{
    if (PyArray_NDIM(target_pressures) != 2) {
        PyErr_SetString(PyExc_ValueError, "target_pressures must be a 2D float64 array with shape (ncol, nout)");
        return -1;
    }
    if (PyArray_DIM(target_pressures, 0) != PyArray_DIM(values, 0)) {
        PyErr_SetString(PyExc_ValueError, "target_pressures.shape[0] must match values.shape[0]");
        return -1;
    }
    if (PyArray_DIM(target_pressures, 1) > INT_MAX) {
        PyErr_SetString(PyExc_ValueError, "target pressure dimension exceeds the native int interface");
        return -1;
    }
    return 0;
}

static int check_same_shape_2d(PyArrayObject *reference, PyArrayObject *candidate, const char *name)
{
    if (PyArray_NDIM(candidate) != 2 ||
        PyArray_DIM(candidate, 0) != PyArray_DIM(reference, 0) ||
        PyArray_DIM(candidate, 1) != PyArray_DIM(reference, 1)) {
        PyErr_Format(PyExc_ValueError, "%s must have shape (%ld, %ld)",
                     name,
                     (long)PyArray_DIM(reference, 0),
                     (long)PyArray_DIM(reference, 1));
        return -1;
    }
    return 0;
}

static int check_target_hybrid(PyArrayObject *target_ak, PyArrayObject *target_bk)
{
    if (PyArray_NDIM(target_ak) != 1 || PyArray_NDIM(target_bk) != 1) {
        PyErr_SetString(PyExc_ValueError, "target_ak and target_bk must be 1D float64 arrays");
        return -1;
    }
    if (PyArray_DIM(target_ak, 0) != PyArray_DIM(target_bk, 0)) {
        PyErr_SetString(PyExc_ValueError, "target_ak and target_bk must have matching lengths");
        return -1;
    }
    if (PyArray_DIM(target_ak, 0) < 2) {
        PyErr_SetString(PyExc_ValueError, "target hybrid coefficients must contain at least two half levels");
        return -1;
    }
    if (PyArray_DIM(target_ak, 0) - 1 > INT_MAX) {
        PyErr_SetString(PyExc_ValueError, "target hybrid dimension exceeds the native int interface");
        return -1;
    }
    return 0;
}

static PyObject *pressure_ppq(PyObject *self, PyObject *args)
{
    PyObject *values_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *levels_obj = NULL;
    PyArrayObject *values = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *levels = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOO", &values_obj, &ak_obj, &bk_obj, &ps_obj, &levels_obj)) {
        return NULL;
    }
    values = as_fortran_array(values_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    levels = as_c_array(levels_obj, NPY_DOUBLE, 1, 1);
    if (!values || !ak || !bk || !ps || !levels || check_common(values, ak, bk, ps, levels) != 0) {
        goto fail;
    }
    ncol = (int)PyArray_DIM(values, 0);
    nlev = (int)PyArray_DIM(values, 1);
    nout = (int)PyArray_DIM(levels, 0);
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_pressure_ppq_c(&ncol, &nlev, &nout,
                           (const double *)PyArray_DATA(values),
                           (const double *)PyArray_DATA(ak),
                           (const double *)PyArray_DATA(bk),
                           (const double *)PyArray_DATA(ps),
                           (const double *)PyArray_DATA(levels),
                           (double *)PyArray_DATA(out),
                           &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPQ pressure interpolation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    return (PyObject *)out;

fail:
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *pressure_ppt(PyObject *self, PyObject *args)
{
    PyObject *values_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *levels_obj = NULL;
    PyArrayObject *values = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *levels = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOO", &values_obj, &ak_obj, &bk_obj, &ps_obj, &levels_obj)) {
        return NULL;
    }
    values = as_fortran_array(values_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    levels = as_c_array(levels_obj, NPY_DOUBLE, 1, 1);
    if (!values || !ak || !bk || !ps || !levels || check_common(values, ak, bk, ps, levels) != 0) {
        goto fail;
    }
    ncol = (int)PyArray_DIM(values, 0);
    nlev = (int)PyArray_DIM(values, 1);
    nout = (int)PyArray_DIM(levels, 0);
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_pressure_ppt_c(&ncol, &nlev, &nout,
                           (const double *)PyArray_DATA(values),
                           (const double *)PyArray_DATA(ak),
                           (const double *)PyArray_DATA(bk),
                           (const double *)PyArray_DATA(ps),
                           (const double *)PyArray_DATA(levels),
                           (double *)PyArray_DATA(out),
                           &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPT pressure interpolation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    return (PyObject *)out;

fail:
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *pressure_ppuv(PyObject *self, PyObject *args)
{
    PyObject *u_obj = NULL, *v_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *levels_obj = NULL;
    PyArrayObject *u = NULL, *v = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *levels = NULL;
    PyArrayObject *u_out = NULL, *v_out = NULL;
    PyObject *tuple = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOO", &u_obj, &v_obj, &ak_obj, &bk_obj, &ps_obj, &levels_obj)) {
        return NULL;
    }
    u = as_fortran_array(u_obj, NPY_DOUBLE, 2, 2);
    v = as_fortran_array(v_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    levels = as_c_array(levels_obj, NPY_DOUBLE, 1, 1);
    if (!u || !v || !ak || !bk || !ps || !levels || check_common(u, ak, bk, ps, levels) != 0) {
        goto fail;
    }
    if (PyArray_NDIM(v) != 2 || PyArray_DIM(v, 0) != PyArray_DIM(u, 0) || PyArray_DIM(v, 1) != PyArray_DIM(u, 1)) {
        PyErr_SetString(PyExc_ValueError, "u and v must have matching 2D shapes");
        goto fail;
    }
    ncol = (int)PyArray_DIM(u, 0);
    nlev = (int)PyArray_DIM(u, 1);
    nout = (int)PyArray_DIM(levels, 0);
    dims[0] = ncol;
    dims[1] = nout;
    u_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    v_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!u_out || !v_out) {
        goto fail;
    }
    fullpos_pressure_ppuv_c(&ncol, &nlev, &nout,
                            (const double *)PyArray_DATA(u),
                            (const double *)PyArray_DATA(v),
                            (const double *)PyArray_DATA(ak),
                            (const double *)PyArray_DATA(bk),
                            (const double *)PyArray_DATA(ps),
                            (const double *)PyArray_DATA(levels),
                            (double *)PyArray_DATA(u_out),
                            (double *)PyArray_DATA(v_out),
                            &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPUV pressure interpolation failed with ierr=%d", ierr);
        goto fail;
    }
    tuple = PyTuple_New(2);
    if (!tuple) {
        goto fail;
    }
    PyTuple_SET_ITEM(tuple, 0, (PyObject *)u_out);
    PyTuple_SET_ITEM(tuple, 1, (PyObject *)v_out);
    u_out = NULL;
    v_out = NULL;
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    return tuple;

fail:
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    Py_XDECREF(u_out);
    Py_XDECREF(v_out);
    Py_XDECREF(tuple);
    return NULL;
}

static PyObject *column_pressure_ppq(PyObject *self, PyObject *args)
{
    PyObject *values_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *targets_obj = NULL;
    PyArrayObject *values = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *targets = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOO", &values_obj, &ak_obj, &bk_obj, &ps_obj, &targets_obj)) {
        return NULL;
    }
    values = as_fortran_array(values_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    targets = as_fortran_array(targets_obj, NPY_DOUBLE, 2, 2);
    if (!values || !ak || !bk || !ps || !targets || check_common(values, ak, bk, ps, ps) != 0 || check_column_pressures(values, targets) != 0) {
        goto fail;
    }
    ncol = (int)PyArray_DIM(values, 0);
    nlev = (int)PyArray_DIM(values, 1);
    nout = (int)PyArray_DIM(targets, 1);
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_column_pressure_ppq_c(&ncol, &nlev, &nout,
                                  (const double *)PyArray_DATA(values),
                                  (const double *)PyArray_DATA(ak),
                                  (const double *)PyArray_DATA(bk),
                                  (const double *)PyArray_DATA(ps),
                                  (const double *)PyArray_DATA(targets),
                                  (double *)PyArray_DATA(out),
                                  &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPQ column-pressure interpolation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(targets);
    return (PyObject *)out;

fail:
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(targets);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *column_pressure_ppt(PyObject *self, PyObject *args)
{
    PyObject *values_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *targets_obj = NULL;
    PyArrayObject *values = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *targets = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOO", &values_obj, &ak_obj, &bk_obj, &ps_obj, &targets_obj)) {
        return NULL;
    }
    values = as_fortran_array(values_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    targets = as_fortran_array(targets_obj, NPY_DOUBLE, 2, 2);
    if (!values || !ak || !bk || !ps || !targets || check_common(values, ak, bk, ps, ps) != 0 || check_column_pressures(values, targets) != 0) {
        goto fail;
    }
    ncol = (int)PyArray_DIM(values, 0);
    nlev = (int)PyArray_DIM(values, 1);
    nout = (int)PyArray_DIM(targets, 1);
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_column_pressure_ppt_c(&ncol, &nlev, &nout,
                                  (const double *)PyArray_DATA(values),
                                  (const double *)PyArray_DATA(ak),
                                  (const double *)PyArray_DATA(bk),
                                  (const double *)PyArray_DATA(ps),
                                  (const double *)PyArray_DATA(targets),
                                  (double *)PyArray_DATA(out),
                                  &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPT column-pressure interpolation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(targets);
    return (PyObject *)out;

fail:
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(targets);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *apache_column_pressure_tuvq(PyObject *self, PyObject *args)
{
    PyObject *t_obj = NULL, *u_obj = NULL, *v_obj = NULL, *q_obj = NULL;
    PyObject *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *targets_obj = NULL;
    PyObject *target_ps_obj = NULL;
    PyArrayObject *t = NULL, *u = NULL, *v = NULL, *q = NULL;
    PyArrayObject *ak = NULL, *bk = NULL, *ps = NULL, *targets = NULL, *target_ps = NULL;
    PyArrayObject *t_out = NULL, *u_out = NULL, *v_out = NULL, *q_out = NULL;
    int ncol, nlev, nout, ierr = 0, enable_lescale = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOOOO|Op", &t_obj, &u_obj, &v_obj, &q_obj, &ak_obj, &bk_obj, &ps_obj, &targets_obj, &target_ps_obj, &enable_lescale)) {
        return NULL;
    }
    if (target_ps_obj && PyBool_Check(target_ps_obj)) {
        enable_lescale = PyObject_IsTrue(target_ps_obj);
        target_ps_obj = NULL;
    }
    t = as_fortran_array(t_obj, NPY_DOUBLE, 2, 2);
    u = as_fortran_array(u_obj, NPY_DOUBLE, 2, 2);
    v = as_fortran_array(v_obj, NPY_DOUBLE, 2, 2);
    q = as_fortran_array(q_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    targets = as_fortran_array(targets_obj, NPY_DOUBLE, 2, 2);
    if (target_ps_obj && target_ps_obj != Py_None) {
        target_ps = as_c_array(target_ps_obj, NPY_DOUBLE, 1, 1);
    } else {
        target_ps = (PyArrayObject *)ps;
        Py_INCREF(target_ps);
    }
    if (!t || !u || !v || !q || !ak || !bk || !ps || !targets ||
        !target_ps ||
        check_same_shape_2d(t, u, "u") != 0 ||
        check_same_shape_2d(t, v, "v") != 0 ||
        check_same_shape_2d(t, q, "q") != 0 ||
        check_common(t, ak, bk, ps, ps) != 0 ||
        check_column_pressures(t, targets) != 0 ||
        PyArray_DIM(target_ps, 0) != PyArray_DIM(t, 0)) {
        if (target_ps && PyArray_DIM(target_ps, 0) != PyArray_DIM(t, 0)) {
            PyErr_SetString(PyExc_ValueError, "target_surface_pressure length must match values.shape[0]");
        }
        goto fail;
    }
    ncol = (int)PyArray_DIM(t, 0);
    nlev = (int)PyArray_DIM(t, 1);
    nout = (int)PyArray_DIM(targets, 1);
    dims[0] = ncol;
    dims[1] = nout;
    t_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    u_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    v_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    q_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!t_out || !u_out || !v_out || !q_out) {
        goto fail;
    }
    fullpos_apache_column_pressure_tuvq_c(&ncol, &nlev, &nout,
                                          (const double *)PyArray_DATA(t),
                                          (const double *)PyArray_DATA(u),
                                          (const double *)PyArray_DATA(v),
                                          (const double *)PyArray_DATA(q),
                                          (const double *)PyArray_DATA(ak),
                                          (const double *)PyArray_DATA(bk),
                                          (const double *)PyArray_DATA(ps),
                                          (const double *)PyArray_DATA(targets),
                                          (const double *)PyArray_DATA(target_ps),
                                          &enable_lescale,
                                          (double *)PyArray_DATA(t_out),
                                          (double *)PyArray_DATA(u_out),
                                          (double *)PyArray_DATA(v_out),
                                          (double *)PyArray_DATA(q_out),
                                          &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS APACHE column-pressure interpolation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(t);
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(q);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(targets);
    Py_XDECREF(target_ps);
    return Py_BuildValue("NNNN", (PyObject *)t_out, (PyObject *)u_out, (PyObject *)v_out, (PyObject *)q_out);

fail:
    Py_XDECREF(t);
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(q);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(targets);
    Py_XDECREF(target_ps);
    Py_XDECREF(t_out);
    Py_XDECREF(u_out);
    Py_XDECREF(v_out);
    Py_XDECREF(q_out);
    return NULL;
}

static PyObject *column_pressure_ppuv(PyObject *self, PyObject *args)
{
    PyObject *u_obj = NULL, *v_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *targets_obj = NULL;
    PyArrayObject *u = NULL, *v = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *targets = NULL;
    PyArrayObject *u_out = NULL, *v_out = NULL;
    PyObject *tuple = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOO", &u_obj, &v_obj, &ak_obj, &bk_obj, &ps_obj, &targets_obj)) {
        return NULL;
    }
    u = as_fortran_array(u_obj, NPY_DOUBLE, 2, 2);
    v = as_fortran_array(v_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    targets = as_fortran_array(targets_obj, NPY_DOUBLE, 2, 2);
    if (!u || !v || !ak || !bk || !ps || !targets || check_common(u, ak, bk, ps, ps) != 0 || check_column_pressures(u, targets) != 0) {
        goto fail;
    }
    if (PyArray_NDIM(v) != 2 || PyArray_DIM(v, 0) != PyArray_DIM(u, 0) || PyArray_DIM(v, 1) != PyArray_DIM(u, 1)) {
        PyErr_SetString(PyExc_ValueError, "u and v must have matching 2D shapes");
        goto fail;
    }
    ncol = (int)PyArray_DIM(u, 0);
    nlev = (int)PyArray_DIM(u, 1);
    nout = (int)PyArray_DIM(targets, 1);
    dims[0] = ncol;
    dims[1] = nout;
    u_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    v_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!u_out || !v_out) {
        goto fail;
    }
    fullpos_column_pressure_ppuv_c(&ncol, &nlev, &nout,
                                   (const double *)PyArray_DATA(u),
                                   (const double *)PyArray_DATA(v),
                                   (const double *)PyArray_DATA(ak),
                                   (const double *)PyArray_DATA(bk),
                                   (const double *)PyArray_DATA(ps),
                                   (const double *)PyArray_DATA(targets),
                                   (double *)PyArray_DATA(u_out),
                                   (double *)PyArray_DATA(v_out),
                                   &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPUV column-pressure interpolation failed with ierr=%d", ierr);
        goto fail;
    }
    tuple = PyTuple_New(2);
    if (!tuple) {
        goto fail;
    }
    PyTuple_SET_ITEM(tuple, 0, (PyObject *)u_out);
    PyTuple_SET_ITEM(tuple, 1, (PyObject *)v_out);
    u_out = NULL;
    v_out = NULL;
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(targets);
    return tuple;

fail:
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(targets);
    Py_XDECREF(u_out);
    Py_XDECREF(v_out);
    Py_XDECREF(tuple);
    return NULL;
}

static PyObject *hybrid_pressure_ppq(PyObject *self, PyObject *args)
{
    PyObject *values_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *target_ak_obj = NULL, *target_bk_obj = NULL;
    PyArrayObject *values = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *target_ak = NULL, *target_bk = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOO", &values_obj, &ak_obj, &bk_obj, &ps_obj, &target_ak_obj, &target_bk_obj)) {
        return NULL;
    }
    values = as_fortran_array(values_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    target_ak = as_c_array(target_ak_obj, NPY_DOUBLE, 1, 1);
    target_bk = as_c_array(target_bk_obj, NPY_DOUBLE, 1, 1);
    if (!values || !ak || !bk || !ps || !target_ak || !target_bk || check_common(values, ak, bk, ps, ps) != 0 || check_target_hybrid(target_ak, target_bk) != 0) {
        goto fail;
    }
    ncol = (int)PyArray_DIM(values, 0);
    nlev = (int)PyArray_DIM(values, 1);
    nout = (int)PyArray_DIM(target_ak, 0) - 1;
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_hybrid_pressure_ppq_c(&ncol, &nlev, &nout,
                                  (const double *)PyArray_DATA(values),
                                  (const double *)PyArray_DATA(ak),
                                  (const double *)PyArray_DATA(bk),
                                  (const double *)PyArray_DATA(ps),
                                  (const double *)PyArray_DATA(target_ak),
                                  (const double *)PyArray_DATA(target_bk),
                                  (double *)PyArray_DATA(out),
                                  &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPQ hybrid-target interpolation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(target_ak);
    Py_XDECREF(target_bk);
    return (PyObject *)out;

fail:
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(target_ak);
    Py_XDECREF(target_bk);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *hybrid_pressure_ppt(PyObject *self, PyObject *args)
{
    PyObject *values_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *target_ak_obj = NULL, *target_bk_obj = NULL;
    PyArrayObject *values = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *target_ak = NULL, *target_bk = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOO", &values_obj, &ak_obj, &bk_obj, &ps_obj, &target_ak_obj, &target_bk_obj)) {
        return NULL;
    }
    values = as_fortran_array(values_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    target_ak = as_c_array(target_ak_obj, NPY_DOUBLE, 1, 1);
    target_bk = as_c_array(target_bk_obj, NPY_DOUBLE, 1, 1);
    if (!values || !ak || !bk || !ps || !target_ak || !target_bk || check_common(values, ak, bk, ps, ps) != 0 || check_target_hybrid(target_ak, target_bk) != 0) {
        goto fail;
    }
    ncol = (int)PyArray_DIM(values, 0);
    nlev = (int)PyArray_DIM(values, 1);
    nout = (int)PyArray_DIM(target_ak, 0) - 1;
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_hybrid_pressure_ppt_c(&ncol, &nlev, &nout,
                                  (const double *)PyArray_DATA(values),
                                  (const double *)PyArray_DATA(ak),
                                  (const double *)PyArray_DATA(bk),
                                  (const double *)PyArray_DATA(ps),
                                  (const double *)PyArray_DATA(target_ak),
                                  (const double *)PyArray_DATA(target_bk),
                                  (double *)PyArray_DATA(out),
                                  &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPT hybrid-target interpolation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(target_ak);
    Py_XDECREF(target_bk);
    return (PyObject *)out;

fail:
    Py_XDECREF(values);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(target_ak);
    Py_XDECREF(target_bk);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *hybrid_pressure_ppuv(PyObject *self, PyObject *args)
{
    PyObject *u_obj = NULL, *v_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *target_ak_obj = NULL, *target_bk_obj = NULL;
    PyArrayObject *u = NULL, *v = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *target_ak = NULL, *target_bk = NULL;
    PyArrayObject *u_out = NULL, *v_out = NULL;
    PyObject *tuple = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOOO", &u_obj, &v_obj, &ak_obj, &bk_obj, &ps_obj, &target_ak_obj, &target_bk_obj)) {
        return NULL;
    }
    u = as_fortran_array(u_obj, NPY_DOUBLE, 2, 2);
    v = as_fortran_array(v_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    target_ak = as_c_array(target_ak_obj, NPY_DOUBLE, 1, 1);
    target_bk = as_c_array(target_bk_obj, NPY_DOUBLE, 1, 1);
    if (!u || !v || !ak || !bk || !ps || !target_ak || !target_bk || check_common(u, ak, bk, ps, ps) != 0 || check_target_hybrid(target_ak, target_bk) != 0) {
        goto fail;
    }
    if (PyArray_NDIM(v) != 2 || PyArray_DIM(v, 0) != PyArray_DIM(u, 0) || PyArray_DIM(v, 1) != PyArray_DIM(u, 1)) {
        PyErr_SetString(PyExc_ValueError, "u and v must have matching 2D shapes");
        goto fail;
    }
    ncol = (int)PyArray_DIM(u, 0);
    nlev = (int)PyArray_DIM(u, 1);
    nout = (int)PyArray_DIM(target_ak, 0) - 1;
    dims[0] = ncol;
    dims[1] = nout;
    u_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    v_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!u_out || !v_out) {
        goto fail;
    }
    fullpos_hybrid_pressure_ppuv_c(&ncol, &nlev, &nout,
                                   (const double *)PyArray_DATA(u),
                                   (const double *)PyArray_DATA(v),
                                   (const double *)PyArray_DATA(ak),
                                   (const double *)PyArray_DATA(bk),
                                   (const double *)PyArray_DATA(ps),
                                   (const double *)PyArray_DATA(target_ak),
                                   (const double *)PyArray_DATA(target_bk),
                                   (double *)PyArray_DATA(u_out),
                                   (double *)PyArray_DATA(v_out),
                                   &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPUV hybrid-target interpolation failed with ierr=%d", ierr);
        goto fail;
    }
    tuple = PyTuple_New(2);
    if (!tuple) {
        goto fail;
    }
    PyTuple_SET_ITEM(tuple, 0, (PyObject *)u_out);
    PyTuple_SET_ITEM(tuple, 1, (PyObject *)v_out);
    u_out = NULL;
    v_out = NULL;
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(target_ak);
    Py_XDECREF(target_bk);
    return tuple;

fail:
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(target_ak);
    Py_XDECREF(target_bk);
    Py_XDECREF(u_out);
    Py_XDECREF(v_out);
    Py_XDECREF(tuple);
    return NULL;
}

static PyObject *eta_pressures(PyObject *self, PyObject *args)
{
    PyObject *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *levels_obj = NULL;
    PyArrayObject *ak = NULL, *bk = NULL, *ps = NULL, *levels = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOO", &ak_obj, &bk_obj, &ps_obj, &levels_obj)) {
        return NULL;
    }
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    levels = as_c_array(levels_obj, NPY_DOUBLE, 1, 1);
    if (!ak || !bk || !ps || !levels) {
        goto fail;
    }
    if (PyArray_DIM(ak, 0) != PyArray_DIM(bk, 0)) {
        PyErr_SetString(PyExc_ValueError, "ak and bk must have matching lengths");
        goto fail;
    }
    if (PyArray_DIM(ak, 0) < 3) {
        PyErr_SetString(PyExc_ValueError, "ak/bk must contain at least three half levels");
        goto fail;
    }
    if (PyArray_DIM(ps, 0) <= 0) {
        PyErr_SetString(PyExc_ValueError, "surface pressure must contain at least one column");
        goto fail;
    }
    if (PyArray_DIM(levels, 0) <= 0) {
        PyErr_SetString(PyExc_ValueError, "eta levels must contain at least one output level");
        goto fail;
    }
    if (PyArray_DIM(ak, 0) - 1 > INT_MAX || PyArray_DIM(ps, 0) > INT_MAX || PyArray_DIM(levels, 0) > INT_MAX) {
        PyErr_SetString(PyExc_OverflowError, "input dimensions exceed native int range");
        goto fail;
    }
    nlev = (int)PyArray_DIM(ak, 0) - 1;
    ncol = (int)PyArray_DIM(ps, 0);
    nout = (int)PyArray_DIM(levels, 0);
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_eta_pressures_c(&ncol, &nlev, &nout,
                            (const double *)PyArray_DATA(ak),
                            (const double *)PyArray_DATA(bk),
                            (const double *)PyArray_DATA(ps),
                            (const double *)PyArray_DATA(levels),
                            (double *)PyArray_DATA(out),
                            &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPLETA eta target-pressure computation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    return (PyObject *)out;

fail:
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *theta_pressures(PyObject *self, PyObject *args)
{
    PyObject *temperature_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *levels_obj = NULL;
    PyArrayObject *temperature = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *levels = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOO", &temperature_obj, &ak_obj, &bk_obj, &ps_obj, &levels_obj)) {
        return NULL;
    }
    temperature = as_fortran_array(temperature_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    levels = as_c_array(levels_obj, NPY_DOUBLE, 1, 1);
    if (!temperature || !ak || !bk || !ps || !levels || check_common(temperature, ak, bk, ps, levels) != 0) {
        goto fail;
    }
    ncol = (int)PyArray_DIM(temperature, 0);
    nlev = (int)PyArray_DIM(temperature, 1);
    nout = (int)PyArray_DIM(levels, 0);
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_theta_pressures_c(&ncol, &nlev, &nout,
                              (const double *)PyArray_DATA(temperature),
                              (const double *)PyArray_DATA(ak),
                              (const double *)PyArray_DATA(bk),
                              (const double *)PyArray_DATA(ps),
                              (const double *)PyArray_DATA(levels),
                              (double *)PyArray_DATA(out),
                              &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPLTETA target-pressure computation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(temperature);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    return (PyObject *)out;

fail:
    Py_XDECREF(temperature);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *temperature_pressures(PyObject *self, PyObject *args)
{
    PyObject *temperature_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *levels_obj = NULL;
    PyArrayObject *temperature = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *levels = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOO", &temperature_obj, &ak_obj, &bk_obj, &ps_obj, &levels_obj)) {
        return NULL;
    }
    temperature = as_fortran_array(temperature_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    levels = as_c_array(levels_obj, NPY_DOUBLE, 1, 1);
    if (!temperature || !ak || !bk || !ps || !levels || check_common(temperature, ak, bk, ps, levels) != 0) {
        goto fail;
    }
    ncol = (int)PyArray_DIM(temperature, 0);
    nlev = (int)PyArray_DIM(temperature, 1);
    nout = (int)PyArray_DIM(levels, 0);
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_temperature_pressures_c(&ncol, &nlev, &nout,
                                    (const double *)PyArray_DATA(temperature),
                                    (const double *)PyArray_DATA(ak),
                                    (const double *)PyArray_DATA(bk),
                                    (const double *)PyArray_DATA(ps),
                                    (const double *)PyArray_DATA(levels),
                                    (double *)PyArray_DATA(out),
                                    &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS temperature target-pressure computation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(temperature);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    return (PyObject *)out;

fail:
    Py_XDECREF(temperature);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(levels);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *height_above_orography_pressures(PyObject *self, PyObject *args)
{
    PyObject *temperature_obj = NULL, *humidity_obj = NULL, *ak_obj = NULL, *bk_obj = NULL;
    PyObject *ps_obj = NULL, *orog_obj = NULL, *levels_obj = NULL;
    PyArrayObject *temperature = NULL, *humidity = NULL, *ak = NULL, *bk = NULL;
    PyArrayObject *ps = NULL, *orog = NULL, *levels = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0, has_humidity = 0;
    const double *humidity_data = NULL;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOOO", &temperature_obj, &humidity_obj, &ak_obj, &bk_obj, &ps_obj, &orog_obj, &levels_obj)) {
        return NULL;
    }
    temperature = as_fortran_array(temperature_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    orog = as_c_array(orog_obj, NPY_DOUBLE, 1, 1);
    levels = as_c_array(levels_obj, NPY_DOUBLE, 1, 1);
    if (!temperature || !ak || !bk || !ps || !orog || !levels || check_common(temperature, ak, bk, ps, levels) != 0) {
        goto fail;
    }
    if (humidity_obj != Py_None) {
        humidity = as_fortran_array(humidity_obj, NPY_DOUBLE, 2, 2);
        if (!humidity || check_same_shape_2d(temperature, humidity, "specific_humidity") != 0) {
            goto fail;
        }
        has_humidity = 1;
        humidity_data = (const double *)PyArray_DATA(humidity);
    } else {
        humidity_data = (const double *)PyArray_DATA(temperature);
    }
    if (PyArray_DIM(orog, 0) != PyArray_DIM(temperature, 0)) {
        PyErr_SetString(PyExc_ValueError, "orography_geopotential length must match temperature.shape[0]");
        goto fail;
    }

    ncol = (int)PyArray_DIM(temperature, 0);
    nlev = (int)PyArray_DIM(temperature, 1);
    nout = (int)PyArray_DIM(levels, 0);
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_height_above_orography_pressures_c(&ncol, &nlev, &nout,
                                               (const double *)PyArray_DATA(temperature),
                                               humidity_data,
                                               &has_humidity,
                                               (const double *)PyArray_DATA(ak),
                                               (const double *)PyArray_DATA(bk),
                                               (const double *)PyArray_DATA(ps),
                                               (const double *)PyArray_DATA(orog),
                                               (const double *)PyArray_DATA(levels),
                                               (double *)PyArray_DATA(out),
                                               &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS height-above-orography target-pressure computation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(temperature);
    Py_XDECREF(humidity);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(orog);
    Py_XDECREF(levels);
    return (PyObject *)out;

fail:
    Py_XDECREF(temperature);
    Py_XDECREF(humidity);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(orog);
    Py_XDECREF(levels);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *height_above_sea_pressures(PyObject *self, PyObject *args)
{
    PyObject *temperature_obj = NULL, *humidity_obj = NULL, *ak_obj = NULL, *bk_obj = NULL;
    PyObject *ps_obj = NULL, *orog_obj = NULL, *levels_obj = NULL;
    PyArrayObject *temperature = NULL, *humidity = NULL, *ak = NULL, *bk = NULL;
    PyArrayObject *ps = NULL, *orog = NULL, *levels = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0, has_humidity = 0;
    const double *humidity_data = NULL;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOOO", &temperature_obj, &humidity_obj, &ak_obj, &bk_obj, &ps_obj, &orog_obj, &levels_obj)) {
        return NULL;
    }
    temperature = as_fortran_array(temperature_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    orog = as_c_array(orog_obj, NPY_DOUBLE, 1, 1);
    levels = as_c_array(levels_obj, NPY_DOUBLE, 1, 1);
    if (!temperature || !ak || !bk || !ps || !orog || !levels || check_common(temperature, ak, bk, ps, levels) != 0) {
        goto fail;
    }
    if (humidity_obj != Py_None) {
        humidity = as_fortran_array(humidity_obj, NPY_DOUBLE, 2, 2);
        if (!humidity || check_same_shape_2d(temperature, humidity, "specific_humidity") != 0) {
            goto fail;
        }
        has_humidity = 1;
        humidity_data = (const double *)PyArray_DATA(humidity);
    } else {
        humidity_data = (const double *)PyArray_DATA(temperature);
    }
    if (PyArray_DIM(orog, 0) != PyArray_DIM(temperature, 0)) {
        PyErr_SetString(PyExc_ValueError, "orography_geopotential length must match temperature.shape[0]");
        goto fail;
    }

    ncol = (int)PyArray_DIM(temperature, 0);
    nlev = (int)PyArray_DIM(temperature, 1);
    nout = (int)PyArray_DIM(levels, 0);
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_height_above_sea_pressures_c(&ncol, &nlev, &nout,
                                         (const double *)PyArray_DATA(temperature),
                                         humidity_data,
                                         &has_humidity,
                                         (const double *)PyArray_DATA(ak),
                                         (const double *)PyArray_DATA(bk),
                                         (const double *)PyArray_DATA(ps),
                                         (const double *)PyArray_DATA(orog),
                                         (const double *)PyArray_DATA(levels),
                                         (double *)PyArray_DATA(out),
                                         &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS height-above-sea target-pressure computation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(temperature);
    Py_XDECREF(humidity);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(orog);
    Py_XDECREF(levels);
    return (PyObject *)out;

fail:
    Py_XDECREF(temperature);
    Py_XDECREF(humidity);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(orog);
    Py_XDECREF(levels);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *potential_vorticity_pressures(PyObject *self, PyObject *args)
{
    PyObject *pv_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *coriolis_obj = NULL, *levels_obj = NULL;
    PyArrayObject *pv = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *coriolis = NULL, *levels = NULL, *out = NULL;
    int ncol, nlev, nout, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOO", &pv_obj, &ak_obj, &bk_obj, &ps_obj, &coriolis_obj, &levels_obj)) {
        return NULL;
    }
    pv = as_fortran_array(pv_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    coriolis = as_c_array(coriolis_obj, NPY_DOUBLE, 1, 1);
    levels = as_c_array(levels_obj, NPY_DOUBLE, 1, 1);
    if (!pv || !ak || !bk || !ps || !coriolis || !levels || check_common(pv, ak, bk, ps, levels) != 0) {
        goto fail;
    }
    if (PyArray_DIM(coriolis, 0) != PyArray_DIM(pv, 0)) {
        PyErr_SetString(PyExc_ValueError, "coriolis length must match pv.shape[0]");
        goto fail;
    }
    ncol = (int)PyArray_DIM(pv, 0);
    nlev = (int)PyArray_DIM(pv, 1);
    nout = (int)PyArray_DIM(levels, 0);
    dims[0] = ncol;
    dims[1] = nout;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }
    fullpos_potential_vorticity_pressures_c(&ncol, &nlev, &nout,
                                            (const double *)PyArray_DATA(pv),
                                            (const double *)PyArray_DATA(ak),
                                            (const double *)PyArray_DATA(bk),
                                            (const double *)PyArray_DATA(ps),
                                            (const double *)PyArray_DATA(coriolis),
                                            (const double *)PyArray_DATA(levels),
                                            (double *)PyArray_DATA(out),
                                            &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS PPLTP potential-vorticity pressure computation failed with ierr=%d", ierr);
        goto fail;
    }
    Py_XDECREF(pv);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(coriolis);
    Py_XDECREF(levels);
    return (PyObject *)out;

fail:
    Py_XDECREF(pv);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(coriolis);
    Py_XDECREF(levels);
    Py_XDECREF(out);
    return NULL;
}

static PyObject *diagnose_potential_vorticity(PyObject *self, PyObject *args)
{
    PyObject *u_obj = NULL, *v_obj = NULL, *temperature_obj = NULL, *vort_obj = NULL;
    PyObject *tm_obj = NULL, *tl_obj = NULL, *spm_obj = NULL, *spl_obj = NULL;
    PyObject *kappa_obj = NULL, *ak_obj = NULL, *bk_obj = NULL, *ps_obj = NULL, *coriolis_obj = NULL;
    PyArrayObject *u = NULL, *v = NULL, *temperature = NULL, *vort = NULL;
    PyArrayObject *tm = NULL, *tl = NULL, *spm = NULL, *spl = NULL;
    PyArrayObject *kappa = NULL, *ak = NULL, *bk = NULL, *ps = NULL, *coriolis = NULL;
    PyArrayObject *pv_out = NULL, *theta_out = NULL;
    int ncol, nlev, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "OOOOOOOOOOOOO",
                          &u_obj,
                          &v_obj,
                          &temperature_obj,
                          &vort_obj,
                          &tm_obj,
                          &tl_obj,
                          &spm_obj,
                          &spl_obj,
                          &kappa_obj,
                          &ak_obj,
                          &bk_obj,
                          &ps_obj,
                          &coriolis_obj)) {
        return NULL;
    }

    u = as_fortran_array(u_obj, NPY_DOUBLE, 2, 2);
    v = as_fortran_array(v_obj, NPY_DOUBLE, 2, 2);
    temperature = as_fortran_array(temperature_obj, NPY_DOUBLE, 2, 2);
    vort = as_fortran_array(vort_obj, NPY_DOUBLE, 2, 2);
    tm = as_fortran_array(tm_obj, NPY_DOUBLE, 2, 2);
    tl = as_fortran_array(tl_obj, NPY_DOUBLE, 2, 2);
    spm = as_c_array(spm_obj, NPY_DOUBLE, 1, 1);
    spl = as_c_array(spl_obj, NPY_DOUBLE, 1, 1);
    kappa = as_fortran_array(kappa_obj, NPY_DOUBLE, 2, 2);
    ak = as_c_array(ak_obj, NPY_DOUBLE, 1, 1);
    bk = as_c_array(bk_obj, NPY_DOUBLE, 1, 1);
    ps = as_c_array(ps_obj, NPY_DOUBLE, 1, 1);
    coriolis = as_c_array(coriolis_obj, NPY_DOUBLE, 1, 1);
    if (!u || !v || !temperature || !vort || !tm || !tl || !spm || !spl || !kappa || !ak || !bk || !ps || !coriolis) {
        goto fail;
    }
    if (PyArray_DIM(u, 0) > INT_MAX || PyArray_DIM(u, 1) > INT_MAX) {
        PyErr_SetString(PyExc_ValueError, "input dimensions exceed the native int interface");
        goto fail;
    }
    if (PyArray_DIM(ak, 0) != PyArray_DIM(u, 1) + 1 || PyArray_DIM(bk, 0) != PyArray_DIM(u, 1) + 1) {
        PyErr_SetString(PyExc_ValueError, "ak and bk must have nlev + 1 half-level coefficients");
        goto fail;
    }
    if (PyArray_DIM(ps, 0) != PyArray_DIM(u, 0) || PyArray_DIM(coriolis, 0) != PyArray_DIM(u, 0) ||
        PyArray_DIM(spm, 0) != PyArray_DIM(u, 0) || PyArray_DIM(spl, 0) != PyArray_DIM(u, 0)) {
        PyErr_SetString(PyExc_ValueError, "surface vectors must have length ncol");
        goto fail;
    }
    if (check_same_shape_2d(u, v, "v") != 0 ||
        check_same_shape_2d(u, temperature, "temperature") != 0 ||
        check_same_shape_2d(u, vort, "relative_vorticity") != 0 ||
        check_same_shape_2d(u, tm, "temperature_meridional_gradient") != 0 ||
        check_same_shape_2d(u, tl, "temperature_zonal_gradient") != 0 ||
        check_same_shape_2d(u, kappa, "kappa") != 0) {
        goto fail;
    }

    ncol = (int)PyArray_DIM(u, 0);
    nlev = (int)PyArray_DIM(u, 1);
    dims[0] = ncol;
    dims[1] = nlev;
    pv_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    theta_out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!pv_out || !theta_out) {
        goto fail;
    }

    fullpos_diagnose_potential_vorticity_c(&ncol,
                                           &nlev,
                                           (const double *)PyArray_DATA(u),
                                           (const double *)PyArray_DATA(v),
                                           (const double *)PyArray_DATA(temperature),
                                           (const double *)PyArray_DATA(vort),
                                           (const double *)PyArray_DATA(tm),
                                           (const double *)PyArray_DATA(tl),
                                           (const double *)PyArray_DATA(spm),
                                           (const double *)PyArray_DATA(spl),
                                           (const double *)PyArray_DATA(kappa),
                                           (const double *)PyArray_DATA(ak),
                                           (const double *)PyArray_DATA(bk),
                                           (const double *)PyArray_DATA(ps),
                                           (const double *)PyArray_DATA(coriolis),
                                           (double *)PyArray_DATA(pv_out),
                                           (double *)PyArray_DATA(theta_out),
                                           &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS GPPVO potential-vorticity diagnostic failed with ierr=%d", ierr);
        goto fail;
    }

    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(temperature);
    Py_XDECREF(vort);
    Py_XDECREF(tm);
    Py_XDECREF(tl);
    Py_XDECREF(spm);
    Py_XDECREF(spl);
    Py_XDECREF(kappa);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(coriolis);
    return Py_BuildValue("NN", (PyObject *)pv_out, (PyObject *)theta_out);

fail:
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(temperature);
    Py_XDECREF(vort);
    Py_XDECREF(tm);
    Py_XDECREF(tl);
    Py_XDECREF(spm);
    Py_XDECREF(spl);
    Py_XDECREF(kappa);
    Py_XDECREF(ak);
    Py_XDECREF(bk);
    Py_XDECREF(ps);
    Py_XDECREF(coriolis);
    Py_XDECREF(pv_out);
    Py_XDECREF(theta_out);
    return NULL;
}

static PyObject *gprcp_kappa(PyObject *self, PyObject *args)
{
    PyObject *q_obj = NULL;
    PyArrayObject *q = NULL;
    PyArrayObject *out = NULL;
    int ncol, nlev, ierr = 0;
    npy_intp dims[2];
    (void)self;

    if (!PyArg_ParseTuple(args, "O", &q_obj)) {
        return NULL;
    }

    q = as_fortran_array(q_obj, NPY_DOUBLE, 2, 2);
    if (!q) {
        goto fail;
    }
    if (PyArray_DIM(q, 0) > INT_MAX || PyArray_DIM(q, 1) > INT_MAX) {
        PyErr_SetString(PyExc_ValueError, "input dimensions exceed the native int interface");
        goto fail;
    }
    ncol = (int)PyArray_DIM(q, 0);
    nlev = (int)PyArray_DIM(q, 1);
    dims[0] = ncol;
    dims[1] = nlev;
    out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
    if (!out) {
        goto fail;
    }

    fullpos_gprcp_kappa_c(&ncol,
                          &nlev,
                          (const double *)PyArray_DATA(q),
                          (double *)PyArray_DATA(out),
                          &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_RuntimeError, "FULLPOS GPRCP kappa calculation failed with ierr=%d", ierr);
        goto fail;
    }

    Py_XDECREF(q);
    return (PyObject *)out;

fail:
    Py_XDECREF(q);
    Py_XDECREF(out);
    return NULL;
}

static PyMethodDef methods[] = {
    {"pressure_ppq", pressure_ppq, METH_VARARGS, "Interpolate scalar fields with native FULLPOS PPQ."},
    {"pressure_ppt", pressure_ppt, METH_VARARGS, "Interpolate temperature with native FULLPOS PPT."},
    {"pressure_ppuv", pressure_ppuv, METH_VARARGS, "Interpolate wind components with native FULLPOS PPUV."},
    {"apache_column_pressure_tuvq", apache_column_pressure_tuvq, METH_VARARGS, "Interpolate T/U/V/Q to per-column pressures with native FULLPOS APACHE."},
    {"column_pressure_ppq", column_pressure_ppq, METH_VARARGS, "Interpolate scalar fields to per-column pressures with native FULLPOS PPQ."},
    {"column_pressure_ppt", column_pressure_ppt, METH_VARARGS, "Interpolate temperature to per-column pressures with native FULLPOS PPT."},
    {"column_pressure_ppuv", column_pressure_ppuv, METH_VARARGS, "Interpolate wind components to per-column pressures with native FULLPOS PPUV."},
    {"hybrid_pressure_ppq", hybrid_pressure_ppq, METH_VARARGS, "Interpolate scalar fields to target hybrid levels with native FULLPOS PPQ."},
    {"hybrid_pressure_ppt", hybrid_pressure_ppt, METH_VARARGS, "Interpolate temperature to target hybrid levels with native FULLPOS PPT."},
    {"hybrid_pressure_ppuv", hybrid_pressure_ppuv, METH_VARARGS, "Interpolate wind components to target hybrid levels with native FULLPOS PPUV."},
    {"eta_pressures", eta_pressures, METH_VARARGS, "Compute eta/model-level-index target pressures with native FULLPOS PPLETA/GPHPRE."},
    {"theta_pressures", theta_pressures, METH_VARARGS, "Compute potential-temperature target pressures with native FULLPOS GPTET/PPLTETA."},
    {"temperature_pressures", temperature_pressures, METH_VARARGS, "Compute temperature target pressures with native FULLPOS PPLTW/FPPS."},
    {"height_above_orography_pressures", height_above_orography_pressures, METH_VARARGS, "Compute height-above-orography target pressures with native FULLPOS GPHPRE/GPGEO/FPPS."},
    {"height_above_sea_pressures", height_above_sea_pressures, METH_VARARGS, "Compute height-above-sea target pressures with native FULLPOS GPHPRE/GPGEO/FPPS."},
    {"potential_vorticity_pressures", potential_vorticity_pressures, METH_VARARGS, "Compute potential-vorticity target pressures with native FULLPOS PPLTP."},
    {"gprcp_kappa", gprcp_kappa, METH_VARARGS, "Compute R/Cp from specific humidity with native FULLPOS GPRCP."},
    {"diagnose_potential_vorticity", diagnose_potential_vorticity, METH_VARARGS, "Diagnose model-level potential vorticity and potential temperature with native FULLPOS GPPVO."},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_vertical_native",
    "Native FULLPOS vertical interpolation wrappers.",
    -1,
    methods,
};

PyMODINIT_FUNC PyInit__vertical_native(void)
{
    import_array();
    return PyModule_Create(&moduledef);
}
