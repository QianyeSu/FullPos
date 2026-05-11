#define PY_SSIZE_T_CLEAN
#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_1_20_API_VERSION
#include <numpy/arrayobject.h>

#include <stdlib.h>
#include <string.h>
#include <limits.h>

#include "ectrans/transi.h"

static int native_initialized = 0;

extern void fullpos_fpint4_c(
    const int *kaslb1,
    const int *kfprow,
    const int *kfields,
    const int *kgpst,
    const int *kgpend,
    const int *kfproma,
    const int *kfldbuf,
    const int *kder,
    const int *kdmp,
    const int *ldml_i,
    const int *ldnrst_i,
    const int *kbox,
    const int *kl0,
    const double *pwxx,
    const double *pwxy,
    const int *ldmask_i,
    const double *pbuf,
    double *prow,
    const double *pundef);

extern void fullpos_fpint12_c(
    const int *kaslb1,
    const int *kfprow,
    const int *kfields,
    const int *kgpst,
    const int *kgpend,
    const int *kfproma,
    const int *kfldbuf,
    const int *kder,
    const int *kdmp,
    const int *ldml_i,
    const int *ldnrst_i,
    const int *ldmono_i,
    const int *kbox,
    const int *kl0,
    const double *pwxx,
    const double *pwxy,
    const int *ldmask_i,
    const double *pbuf,
    double *prow,
    const double *pundef);

extern void fullpos_fpavg_c(
    const int *kaslb1,
    const int *kslwide,
    const int *kfields,
    const int *kmasks,
    const int *kgpst,
    const int *kgpend,
    const int *kfproma,
    const int *kfldbuf,
    const int *kmask,
    const int *ks0,
    const int *ldmask_i,
    const double *pbuf,
    const double *pmask,
    const double *pundef,
    double *prow);

extern void fullpos_fpnear_c(
    const int *kaslb1,
    const int *kslwide,
    const int *kfields,
    const int *kmasks,
    const int *kgpst,
    const int *kgpend,
    const int *kfproma,
    const int *kfldbuf,
    const int *kmask,
    const int *ks0,
    const int *ldmask_i,
    const double *pbuf,
    const double *pmask,
    const double *pundef,
    double *prow);

extern void fullpos_horizontal_regular_c(
    const int *nlat,
    const int *nsrc_points,
    const int *ntarget_points,
    const int *kbinl,
    const double *values,
    const int *nloen,
    const double *source_lats,
    const double *target_lats,
    const double *target_lons,
    double *output,
    int *ierr);

extern void fullpos_horizontal_halo_c(
    const int *nlat,
    const int *nsrc_points,
    const int *ntarget_points,
    const int *kslwide,
    const int *use_near,
    const double *values,
    const int *nloen,
    const double *source_lats,
    const double *target_lats,
    const double *target_lons,
    double *output,
    int *ierr);

static int set_trans_error(const char *call, int errcode)
{
    PyErr_Format(PyExc_RuntimeError, "%s failed: %s", call, trans_error_msg(errcode));
    return -1;
}

#define CHECK_TRANS(call)                           \
    do {                                           \
        int _err = (call);                         \
        if (_err != TRANS_SUCCESS) {               \
            set_trans_error(#call, _err);          \
            goto fail;                             \
        }                                          \
    } while (0)

static int setup_transform(struct Trans_t *trans, PyArrayObject *pl, int trunc)
{
    int ndgl = (int)PyArray_SIZE(pl);
    int *nloen = (int *)PyArray_DATA(pl);

    CHECK_TRANS(trans_new(trans));
    CHECK_TRANS(trans_set_resol(trans, ndgl, nloen));
    CHECK_TRANS(trans_set_trunc(trans, trunc));
    CHECK_TRANS(trans_setup(trans));
    return 0;

fail:
    return -1;
}

static int ensure_native_initialized(void)
{
    if (native_initialized) {
        return 0;
    }
    int err = trans_use_mpi(0);
    if (err != TRANS_SUCCESS) {
        return set_trans_error("trans_use_mpi(0)", err);
    }
    err = trans_init();
    if (err != TRANS_SUCCESS) {
        return set_trans_error("trans_init()", err);
    }
    native_initialized = 1;
    return 0;
}

static PyObject *ectrans_regrid(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"values", "source_pl", "target_pl", "trunc", NULL};
    PyObject *values_obj = NULL;
    PyObject *source_pl_obj = NULL;
    PyObject *target_pl_obj = NULL;
    int trunc = -1;

    PyArrayObject *values = NULL;
    PyArrayObject *source_pl = NULL;
    PyArrayObject *target_pl = NULL;
    PyArrayObject *out = NULL;

    struct Trans_t src;
    struct Trans_t dst;
    int src_ready = 0;
    int dst_ready = 0;
    int *nfrom = NULL;
    int *nto = NULL;

    int input_ndim = 0;
    int nfld = 1;
    npy_intp input_points = 0;

    double *rgp_src = NULL;
    double *rspec_src = NULL;
    double *rspec_global = NULL;
    double *rspec_dst = NULL;
    double *rgp_dst = NULL;

    (void)self;
    memset(&src, 0, sizeof(src));
    memset(&dst, 0, sizeof(dst));

    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "OOOi", kwlist, &values_obj, &source_pl_obj, &target_pl_obj, &trunc)) {
        return NULL;
    }
    if (trunc < 0) {
        PyErr_SetString(PyExc_ValueError, "trunc must be non-negative");
        return NULL;
    }

    values = (PyArrayObject *)PyArray_FROM_OTF(values_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    source_pl = (PyArrayObject *)PyArray_FROM_OTF(source_pl_obj, NPY_INT, NPY_ARRAY_IN_ARRAY);
    target_pl = (PyArrayObject *)PyArray_FROM_OTF(target_pl_obj, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (!values || !source_pl || !target_pl) {
        goto fail;
    }
    input_ndim = PyArray_NDIM(values);
    if ((input_ndim != 1 && input_ndim != 2) ||
        PyArray_NDIM(source_pl) != 1 ||
        PyArray_NDIM(target_pl) != 1) {
        PyErr_SetString(PyExc_ValueError, "values must be 1D or 2D; source_pl and target_pl must be 1D arrays");
        goto fail;
    }
    if (input_ndim == 1) {
        nfld = 1;
        input_points = PyArray_DIM(values, 0);
    } else {
        if (PyArray_DIM(values, 0) > INT_MAX) {
            PyErr_SetString(PyExc_ValueError, "too many fields for TRANSI int interface");
            goto fail;
        }
        nfld = (int)PyArray_DIM(values, 0);
        input_points = PyArray_DIM(values, 1);
    }

    if (ensure_native_initialized() != 0) {
        goto fail;
    }

    if (setup_transform(&src, source_pl, trunc) != 0) {
        goto fail;
    }
    src_ready = 1;
    if (setup_transform(&dst, target_pl, trunc) != 0) {
        goto fail;
    }
    dst_ready = 1;

    if ((size_t)input_points != (size_t)src.ngptotg) {
        PyErr_Format(
            PyExc_ValueError,
            "values has %zd points per field, source grid expects %d",
            input_points,
            src.ngptotg);
        goto fail;
    }

    rgp_src = (double *)calloc((size_t)nfld * (size_t)src.ngptot, sizeof(double));
    rspec_src = (double *)calloc((size_t)nfld * (size_t)src.nspec2, sizeof(double));
    rspec_global = (double *)calloc((size_t)nfld * (size_t)src.nspec2g, sizeof(double));
    rspec_dst = (double *)calloc((size_t)nfld * (size_t)dst.nspec2, sizeof(double));
    rgp_dst = (double *)calloc((size_t)nfld * (size_t)dst.ngptot, sizeof(double));
    nfrom = (int *)malloc((size_t)nfld * sizeof(int));
    nto = (int *)malloc((size_t)nfld * sizeof(int));
    if (!rgp_src || !rspec_src || !rspec_global || !rspec_dst || !rgp_dst || !nfrom || !nto) {
        PyErr_NoMemory();
        goto fail;
    }
    for (int i = 0; i < nfld; ++i) {
        nfrom[i] = 1;
        nto[i] = 1;
    }

    {
        struct DistGrid_t distgrid = new_distgrid(&src);
        distgrid.nfrom = nfrom;
        distgrid.rgpg = (const double *)PyArray_DATA(values);
        distgrid.rgp = rgp_src;
        distgrid.nfld = nfld;
        CHECK_TRANS(trans_distgrid(&distgrid));
    }

    {
        struct DirTrans_t dirtrans = new_dirtrans(&src);
        dirtrans.nscalar = nfld;
        dirtrans.rgp = rgp_src;
        dirtrans.rspscalar = rspec_src;
        CHECK_TRANS(trans_dirtrans(&dirtrans));
    }

    {
        struct GathSpec_t gathspec = new_gathspec(&src);
        gathspec.rspec = rspec_src;
        gathspec.rspecg = rspec_global;
        gathspec.nfld = nfld;
        gathspec.nto = nto;
        CHECK_TRANS(trans_gathspec(&gathspec));
    }

    {
        struct DistSpec_t distspec = new_distspec(&dst);
        distspec.nfrom = nfrom;
        distspec.rspecg = rspec_global;
        distspec.rspec = rspec_dst;
        distspec.nfld = nfld;
        CHECK_TRANS(trans_distspec(&distspec));
    }

    {
        struct InvTrans_t invtrans = new_invtrans(&dst);
        invtrans.nscalar = nfld;
        invtrans.rspscalar = rspec_dst;
        invtrans.rgp = rgp_dst;
        CHECK_TRANS(trans_invtrans(&invtrans));
    }

    {
        npy_intp dims[2] = {(npy_intp)nfld, (npy_intp)dst.ngptotg};
        if (input_ndim == 1) {
            out = (PyArrayObject *)PyArray_SimpleNew(1, &dims[1], NPY_DOUBLE);
        } else {
            out = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        }
        if (!out) {
            goto fail;
        }

        struct GathGrid_t gathgrid = new_gathgrid(&dst);
        gathgrid.rgp = rgp_dst;
        gathgrid.rgpg = (double *)PyArray_DATA(out);
        gathgrid.nfld = nfld;
        gathgrid.nto = nto;
        CHECK_TRANS(trans_gathgrid(&gathgrid));
    }

    CHECK_TRANS(trans_delete(&src));
    src_ready = 0;
    CHECK_TRANS(trans_delete(&dst));
    dst_ready = 0;

    free(rgp_src);
    free(rspec_src);
    free(rspec_global);
    free(rspec_dst);
    free(rgp_dst);
    free(nfrom);
    free(nto);
    Py_XDECREF(values);
    Py_XDECREF(source_pl);
    Py_XDECREF(target_pl);
    return (PyObject *)out;

fail:
    if (src_ready) {
        trans_delete(&src);
    }
    if (dst_ready) {
        trans_delete(&dst);
    }
    free(rgp_src);
    free(rspec_src);
    free(rspec_global);
    free(rspec_dst);
    free(rgp_dst);
    free(nfrom);
    free(nto);
    Py_XDECREF(out);
    Py_XDECREF(values);
    Py_XDECREF(source_pl);
    Py_XDECREF(target_pl);
    return NULL;
}

static PyObject *ectrans_fit(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"values", "pl", "trunc", NULL};
    PyObject *values_obj = NULL;
    PyObject *pl_obj = NULL;
    int trunc = -1;

    PyArrayObject *values = NULL;
    PyArrayObject *pl = NULL;
    PyArrayObject *out = NULL;

    struct Trans_t trans;
    int trans_ready = 0;
    int *nfrom = NULL;
    int *nto = NULL;

    int input_ndim = 0;
    int nfld = 1;
    npy_intp input_points = 0;

    double *rgp = NULL;
    double *rspec = NULL;
    double *rspec_global = NULL;

    (void)self;
    memset(&trans, 0, sizeof(trans));

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOi", kwlist, &values_obj, &pl_obj, &trunc)) {
        return NULL;
    }
    if (trunc < 0) {
        PyErr_SetString(PyExc_ValueError, "trunc must be non-negative");
        return NULL;
    }

    values = (PyArrayObject *)PyArray_FROM_OTF(values_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    pl = (PyArrayObject *)PyArray_FROM_OTF(pl_obj, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (!values || !pl) {
        goto fail;
    }
    input_ndim = PyArray_NDIM(values);
    if ((input_ndim != 1 && input_ndim != 2) || PyArray_NDIM(pl) != 1) {
        PyErr_SetString(PyExc_ValueError, "values must be 1D or 2D; pl must be a 1D array");
        goto fail;
    }
    if (input_ndim == 1) {
        nfld = 1;
        input_points = PyArray_DIM(values, 0);
    } else {
        if (PyArray_DIM(values, 0) > INT_MAX) {
            PyErr_SetString(PyExc_ValueError, "too many fields for TRANSI int interface");
            goto fail;
        }
        nfld = (int)PyArray_DIM(values, 0);
        input_points = PyArray_DIM(values, 1);
    }

    if (ensure_native_initialized() != 0) {
        goto fail;
    }
    if (setup_transform(&trans, pl, trunc) != 0) {
        goto fail;
    }
    trans_ready = 1;

    if ((size_t)input_points != (size_t)trans.ngptotg) {
        PyErr_Format(
            PyExc_ValueError,
            "values has %zd points per field, grid expects %d",
            input_points,
            trans.ngptotg);
        goto fail;
    }

    rgp = (double *)calloc((size_t)nfld * (size_t)trans.ngptot, sizeof(double));
    rspec = (double *)calloc((size_t)nfld * (size_t)trans.nspec2, sizeof(double));
    rspec_global = (double *)calloc((size_t)nfld * (size_t)trans.nspec2g, sizeof(double));
    nfrom = (int *)malloc((size_t)nfld * sizeof(int));
    nto = (int *)malloc((size_t)nfld * sizeof(int));
    if (!rgp || !rspec || !rspec_global || !nfrom || !nto) {
        PyErr_NoMemory();
        goto fail;
    }
    for (int i = 0; i < nfld; ++i) {
        nfrom[i] = 1;
        nto[i] = 1;
    }

    {
        struct DistGrid_t distgrid = new_distgrid(&trans);
        distgrid.nfrom = nfrom;
        distgrid.rgpg = (const double *)PyArray_DATA(values);
        distgrid.rgp = rgp;
        distgrid.nfld = nfld;
        CHECK_TRANS(trans_distgrid(&distgrid));
    }

    {
        struct DirTrans_t dirtrans = new_dirtrans(&trans);
        dirtrans.nscalar = nfld;
        dirtrans.rgp = rgp;
        dirtrans.rspscalar = rspec;
        CHECK_TRANS(trans_dirtrans(&dirtrans));
    }

    {
        struct GathSpec_t gathspec = new_gathspec(&trans);
        gathspec.rspec = rspec;
        gathspec.rspecg = rspec_global;
        gathspec.nfld = nfld;
        gathspec.nto = nto;
        CHECK_TRANS(trans_gathspec(&gathspec));
    }

    {
        npy_intp dims[2] = {(npy_intp)nfld, (npy_intp)trans.nspec2g};
        if (input_ndim == 1) {
            out = (PyArrayObject *)PyArray_SimpleNew(1, &dims[1], NPY_DOUBLE);
        } else {
            out = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        }
        if (!out) {
            goto fail;
        }
        double *out_data = (double *)PyArray_DATA(out);
        if (nfld == 1) {
            memcpy(out_data, rspec_global, (size_t)trans.nspec2g * sizeof(double));
        } else {
            for (int f = 0; f < nfld; ++f) {
                for (int s = 0; s < trans.nspec2g; ++s) {
                    out_data[(size_t)f * (size_t)trans.nspec2g + (size_t)s] =
                        rspec_global[(size_t)s * (size_t)nfld + (size_t)f];
                }
            }
        }
    }

    CHECK_TRANS(trans_delete(&trans));
    trans_ready = 0;

    free(rgp);
    free(rspec);
    free(rspec_global);
    free(nfrom);
    free(nto);
    Py_XDECREF(values);
    Py_XDECREF(pl);
    return (PyObject *)out;

fail:
    if (trans_ready) {
        trans_delete(&trans);
    }
    free(rgp);
    free(rspec);
    free(rspec_global);
    free(nfrom);
    free(nto);
    Py_XDECREF(out);
    Py_XDECREF(values);
    Py_XDECREF(pl);
    return NULL;
}

static PyObject *ectrans_synthesis(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"coefficients", "pl", "trunc", NULL};
    PyObject *coefficients_obj = NULL;
    PyObject *pl_obj = NULL;
    int trunc = -1;

    PyArrayObject *coefficients = NULL;
    PyArrayObject *pl = NULL;
    PyArrayObject *out = NULL;

    struct Trans_t trans;
    int trans_ready = 0;
    int *nfrom = NULL;
    int *nto = NULL;

    int input_ndim = 0;
    int nfld = 1;
    npy_intp input_coeffs = 0;

    double *rspec_global = NULL;
    double *rspec = NULL;
    double *rgp = NULL;

    (void)self;
    memset(&trans, 0, sizeof(trans));

    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "OOi", kwlist, &coefficients_obj, &pl_obj, &trunc)) {
        return NULL;
    }
    if (trunc < 0) {
        PyErr_SetString(PyExc_ValueError, "trunc must be non-negative");
        return NULL;
    }

    coefficients =
        (PyArrayObject *)PyArray_FROM_OTF(coefficients_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    pl = (PyArrayObject *)PyArray_FROM_OTF(pl_obj, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (!coefficients || !pl) {
        goto fail;
    }
    input_ndim = PyArray_NDIM(coefficients);
    if ((input_ndim != 1 && input_ndim != 2) || PyArray_NDIM(pl) != 1) {
        PyErr_SetString(
            PyExc_ValueError,
            "coefficients must be 1D or 2D with shape (nfield, nspec2); pl must be a 1D array");
        goto fail;
    }
    if (input_ndim == 1) {
        nfld = 1;
        input_coeffs = PyArray_DIM(coefficients, 0);
    } else {
        if (PyArray_DIM(coefficients, 0) > INT_MAX) {
            PyErr_SetString(PyExc_ValueError, "too many fields for TRANSI int interface");
            goto fail;
        }
        nfld = (int)PyArray_DIM(coefficients, 0);
        input_coeffs = PyArray_DIM(coefficients, 1);
    }

    if (ensure_native_initialized() != 0) {
        goto fail;
    }
    if (setup_transform(&trans, pl, trunc) != 0) {
        goto fail;
    }
    trans_ready = 1;

    if ((size_t)input_coeffs != (size_t)trans.nspec2g) {
        PyErr_Format(
            PyExc_ValueError,
            "coefficients has %zd values per field, truncation T%d expects %d",
            input_coeffs,
            trunc,
            trans.nspec2g);
        goto fail;
    }

    rspec_global = (double *)calloc((size_t)nfld * (size_t)trans.nspec2g, sizeof(double));
    rspec = (double *)calloc((size_t)nfld * (size_t)trans.nspec2, sizeof(double));
    rgp = (double *)calloc((size_t)nfld * (size_t)trans.ngptot, sizeof(double));
    nfrom = (int *)malloc((size_t)nfld * sizeof(int));
    nto = (int *)malloc((size_t)nfld * sizeof(int));
    if (!rspec_global || !rspec || !rgp || !nfrom || !nto) {
        PyErr_NoMemory();
        goto fail;
    }
    for (int i = 0; i < nfld; ++i) {
        nfrom[i] = 1;
        nto[i] = 1;
    }

    {
        const double *coeff_data = (const double *)PyArray_DATA(coefficients);
        if (nfld == 1) {
            memcpy(rspec_global, coeff_data, (size_t)trans.nspec2g * sizeof(double));
        } else {
            for (int f = 0; f < nfld; ++f) {
                for (int s = 0; s < trans.nspec2g; ++s) {
                    rspec_global[(size_t)s * (size_t)nfld + (size_t)f] =
                        coeff_data[(size_t)f * (size_t)trans.nspec2g + (size_t)s];
                }
            }
        }
    }

    {
        struct DistSpec_t distspec = new_distspec(&trans);
        distspec.nfrom = nfrom;
        distspec.rspecg = rspec_global;
        distspec.rspec = rspec;
        distspec.nfld = nfld;
        CHECK_TRANS(trans_distspec(&distspec));
    }

    {
        struct InvTrans_t invtrans = new_invtrans(&trans);
        invtrans.nscalar = nfld;
        invtrans.rspscalar = rspec;
        invtrans.rgp = rgp;
        CHECK_TRANS(trans_invtrans(&invtrans));
    }

    {
        npy_intp dims[2] = {(npy_intp)nfld, (npy_intp)trans.ngptotg};
        if (input_ndim == 1) {
            out = (PyArrayObject *)PyArray_SimpleNew(1, &dims[1], NPY_DOUBLE);
        } else {
            out = (PyArrayObject *)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
        }
        if (!out) {
            goto fail;
        }

        struct GathGrid_t gathgrid = new_gathgrid(&trans);
        gathgrid.rgp = rgp;
        gathgrid.rgpg = (double *)PyArray_DATA(out);
        gathgrid.nfld = nfld;
        gathgrid.nto = nto;
        CHECK_TRANS(trans_gathgrid(&gathgrid));
    }

    CHECK_TRANS(trans_delete(&trans));
    trans_ready = 0;

    free(rspec_global);
    free(rspec);
    free(rgp);
    free(nfrom);
    free(nto);
    Py_XDECREF(coefficients);
    Py_XDECREF(pl);
    return (PyObject *)out;

fail:
    if (trans_ready) {
        trans_delete(&trans);
    }
    free(rspec_global);
    free(rspec);
    free(rgp);
    free(nfrom);
    free(nto);
    Py_XDECREF(out);
    Py_XDECREF(coefficients);
    Py_XDECREF(pl);
    return NULL;
}

static PyObject *ectrans_vector_scalar_diagnostics(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"u", "v", "scalars", "pl", "trunc", NULL};
    PyObject *u_obj = NULL;
    PyObject *v_obj = NULL;
    PyObject *scalars_obj = NULL;
    PyObject *pl_obj = NULL;
    int trunc = -1;

    PyArrayObject *u = NULL;
    PyArrayObject *v = NULL;
    PyArrayObject *scalars = NULL;
    PyArrayObject *pl = NULL;
    PyArrayObject *vort_out = NULL;
    PyArrayObject *div_out = NULL;
    PyArrayObject *scalar_ns_out = NULL;
    PyArrayObject *scalar_ew_out = NULL;

    struct Trans_t trans;
    int trans_ready = 0;
    int *nfrom = NULL;
    int *nto = NULL;

    int u_ndim = 0;
    int scalar_ndim = 0;
    int nvordiv = 0;
    int nscalar = 0;
    int input_nfld = 0;
    int output_nfld = 0;
    npy_intp input_points = 0;
    npy_intp scalar_points = 0;

    double *rgpg_in = NULL;
    double *rgp_in = NULL;
    double *rspvor = NULL;
    double *rspdiv = NULL;
    double *rspscalar = NULL;
    double *rgp_out = NULL;
    double *rgpg_out = NULL;

    (void)self;
    memset(&trans, 0, sizeof(trans));

    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs, "OOOOi", kwlist, &u_obj, &v_obj, &scalars_obj, &pl_obj, &trunc)) {
        return NULL;
    }
    if (trunc < 0) {
        PyErr_SetString(PyExc_ValueError, "trunc must be non-negative");
        return NULL;
    }

    u = (PyArrayObject *)PyArray_FROM_OTF(u_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    v = (PyArrayObject *)PyArray_FROM_OTF(v_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    scalars = (PyArrayObject *)PyArray_FROM_OTF(scalars_obj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
    pl = (PyArrayObject *)PyArray_FROM_OTF(pl_obj, NPY_INT, NPY_ARRAY_IN_ARRAY);
    if (!u || !v || !scalars || !pl) {
        goto fail;
    }

    u_ndim = PyArray_NDIM(u);
    scalar_ndim = PyArray_NDIM(scalars);
    if ((u_ndim != 1 && u_ndim != 2) ||
        PyArray_NDIM(v) != u_ndim ||
        (scalar_ndim != 1 && scalar_ndim != 2) ||
        PyArray_NDIM(pl) != 1) {
        PyErr_SetString(PyExc_ValueError, "u/v/scalars must be 1D or 2D arrays; pl must be a 1D array");
        goto fail;
    }
    if (u_ndim == 1) {
        nvordiv = 1;
        input_points = PyArray_DIM(u, 0);
        if (PyArray_DIM(v, 0) != input_points) {
            PyErr_SetString(PyExc_ValueError, "v must have the same shape as u");
            goto fail;
        }
    } else {
        if (PyArray_DIM(u, 0) > INT_MAX) {
            PyErr_SetString(PyExc_ValueError, "too many vector fields for TRANSI int interface");
            goto fail;
        }
        nvordiv = (int)PyArray_DIM(u, 0);
        input_points = PyArray_DIM(u, 1);
        if (PyArray_DIM(v, 0) != PyArray_DIM(u, 0) || PyArray_DIM(v, 1) != input_points) {
            PyErr_SetString(PyExc_ValueError, "v must have the same shape as u");
            goto fail;
        }
    }
    if (scalar_ndim == 1) {
        nscalar = 1;
        scalar_points = PyArray_DIM(scalars, 0);
    } else {
        if (PyArray_DIM(scalars, 0) > INT_MAX) {
            PyErr_SetString(PyExc_ValueError, "too many scalar fields for TRANSI int interface");
            goto fail;
        }
        nscalar = (int)PyArray_DIM(scalars, 0);
        scalar_points = PyArray_DIM(scalars, 1);
    }
    if (nvordiv <= 0 || nscalar <= 0) {
        PyErr_SetString(PyExc_ValueError, "at least one vector and one scalar field are required");
        goto fail;
    }
    if (scalar_points != input_points) {
        PyErr_SetString(PyExc_ValueError, "scalars must have the same point count as u/v");
        goto fail;
    }
    if (nvordiv > (INT_MAX - nscalar) / 4) {
        PyErr_SetString(PyExc_ValueError, "too many diagnostic fields for TRANSI int interface");
        goto fail;
    }
    input_nfld = 2 * nvordiv + nscalar;
    output_nfld = 4 * nvordiv + 3 * nscalar;

    if (ensure_native_initialized() != 0) {
        goto fail;
    }
    if (setup_transform(&trans, pl, trunc) != 0) {
        goto fail;
    }
    trans_ready = 1;
    if ((size_t)input_points != (size_t)trans.ngptotg) {
        PyErr_Format(
            PyExc_ValueError,
            "u/v/scalars have %zd points per field, grid expects %d",
            input_points,
            trans.ngptotg);
        goto fail;
    }

    rgpg_in = (double *)calloc((size_t)input_nfld * (size_t)trans.ngptotg, sizeof(double));
    rgp_in = (double *)calloc((size_t)input_nfld * (size_t)trans.ngptot, sizeof(double));
    rspvor = (double *)calloc((size_t)nvordiv * (size_t)trans.nspec2, sizeof(double));
    rspdiv = (double *)calloc((size_t)nvordiv * (size_t)trans.nspec2, sizeof(double));
    rspscalar = (double *)calloc((size_t)nscalar * (size_t)trans.nspec2, sizeof(double));
    rgp_out = (double *)calloc((size_t)output_nfld * (size_t)trans.ngptot, sizeof(double));
    rgpg_out = (double *)calloc((size_t)output_nfld * (size_t)trans.ngptotg, sizeof(double));
    nfrom = (int *)malloc((size_t)input_nfld * sizeof(int));
    nto = (int *)malloc((size_t)output_nfld * sizeof(int));
    if (!rgpg_in || !rgp_in || !rspvor || !rspdiv || !rspscalar || !rgp_out || !rgpg_out || !nfrom || !nto) {
        PyErr_NoMemory();
        goto fail;
    }
    for (int i = 0; i < input_nfld; ++i) {
        nfrom[i] = 1;
    }
    for (int i = 0; i < output_nfld; ++i) {
        nto[i] = 1;
    }

    {
        const double *u_data = (const double *)PyArray_DATA(u);
        const double *v_data = (const double *)PyArray_DATA(v);
        const double *scalar_data = (const double *)PyArray_DATA(scalars);
        size_t points = (size_t)input_points;
        for (int f = 0; f < nvordiv; ++f) {
            memcpy(rgpg_in + (size_t)f * points, u_data + (size_t)f * points, points * sizeof(double));
            memcpy(rgpg_in + (size_t)(nvordiv + f) * points, v_data + (size_t)f * points, points * sizeof(double));
        }
        for (int f = 0; f < nscalar; ++f) {
            memcpy(rgpg_in + (size_t)(2 * nvordiv + f) * points, scalar_data + (size_t)f * points, points * sizeof(double));
        }
    }

    {
        struct DistGrid_t distgrid = new_distgrid(&trans);
        distgrid.nfrom = nfrom;
        distgrid.rgpg = rgpg_in;
        distgrid.rgp = rgp_in;
        distgrid.nfld = input_nfld;
        CHECK_TRANS(trans_distgrid(&distgrid));
    }

    {
        struct DirTrans_t dirtrans = new_dirtrans(&trans);
        dirtrans.nvordiv = nvordiv;
        dirtrans.nscalar = nscalar;
        dirtrans.rgp = rgp_in;
        dirtrans.rspvor = rspvor;
        dirtrans.rspdiv = rspdiv;
        dirtrans.rspscalar = rspscalar;
        CHECK_TRANS(trans_dirtrans(&dirtrans));
    }

    {
        struct InvTrans_t invtrans = new_invtrans(&trans);
        invtrans.nvordiv = nvordiv;
        invtrans.nscalar = nscalar;
        invtrans.lvordivgp = 1;
        invtrans.lscalarders = 1;
        invtrans.luvder_EW = 0;
        invtrans.rspvor = rspvor;
        invtrans.rspdiv = rspdiv;
        invtrans.rspscalar = rspscalar;
        invtrans.rgp = rgp_out;
        CHECK_TRANS(trans_invtrans(&invtrans));
    }

    {
        struct GathGrid_t gathgrid = new_gathgrid(&trans);
        gathgrid.rgp = rgp_out;
        gathgrid.rgpg = rgpg_out;
        gathgrid.nfld = output_nfld;
        gathgrid.nto = nto;
        CHECK_TRANS(trans_gathgrid(&gathgrid));
    }

    {
        npy_intp vector_dims[2] = {(npy_intp)nvordiv, input_points};
        npy_intp scalar_dims[2] = {(npy_intp)nscalar, input_points};
        size_t points = (size_t)input_points;
        int scalar_ns_offset = 4 * nvordiv + nscalar;
        int scalar_ew_offset = 4 * nvordiv + 2 * nscalar;

        vort_out = (PyArrayObject *)PyArray_SimpleNew(2, vector_dims, NPY_DOUBLE);
        div_out = (PyArrayObject *)PyArray_SimpleNew(2, vector_dims, NPY_DOUBLE);
        scalar_ns_out = (PyArrayObject *)PyArray_SimpleNew(2, scalar_dims, NPY_DOUBLE);
        scalar_ew_out = (PyArrayObject *)PyArray_SimpleNew(2, scalar_dims, NPY_DOUBLE);
        if (!vort_out || !div_out || !scalar_ns_out || !scalar_ew_out) {
            goto fail;
        }
        memcpy(PyArray_DATA(vort_out), rgpg_out, (size_t)nvordiv * points * sizeof(double));
        memcpy(PyArray_DATA(div_out), rgpg_out + (size_t)nvordiv * points, (size_t)nvordiv * points * sizeof(double));
        memcpy(PyArray_DATA(scalar_ns_out), rgpg_out + (size_t)scalar_ns_offset * points, (size_t)nscalar * points * sizeof(double));
        memcpy(PyArray_DATA(scalar_ew_out), rgpg_out + (size_t)scalar_ew_offset * points, (size_t)nscalar * points * sizeof(double));
    }

    CHECK_TRANS(trans_delete(&trans));
    trans_ready = 0;

    free(rgpg_in);
    free(rgp_in);
    free(rspvor);
    free(rspdiv);
    free(rspscalar);
    free(rgp_out);
    free(rgpg_out);
    free(nfrom);
    free(nto);
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(scalars);
    Py_XDECREF(pl);
    return Py_BuildValue("NNNN", (PyObject *)vort_out, (PyObject *)div_out, (PyObject *)scalar_ns_out, (PyObject *)scalar_ew_out);

fail:
    if (trans_ready) {
        trans_delete(&trans);
    }
    free(rgpg_in);
    free(rgp_in);
    free(rspvor);
    free(rspdiv);
    free(rspscalar);
    free(rgp_out);
    free(rgpg_out);
    free(nfrom);
    free(nto);
    Py_XDECREF(vort_out);
    Py_XDECREF(div_out);
    Py_XDECREF(scalar_ns_out);
    Py_XDECREF(scalar_ew_out);
    Py_XDECREF(u);
    Py_XDECREF(v);
    Py_XDECREF(scalars);
    Py_XDECREF(pl);
    return NULL;
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

static PyObject *ectrans_fpint4_kernel(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {
        "pbuf",
        "kl0",
        "pwxx",
        "pwxy",
        "kbox",
        "kder",
        "kdmp",
        "ldnrst",
        "ldmask",
        "ldml",
        "pundef",
        NULL};
    PyObject *pbuf_obj = NULL;
    PyObject *kl0_obj = NULL;
    PyObject *pwxx_obj = NULL;
    PyObject *pwxy_obj = NULL;
    PyObject *kbox_obj = NULL;
    PyObject *kder_obj = NULL;
    PyObject *kdmp_obj = NULL;
    PyObject *ldnrst_obj = NULL;
    PyObject *ldmask_obj = NULL;
    int ldml = 0;
    double pundef = 1.0e20;

    PyArrayObject *pbuf = NULL;
    PyArrayObject *kl0 = NULL;
    PyArrayObject *pwxx = NULL;
    PyArrayObject *pwxy = NULL;
    PyArrayObject *kbox = NULL;
    PyArrayObject *kder = NULL;
    PyArrayObject *kdmp = NULL;
    PyArrayObject *ldnrst = NULL;
    PyArrayObject *ldmask = NULL;
    PyArrayObject *out = NULL;

    (void)self;

    if (!PyArg_ParseTupleAndKeywords(
            args,
            kwargs,
            "OOOOOOOOO|id",
            kwlist,
            &pbuf_obj,
            &kl0_obj,
            &pwxx_obj,
            &pwxy_obj,
            &kbox_obj,
            &kder_obj,
            &kdmp_obj,
            &ldnrst_obj,
            &ldmask_obj,
            &ldml,
            &pundef)) {
        return NULL;
    }

    pbuf = as_fortran_array(pbuf_obj, NPY_DOUBLE, 1, 1);
    kl0 = as_fortran_array(kl0_obj, NPY_INT, 2, 2);
    pwxx = as_fortran_array(pwxx_obj, NPY_DOUBLE, 2, 2);
    pwxy = as_fortran_array(pwxy_obj, NPY_DOUBLE, 2, 2);
    kbox = as_fortran_array(kbox_obj, NPY_INT, 1, 1);
    kder = as_fortran_array(kder_obj, NPY_INT, 1, 1);
    kdmp = as_fortran_array(kdmp_obj, NPY_INT, 1, 1);
    ldnrst = as_fortran_array(ldnrst_obj, NPY_INT, 1, 1);
    ldmask = as_fortran_array(ldmask_obj, NPY_INT, 1, 1);
    if (!pbuf || !kl0 || !pwxx || !pwxy || !kbox || !kder || !kdmp || !ldnrst || !ldmask) {
        goto fail;
    }

    int kfproma = (int)PyArray_DIM(kl0, 0);
    int kfprow = (int)PyArray_DIM(kl0, 1);
    int kfields = (int)PyArray_DIM(kder, 0);
    int kfldbuf = (int)PyArray_DIM(kdmp, 0);
    int kaslb1 = (int)(PyArray_SIZE(pbuf) / kfldbuf);
    int kgpst = 1;
    int kgpend = kfproma;

    if (kfprow < 4) {
        PyErr_SetString(PyExc_ValueError, "kl0 must have at least 4 raw addresses for FPINT4");
        goto fail;
    }
    if (PyArray_DIM(pwxx, 0) != kfproma || PyArray_DIM(pwxx, 1) != 4) {
        PyErr_SetString(PyExc_ValueError, "pwxx must have shape (kfproma, 4)");
        goto fail;
    }
    if (PyArray_DIM(pwxy, 0) != kfproma || PyArray_DIM(pwxy, 1) != 12) {
        PyErr_SetString(PyExc_ValueError, "pwxy must have shape (kfproma, 12)");
        goto fail;
    }
    if (PyArray_DIM(kbox, 0) != kfproma ||
        PyArray_DIM(ldnrst, 0) != kfields ||
        PyArray_DIM(ldmask, 0) != kfields) {
        PyErr_SetString(PyExc_ValueError, "kbox, ldnrst, and ldmask dimensions are inconsistent");
        goto fail;
    }
    if (kfldbuf <= 0 || (PyArray_SIZE(pbuf) % kfldbuf) != 0) {
        PyErr_SetString(PyExc_ValueError, "pbuf size must be divisible by kdmp size");
        goto fail;
    }

    {
        npy_intp dims[2] = {(npy_intp)kfproma, (npy_intp)kfields};
        out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
        if (!out) {
            goto fail;
        }
    }

    fullpos_fpint4_c(
        &kaslb1,
        &kfprow,
        &kfields,
        &kgpst,
        &kgpend,
        &kfproma,
        &kfldbuf,
        (const int *)PyArray_DATA(kder),
        (const int *)PyArray_DATA(kdmp),
        &ldml,
        (const int *)PyArray_DATA(ldnrst),
        (const int *)PyArray_DATA(kbox),
        (const int *)PyArray_DATA(kl0),
        (const double *)PyArray_DATA(pwxx),
        (const double *)PyArray_DATA(pwxy),
        (const int *)PyArray_DATA(ldmask),
        (const double *)PyArray_DATA(pbuf),
        (double *)PyArray_DATA(out),
        &pundef);

    Py_XDECREF(pbuf);
    Py_XDECREF(kl0);
    Py_XDECREF(pwxx);
    Py_XDECREF(pwxy);
    Py_XDECREF(kbox);
    Py_XDECREF(kder);
    Py_XDECREF(kdmp);
    Py_XDECREF(ldnrst);
    Py_XDECREF(ldmask);
    return (PyObject *)out;

fail:
    Py_XDECREF(out);
    Py_XDECREF(pbuf);
    Py_XDECREF(kl0);
    Py_XDECREF(pwxx);
    Py_XDECREF(pwxy);
    Py_XDECREF(kbox);
    Py_XDECREF(kder);
    Py_XDECREF(kdmp);
    Py_XDECREF(ldnrst);
    Py_XDECREF(ldmask);
    return NULL;
}

static PyObject *ectrans_fpint12_kernel(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {
        "pbuf", "kl0", "pwxx", "pwxy", "kbox", "kder", "kdmp",
        "ldnrst", "ldmono", "ldmask", "ldml", "pundef", NULL};
    PyObject *pbuf_obj = NULL, *kl0_obj = NULL, *pwxx_obj = NULL, *pwxy_obj = NULL;
    PyObject *kbox_obj = NULL, *kder_obj = NULL, *kdmp_obj = NULL;
    PyObject *ldnrst_obj = NULL, *ldmono_obj = NULL, *ldmask_obj = NULL;
    int ldml = 0;
    double pundef = 1.0e20;
    PyArrayObject *pbuf = NULL, *kl0 = NULL, *pwxx = NULL, *pwxy = NULL;
    PyArrayObject *kbox = NULL, *kder = NULL, *kdmp = NULL;
    PyArrayObject *ldnrst = NULL, *ldmono = NULL, *ldmask = NULL, *out = NULL;
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOOOOOOO|id", kwlist,
                                     &pbuf_obj, &kl0_obj, &pwxx_obj, &pwxy_obj,
                                     &kbox_obj, &kder_obj, &kdmp_obj,
                                     &ldnrst_obj, &ldmono_obj, &ldmask_obj,
                                     &ldml, &pundef)) {
        return NULL;
    }
    pbuf = as_fortran_array(pbuf_obj, NPY_DOUBLE, 1, 1);
    kl0 = as_fortran_array(kl0_obj, NPY_INT, 2, 2);
    pwxx = as_fortran_array(pwxx_obj, NPY_DOUBLE, 2, 2);
    pwxy = as_fortran_array(pwxy_obj, NPY_DOUBLE, 2, 2);
    kbox = as_fortran_array(kbox_obj, NPY_INT, 1, 1);
    kder = as_fortran_array(kder_obj, NPY_INT, 1, 1);
    kdmp = as_fortran_array(kdmp_obj, NPY_INT, 1, 1);
    ldnrst = as_fortran_array(ldnrst_obj, NPY_INT, 1, 1);
    ldmono = as_fortran_array(ldmono_obj, NPY_INT, 1, 1);
    ldmask = as_fortran_array(ldmask_obj, NPY_INT, 1, 1);
    if (!pbuf || !kl0 || !pwxx || !pwxy || !kbox || !kder || !kdmp || !ldnrst || !ldmono || !ldmask) {
        goto fail;
    }

    int kfproma = (int)PyArray_DIM(kl0, 0);
    int kfprow = (int)PyArray_DIM(kl0, 1);
    int kfields = (int)PyArray_DIM(kder, 0);
    int kfldbuf = (int)PyArray_DIM(kdmp, 0);
    int kaslb1 = (int)(PyArray_SIZE(pbuf) / kfldbuf);
    int kgpst = 1;
    int kgpend = kfproma;

    if (kfprow < 4 || PyArray_DIM(pwxx, 0) != kfproma || PyArray_DIM(pwxx, 1) != 12 ||
        PyArray_DIM(pwxy, 0) != kfproma || PyArray_DIM(pwxy, 1) != 12 ||
        PyArray_DIM(kbox, 0) != kfproma || PyArray_DIM(ldnrst, 0) != kfields ||
        PyArray_DIM(ldmono, 0) != kfields || PyArray_DIM(ldmask, 0) != kfields ||
        kfldbuf <= 0 || (PyArray_SIZE(pbuf) % kfldbuf) != 0) {
        PyErr_SetString(PyExc_ValueError, "inconsistent FPINT12 kernel array dimensions");
        goto fail;
    }

    {
        npy_intp dims[2] = {(npy_intp)kfproma, (npy_intp)kfields};
        out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
        if (!out) {
            goto fail;
        }
    }
    fullpos_fpint12_c(&kaslb1, &kfprow, &kfields, &kgpst, &kgpend, &kfproma, &kfldbuf,
                      (const int *)PyArray_DATA(kder), (const int *)PyArray_DATA(kdmp),
                      &ldml, (const int *)PyArray_DATA(ldnrst), (const int *)PyArray_DATA(ldmono),
                      (const int *)PyArray_DATA(kbox), (const int *)PyArray_DATA(kl0),
                      (const double *)PyArray_DATA(pwxx), (const double *)PyArray_DATA(pwxy),
                      (const int *)PyArray_DATA(ldmask), (const double *)PyArray_DATA(pbuf),
                      (double *)PyArray_DATA(out), &pundef);

    Py_XDECREF(pbuf); Py_XDECREF(kl0); Py_XDECREF(pwxx); Py_XDECREF(pwxy);
    Py_XDECREF(kbox); Py_XDECREF(kder); Py_XDECREF(kdmp);
    Py_XDECREF(ldnrst); Py_XDECREF(ldmono); Py_XDECREF(ldmask);
    return (PyObject *)out;
fail:
    Py_XDECREF(out); Py_XDECREF(pbuf); Py_XDECREF(kl0); Py_XDECREF(pwxx); Py_XDECREF(pwxy);
    Py_XDECREF(kbox); Py_XDECREF(kder); Py_XDECREF(kdmp);
    Py_XDECREF(ldnrst); Py_XDECREF(ldmono); Py_XDECREF(ldmask);
    return NULL;
}

static PyObject *fpavg_or_near_kernel(PyObject *args, PyObject *kwargs, int use_near)
{
    static char *kwlist[] = {"pbuf", "ks0", "pmask", "kmask", "ldmask", "pundef", NULL};
    PyObject *pbuf_obj = NULL, *ks0_obj = NULL, *pmask_obj = NULL, *kmask_obj = NULL, *ldmask_obj = NULL;
    double pundef = 1.0e20;
    PyArrayObject *pbuf = NULL, *ks0 = NULL, *pmask = NULL, *kmask = NULL, *ldmask = NULL, *out = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOO|d", kwlist,
                                     &pbuf_obj, &ks0_obj, &pmask_obj, &kmask_obj, &ldmask_obj, &pundef)) {
        return NULL;
    }
    pbuf = as_fortran_array(pbuf_obj, NPY_DOUBLE, 1, 1);
    ks0 = as_fortran_array(ks0_obj, NPY_INT, 2, 2);
    pmask = as_fortran_array(pmask_obj, NPY_DOUBLE, 2, 2);
    kmask = as_fortran_array(kmask_obj, NPY_INT, 1, 1);
    ldmask = as_fortran_array(ldmask_obj, NPY_INT, 1, 1);
    if (!pbuf || !ks0 || !pmask || !kmask || !ldmask) {
        goto fail;
    }

    int kfproma = (int)PyArray_DIM(ks0, 0);
    int kslwide2 = (int)PyArray_DIM(ks0, 1);
    int kslwide = kslwide2 / 2;
    int kmasks = (int)PyArray_DIM(pmask, 1);
    int kfields = (int)PyArray_DIM(kmask, 0);
    int kfldbuf = kfields;
    int kaslb1 = (int)(PyArray_SIZE(pbuf) / kfldbuf);
    int kgpst = 1;
    int kgpend = kfproma;
    if (kslwide <= 0 || kslwide2 != 2 * kslwide ||
        PyArray_DIM(pmask, 0) != kfproma || PyArray_DIM(ldmask, 0) != kfields ||
        kfldbuf <= 0 || (PyArray_SIZE(pbuf) % kfldbuf) != 0) {
        PyErr_SetString(PyExc_ValueError, "inconsistent FPAVG/FPNEAR kernel array dimensions");
        goto fail;
    }
    {
        npy_intp dims[2] = {(npy_intp)kfproma, (npy_intp)kfields};
        out = (PyArrayObject *)PyArray_EMPTY(2, dims, NPY_DOUBLE, 1);
        if (!out) {
            goto fail;
        }
    }
    if (use_near) {
        fullpos_fpnear_c(&kaslb1, &kslwide, &kfields, &kmasks, &kgpst, &kgpend,
                         &kfproma, &kfldbuf, (const int *)PyArray_DATA(kmask),
                         (const int *)PyArray_DATA(ks0), (const int *)PyArray_DATA(ldmask),
                         (const double *)PyArray_DATA(pbuf), (const double *)PyArray_DATA(pmask),
                         &pundef, (double *)PyArray_DATA(out));
    } else {
        fullpos_fpavg_c(&kaslb1, &kslwide, &kfields, &kmasks, &kgpst, &kgpend,
                        &kfproma, &kfldbuf, (const int *)PyArray_DATA(kmask),
                        (const int *)PyArray_DATA(ks0), (const int *)PyArray_DATA(ldmask),
                        (const double *)PyArray_DATA(pbuf), (const double *)PyArray_DATA(pmask),
                        &pundef, (double *)PyArray_DATA(out));
    }
    Py_XDECREF(pbuf); Py_XDECREF(ks0); Py_XDECREF(pmask); Py_XDECREF(kmask); Py_XDECREF(ldmask);
    return (PyObject *)out;
fail:
    Py_XDECREF(out); Py_XDECREF(pbuf); Py_XDECREF(ks0); Py_XDECREF(pmask); Py_XDECREF(kmask); Py_XDECREF(ldmask);
    return NULL;
}

static PyObject *ectrans_fpavg_kernel(PyObject *self, PyObject *args, PyObject *kwargs)
{
    (void)self;
    return fpavg_or_near_kernel(args, kwargs, 0);
}

static PyObject *ectrans_fpnear_kernel(PyObject *self, PyObject *args, PyObject *kwargs)
{
    (void)self;
    return fpavg_or_near_kernel(args, kwargs, 1);
}

static PyObject *ectrans_horizontal_regular_kernel(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"values", "nloen", "source_lats", "target_lats", "target_lons", "method", NULL};
    PyObject *values_obj = NULL, *nloen_obj = NULL, *source_lats_obj = NULL;
    PyObject *target_lats_obj = NULL, *target_lons_obj = NULL;
    const char *method = "bilinear";
    PyArrayObject *values = NULL, *nloen = NULL, *source_lats = NULL, *target_lats = NULL, *target_lons = NULL;
    PyArrayObject *out = NULL;
    int kbinl = 4;
    int ierr = 0;
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOO|s", kwlist,
                                     &values_obj, &nloen_obj, &source_lats_obj,
                                     &target_lats_obj, &target_lons_obj, &method)) {
        return NULL;
    }
    if (strcmp(method, "bilinear") == 0 || strcmp(method, "fpint4") == 0) {
        kbinl = 4;
    } else if (strcmp(method, "quadratic12") == 0 || strcmp(method, "fpint12") == 0) {
        kbinl = 12;
    } else {
        PyErr_SetString(PyExc_ValueError, "method must be 'bilinear'/'fpint4' or 'quadratic12'/'fpint12'");
        return NULL;
    }

    values = as_fortran_array(values_obj, NPY_DOUBLE, 1, 1);
    nloen = as_fortran_array(nloen_obj, NPY_INT, 1, 1);
    source_lats = as_fortran_array(source_lats_obj, NPY_DOUBLE, 1, 1);
    target_lats = as_fortran_array(target_lats_obj, NPY_DOUBLE, 1, 1);
    target_lons = as_fortran_array(target_lons_obj, NPY_DOUBLE, 1, 1);
    if (!values || !nloen || !source_lats || !target_lats || !target_lons) {
        goto fail;
    }
    int nlat = (int)PyArray_DIM(nloen, 0);
    int nsrc_points = (int)PyArray_DIM(values, 0);
    int ntarget_points = (int)PyArray_DIM(target_lats, 0);
    if (PyArray_DIM(source_lats, 0) != nlat ||
        PyArray_DIM(target_lons, 0) != ntarget_points) {
        PyErr_SetString(PyExc_ValueError, "source_lats/nloen or target_lats/target_lons dimensions are inconsistent");
        goto fail;
    }
    {
        npy_intp dims[1] = {(npy_intp)ntarget_points};
        out = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_DOUBLE, 1);
        if (!out) {
            goto fail;
        }
    }

    fullpos_horizontal_regular_c(
        &nlat,
        &nsrc_points,
        &ntarget_points,
        &kbinl,
        (const double *)PyArray_DATA(values),
        (const int *)PyArray_DATA(nloen),
        (const double *)PyArray_DATA(source_lats),
        (const double *)PyArray_DATA(target_lats),
        (const double *)PyArray_DATA(target_lons),
        (double *)PyArray_DATA(out),
        &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_ValueError, "FULLPOS horizontal regular kernel failed with ierr=%d", ierr);
        goto fail;
    }

    Py_XDECREF(values); Py_XDECREF(nloen); Py_XDECREF(source_lats); Py_XDECREF(target_lats); Py_XDECREF(target_lons);
    return (PyObject *)out;
fail:
    Py_XDECREF(out);
    Py_XDECREF(values); Py_XDECREF(nloen); Py_XDECREF(source_lats); Py_XDECREF(target_lats); Py_XDECREF(target_lons);
    return NULL;
}

static PyObject *ectrans_horizontal_halo_kernel(PyObject *self, PyObject *args, PyObject *kwargs)
{
    static char *kwlist[] = {"values", "nloen", "source_lats", "target_lats", "target_lons", "method", "kslwide", NULL};
    PyObject *values_obj = NULL, *nloen_obj = NULL, *source_lats_obj = NULL;
    PyObject *target_lats_obj = NULL, *target_lons_obj = NULL;
    const char *method = "nearest";
    int kslwide = 1;
    int use_near = 1;
    int ierr = 0;
    PyArrayObject *values = NULL, *nloen = NULL, *source_lats = NULL, *target_lats = NULL, *target_lons = NULL;
    PyArrayObject *out = NULL;
    (void)self;

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OOOOO|si", kwlist,
                                     &values_obj, &nloen_obj, &source_lats_obj,
                                     &target_lats_obj, &target_lons_obj, &method, &kslwide)) {
        return NULL;
    }
    if (strcmp(method, "nearest") == 0 || strcmp(method, "fpnear") == 0) {
        use_near = 1;
    } else if (strcmp(method, "average") == 0 || strcmp(method, "fpavg") == 0) {
        use_near = 0;
    } else {
        PyErr_SetString(PyExc_ValueError, "method must be 'nearest'/'fpnear' or 'average'/'fpavg'");
        return NULL;
    }
    if (kslwide <= 0) {
        PyErr_SetString(PyExc_ValueError, "kslwide must be positive");
        return NULL;
    }

    values = as_fortran_array(values_obj, NPY_DOUBLE, 1, 1);
    nloen = as_fortran_array(nloen_obj, NPY_INT, 1, 1);
    source_lats = as_fortran_array(source_lats_obj, NPY_DOUBLE, 1, 1);
    target_lats = as_fortran_array(target_lats_obj, NPY_DOUBLE, 1, 1);
    target_lons = as_fortran_array(target_lons_obj, NPY_DOUBLE, 1, 1);
    if (!values || !nloen || !source_lats || !target_lats || !target_lons) {
        goto fail;
    }
    int nlat = (int)PyArray_DIM(nloen, 0);
    int nsrc_points = (int)PyArray_DIM(values, 0);
    int ntarget_points = (int)PyArray_DIM(target_lats, 0);
    if (PyArray_DIM(source_lats, 0) != nlat ||
        PyArray_DIM(target_lons, 0) != ntarget_points) {
        PyErr_SetString(PyExc_ValueError, "source_lats/nloen or target_lats/target_lons dimensions are inconsistent");
        goto fail;
    }
    {
        npy_intp dims[1] = {(npy_intp)ntarget_points};
        out = (PyArrayObject *)PyArray_EMPTY(1, dims, NPY_DOUBLE, 1);
        if (!out) {
            goto fail;
        }
    }

    fullpos_horizontal_halo_c(
        &nlat,
        &nsrc_points,
        &ntarget_points,
        &kslwide,
        &use_near,
        (const double *)PyArray_DATA(values),
        (const int *)PyArray_DATA(nloen),
        (const double *)PyArray_DATA(source_lats),
        (const double *)PyArray_DATA(target_lats),
        (const double *)PyArray_DATA(target_lons),
        (double *)PyArray_DATA(out),
        &ierr);
    if (ierr != 0) {
        PyErr_Format(PyExc_ValueError, "FULLPOS horizontal halo kernel failed with ierr=%d", ierr);
        goto fail;
    }

    Py_XDECREF(values); Py_XDECREF(nloen); Py_XDECREF(source_lats); Py_XDECREF(target_lats); Py_XDECREF(target_lons);
    return (PyObject *)out;
fail:
    Py_XDECREF(out);
    Py_XDECREF(values); Py_XDECREF(nloen); Py_XDECREF(source_lats); Py_XDECREF(target_lats); Py_XDECREF(target_lons);
    return NULL;
}

static PyMethodDef EctransMethods[] = {
    {"fit", (PyCFunction)ectrans_fit, METH_VARARGS | METH_KEYWORDS,
     "Transform global Gaussian scalar field(s) to ECTRANS spectral coefficients."},
    {"fpint4_kernel", (PyCFunction)ectrans_fpint4_kernel, METH_VARARGS | METH_KEYWORDS,
     "Call the native OpenIFS/FULLPOS FPINT4 interpolation kernel with precomputed addresses and weights."},
    {"fpint12_kernel", (PyCFunction)ectrans_fpint12_kernel, METH_VARARGS | METH_KEYWORDS,
     "Call the native OpenIFS/FULLPOS FPINT12 interpolation kernel with precomputed addresses and weights."},
    {"fpavg_kernel", (PyCFunction)ectrans_fpavg_kernel, METH_VARARGS | METH_KEYWORDS,
     "Call the native OpenIFS/FULLPOS FPAVG interpolation kernel with precomputed halo addresses and masks."},
    {"fpnear_kernel", (PyCFunction)ectrans_fpnear_kernel, METH_VARARGS | METH_KEYWORDS,
     "Call the native OpenIFS/FULLPOS FPNEAR interpolation kernel with precomputed halo addresses and masks."},
    {"horizontal_regular_kernel", (PyCFunction)ectrans_horizontal_regular_kernel, METH_VARARGS | METH_KEYWORDS,
     "Call native OpenIFS/FULLPOS SUHOW1/SUHOW2 plus FPINT4/FPINT12 for regular global rows."},
    {"horizontal_halo_kernel", (PyCFunction)ectrans_horizontal_halo_kernel, METH_VARARGS | METH_KEYWORDS,
     "Call native OpenIFS/FULLPOS SUHOX1 plus FPNEAR/FPAVG for nearest/average halo interpolation."},
    {"regrid", (PyCFunction)ectrans_regrid, METH_VARARGS | METH_KEYWORDS,
     "Regrid one global Gaussian scalar field through ECTRANS spectral coefficients."},
    {"synthesis", (PyCFunction)ectrans_synthesis, METH_VARARGS | METH_KEYWORDS,
     "Transform ECTRANS spectral coefficients to global Gaussian scalar field(s)."},
    {"vector_scalar_diagnostics", (PyCFunction)ectrans_vector_scalar_diagnostics, METH_VARARGS | METH_KEYWORDS,
     "Diagnose vorticity and scalar horizontal derivatives with native ECTRANS transforms."},
    {NULL, NULL, 0, NULL},
};

static struct PyModuleDef EctransModule = {
    PyModuleDef_HEAD_INIT,
    "_ectrans",
    NULL,
    -1,
    EctransMethods,
};

PyMODINIT_FUNC PyInit__ectrans(void)
{
    PyObject *module = PyModule_Create(&EctransModule);
    if (!module) {
        return NULL;
    }
    import_array();
    return module;
}
