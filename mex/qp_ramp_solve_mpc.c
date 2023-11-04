#include "mex.h"
#include "../qp_ramp.h"
#include <stdlib.h>

static double *invq;
static double *y;
static double *v;
static uint8_t *setarr1;
static ssize_t *setarr2;
static ssize_t *setarr3;
static iterable_set_t a_set;

static void allocate(size_t c) {
	invq = (double*)mxMalloc((mwSize)(c*c*sizeof(double)));
	y = (double*)mxMalloc((mwSize)(c*sizeof(double)));
	v = (double*)mxMalloc((mwSize)(c*sizeof(double)));
	setarr1 = (uint8_t*)mxMalloc((mwSize)(c*sizeof(uint8_t))); 
	setarr2 = (ssize_t*)mxMalloc((mwSize)(c*sizeof(ssize_t))); 
	setarr3 = (ssize_t*)mxMalloc((mwSize)(c*sizeof(ssize_t))); 
    a_set = (iterable_set_t) {
        .capacity = c,
        .elements = setarr1,
        .next = setarr2,
        .prev = setarr3,
    };
}

static void deallocate() {
	mxFree(invq);
	mxFree(y);
	mxFree(v);
	mxFree(setarr1);
	mxFree(setarr2);
	mxFree(setarr3);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("MyToolbox:qp_ramp_solve_mpc:nrhs","Seven inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:qp_ramp_solve_mpc:nlhs","One output required.");
    }

    if (!mxIsScalar(prhs[0])) {
        mexErrMsgTxt("Input p must be scalar.");
    }
    for (int i = 1; i < nrhs; ++i) {
        if( !mxIsDouble(prhs[i]) || 
             mxIsComplex(prhs[i])) {
            mexErrMsgIdAndTxt("MyToolbox:qp_ramp_solve_mpc:notDouble","Input matrix must be type double.");
        }
    }

    // Note that transposed matrices are needed because of MATLABs column-major memory layout.
    size_t c = mxGetM(prhs[1]);
    size_t n = mxGetM(prhs[2]);
    size_t m = mxGetM(prhs[5]);

    if(mxGetM(prhs[1])!=c || mxGetN(prhs[1])!=c) {
        mexErrMsgTxt("Input neg_g_invh_gt must be c by c.");
    }
    if(mxGetM(prhs[2])!=n || mxGetN(prhs[2])!=c) {
        mexErrMsgTxt("Input neg_s must be n by c.");
    }
    if(mxGetM(prhs[3])!=1 || mxGetN(prhs[3])!=c) {
        mexErrMsgTxt("Input neg_w must be 1 by c.");
    }
    if(mxGetM(prhs[4])!=n || mxGetN(prhs[4])!=m) {
        mexErrMsgTxt("Input neg_invh_f must be n by m.");
    }
    if(mxGetM(prhs[5])!=m || mxGetN(prhs[5])!=c) {
        mexErrMsgTxt("Input neg_g_invh must be m by c.");
    }
    if(mxGetM(prhs[6])!=1 || mxGetN(prhs[6])!=n) {
        mexErrMsgTxt("Input x must be 1 by n.");
    }

    size_t p = (size_t)mxGetScalar(prhs[0]);
    double *neg_g_invh_gt = mxGetPr(prhs[1]);
    double *neg_s = mxGetPr(prhs[2]);
    double *neg_w = mxGetPr(prhs[3]);
    double *neg_invh_f = mxGetPr(prhs[4]);
    double *neg_g_invh = mxGetPr(prhs[5]);
    double *x = mxGetPr(prhs[6]);

    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)m,mxREAL);
    double *u = mxGetPr(plhs[0]);

    allocate(c);
    set_init(&a_set);
    qp_ramp_solve_mpc(c, n, m, p, (double(*)[])neg_g_invh_gt, (double(*)[])neg_s, neg_w, (double(*)[])neg_invh_f, (double(*)[])neg_g_invh, x, (double(*)[])invq, &a_set, y, v, u);
    deallocate();
}
 
