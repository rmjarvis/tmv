#ifndef FBLAS_H
#define FBLAS_H

// Converted from cblas.h to use the function names that appear in a 
// typical fortran library.
// Specifically, this was done to work with Goto BLAS, which does not
// have a C interface.

//
// ============
// Level 1 BLAS
// ============
//

struct cdouble { double r; double i;};
struct cfloat { float r; float i;};

double dsdot_(const int& N,
    const float* X, const int& incX, const float* Y, const int& incY);
double ddot_(const int& N, const double* X, const int& incX,
    const double* Y, const int& incY);
float sdsdot_(const int& N, const float& alpha, 
    const float* X, const int& incX, const float* Y, const int& incY);
float sdot_(const int& N, const float* X, const int& incX,
    const float* Y, const int& incY);
void zdotu_(cdouble* ret,
    const int& N, const cdouble* X, const int& incX,
    const cdouble* Y, const int& incY);
void zdotc_(cdouble* ret, 
    const int& N, const cdouble* X, const int& incX,
    const cdouble* Y, const int& incY);
void cdotu_(cfloat* ret,
    const int& N, const cfloat* X, const int& incX,
    const cfloat* Y, const int& incY);
void cdotc_(cfloat* ret,
    const int& N, const cfloat* X, const int& incX,
    const cfloat* Y, const int& incY);

double dnrm2_(const int& N, const double* X, const int& incX);
float snrm2_(const int& N, const float* X, const int& incX);
double dznrm2_(const int& N, const cdouble* X, const int& incX);
float scnrm2_(const int& N, const cfloat* X, const int& incX);

double dasum_(const int& N, const double* X, const int& incX);
float sasum_(const int& N, const float* X, const int& incX);
double dzasum_(const int& N, const cdouble* X, const int& incX);
float scasum_(const int& N, const cfloat* X, const int& incX);

int idamax_(const int& N, const double* X, const int& incX);
int isamax_(const int& N, const float* X, const int& incX);
int izamax_(const int& N, const cdouble* X, const int& incX);
int icamax_(const int& N, const cfloat* X, const int& incX);

void dswap_(const int& N, double* X, const int& incX,
    double* Y, const int& incY);
void sswap_(const int& N, float* X, const int& incX,
    float* Y, const int& incY);
void zswap_(const int& N, cdouble* X, const int& incX,
    cdouble* Y, const int& incY);
void cswap_(const int& N, cfloat* X, const int& incX,
    cfloat* Y, const int& incY);

void dcopy_(const int& N, const double* X, const int& incX,
    double* Y, const int& incY);
void scopy_(const int& N, const float* X, const int& incX,
    float* Y, const int& incY);
void zcopy_(const int& N, const cdouble* X, const int& incX,
    cdouble* Y, const int& incY);
void ccopy_(const int& N, const cfloat* X, const int& incX,
    cfloat* Y, const int& incY);

void daxpy_(const int& N, const double& alpha,
    const double* X, const int& incX, double* Y, const int& incY);
void saxpy_(const int& N, const float& alpha,
    const float* X, const int& incX, float* Y, const int& incY);
void zaxpy_(const int& N, const cdouble* alpha,
    const cdouble* X, const int& incX, cdouble* Y, const int& incY);
void caxpy_(const int& N, const cfloat* alpha, 
    const cfloat* X, const int& incX, cfloat* Y, const int& incY);

void drotg_(double& a, double& b, double& c, double& s);
void srotg_(float& a, float& b, float& c, float& s);
void zrotg_(cdouble* a, cdouble* b, cdouble* c, cdouble* s);
void crotg_(cfloat* a, cfloat* b, cfloat* c, cfloat* s);

void drot_(const int& N, double* X, const int& incX,
    double* Y, const int& incY, const double& c, const double& s);
void srot_(const int& N, float* X, const int& incX,
    float* Y, const int& incY, const float& c, const float& s);
//void zdrot_(const int& N, cdouble* X, const int& incX,
    //cdouble* Y, const int& incY, const double& c, const double& s);
//void csrot_(const int& N, cfloat* X, const int& incX,
    //cfloat* Y, const int& incY, const float& c, const float& s);

void drotmg_(double& d1, double& d2, double& b1, const double& b2, double* P);
void srotmg_(float& d1, float& d2, float& b1, const float& b2, float* P);

void drotm_(const int& N, double* X, const int& incX,
    double* Y, const int& incY, const double* P);
void srotm_(const int& N, float* X, const int& incX,
    float* Y, const int& incY, const float* P);

void dscal_(const int& N, const double& alpha, double* X, const int& incX);
void sscal_(const int& N, const float& alpha, float* X, const int& incX);
void zscal_(const int& N, const cdouble* alpha, cdouble* X, const int& incX);
void zdscal_(const int& N, const double& alpha, cdouble* X, const int& incX);
void cscal_(const int& N, const cfloat* alpha, cfloat* X, const int& incX);
void csscal_(const int& N, const float& alpha, cfloat* X, const int& incX);


// ============
// Level 2 BLAS
// ============

void dgemv_(const char& TransA, const int& M, const int& N,
    const double& alpha, const double* A, const int& lda,
    const double* X, const int& incX, const double& beta,
    double* Y, const int& incY, const int tlen);
void sgemv_(const char& TransA, const int& M, const int& N,
    const float& alpha, const float* A, const int& lda,
    const float* X, const int& incX, const float& beta,
    float* Y, const int& incY, const int tlen);
void zgemv_(const char& TransA, const int& M, const int& N,
    const cdouble* alpha, const cdouble* A, const int& lda,
    const cdouble* X, const int& incX, const cdouble* beta,
    cdouble* Y, const int& incY, const int tlen);
void cgemv_(const char& TransA, const int& M, const int& N,
    const cfloat* alpha, const cfloat* A, const int& lda,
    const cfloat* X, const int& incX, const cfloat* beta,
    cfloat* Y, const int& incY, const int tlen);

void dgbmv_(const char& TransA, const int& M, const int& N,
    const int& KL, const int& KU, const double& alpha,
    const double* A, const int& lda, const double* X, const int& incX,
    const double& beta, double* Y, const int& incY, const int tlen);
void sgbmv_(const char& TransA, const int& M, const int& N,
    const int& KL, const int& KU, const float& alpha,
    const float* A, const int& lda, const float* X, const int& incX,
    const float& beta, float* Y, const int& incY, const int tlen);
void zgbmv_(const char& TransA, const int& M, const int& N,
    const int& KL, const int& KU, const cdouble* alpha,
    const cdouble* A, const int& lda, const cdouble* X, const int& incX,
    const cdouble* beta, cdouble* Y, const int& incY, const int tlen);
void cgbmv_(const char& TransA, const int& M, const int& N,
    const int& KL, const int& KU, const cfloat* alpha,
    const cfloat* A, const int& lda, const cfloat* X, const int& incX, 
    const cfloat* beta, cfloat* Y, const int& incY, const int tlen);

void dtrmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const double* A, const int& lda,
    double* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void strmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const float* A, const int& lda,
    float* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ztrmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const cdouble* A, const int& lda,
    cdouble* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ctrmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const cfloat* A, const int& lda,
    cfloat* X, const int& incX,
    const int ulen, const int tlen, const int dlen);

void dtbmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const int& K, const double* A, const int& lda,
    double* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void stbmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const int& K, const float* A, const int& lda,
    float* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ztbmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const int& K, const cdouble* A, const int& lda,
    cdouble* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ctbmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const int& K, const cfloat* A, const int& lda,
    cfloat* X, const int& incX,
    const int ulen, const int tlen, const int dlen);

void dtpmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const double* Ap, double* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void stpmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const float* Ap, float* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ztpmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const cdouble* Ap, cdouble* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ctpmv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const cfloat* Ap, cfloat* X, const int& incX,
    const int ulen, const int tlen, const int dlen);

void dtrsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const double* A, const int& lda,
    double* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void strsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const float* A, const int& lda,
    float* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ztrsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const cdouble* A, const int& lda,
    cdouble* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ctrsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const cfloat* A, const int& lda,
    cfloat* X, const int& incX,
    const int ulen, const int tlen, const int dlen);

void dtbsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const int& K, const double* A, const int& lda,
    double* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void stbsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const int& K, const float* A, const int& lda,
    float* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ztbsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const int& K, const cdouble* A, const int& lda,
    cdouble* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ctbsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const int& K, const cfloat* A, const int& lda,
    cfloat* X, const int& incX,
    const int ulen, const int tlen, const int dlen);

void dtpsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const double* Ap, double* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void stpsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const float* Ap, float* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ztpsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const cdouble* Ap, cdouble* X, const int& incX,
    const int ulen, const int tlen, const int dlen);
void ctpsv_(const char& Uplo, const char& TransA, const char& Diag,
    const int& N, const cfloat* Ap, cfloat* X, const int& incX,
    const int ulen, const int tlen, const int dlen);

void dsymv_(const char& Uplo, const int& N, const double& alpha,
    const double* A, const int& lda, const double* X, const int& incX,
    const double& beta, double* Y, const int& incY, const int ulen);
void ssymv_(const char& Uplo, const int& N, const float& alpha,
    const float* A, const int& lda, const float* X, const int& incX,
    const float& beta, float* Y, const int& incY, const int ulen);
void zhemv_(const char& Uplo, const int& N, const cdouble* alpha,
    const cdouble* A, const int& lda, const cdouble* X, const int& incX,
    const cdouble* beta, cdouble* Y, const int& incY, const int ulen);
void chemv_(const char& Uplo, const int& N, const cfloat* alpha, 
    const cfloat* A, const int& lda, const cfloat* X, const int& incX,
    const cfloat* beta, cfloat* Y, const int& incY, const int ulen);

void dsbmv_(const char& Uplo, const int& N, const int& K,
    const double& alpha, const double* A, const int& lda,
    const double* X, const int& incX,
    const double& beta, double* Y, const int& incY, const int ulen);
void ssbmv_(const char& Uplo, const int& N, const int& K,
    const float& alpha, const float* A, const int& lda,
    const float* X, const int& incX,
    const float& beta, float* Y, const int& incY, const int ulen);
void zhbmv_(const char& Uplo, const int& N, const int& K,
    const cdouble* alpha, const cdouble* A, const int& lda,
    const cdouble* X, const int& incX, const cdouble* beta,
    cdouble* Y, const int& incY, const int ulen);
void chbmv_(const char& Uplo, const int& N, const int& K,
    const cfloat* alpha, const cfloat* A, const int& lda,
    const cfloat* X, const int& incX,
    const cfloat* beta, cfloat* Y, const int& incY, const int ulen);

void dspmv_(const char& Uplo, const int& N, const double& alpha,
    const double* Ap, const double* X, const int& incX,
    const double& beta, double* Y, const int& incY, const int ulen);
void sspmv_(const char& Uplo, const int& N, const float& alpha,
    const float* Ap, const float* X, const int& incX,
    const float& beta, float* Y, const int& incY, const int ulen);
void zhpmv_(const char& Uplo, const int& N, const cdouble* alpha,
    const cdouble* Ap, const cdouble* X, const int& incX,
    const cdouble* beta, cdouble* Y, const int& incY, const int ulen);
void chpmv_(const char& Uplo, const int& N, const cfloat* alpha,
    const cfloat* Ap, const cfloat* X, const int& incX,
    const cfloat* beta, cfloat* Y, const int& incY, const int ulen);

void dger_(const int& M, const int& N, const double& alpha,
    const double* X, const int& incX,
    const double* Y, const int& incY, double* A, const int& lda);
void sger_(const int& M, const int& N,
    const float& alpha, const float* X, const int& incX,
    const float* Y, const int& incY, float* A, const int& lda);
void zgeru_(const int& M, const int& N, const cdouble* alpha,
    const cdouble* X, const int& incX, const cdouble* Y, const int& incY,
    cdouble* A, const int& lda);
void zgerc_(const int& M, const int& N, const cdouble* alpha,
    const cdouble* X, const int& incX, const cdouble* Y, const int& incY,
    cdouble* A, const int& lda);
void cgeru_(const int& M, const int& N, const cfloat* alpha,
    const cfloat* X, const int& incX, const cfloat* Y, const int& incY,
    cfloat* A, const int& lda);
void cgerc_(const int& M, const int& N, const cfloat* alpha,
    const cfloat* X, const int& incX, const cfloat* Y, const int& incY,
    cfloat* A, const int& lda);

void dsyr_(const char& Uplo, const int& N, const double& alpha,
    const double* X, const int& incX,
    double* A, const int& lda, const int ulen);
void ssyr_(const char& Uplo, const int& N, const float& alpha,
    const float* X, const int& incX,
    float* A, const int& lda, const int ulen);
void zher_(const char& Uplo, const int& N, const double& alpha,
    const cdouble* X, const int& incX,
    cdouble* A, const int& lda, const int ulen);
void cher_(const char& Uplo, const int& N, const float& alpha,
    const cfloat* X, const int& incX,
    cfloat* A, const int& lda, const int ulen);

void dspr_(const char& Uplo, const int& N, const double& alpha,
    const double* X, const int& incX, double* Ap, const int ulen);
void sspr_(const char& Uplo, const int& N, const float& alpha,
    const float* X, const int& incX, float* Ap, const int ulen);
void zhpr_(const char& Uplo, const int& N, const double& alpha,
    const cdouble* X, const int& incX, cdouble* A, const int ulen);
void chpr_(const char& Uplo, const int& N, const float& alpha,
    const cfloat* X, const int& incX, cfloat* A, const int ulen);

void dsyr2_(const char& Uplo, const int& N, const double& alpha,
    const double* X, const int& incX, const double* Y, const int& incY,
    double* A, const int& lda, const int ulen);
void ssyr2_(const char& Uplo, const int& N, const float& alpha,
    const float* X, const int& incX, const float* Y, const int& incY,
    float* A, const int& lda, const int ulen);
void zher2_(const char& Uplo, const int& N,
    const cdouble* alpha, const cdouble* X, const int& incX,
    const cdouble* Y, const int& incY,
    cdouble* A, const int& lda, const int ulen);
void cher2_(const char& Uplo, const int& N, const cfloat* alpha,
    const cfloat* X, const int& incX, const cfloat* Y, const int& incY,
    cfloat* A, const int& lda, const int ulen);

void dspr2_(const char& Uplo, const int& N, const double& alpha,
    const double* X, const int& incX,
    const double* Y, const int& incY, double* A, const int ulen);
void sspr2_(const char& Uplo, const int& N, const float& alpha,
    const float* X, const int& incX,
    const float* Y, const int& incY, float* A, const int ulen);
void zhpr2_(const char& Uplo, const int& N, const cdouble* alpha,
    const cdouble* X, const int& incX, const cdouble* Y, const int& incY,
    cdouble* Ap, const int ulen);
void chpr2_(const char& Uplo, const int& N, const cfloat* alpha,
    const cfloat* X, const int& incX, const cfloat* Y, const int& incY,
    cfloat* Ap, const int ulen);


// ============
// Level 3 BLAS
// ============

void dgemm_(const char& TransA, const char& TransB,
    const int& M, const int& N, const int& K, const double& alpha, 
    const double* A, const int& lda, const double* B, const int& ldb,
    const double& beta, double* C, const int& ldc,
    const int talen, const int tblen);
void sgemm_(const char& TransA, const char& TransB,
    const int& M, const int& N, const int& K,
    const float& alpha, const float* A, const int& lda,
    const float* B, const int& ldb, const float& beta,
    float* C, const int& ldc, const int talen, const int tblen);
void zgemm_(const char& TransA, const char& TransB,
    const int& M, const int& N, const int& K,
    const cdouble* alpha, const cdouble* A, const int& lda,
    const cdouble* B, const int& ldb, const cdouble* beta,
    cdouble* C, const int& ldc, const int talen, const int tblen);
void cgemm_(const char& TransA, const char& TransB,
    const int& M, const int& N, const int& K,
    const cfloat* alpha, const cfloat* A, const int& lda,
    const cfloat* B, const int& ldb, const cfloat* beta,
    cfloat* C, const int& ldc, const int talen, const int tblen);

void dsymm_(const char& Side, const char& Uplo, const int& M, const int& N,
    const double& alpha, const double* A, const int& lda,
    const double* B, const int& ldb, const double& beta,
    double* C, const int& ldc, const int slen, const int ulen);
void ssymm_(const char& Side, const char& Uplo, const int& M, const int& N,
    const float& alpha, const float* A, const int& lda,
    const float* B, const int& ldb, const float& beta,
    float* C, const int& ldc, const int slen, const int ulen);
void zsymm_(const char& Side, const char& Uplo, const int& M, const int& N,
    const cdouble* alpha, const cdouble* A, const int& lda,
    const cdouble* B, const int& ldb, const cdouble* beta,
    cdouble* C, const int& ldc, const int slen, const int ulen);
void zhemm_(const char& Side, const char& Uplo, const int& M, const int& N,
    const cdouble* alpha, const cdouble* A, const int& lda,
    const cdouble* B, const int& ldb, const cdouble* beta,
    cdouble* C, const int& ldc, const int slen, const int ulen);
void csymm_(const char& Side, const char& Uplo, const int& M, const int& N,
    const cfloat* alpha, const cfloat* A, const int& lda,
    const cfloat* B, const int& ldb, const cfloat* beta,
    cfloat* C, const int& ldc, const int slen, const int ulen);
void chemm_(const char& Side, const char& Uplo, const int& M, const int& N,
    const cfloat* alpha, const cfloat* A, const int& lda,
    const cfloat* B, const int& ldb, const cfloat* beta,
    cfloat* C, const int& ldc, const int slen, const int ulen);

void dsyrk_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const double& alpha, const double* A, const int& lda,
    const double& beta, double* C, const int& ldc,
    const int ulen, const int tlen);
void ssyrk_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const float& alpha, const float* A, const int& lda,
    const float& beta, float* C, const int& ldc,
    const int ulen, const int tlen);
void zsyrk_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const cdouble* alpha, const cdouble* A, const int& lda,
    const cdouble* beta, cdouble* C, const int& ldc,
    const int ulen, const int tlen);
void zherk_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const double& alpha, const cdouble* A, const int& lda,
    const double& beta, cdouble* C, const int& ldc,
    const int ulen, const int tlen);
void csyrk_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const cfloat* alpha, const cfloat* A, const int& lda,
    const cfloat* beta, cfloat* C, const int& ldc,
    const int ulen, const int tlen);
void cherk_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const float& alpha, const cfloat* A, const int& lda,
    const float& beta, cfloat* C, const int& ldc,
    const int ulen, const int tlen);

void dsyr2k_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const double& alpha, const double* A, const int& lda,
    const double* B, const int& ldb, const double& beta,
    double* C, const int& ldc, const int ulen, const int tlen);
void ssyr2k_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const float& alpha, const float* A, const int& lda,
    const float* B, const int& ldb, const float& beta,
    float* C, const int& ldc, const int ulen, const int tlen);
void zsyr2k_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const cdouble* alpha, const cdouble* A, const int& lda,
    const cdouble* B, const int& ldb, const cdouble* beta,
    cdouble* C, const int& ldc, const int ulen, const int tlen);
void zher2k_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const cdouble* alpha, const cdouble* A, const int& lda,
    const cdouble* B, const int& ldb, const double& beta,
    cdouble* C, const int& ldc, const int ulen, const int tlen);
void csyr2k_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const cfloat* alpha, const cfloat* A, const int& lda,
    const cfloat* B, const int& ldb, const cfloat* beta,
    cfloat* C, const int& ldc, const int ulen, const int tlen);
void cher2k_(const char& Uplo, const char& Trans, const int& N, const int& K,
    const cfloat* alpha, const cfloat* A, const int& lda,
    const cfloat* B, const int& ldb, const float& beta,
    cfloat* C, const int& ldc, const int ulen, const int tlen);

void dtrmm_(const char& Side, const char& Uplo, const char& TransA,
    const char& Diag, const int& M, const int& N,
    const double& alpha, const double* A, const int& lda,
    double* B, const int& ldb,
    const int slen, const int ulen, const int tlen, const int dlen);
void strmm_(const char& Side, const char& Uplo, const char& TransA,
    const char& Diag, const int& M, const int& N,
    const float& alpha, const float* A, const int& lda,
    float* B, const int& ldb,
    const int slen, const int ulen, const int tlen, const int dlen);
void ztrmm_(const char& Side, const char& Uplo, const char& TransA,
    const char& Diag, const int& M, const int& N,
    const cdouble* alpha, const cdouble* A, const int& lda,
    cdouble* B, const int& ldb,
    const int slen, const int ulen, const int tlen, const int dlen);
void ctrmm_(const char& Side, const char& Uplo, const char& TransA,
    const char& Diag, const int& M, const int& N,
    const cfloat* alpha, const cfloat* A, const int& lda,
    cfloat* B, const int& ldb,
    const int slen, const int ulen, const int tlen, const int dlen);

void dtrsm_(const char& Side, const char& Uplo, const char& TransA,
    const char& Diag, const int& M, const int& N,
    const double& alpha, const double* A, const int& lda,
    double* B, const int& ldb,
    const int slen, const int ulen, const int tlen, const int dlen);
void strsm_(const char& Side, const char& Uplo, const char& TransA,
    const char& Diag, const int& M, const int& N,
    const float& alpha, const float* A, const int& lda,
    float* B, const int& ldb,
    const int slen, const int ulen, const int tlen, const int dlen);
void ztrsm_(const char& Side, const char& Uplo, const char& TransA,
    const char& Diag, const int& M, const int& N,
    const cdouble* alpha, const cdouble* A, const int& lda,
    cdouble* B, const int& ldb,
    const int slen, const int ulen, const int tlen, const int dlen);
void ctrsm_(const char& Side, const char& Uplo, const char& TransA,
    const char& Diag, const int& M, const int& N,
    const cfloat* alpha, const cfloat* A, const int& lda,
    cfloat* B, const int& ldb,
    const int slen, const int ulen, const int tlen, const int dlen);

#endif 
