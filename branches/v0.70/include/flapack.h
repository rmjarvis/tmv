#ifndef FLAPACK_H
#define FLAPACK_H

// Converted from the standard CLAPACK file clapack.h
// to work with the standard Fortran LAPACK library

#ifndef FBLAS_H
struct cdouble { double r; double i; };
struct cfloat { float r; float i; };
#endif

int dgbtrf_(const int& m, const int& n, const int& kl, const int& ku,
    double* ab, const int& ldab, int* ipiv, int* info);
int sgbtrf_(const int& m, const int& n, const int& kl, const int& ku,
    float* ab, const int& ldab, int* ipiv, int* info);
int zgbtrf_(const int& m, const int& n, const int& kl, const int& ku,
    cdouble* ab, const int& ldab, int* ipiv, int* info);
int cgbtrf_(const int& m, const int& n, const int& kl, const int& ku,
    cfloat* ab, const int& ldab, int* ipiv, int* info);

int dgttrf_(const int& n, double* dl, double* d, 
    double* du, double* du2, int* ipiv, int* info);
int sgttrf_(const int& n, float* dl, float* d, 
    float* du, float* du2, int* ipiv, int* info);
int zgttrf_(const int& n, cdouble* dl, cdouble*  d, 
    cdouble* du, cdouble* du2, int* ipiv, int*  info);
int cgttrf_(const int& n, cfloat* dl, cfloat* d, 
    cfloat* du, cfloat* du2, int* ipiv, int* info);

int dgbtrs_(const char& trans, const int& n, const int& kl, const int& ku,
    const int& nrhs, const double* ab, const int& ldab, int* ipiv, 
    double* b, const int& ldb, int* info);
int sgbtrs_(const char& trans, const int& n, const int& kl, const int& ku, 
    const int& nrhs, const float* ab, const int& ldab, int* ipiv, 
    float* b, const int& ldb, int* info);
int zgbtrs_(const char& trans, const int& n, const int& kl, const int& ku,
    const int& nrhs, const cdouble* ab, const int& ldab, int* ipiv, 
    cdouble* b, const int& ldb, int* info);
int cgbtrs_(const char& trans, const int& n, const int& kl, const int& ku,
    const int& nrhs, const cfloat* ab, const int& ldab, int* ipiv, 
    cfloat* b, const int& ldb, int* info);

int dgttrs_(const char& trans, const int& n, const int& nrhs, 
    const double* dl, const double* d, const double* du, const double* du2,
    int* ipiv, double* b, const int& ldb, int* info);
int sgttrs_(const char& trans, const int& n, const int& nrhs, 
    const float* dl, const float* d, const float* du, const float* du2, 
    const int* ipiv, float* b, const int& ldb, int* info);
int zgttrs_(const char& trans, const int& n, const int& nrhs, 
    const cdouble* dl, const cdouble* d, const cdouble* du, const cdouble* du2, 
    const int* ipiv, cdouble* b, const int& ldb, int* info);
int cgttrs_(const char& trans, const int& n, const int& nrhs, 
    const cfloat* dl, const cfloat* d, const cfloat* du, const cfloat* du2,
    const int* ipiv, cfloat* b, const int& ldb, int* info);

int dgbbrd_(const char& vect, const int& m, const int& n, const int& ncc,
    const int& kl, const int& ku, const double* ab, const int& ldab, 
    const double* d, const double* e, double* q, const int& ldq,
    double* pt, const int& ldpt, double* c, const int& ldc,
    double* work, int* info);
int sgbbrd_(const char& vect, const int& m, const int& n, const int& ncc,
    const int& kl, const int& ku, const float* ab, const int& ldab, 
    const float* d, const float* e, float* q, const int& ldq,
    float* pt, const int& ldpt, float* c, const int& ldc,
    float* work, int* info);
int zgbbrd_(const char& vect, const int& m, const int& n, const int& ncc,
    const int& kl, const int& ku, const cdouble* ab, const int& ldab, 
    const double* d, const double* e, cdouble* q, const int& ldq,
    cdouble* pt, const int& ldpt, cdouble* c, const int& ldc,
    cdouble* work, double* rwork, int* info);
int cgbbrd_(const char& vect, const int& m, const int& n, const int& ncc,
    const int& kl, const int& ku, const cfloat* ab, const int& ldab,
    const float* d,const  float* e, cfloat* q, const int& ldq,
    cfloat* pt, const int& ldpt, cfloat* c, const int& ldc, 
    cfloat* work, float* rwork, int* info);

int dtbtrs_(const char& uplo, const char& trans, const char& diag,
    const int& n, const int& kd, const int& nrhs, 
    const double* ab, const int& ldab, double* b, const int& ldb, int* info);
int stbtrs_(const char& uplo, const char& trans, const char& diag,
    const int& n, const int& kd, const int& nrhs,
    const float* ab, const int& ldab, float* b, const int& ldb, int* info);
int ztbtrs_(const char& uplo, const char& trans, const char& diag,
    const int& n, const int& kd, const int& nrhs,
    const cdouble* ab, const int& ldab, cdouble* b, const int& ldb, int* info);
int ctbtrs_(const char& uplo, const char& trans, const char& diag,
    const int& n, const int& kd, const int& nrhs,
    const cfloat* ab, const int& ldab, cfloat* b, const int& ldb, int* info);

int dorglq_(const int& m, const int& n, const int& k,
    double* a, const int& lda, const double* tau,
    double* work, const int& lwork, int* info);
int sorglq_(const int& m, const int& n, const int& k,
    float* a, const int& lda, const float* tau,
    float* work, const int& lwork, int* info);
int dorgqr_(const int& m, const int& n, const int& k,
    double* a, const int& lda, const double* tau,
    double* work, const int& lwork, int* info);
int sorgqr_(const int& m, const int& n, const int& k,
    float* a, const int& lda, const float* tau,
    float* work, const int& lwork, int* info);
int zunglq_(const int& m, const int& n, const int& k, 
    cdouble* a, const int& lda, const cdouble* tau,
    cdouble* work, const int& lwork, int* info);
int zungqr_(const int& m, const int& n, const int& k, 
    cdouble* a, const int& lda, const cdouble* tau,
    cdouble* work, const int& lwork, int* info);
int cunglq_(const int& m, const int& n, const int& k,
    cfloat* a, const int& lda, const cfloat* tau, 
    cfloat* work, const int& lwork, int* info);
int cungqr_(const int& m, const int& n, const int& k,
    cfloat* a, const int& lda, const cfloat* tau,
    cfloat* work, const int& lwork, int* info);

int zrot_(const int& n, cdouble* cx, const int& incx,
    cdouble* cy, const int& incy, const double& c, const cdouble* s);
int crot_(const int& n, cfloat* cx, const int& incx,
    cfloat* cy, const int& incy, const float& c, const cfloat* s);

int dgetrf_(const int& m, const int& n, double* a, const int& lda,
    int* ipiv, int* info);
int sgetrf_(const int& m, const int& n, float* a, const int& lda, 
    int* ipiv, int* info);
int zgetrf_(const int& m, const int& n, cdouble* a, const int& lda,
    int* ipiv, int* info);
int cgetrf_(const int& m, const int& n, cfloat* a, const int& lda,
    int* ipiv, int* info);

int dgetrs_(const char& trans, const int& n, const int& nrhs, const double* a,
    const int& lda, const int* ipiv, double* b, const int& ldb, int* info);
int sgetrs_(const char& trans, const int& n, const int& nrhs, const float* a, 
    const int& lda, const int* ipiv, float* b, const int& ldb, int* info);
int zgetrs_(const char& trans, const int& n, const int& nrhs, const cdouble* a,
    const int& lda, const int* ipiv, cdouble* b, const int& ldb, int* info);
int cgetrs_(const char& trans, const int& n, const int& nrhs, const cfloat* a,
    const int& lda, const int* ipiv, cfloat* b, const int& ldb, int* info);

int dgetri_(const int& n, double* a, const int& lda, const int* ipiv,
    double* work, const int& lwork, int* info);
int sgetri_(const int& n, float* a, const int& lda, const int* ipiv,
    float* work, const int& lwork, int* info);
int zgetri_(const int& n, cdouble* a, const int& lda, const int* ipiv,
    cdouble* work, const int& lwork, int* info);
int cgetri_(const int& n, cfloat* a, const int& lda, const int* ipiv,
    cfloat* work, const int& lwork, int* info);

int dlacpy_(const char& uplo, const int& m, const int& n, 
    const double* a, const int& lda, double* b, const int& ldb);
int slacpy_(const char& uplo, const int& m, const int& n,
    const float* a, const int& lda, float* b, const int& ldb);
int zlacpy_(const char& uplo, const int& m, const int& n,
    const cdouble* a, const int& lda, cdouble* b, const int& ldb);
int clacpy_(const char& uplo, const int& m, const int& n,
    const cfloat* a, const int& lda, cfloat* b, const int& ldb);

int dlagtm_(const char& trans, const int& n, const int& nrhs, 
    const double& alpha, const double* dl, const double* d,
    const double* du, const double* x, const int& ldx,
    const double& beta, double* b, const int& ldb);
int slagtm_(const char& trans, const int& n, const int& nrhs,
    const float& alpha, const float* dl, const float* d,
    const float* du, const float* x, const int& ldx,
    const float& beta, float* b, const int& ldb);
int zlagtm_(const char& trans, const int& n, const int& nrhs, 
    const double& alpha, const cdouble* dl, const cdouble* d,
    const cdouble* du, const cdouble* x, const int& ldx,
    const double& beta, cdouble* b, const int& ldb);
int clagtm_(const char& trans, const int& n, const int& nrhs, 
    const float& alpha, const cfloat* dl, const cfloat* d,
    const cfloat* du, const cfloat* x, const int& ldx,
    const float& beta, cfloat* b, const int& ldb);

int zsymv_(const char& uplo, const int& n, const cdouble* alpha, 
    const cdouble* a, const int& lda, const cdouble* x, const int& incx, 
    const cdouble* beta, cdouble* y, const int& incy);
int csymv_(const char& uplo, const int& n, const cfloat* alpha, 
    const cfloat* a, const int& lda, const cfloat* x, const int& incx,
    const cfloat* beta, cfloat* y, const int& incy);

int dgelqf_(const int& m, const int& n, double* a, const int& lda,
    double* tau, double* work, const int& lwork, int* info);
int sgelqf_(const int& m, const int& n, float* a, const int& lda, 
    float* tau, float* work, const int& lwork, int* info);
int zgelqf_(const int& m, const int& n, cdouble* a, const int& lda,
    cdouble* tau, cdouble* work, const int& lwork, int* info);
int cgelqf_(const int& m, const int& n, cfloat* a, const int& lda,
    cfloat* tau, cfloat* work, const int& lwork, int* info);
int dgeqrf_(const int& m, const int& n, double* a, const int& lda,
    double* tau, double* work, const int& lwork, int* info);
int sgeqrf_(const int& m, const int& n, float* a, const int& lda, 
    float* tau, float* work, const int& lwork, int* info);
int zgeqrf_(const int& m, const int& n, cdouble* a, const int& lda,
    cdouble* tau, cdouble* work, const int& lwork, int* info);
int cgeqrf_(const int& m, const int& n, cfloat* a, const int& lda,
    cfloat* tau, cfloat* work, const int& lwork, int* info);

int dormlq_(const char& side, const char& trans, const int& m, const int& n,
    const int& k, const double* a, const int& lda, const double* tau,
    double* c, const int& ldc, double* work, const int& lwork, int* info);
int sormlq_(const char& side, const char& trans, const int& m, const int& n, 
    const int& k, const float* a, const int& lda, const float* tau,
    float* c, const int& ldc, float* work, const int& lwork, int* info);
int dormqr_(const char& side, const char& trans, const int& m, const int& n, 
    const int& k, const double* a, const int& lda, const double* tau,
    double* c, const int& ldc, double* work, const int& lwork, int* info);
int sormqr_(const char& side, const char& trans, const int& m, const int& n, 
    const int& k, const float* a, const int& lda, const float* tau,
    float* c, const int& ldc, float* work, const int& lwork, int* info);
int zunmlq_(const char& side, const char& trans, const int& m, const int& n, 
    const int& k, const cdouble* a, const int& lda, const cdouble* tau, 
    cdouble* c, const int& ldc, cdouble* work, const int& lwork, int* info);
int cunmlq_(const char& side, const char& trans, const int& m, const int& n, 
    const int& k, const cfloat* a, const int& lda, const cfloat* tau, 
    cfloat* c, const int& ldc, cfloat* work, const int& lwork, int* info);
int zunmqr_(const char& side, const char& trans, const int& m, const int& n, 
    const int& k, const cdouble* a, const int& lda, const cdouble* tau,
    cdouble* c, const int& ldc, cdouble* work, const int& lwork, int* info);
int cunmqr_(const char& side, const char& trans, const int& m, const int& n, 
    const int& k, const cfloat* a, const int& lda, const cfloat* tau, 
    cfloat* c, const int& ldc, cfloat* work, const int& lwork, int* info);

int dgeqp3_(const int& m, const int& n, double* a, const int& lda, int* jpvt,
    double* tau, double* work, const int& lwork, int* info);
int sgeqp3_(const int& m, const int& n, float* a, const int& lda, int* jpvt,
    float* tau, float* work, const int& lwork, int* info);
int zgeqp3_(const int& m, const int& n, cdouble* a, const int& lda, int* jpvt,
    cdouble* tau, cdouble* work, const int& lwork, double* rwork, int* info);
int cgeqp3_(const int& m, const int& n, cfloat* a, const int& lda, int* jpvt,
    cfloat* tau, cfloat* work, const int& lwork, float* rwork, int* info);

int zsyr_(const char& uplo, const int& n, const cdouble* alpha, 
    const cdouble* x, const int& incx, cdouble* a, const int& lda);
int csyr_(const char& uplo, const int& n, const cfloat* alpha, 
    const cfloat* x, const int& incx, cfloat* a, const int& lda);

int dlauum_(const char& uplo, const int& n, double* a, const int& lda,
    int* info);
int slauum_(const char& uplo, const int& n, float* a, const int& lda, 
    int* info);
int zlauum_(const char& uplo, const int& n, cdouble* a, const int& lda,
    int* info);
int clauum_(const char& uplo, const int& n, cfloat* a, const int& lda,
    int* info);

int dbdsdc_(const char& uplo, const char& compq, const int& n, double* d,
    double* e, double* u, const int& ldu, double* vt, const int& ldvt,
    double* q, int* iq, double* work, int* iwork, int* info);
int sbdsdc_(const char& uplo, const char& compq, const int& n, float* d, 
    float* e, float* u, const int& ldu, float* vt, const int& ldvt,
    float* q, int* iq, float* work, int* iwork, int* info);

int dgebrd_(const int& m, const int& n, double* a, const int& lda,
    double* d, double* e, double* tauq, double* taup,
    double* work, const int& lwork, int* info);
int sgebrd_(const int& m, const int& n, float* a, const int& lda, 
    float* d, float* e, float* tauq, float* taup, 
    float* work, const int& lwork, int* info);
int zgebrd_(const int& m, const int& n, cdouble* a, const int& lda,
    double* d, double* e, cdouble* tauq, cdouble* taup,
    cdouble* work, const int& lwork, int* info);
int cgebrd_(const int& m, const int& n, cfloat* a, const int& lda,
    float* d, float* e, cfloat* tauq, cfloat* taup,
    cfloat* work, const int& lwork, int* info);

int dpbtrf_(const char& uplo, const int& n, const int& kd, double* ab,
    const int& ldab, int* info);
int spbtrf_(const char& uplo, const int& n, const int& kd, float* ab, 
    const int& ldab, int* info);
int zpbtrf_(const char& uplo, const int& n, const int& kd, cdouble* ab,
    const int& ldab, int* info);
int cpbtrf_(const char& uplo, const int& n, const int& kd, cfloat* ab,
    const int& ldab, int* info);

int dpttrf_(const int& n, double* d, double* e, int* info);
int spttrf_(const int& n, float* d, float* e, int* info);
int zpttrf_(const int& n, double* d, cdouble* e, int* info);
int cpttrf_(const int& n, float* d, cfloat* e, int* info);

int dsbtrd_(const char& vect, const char& uplo, const int& n, const int& kd, 
    const double* ab, const int& ldab, double* d, double* e, 
    double* q, const int& ldq, double* work, int* info);
int ssbtrd_(const char& vect, const char& uplo, const int& n, const int& kd, 
    const float* ab, const int& ldab, float* d, float* e,
    float* q, const int& ldq, float* work, int* info);
int zhbtrd_(const char& vect, const char& uplo, const int& n, const int& kd, 
    const cdouble* ab, const int& ldab, double* d, double* e, 
    cdouble* q, const int& ldq, cdouble* work, int* info);
int chbtrd_(const char& vect, const char& uplo, const int& n, const int& kd, 
    const cfloat* ab, const int& ldab, float* d, float* e,
    cfloat* q, const int& ldq, cfloat* work, int* info);

int dpotrf_(const char& uplo, const int& n, double* a, const int& lda,
    int* info);
int spotrf_(const char& uplo, const int& n, float* a, const int& lda, 
    int* info);
int zpotrf_(const char& uplo, const int& n, cdouble* a, const int& lda,
    int* info);
int cpotrf_(const char& uplo, const int& n, cfloat* a, const int& lda,
    int* info);

int dpotri_(const char& uplo, const int& n, double* a, const int& lda,
    int* info);
int spotri_(const char& uplo, const int& n, float* a, const int& lda, 
    int* info);
int zpotri_(const char& uplo, const int& n, cdouble* a, const int& lda,
    int* info);
int cpotri_(const char& uplo, const int& n, cfloat* a, const int& lda,
    int* info);

int dsytrf_(const char& uplo, const int& n, double* a, const int& lda,
    int* ipiv, double* work, const int& lwork, int* info);
int ssytrf_(const char& uplo, const int& n, float* a, const int& lda, 
    int* ipiv, float* work, const int& lwork, int* info);
int zsytrf_(const char& uplo, const int& n, cdouble* a, const int& lda,
    int* ipiv, cdouble* work, const int& lwork, int* info);
int zhetrf_(const char& uplo, const int& n, cdouble* a, const int& lda,
    int* ipiv, cdouble* work, const int& lwork, int* info);
int csytrf_(const char& uplo, const int& n, cfloat* a, const int& lda,
    int* ipiv, cfloat* work, const int& lwork, int* info);
int chetrf_(const char& uplo, const int& n, cfloat* a, const int& lda,
    int* ipiv, cfloat* work, const int& lwork, int* info);

int dstedc_(const char& compz, const int& n, double* d, double* e,
    double* z, const int& ldz, double* work, const int& lwork, 
    int* iwork, const int& liwork, int* info);
int sstedc_(const char& compz, const int& n, float* d, float* e, 
    float* z, const int& ldz, float* work, const int& lwork,
    int* iwork, const int& liwork, int* info);
int zstedc_(const char& compz, const int& n, double* d, double* e,
    cdouble* z, const int& ldz, cdouble* work, const int& lwork,
    double* rwork, const int& lrwork, int* iwork, const int& liwork, int* info);
int cstedc_(const char& compz, const int& n, float* d, float* e, 
    cfloat* z, const int& ldz, cfloat* work, const int& lwork, 
    float* rwork, const int& lrwork, int* iwork, const int& liwork, int* info);

int dstegr_(const char& jobz, const char& range, const int& n, 
    double* d, double* e, const double& vl, const double& vu, 
    const int& il, const int& iu, const double& abstol, int* m,
    double* w, double* z, const int& ldz, int* isuppz, double* work,
    const int& lwork, int* iwork, const int& liwork, int* info);
int sstegr_(const char& jobz, const char& range, const int& n,
    float* d, float* e, const float& vl, const float& vu,
    const int& il, const int& iu, const float& abstol, int* m, 
    float* w, float* z, const int& ldz, int* isuppz, float* work, 
    const int& lwork, int* iwork, const int& liwork, int* info);

int dsytrd_(const char& uplo, const int& n, double* a, const int& lda,
    double* d, double* e, double* tau, double* work, const int& lwork,
    int* info);
int ssytrd_(const char& uplo, const int& n, float* a, const int& lda, 
    float* d, float* e, float* tau, float* work, const int& lwork,
    int* info);
int zhetrd_(const char& uplo, const int& n, cdouble* a, const int& lda,
    double* d, double* e, cdouble* tau, cdouble* work, const int& lwork,
    int* info);
int chetrd_(const char& uplo, const int& n, cfloat* a, const int& lda,
    float* d, float* e, cfloat* tau, cfloat* work, const int& lwork, 
    int* info);

int dtrtri_(const char& uplo, const char& diag, const int& n, 
    double* a, const int& lda, int* info);
int strtri_(const char& uplo, const char& diag, const int& n,
    float* a, const int& lda, int* info);
int ztrtri_(const char& uplo, const char& diag, const int& n, 
    cdouble* a, const int& lda, int* info);
int ctrtri_(const char& uplo, const char& diag, const int& n,
    cfloat* a, const int& lda, int* info);

int zlacgv_(const int& n, cdouble* x, const int& incx);
int clacgv_(const int& n, cfloat* x, const int& incx);

// The rest are not used by TMV, so I haven't bothered to convert them
// to the correct const-ness and such.
// But I leave them here for completeness and in case I eventaully
// want to promote some of them to active duty.

int cbdsqr_(char *uplo, int *n, int *ncvt, int *
    nru, int *ncc, float *d__, float *e, void *vt, int *ldvt, 
    void *u, int *ldu, void *c__, int *ldc, float *rwork, 
    int *info);


int cgbcon_(char *norm, int *n, int *kl, int *ku,
    void *ab, int *ldab, int *ipiv, float *anorm, float *rcond, 
    void *work, float *rwork, int *info);

int cgbequ_(int *m, int *n, int *kl, int *ku,
    void *ab, int *ldab, float *r__, float *c__, float *rowcnd, float 
    *colcnd, float *amax, int *info);

int cgbrfs_(char *trans, int *n, int *kl, int *
    ku, int *nrhs, void *ab, int *ldab, void *afb, int *
    ldafb, int *ipiv, void *b, int *ldb, void *x, int *
    ldx, float *ferr, float *berr, void *work, float *rwork, int *
    info);

int cgbsv_(int *n, int *kl, int *ku, int *
    nrhs, void *ab, int *ldab, int *ipiv, void *b, int *
    ldb, int *info);

int cgbsvx_(char *fact, char *trans, int *n, int *kl,
    int *ku, int *nrhs, void *ab, int *ldab, void *afb,
    int *ldafb, int *ipiv, char *equed, float *r__, float *c__, 
    void *b, int *ldb, void *x, int *ldx, float *rcond, float 
    *ferr, float *berr, void *work, float *rwork, int *info);

int cgbtf2_(int *m, int *n, int *kl, int *ku,
    void *ab, int *ldab, int *ipiv, int *info);



int cgebak_(char *job, char *side, int *n, int *ilo, 
    int *ihi, float *scale, int *m, void *v, int *ldv, 
    int *info);

int cgebal_(char *job, int *n, void *a, int *lda, 
    int *ilo, int *ihi, float *scale, int *info);

int cgebd2_(int *m, int *n, void *a, int *lda,
    float *d__, float *e, void *tauq, void *taup, void *work, 
    int *info);


int cgecon_(char *norm, int *n, void *a, int *lda,
    float *anorm, float *rcond, void *work, float *rwork, int *info);

int cgeequ_(int *m, int *n, void *a, int *lda,
    float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, 
    int *info);

int cgeev_(char *jobvl, char *jobvr, int *n, void *a, 
    int *lda, void *w, void *vl, int *ldvl, void *vr, 
    int *ldvr, void *work, int *lwork, float *rwork, int *
    info);

int cgeevx_(char *balanc, char *jobvl, char *jobvr, char *
    sense, int *n, void *a, int *lda, void *w, void *vl, 
    int *ldvl, void *vr, int *ldvr, int *ilo, int *ihi,
    float *scale, float *abnrm, float *rconde, float *rcondv, void *work, 
    int *lwork, float *rwork, int *info);

int cgegs_(char *jobvsl, char *jobvsr, int *n, void *
    a, int *lda, void *b, int *ldb, void *alpha, void *
    beta, void *vsl, int *ldvsl, void *vsr, int *ldvsr, 
    void *work, int *lwork, float *rwork, int *info);

int cgegv_(char *jobvl, char *jobvr, int *n, void *a, 
    int *lda, void *b, int *ldb, void *alpha, void *beta,
    void *vl, int *ldvl, void *vr, int *ldvr, void *
    work, int *lwork, float *rwork, int *info);

int cgehd2_(int *n, int *ilo, int *ihi, void *
    a, int *lda, void *tau, void *work, int *info);

int cgehrd_(int *n, int *ilo, int *ihi, void *
    a, int *lda, void *tau, void *work, int *lwork, int 
    *info);

int cgelq2_(int *m, int *n, void *a, int *lda,
    void *tau, void *work, int *info);


int cgels_(char *trans, int *m, int *n, int *
    nrhs, void *a, int *lda, void *b, int *ldb, void *
    work, int *lwork, int *info);

int cgelsx_(int *m, int *n, int *nrhs, void *
    a, int *lda, void *b, int *ldb, int *jpvt, float *rcond,
    int *rank, void *work, float *rwork, int *info);

int cgelsy_(int *m, int *n, int *nrhs, void *
    a, int *lda, void *b, int *ldb, int *jpvt, float *rcond,
    int *rank, void *work, int *lwork, float *rwork, int *
    info);

int cgeql2_(int *m, int *n, void *a, int *lda,
    void *tau, void *work, int *info);

int cgeqlf_(int *m, int *n, void *a, int *lda,
    void *tau, void *work, int *lwork, int *info);


int cgeqpf_(int *m, int *n, void *a, int *lda,
    int *jpvt, void *tau, void *work, float *rwork, int *
    info);

int cgeqr2_(int *m, int *n, void *a, int *lda,
    void *tau, void *work, int *info);


int cgerfs_(char *trans, int *n, int *nrhs, void *
    a, int *lda, void *af, int *ldaf, int *ipiv, void *
    b, int *ldb, void *x, int *ldx, float *ferr, float *berr, 
    void *work, float *rwork, int *info);

int cgerq2_(int *m, int *n, void *a, int *lda,
    void *tau, void *work, int *info);

int cgerqf_(int *m, int *n, void *a, int *lda,
    void *tau, void *work, int *lwork, int *info);

int cgesc2_(int *n, void *a, int *lda, void *
    rhs, int *ipiv, int *jpiv, float *scale);

int cgesv_(int *n, int *nrhs, void *a, int *
    lda, int *ipiv, void *b, int *ldb, int *info);

int cgesvx_(char *fact, char *trans, int *n, int *
    nrhs, void *a, int *lda, void *af, int *ldaf, int *
    ipiv, char *equed, float *r__, float *c__, void *b, int *ldb, 
    void *x, int *ldx, float *rcond, float *ferr, float *berr, 
    void *work, float *rwork, int *info);

int cgetc2_(int *n, void *a, int *lda, int *
    ipiv, int *jpiv, int *info);

int cgetf2_(int *m, int *n, void *a, int *lda,
    int *ipiv, int *info);




int cggbak_(char *job, char *side, int *n, int *ilo, 
    int *ihi, float *lscale, float *rscale, int *m, void *v, 
    int *ldv, int *info);

int cggbal_(char *job, int *n, void *a, int *lda, 
    void *b, int *ldb, int *ilo, int *ihi, float *lscale, 
    float *rscale, float *work, int *info);

int cggev_(char *jobvl, char *jobvr, int *n, void *a, 
    int *lda, void *b, int *ldb, void *alpha, void *beta,
    void *vl, int *ldvl, void *vr, int *ldvr, void *
    work, int *lwork, float *rwork, int *info);

int cggevx_(char *balanc, char *jobvl, char *jobvr, char *
    sense, int *n, void *a, int *lda, void *b, int *ldb,
    void *alpha, void *beta, void *vl, int *ldvl, void *
    vr, int *ldvr, int *ilo, int *ihi, float *lscale, float *
    rscale, float *abnrm, float *bbnrm, float *rconde, float *rcondv, void 
    *work, int *lwork, float *rwork, int *iwork, int *bwork, 
    int *info);

int cggglm_(int *n, int *m, int *p, void *a, 
    int *lda, void *b, int *ldb, void *d__, void *x, 
    void *y, void *work, int *lwork, int *info);

int cgghrd_(char *compq, char *compz, int *n, int *
    ilo, int *ihi, void *a, int *lda, void *b, int *ldb,
    void *q, int *ldq, void *z__, int *ldz, int *info);

int cgglse_(int *m, int *n, int *p, void *a, 
    int *lda, void *b, int *ldb, void *c__, void *d__, 
    void *x, void *work, int *lwork, int *info);

int cggqrf_(int *n, int *m, int *p, void *a, 
    int *lda, void *taua, void *b, int *ldb, void *taub, 
    void *work, int *lwork, int *info);

int cggrqf_(int *m, int *p, int *n, void *a, 
    int *lda, void *taua, void *b, int *ldb, void *taub, 
    void *work, int *lwork, int *info);

int cggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
    int *n, int *p, int *k, int *l, void *a, int *
    lda, void *b, int *ldb, float *alpha, float *beta, void *u, 
    int *ldu, void *v, int *ldv, void *q, int *ldq, 
    void *work, float *rwork, int *iwork, int *info);

int cggsvp_(char *jobu, char *jobv, char *jobq, int *m, 
    int *p, int *n, void *a, int *lda, void *b, int 
    *ldb, float *tola, float *tolb, int *k, int *l, void *u, 
    int *ldu, void *v, int *ldv, void *q, int *ldq, 
    int *iwork, float *rwork, void *tau, void *work, int *
    info);

int cgtcon_(char *norm, int *n, void *dl, void *
    d__, void *du, void *du2, int *ipiv, float *anorm, float *
    rcond, void *work, int *info);

int cgtrfs_(char *trans, int *n, int *nrhs, void *
    dl, void *d__, void *du, void *dlf, void *df, void *
    duf, void *du2, int *ipiv, void *b, int *ldb, void *
    x, int *ldx, float *ferr, float *berr, void *work, float *rwork, 
    int *info);

int cgtsv_(int *n, int *nrhs, void *dl, void *
    d__, void *du, void *b, int *ldb, int *info);

int cgtsvx_(char *fact, char *trans, int *n, int *
    nrhs, void *dl, void *d__, void *du, void *dlf, void *
    df, void *duf, void *du2, int *ipiv, void *b, int *
    ldb, void *x, int *ldx, float *rcond, float *ferr, float *berr, 
    void *work, float *rwork, int *info);



int cgtts2_(int *itrans, int *n, int *nrhs, 
    void *dl, void *d__, void *du, void *du2, int *ipiv, 
    void *b, int *ldb);

int chbev_(char *jobz, char *uplo, int *n, int *kd, 
    void *ab, int *ldab, float *w, void *z__, int *ldz, 
    void *work, float *rwork, int *info);

int chbevd_(char *jobz, char *uplo, int *n, int *kd, 
    void *ab, int *ldab, float *w, void *z__, int *ldz, 
    void *work, int *lwork, float *rwork, int *lrwork, int *
    iwork, int *liwork, int *info);

int chbevx_(char *jobz, char *range, char *uplo, int *n, 
    int *kd, void *ab, int *ldab, void *q, int *ldq, 
    float *vl, float *vu, int *il, int *iu, float *abstol, int *
    m, float *w, void *z__, int *ldz, void *work, float *rwork, 
    int *iwork, int *ifail, int *info);

int chbgst_(char *vect, char *uplo, int *n, int *ka, 
    int *kb, void *ab, int *ldab, void *bb, int *ldbb, 
    void *x, int *ldx, void *work, float *rwork, int *info);

int chbgv_(char *jobz, char *uplo, int *n, int *ka, 
    int *kb, void *ab, int *ldab, void *bb, int *ldbb, 
    float *w, void *z__, int *ldz, void *work, float *rwork, 
    int *info);

int chbgvx_(char *jobz, char *range, char *uplo, int *n, 
    int *ka, int *kb, void *ab, int *ldab, void *bb, 
    int *ldbb, void *q, int *ldq, float *vl, float *vu, int *
    il, int *iu, float *abstol, int *m, float *w, void *z__, 
    int *ldz, void *work, float *rwork, int *iwork, int *
    ifail, int *info);


int checon_(char *uplo, int *n, void *a, int *lda,
    int *ipiv, float *anorm, float *rcond, void *work, int *
    info);

int cheev_(char *jobz, char *uplo, int *n, void *a, 
    int *lda, float *w, void *work, int *lwork, float *rwork, 
    int *info);

int cheevd_(char *jobz, char *uplo, int *n, void *a, 
    int *lda, float *w, void *work, int *lwork, float *rwork, 
    int *lrwork, int *iwork, int *liwork, int *info);

int cheevr_(char *jobz, char *range, char *uplo, int *n, 
    void *a, int *lda, float *vl, float *vu, int *il, int *
    iu, float *abstol, int *m, float *w, void *z__, int *ldz, 
    int *isuppz, void *work, int *lwork, float *rwork, int *
    lrwork, int *iwork, int *liwork, int *info);

int cheevx_(char *jobz, char *range, char *uplo, int *n, 
    void *a, int *lda, float *vl, float *vu, int *il, int *
    iu, float *abstol, int *m, float *w, void *z__, int *ldz, 
    void *work, int *lwork, float *rwork, int *iwork, int *
    ifail, int *info);

int chegs2_(int *itype, char *uplo, int *n, void *
    a, int *lda, void *b, int *ldb, int *info);

int chegst_(int *itype, char *uplo, int *n, void *
    a, int *lda, void *b, int *ldb, int *info);

int chegv_(int *itype, char *jobz, char *uplo, int *
    n, void *a, int *lda, void *b, int *ldb, float *w, 
    void *work, int *lwork, float *rwork, int *info);

int chegvd_(int *itype, char *jobz, char *uplo, int *
    n, void *a, int *lda, void *b, int *ldb, float *w, 
    void *work, int *lwork, float *rwork, int *lrwork, int *
    iwork, int *liwork, int *info);

int chegvx_(int *itype, char *jobz, char *range, char *
    uplo, int *n, void *a, int *lda, void *b, int *ldb, 
    float *vl, float *vu, int *il, int *iu, float *abstol, int *
    m, float *w, void *z__, int *ldz, void *work, int *lwork,
    float *rwork, int *iwork, int *ifail, int *info);

int cherfs_(char *uplo, int *n, int *nrhs, void *
    a, int *lda, void *af, int *ldaf, int *ipiv, void *
    b, int *ldb, void *x, int *ldx, float *ferr, float *berr, 
    void *work, float *rwork, int *info);

int chesv_(char *uplo, int *n, int *nrhs, void *a,
    int *lda, int *ipiv, void *b, int *ldb, void *work,
    int *lwork, int *info);

int chesvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *a, int *lda, void *af, int *ldaf, int *
    ipiv, void *b, int *ldb, void *x, int *ldx, float *rcond,
    float *ferr, float *berr, void *work, int *lwork, float *rwork, 
    int *info);

int chetf2_(char *uplo, int *n, void *a, int *lda,
    int *ipiv, int *info);



int chetri_(char *uplo, int *n, void *a, int *lda,
    int *ipiv, void *work, int *info);

int chetrs_(char *uplo, int *n, int *nrhs, void *
    a, int *lda, int *ipiv, void *b, int *ldb, int *
    info);

int chgeqz_(char *job, char *compq, char *compz, int *n, 
    int *ilo, int *ihi, void *a, int *lda, void *b, 
    int *ldb, void *alpha, void *beta, void *q, int *ldq,
    void *z__, int *ldz, void *work, int *lwork, float *
    rwork, int *info);

int chpcon_(char *uplo, int *n, void *ap, int *
    ipiv, float *anorm, float *rcond, void *work, int *info);

int chpev_(char *jobz, char *uplo, int *n, void *ap, 
    float *w, void *z__, int *ldz, void *work, float *rwork, 
    int *info);

int chpevd_(char *jobz, char *uplo, int *n, void *ap, 
    float *w, void *z__, int *ldz, void *work, int *lwork, 
    float *rwork, int *lrwork, int *iwork, int *liwork, 
    int *info);

int chpevx_(char *jobz, char *range, char *uplo, int *n, 
    void *ap, float *vl, float *vu, int *il, int *iu, float *
    abstol, int *m, float *w, void *z__, int *ldz, void *
    work, float *rwork, int *iwork, int *ifail, int *info);

int chpgst_(int *itype, char *uplo, int *n, void *
    ap, void *bp, int *info);

int chpgv_(int *itype, char *jobz, char *uplo, int *
    n, void *ap, void *bp, float *w, void *z__, int *ldz, 
    void *work, float *rwork, int *info);

int chpgvd_(int *itype, char *jobz, char *uplo, int *
    n, void *ap, void *bp, float *w, void *z__, int *ldz, 
    void *work, int *lwork, float *rwork, int *lrwork, int *
    iwork, int *liwork, int *info);

int chpgvx_(int *itype, char *jobz, char *range, char *
    uplo, int *n, void *ap, void *bp, float *vl, float *vu, 
    int *il, int *iu, float *abstol, int *m, float *w, void *
    z__, int *ldz, void *work, float *rwork, int *iwork, 
    int *ifail, int *info);

int chprfs_(char *uplo, int *n, int *nrhs, void *
    ap, void *afp, int *ipiv, void *b, int *ldb, void *x,
    int *ldx, float *ferr, float *berr, void *work, float *rwork, 
    int *info);

int chpsv_(char *uplo, int *n, int *nrhs, void *
    ap, int *ipiv, void *b, int *ldb, int *info);

int chpsvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *ap, void *afp, int *ipiv, void *b, int *
    ldb, void *x, int *ldx, float *rcond, float *ferr, float *berr, 
    void *work, float *rwork, int *info);

int chptrd_(char *uplo, int *n, void *ap, float *d__, 
    float *e, void *tau, int *info);

int chptrf_(char *uplo, int *n, void *ap, int *
    ipiv, int *info);

int chptri_(char *uplo, int *n, void *ap, int *
    ipiv, void *work, int *info);

int chptrs_(char *uplo, int *n, int *nrhs, void *
    ap, int *ipiv, void *b, int *ldb, int *info);

int chsein_(char *side, char *eigsrc, char *initv, int *
    select, int *n, void *h__, int *ldh, void *w, void *
    vl, int *ldvl, void *vr, int *ldvr, int *mm, int *
    m, void *work, float *rwork, int *ifaill, int *ifailr, 
    int *info);

int chseqr_(char *job, char *compz, int *n, int *ilo,
    int *ihi, void *h__, int *ldh, void *w, void *z__, 
    int *ldz, void *work, int *lwork, int *info);

int clabrd_(int *m, int *n, int *nb, void *a, 
    int *lda, float *d__, float *e, void *tauq, void *taup, 
    void *x, int *ldx, void *y, int *ldy);


int clacon_(int *n, void *v, void *x, float *est, 
    int *kase);

int clacp2_(char *uplo, int *m, int *n, float *a, 
    int *lda, void *b, int *ldb);


int clacrm_(int *m, int *n, void *a, int *lda,
    float *b, int *ldb, void *c__, int *ldc, float *rwork);

int clacrt_(int *n, void *cx, int *incx, void *
    cy, int *incy, void *c__, void *s);

int claed0_(int *qsiz, int *n, float *d__, float *e, 
    void *q, int *ldq, void *qstore, int *ldqs, float *rwork,
    int *iwork, int *info);

int claed7_(int *n, int *cutpnt, int *qsiz, 
    int *tlvls, int *curlvl, int *curpbm, float *d__, void *
    q, int *ldq, float *rho, int *indxq, float *qstore, int *
    qptr, int *prmptr, int *perm, int *givptr, int *
    givcol, float *givnum, void *work, float *rwork, int *iwork, 
    int *info);

int claed8_(int *k, int *n, int *qsiz, void *
    q, int *ldq, float *d__, float *rho, int *cutpnt, float *z__, 
    float *dlamda, void *q2, int *ldq2, float *w, int *indxp, 
    int *indx, int *indxq, int *perm, int *givptr, 
    int *givcol, float *givnum, int *info);

int claein_(int *rightv, int *noinit, int *n, 
    void *h__, int *ldh, void *w, void *v, void *b, 
    int *ldb, float *rwork, float *eps3, float *smlnum, int *info);

int claesy_(void *a, void *b, void *c__, void *
    rt1, void *rt2, void *evscal, void *cs1, void *sn1);

int claev2_(void *a, void *b, void *c__, float *rt1, 
    float *rt2, float *cs1, void *sn1);

int clags2_(int *upper, float *a1, void *a2, float *a3, 
    float *b1, void *b2, float *b3, float *csu, void *snu, float *csv, 
    void *snv, float *csq, void *snq);


int clahef_(char *uplo, int *n, int *nb, int *kb,
    void *a, int *lda, int *ipiv, void *w, int *ldw, 
    int *info);

int clahqr_(int *wantt, int *wantz, int *n, 
    int *ilo, int *ihi, void *h__, int *ldh, void *w, 
    int *iloz, int *ihiz, void *z__, int *ldz, int *
    info);

int clahrd_(int *n, int *k, int *nb, void *a, 
    int *lda, void *tau, void *t, int *ldt, void *y, 
    int *ldy);

int claic1_(int *job, int *j, void *x, float *sest,
    void *w, void *gamma, float *sestpr, void *s, void *c__);

int clals0_(int *icompq, int *nl, int *nr, 
    int *sqre, int *nrhs, void *b, int *ldb, void *bx, 
    int *ldbx, int *perm, int *givptr, int *givcol, 
    int *ldgcol, float *givnum, int *ldgnum, float *poles, float *
    difl, float *difr, float *z__, int *k, float *c__, float *s, float *
    rwork, int *info);

int clalsa_(int *icompq, int *smlsiz, int *n, 
    int *nrhs, void *b, int *ldb, void *bx, int *ldbx, 
    float *u, int *ldu, float *vt, int *k, float *difl, float *difr, 
    float *z__, float *poles, int *givptr, int *givcol, int *
    ldgcol, int *perm, float *givnum, float *c__, float *s, float *rwork, 
    int *iwork, int *info);

int clapll_(int *n, void *x, int *incx, void *
    y, int *incy, float *ssmin);

int clapmt_(int *forwrd, int *m, int *n, void 
    *x, int *ldx, int *k);

int claqgb_(int *m, int *n, int *kl, int *ku,
    void *ab, int *ldab, float *r__, float *c__, float *rowcnd, float 
    *colcnd, float *amax, char *equed);

int claqge_(int *m, int *n, void *a, int *lda,
    float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, char *
    equed);

int claqhb_(char *uplo, int *n, int *kd, void *ab,
    int *ldab, float *s, float *scond, float *amax, char *equed);

int claqhe_(char *uplo, int *n, void *a, int *lda,
    float *s, float *scond, float *amax, char *equed);

int claqhp_(char *uplo, int *n, void *ap, float *s, 
    float *scond, float *amax, char *equed);

int claqp2_(int *m, int *n, int *offset, void 
    *a, int *lda, int *jpvt, void *tau, float *vn1, float *vn2, 
    void *work);

int claqps_(int *m, int *n, int *offset, int 
    *nb, int *kb, void *a, int *lda, int *jpvt, void *
    tau, float *vn1, float *vn2, void *auxv, void *f, int *ldf);

int claqsb_(char *uplo, int *n, int *kd, void *ab,
    int *ldab, float *s, float *scond, float *amax, char *equed);

int claqsp_(char *uplo, int *n, void *ap, float *s, 
    float *scond, float *amax, char *equed);

int claqsy_(char *uplo, int *n, void *a, int *lda,
    float *s, float *scond, float *amax, char *equed);

int clar1v_(int *n, int *b1, int *bn, float *
    sigma, float *d__, float *l, float *ld, float *lld, float *gersch, void 
    *z__, float *ztz, float *mingma, int *r__, int *isuppz, float *
    work);

int clar2v_(int *n, void *x, void *y, void *z__,
    int *incx, float *c__, void *s, int *incc);

int clarcm_(int *m, int *n, float *a, int *lda, 
    void *b, int *ldb, void *c__, int *ldc, float *rwork);

int clarf_(char *side, int *m, int *n, void *v, 
    int *incv, void *tau, void *c__, int *ldc, void *
    work);

int clarfb_(char *side, char *trans, char *direct, char *
    storev, int *m, int *n, int *k, void *v, int *ldv, 
    void *t, int *ldt, void *c__, int *ldc, void *work, 
    int *ldwork);

int clarfg_(int *n, void *alpha, void *x, int *
    incx, void *tau);

int clarft_(char *direct, char *storev, int *n, int *
    k, void *v, int *ldv, void *tau, void *t, int *ldt);

int clarfx_(char *side, int *m, int *n, void *v, 
    void *tau, void *c__, int *ldc, void *work);

int clargv_(int *n, void *x, int *incx, void *
    y, int *incy, float *c__, int *incc);

int clarnv_(int *idist, int *iseed, int *n, 
    void *x);

int clarrv_(int *n, float *d__, float *l, int *isplit, 
    int *m, float *w, int *iblock, float *gersch, float *tol, 
    void *z__, int *ldz, int *isuppz, float *work, int *
    iwork, int *info);

int clartg_(void *f, void *g, float *cs, void *sn, 
    void *r__);

int clartv_(int *n, void *x, int *incx, void *
    y, int *incy, float *c__, void *s, int *incc);

int clarz_(char *side, int *m, int *n, int *l, 
    void *v, int *incv, void *tau, void *c__, int *ldc, 
    void *work);

int clarzb_(char *side, char *trans, char *direct, char *
    storev, int *m, int *n, int *k, int *l, void *v, 
    int *ldv, void *t, int *ldt, void *c__, int *ldc, 
    void *work, int *ldwork);

int clarzt_(char *direct, char *storev, int *n, int *
    k, void *v, int *ldv, void *tau, void *t, int *ldt);

int clascl_(char *type__, int *kl, int *ku, float *
    cfrom, float *cto, int *m, int *n, void *a, int *lda, 
    int *info);

int claset_(char *uplo, int *m, int *n, void *
    alpha, void *beta, void *a, int *lda);

int clasr_(char *side, char *pivot, char *direct, int *m,
    int *n, float *c__, float *s, void *a, int *lda);

int classq_(int *n, void *x, int *incx, float *
    scale, float *sumsq);

int claswp_(int *n, void *a, int *lda, int *
    k1, int *k2, int *ipiv, int *incx);

int clasyf_(char *uplo, int *n, int *nb, int *kb,
    void *a, int *lda, int *ipiv, void *w, int *ldw, 
    int *info);

int clatbs_(char *uplo, char *trans, char *diag, char *
    normin, int *n, int *kd, void *ab, int *ldab, void *
    x, float *scale, float *cnorm, int *info);

int clatdf_(int *ijob, int *n, void *z__, int 
    *ldz, void *rhs, float *rdsum, float *rdscal, int *ipiv, int 
    *jpiv);

int clatps_(char *uplo, char *trans, char *diag, char *
    normin, int *n, void *ap, void *x, float *scale, float *cnorm,
    int *info);

int clatrd_(char *uplo, int *n, int *nb, void *a, 
    int *lda, float *e, void *tau, void *w, int *ldw);

int clatrs_(char *uplo, char *trans, char *diag, char *
    normin, int *n, void *a, int *lda, void *x, float *scale,
    float *cnorm, int *info);

int clatrz_(int *m, int *n, int *l, void *a, 
    int *lda, void *tau, void *work);

int clatzm_(char *side, int *m, int *n, void *v, 
    int *incv, void *tau, void *c1, void *c2, int *ldc, 
    void *work);

int clauu2_(char *uplo, int *n, void *a, int *lda,
    int *info);


int cpbcon_(char *uplo, int *n, int *kd, void *ab,
    int *ldab, float *anorm, float *rcond, void *work, float *rwork, 
    int *info);

int cpbequ_(char *uplo, int *n, int *kd, void *ab,
    int *ldab, float *s, float *scond, float *amax, int *info);

int cpbrfs_(char *uplo, int *n, int *kd, int *
    nrhs, void *ab, int *ldab, void *afb, int *ldafb, 
    void *b, int *ldb, void *x, int *ldx, float *ferr, float *
    berr, void *work, float *rwork, int *info);

int cpbstf_(char *uplo, int *n, int *kd, void *ab,
    int *ldab, int *info);

int cpbsv_(char *uplo, int *n, int *kd, int *
    nrhs, void *ab, int *ldab, void *b, int *ldb, int *
    info);

int cpbsvx_(char *fact, char *uplo, int *n, int *kd, 
    int *nrhs, void *ab, int *ldab, void *afb, int *
    ldafb, char *equed, float *s, void *b, int *ldb, void *x, 
    int *ldx, float *rcond, float *ferr, float *berr, void *work, 
    float *rwork, int *info);

int cpbtf2_(char *uplo, int *n, int *kd, void *ab,
    int *ldab, int *info);


int cpbtrs_(char *uplo, int *n, int *kd, int *
    nrhs, void *ab, int *ldab, void *b, int *ldb, int *
    info);

int cpocon_(char *uplo, int *n, void *a, int *lda,
    float *anorm, float *rcond, void *work, float *rwork, int *info);

int cpoequ_(int *n, void *a, int *lda, float *s, 
    float *scond, float *amax, int *info);

int cporfs_(char *uplo, int *n, int *nrhs, void *
    a, int *lda, void *af, int *ldaf, void *b, int *ldb,
    void *x, int *ldx, float *ferr, float *berr, void *work, 
    float *rwork, int *info);

int cposv_(char *uplo, int *n, int *nrhs, void *a,
    int *lda, void *b, int *ldb, int *info);

int cposvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *a, int *lda, void *af, int *ldaf, char *
    equed, float *s, void *b, int *ldb, void *x, int *ldx, 
    float *rcond, float *ferr, float *berr, void *work, float *rwork, 
    int *info);

int cpotf2_(char *uplo, int *n, void *a, int *lda,
    int *info);



int cpotrs_(char *uplo, int *n, int *nrhs, void *
    a, int *lda, void *b, int *ldb, int *info);

int cppcon_(char *uplo, int *n, void *ap, float *anorm,
    float *rcond, void *work, float *rwork, int *info);

int cppequ_(char *uplo, int *n, void *ap, float *s, 
    float *scond, float *amax, int *info);

int cpprfs_(char *uplo, int *n, int *nrhs, void *
    ap, void *afp, void *b, int *ldb, void *x, int *ldx, 
    float *ferr, float *berr, void *work, float *rwork, int *info);

int cppsv_(char *uplo, int *n, int *nrhs, void *
    ap, void *b, int *ldb, int *info);

int cppsvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *ap, void *afp, char *equed, float *s, void *b, 
    int *ldb, void *x, int *ldx, float *rcond, float *ferr, float 
    *berr, void *work, float *rwork, int *info);

int cpptrf_(char *uplo, int *n, void *ap, int *
    info);

int cpptri_(char *uplo, int *n, void *ap, int *
    info);

int cpptrs_(char *uplo, int *n, int *nrhs, void *
    ap, void *b, int *ldb, int *info);

int cptcon_(int *n, float *d__, void *e, float *anorm, 
    float *rcond, float *rwork, int *info);

int cptrfs_(char *uplo, int *n, int *nrhs, float *d__,
    void *e, float *df, void *ef, void *b, int *ldb, void 
    *x, int *ldx, float *ferr, float *berr, void *work, float *rwork, 
    int *info);

int cptsv_(int *n, int *nrhs, float *d__, void *e, 
    void *b, int *ldb, int *info);

int cptsvx_(char *fact, int *n, int *nrhs, float *d__,
    void *e, float *df, void *ef, void *b, int *ldb, void 
    *x, int *ldx, float *rcond, float *ferr, float *berr, void *work, 
    float *rwork, int *info);


int cpttrs_(char *uplo, int *n, int *nrhs, float *d__,
    void *e, void *b, int *ldb, int *info);

int cptts2_(int *iuplo, int *n, int *nrhs, float *
    d__, void *e, void *b, int *ldb);


int cspcon_(char *uplo, int *n, void *ap, int *
    ipiv, float *anorm, float *rcond, void *work, int *info);

int cspmv_(char *uplo, int *n, void *alpha, void *
    ap, void *x, int *incx, void *beta, void *y, int *
    incy);

int cspr_(char *uplo, int *n, void *alpha, void *x,
    int *incx, void *ap);

int csprfs_(char *uplo, int *n, int *nrhs, void *
    ap, void *afp, int *ipiv, void *b, int *ldb, void *x,
    int *ldx, float *ferr, float *berr, void *work, float *rwork, 
    int *info);

int cspsv_(char *uplo, int *n, int *nrhs, void *
    ap, int *ipiv, void *b, int *ldb, int *info);

int cspsvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *ap, void *afp, int *ipiv, void *b, int *
    ldb, void *x, int *ldx, float *rcond, float *ferr, float *berr, 
    void *work, float *rwork, int *info);

int csptrf_(char *uplo, int *n, void *ap, int *
    ipiv, int *info);

int csptri_(char *uplo, int *n, void *ap, int *
    ipiv, void *work, int *info);

int csptrs_(char *uplo, int *n, int *nrhs, void *
    ap, int *ipiv, void *b, int *ldb, int *info);

//int csrot_(int *n, void *cx, int *incx, void *
    //cy, int *incy, float *c__, float *s);

int csrscl_(int *n, float *sa, void *sx, int *incx);


int cstein_(int *n, float *d__, float *e, int *m, float 
    *w, int *iblock, int *isplit, void *z__, int *ldz, 
    float *work, int *iwork, int *ifail, int *info);

int csteqr_(char *compz, int *n, float *d__, float *e, 
    void *z__, int *ldz, float *work, int *info);

int csycon_(char *uplo, int *n, void *a, int *lda,
    int *ipiv, float *anorm, float *rcond, void *work, int *
    info);



int csyrfs_(char *uplo, int *n, int *nrhs, void *
    a, int *lda, void *af, int *ldaf, int *ipiv, void *
    b, int *ldb, void *x, int *ldx, float *ferr, float *berr, 
    void *work, float *rwork, int *info);

int csysv_(char *uplo, int *n, int *nrhs, void *a,
    int *lda, int *ipiv, void *b, int *ldb, void *work,
    int *lwork, int *info);

int csysvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *a, int *lda, void *af, int *ldaf, int *
    ipiv, void *b, int *ldb, void *x, int *ldx, float *rcond,
    float *ferr, float *berr, void *work, int *lwork, float *rwork, 
    int *info);

int csytf2_(char *uplo, int *n, void *a, int *lda,
    int *ipiv, int *info);


int csytri_(char *uplo, int *n, void *a, int *lda,
    int *ipiv, void *work, int *info);

int csytrs_(char *uplo, int *n, int *nrhs, void *
    a, int *lda, int *ipiv, void *b, int *ldb, int *
    info);

int ctbcon_(char *norm, char *uplo, char *diag, int *n, 
    int *kd, void *ab, int *ldab, float *rcond, void *work, 
    float *rwork, int *info);

int ctbrfs_(char *uplo, char *trans, char *diag, int *n, 
    int *kd, int *nrhs, void *ab, int *ldab, void *b, 
    int *ldb, void *x, int *ldx, float *ferr, float *berr, 
    void *work, float *rwork, int *info);


int ctgevc_(char *side, char *howmny, int *select, 
    int *n, void *a, int *lda, void *b, int *ldb, 
    void *vl, int *ldvl, void *vr, int *ldvr, int *mm, 
    int *m, void *work, float *rwork, int *info);

int ctgex2_(int *wantq, int *wantz, int *n, 
    void *a, int *lda, void *b, int *ldb, void *q, 
    int *ldq, void *z__, int *ldz, int *j1, int *info);

int ctgexc_(int *wantq, int *wantz, int *n, 
    void *a, int *lda, void *b, int *ldb, void *q, 
    int *ldq, void *z__, int *ldz, int *ifst, int *
    ilst, int *info);

int ctgsen_(int *ijob, int *wantq, int *wantz, 
    int *select, int *n, void *a, int *lda, void *b, 
    int *ldb, void *alpha, void *beta, void *q, int *ldq,
    void *z__, int *ldz, int *m, float *pl, float *pr, float *
    dif, void *work, int *lwork, int *iwork, int *liwork, 
    int *info);

int ctgsja_(char *jobu, char *jobv, char *jobq, int *m, 
    int *p, int *n, int *k, int *l, void *a, int *
    lda, void *b, int *ldb, float *tola, float *tolb, float *alpha, 
    float *beta, void *u, int *ldu, void *v, int *ldv, 
    void *q, int *ldq, void *work, int *ncycle, int *
    info);

int ctgsna_(char *job, char *howmny, int *select, 
    int *n, void *a, int *lda, void *b, int *ldb, 
    void *vl, int *ldvl, void *vr, int *ldvr, float *s, float 
    *dif, int *mm, int *m, void *work, int *lwork, int 
    *iwork, int *info);

int ctgsy2_(char *trans, int *ijob, int *m, int *
    n, void *a, int *lda, void *b, int *ldb, void *c__, 
    int *ldc, void *d__, int *ldd, void *e, int *lde, 
    void *f, int *ldf, float *scale, float *rdsum, float *rdscal, 
    int *info);

int ctgsyl_(char *trans, int *ijob, int *m, int *
    n, void *a, int *lda, void *b, int *ldb, void *c__, 
    int *ldc, void *d__, int *ldd, void *e, int *lde, 
    void *f, int *ldf, float *scale, float *dif, void *work, 
    int *lwork, int *iwork, int *info);

int ctpcon_(char *norm, char *uplo, char *diag, int *n, 
    void *ap, float *rcond, void *work, float *rwork, int *info);

int ctprfs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, void *ap, void *b, int *ldb, void *x, 
    int *ldx, float *ferr, float *berr, void *work, float *rwork, 
    int *info);

int ctptri_(char *uplo, char *diag, int *n, void *ap, 
    int *info);

int ctptrs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, void *ap, void *b, int *ldb, int *info);

int ctrcon_(char *norm, char *uplo, char *diag, int *n, 
    void *a, int *lda, float *rcond, void *work, float *rwork, 
    int *info);

int ctrevc_(char *side, char *howmny, int *select, 
    int *n, void *t, int *ldt, void *vl, int *ldvl, 
    void *vr, int *ldvr, int *mm, int *m, void *work, 
    float *rwork, int *info);

int ctrexc_(char *compq, int *n, void *t, int *
    ldt, void *q, int *ldq, int *ifst, int *ilst, int *
    info);

int ctrrfs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, void *a, int *lda, void *b, int *ldb, 
    void *x, int *ldx, float *ferr, float *berr, void *work, float 
    *rwork, int *info);

int ctrsen_(char *job, char *compq, int *select, int 
    *n, void *t, int *ldt, void *q, int *ldq, void *w, 
    int *m, float *s, float *sep, void *work, int *lwork, 
    int *info);

int ctrsna_(char *job, char *howmny, int *select, 
    int *n, void *t, int *ldt, void *vl, int *ldvl, 
    void *vr, int *ldvr, float *s, float *sep, int *mm, int *
    m, void *work, int *ldwork, float *rwork, int *info);

int ctrsyl_(char *trana, char *tranb, int *isgn, int 
    *m, int *n, void *a, int *lda, void *b, int *ldb, 
    void *c__, int *ldc, float *scale, int *info);

int ctrti2_(char *uplo, char *diag, int *n, void *a, 
    int *lda, int *info);


int ctrtrs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, void *a, int *lda, void *b, int *ldb, 
    int *info);

int ctzrqf_(int *m, int *n, void *a, int *lda,
    void *tau, int *info);

int ctzrzf_(int *m, int *n, void *a, int *lda,
    void *tau, void *work, int *lwork, int *info);

int cung2l_(int *m, int *n, int *k, void *a, 
    int *lda, void *tau, void *work, int *info);

int cung2r_(int *m, int *n, int *k, void *a, 
    int *lda, void *tau, void *work, int *info);

int cungbr_(char *vect, int *m, int *n, int *k, 
    void *a, int *lda, void *tau, void *work, int *lwork,
    int *info);

int cunghr_(int *n, int *ilo, int *ihi, void *
    a, int *lda, void *tau, void *work, int *lwork, int 
    *info);

int cungl2_(int *m, int *n, int *k, void *a, 
    int *lda, void *tau, void *work, int *info);


int cungql_(int *m, int *n, int *k, void *a, 
    int *lda, void *tau, void *work, int *lwork, int *
    info);


int cungr2_(int *m, int *n, int *k, void *a, 
    int *lda, void *tau, void *work, int *info);

int cungrq_(int *m, int *n, int *k, void *a, 
    int *lda, void *tau, void *work, int *lwork, int *
    info);

int cungtr_(char *uplo, int *n, void *a, int *lda,
    void *tau, void *work, int *lwork, int *info);

int cunm2l_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, void *c__, 
    int *ldc, void *work, int *info);

int cunm2r_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, void *c__, 
    int *ldc, void *work, int *info);

int cunmbr_(char *vect, char *side, char *trans, int *m, 
    int *n, int *k, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *lwork, int *
    info);

int cunmhr_(char *side, char *trans, int *m, int *n, 
    int *ilo, int *ihi, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *lwork, int *
    info);

int cunml2_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, void *c__, 
    int *ldc, void *work, int *info);


int cunmql_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, void *c__, 
    int *ldc, void *work, int *lwork, int *info);


int cunmr2_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, void *c__, 
    int *ldc, void *work, int *info);

int cunmr3_(char *side, char *trans, int *m, int *n, 
    int *k, int *l, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *info);

int cunmrq_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, void *c__, 
    int *ldc, void *work, int *lwork, int *info);

int cunmrz_(char *side, char *trans, int *m, int *n, 
    int *k, int *l, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *lwork, int *
    info);

int cunmtr_(char *side, char *uplo, char *trans, int *m, 
    int *n, void *a, int *lda, void *tau, void *c__, 
    int *ldc, void *work, int *lwork, int *info);

int cupgtr_(char *uplo, int *n, void *ap, void *
    tau, void *q, int *ldq, void *work, int *info);

int cupmtr_(char *side, char *uplo, char *trans, int *m, 
    int *n, void *ap, void *tau, void *c__, int *ldc, 
    void *work, int *info);


int dbdsqr_(char *uplo, int *n, int *ncvt, int *
    nru, int *ncc, double *d__, double *e, double *vt, 
    int *ldvt, double *u, int *ldu, double *c__, int *
    ldc, double *work, int *info);

int ddisna_(char *job, int *m, int *n, double *
    d__, double *sep, int *info);


int dgbcon_(char *norm, int *n, int *kl, int *ku,
    double *ab, int *ldab, int *ipiv, double *anorm, 
    double *rcond, double *work, int *iwork, int *info);

int dgbequ_(int *m, int *n, int *kl, int *ku,
    double *ab, int *ldab, double *r__, double *c__, 
    double *rowcnd, double *colcnd, double *amax, int *
    info);

int dgbrfs_(char *trans, int *n, int *kl, int *
    ku, int *nrhs, double *ab, int *ldab, double *afb, 
    int *ldafb, int *ipiv, double *b, int *ldb, 
    double *x, int *ldx, double *ferr, double *berr, 
    double *work, int *iwork, int *info);

int dgbsv_(int *n, int *kl, int *ku, int *
    nrhs, double *ab, int *ldab, int *ipiv, double *b, 
    int *ldb, int *info);

int dgbsvx_(char *fact, char *trans, int *n, int *kl,
    int *ku, int *nrhs, double *ab, int *ldab, 
    double *afb, int *ldafb, int *ipiv, char *equed, 
    double *r__, double *c__, double *b, int *ldb, 
    double *x, int *ldx, double *rcond, double *ferr, 
    double *berr, double *work, int *iwork, int *info);

int dgbtf2_(int *m, int *n, int *kl, int *ku,
    double *ab, int *ldab, int *ipiv, int *info);



int dgebak_(char *job, char *side, int *n, int *ilo, 
    int *ihi, double *scale, int *m, double *v, int *
    ldv, int *info);

int dgebal_(char *job, int *n, double *a, int *
    lda, int *ilo, int *ihi, double *scale, int *info);

int dgebd2_(int *m, int *n, double *a, int *
    lda, double *d__, double *e, double *tauq, double *
    taup, double *work, int *info);


int dgecon_(char *norm, int *n, double *a, int *
    lda, double *anorm, double *rcond, double *work, int *
    iwork, int *info);

int dgeequ_(int *m, int *n, double *a, int *
    lda, double *r__, double *c__, double *rowcnd, double 
    *colcnd, double *amax, int *info);

int dgeev_(char *jobvl, char *jobvr, int *n, double *
    a, int *lda, double *wr, double *wi, double *vl, 
    int *ldvl, double *vr, int *ldvr, double *work, 
    int *lwork, int *info);

int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *
    sense, int *n, double *a, int *lda, double *wr, 
    double *wi, double *vl, int *ldvl, double *vr, 
    int *ldvr, int *ilo, int *ihi, double *scale, 
    double *abnrm, double *rconde, double *rcondv, double 
    *work, int *lwork, int *iwork, int *info);

int dgegs_(char *jobvsl, char *jobvsr, int *n, 
    double *a, int *lda, double *b, int *ldb, double *
    alphar, double *alphai, double *beta, double *vsl, 
    int *ldvsl, double *vsr, int *ldvsr, double *work, 
    int *lwork, int *info);

int dgegv_(char *jobvl, char *jobvr, int *n, double *
    a, int *lda, double *b, int *ldb, double *alphar, 
    double *alphai, double *beta, double *vl, int *ldvl, 
    double *vr, int *ldvr, double *work, int *lwork, 
    int *info);

int dgehd2_(int *n, int *ilo, int *ihi, 
    double *a, int *lda, double *tau, double *work, 
    int *info);

int dgehrd_(int *n, int *ilo, int *ihi, 
    double *a, int *lda, double *tau, double *work, 
    int *lwork, int *info);

int dgelq2_(int *m, int *n, double *a, int *
    lda, double *tau, double *work, int *info);


int dgels_(char *trans, int *m, int *n, int *
    nrhs, double *a, int *lda, double *b, int *ldb, 
    double *work, int *lwork, int *info);

int dgelsd_(int *m, int *n, int *nrhs, 
    double *a, int *lda, double *b, int *ldb, double *
    s, double *rcond, int *rank, double *work, int *lwork,
    int *iwork, int *info);

int dgelss_(int *m, int *n, int *nrhs, 
    double *a, int *lda, double *b, int *ldb, double *
    s, double *rcond, int *rank, double *work, int *lwork,
    int *info);

int dgelsx_(int *m, int *n, int *nrhs, 
    double *a, int *lda, double *b, int *ldb, int *
    jpvt, double *rcond, int *rank, double *work, int *
    info);

int dgelsy_(int *m, int *n, int *nrhs, 
    double *a, int *lda, double *b, int *ldb, int *
    jpvt, double *rcond, int *rank, double *work, int *
    lwork, int *info);

int dgeql2_(int *m, int *n, double *a, int *
    lda, double *tau, double *work, int *info);

int dgeqlf_(int *m, int *n, double *a, int *
    lda, double *tau, double *work, int *lwork, int *info);


int dgeqpf_(int *m, int *n, double *a, int *
    lda, int *jpvt, double *tau, double *work, int *info);

int dgeqr2_(int *m, int *n, double *a, int *
    lda, double *tau, double *work, int *info);


int dgerfs_(char *trans, int *n, int *nrhs, 
    double *a, int *lda, double *af, int *ldaf, int *
    ipiv, double *b, int *ldb, double *x, int *ldx, 
    double *ferr, double *berr, double *work, int *iwork, 
    int *info);

int dgerq2_(int *m, int *n, double *a, int *
    lda, double *tau, double *work, int *info);

int dgerqf_(int *m, int *n, double *a, int *
    lda, double *tau, double *work, int *lwork, int *info);

int dgesc2_(int *n, double *a, int *lda, 
    double *rhs, int *ipiv, int *jpiv, double *scale);

int dgesdd_(char *jobz, int *m, int *n, double *
    a, int *lda, double *s, double *u, int *ldu, 
    double *vt, int *ldvt, double *work, int *lwork, 
    int *iwork, int *info);

int dgesv_(int *n, int *nrhs, double *a, int 
    *lda, int *ipiv, double *b, int *ldb, int *info);

int dgesvd_(char *jobu, char *jobvt, int *m, int *n, 
    double *a, int *lda, double *s, double *u, int *
    ldu, double *vt, int *ldvt, double *work, int *lwork, 
    int *info);

int dgesvx_(char *fact, char *trans, int *n, int *
    nrhs, double *a, int *lda, double *af, int *ldaf, 
    int *ipiv, char *equed, double *r__, double *c__, 
    double *b, int *ldb, double *x, int *ldx, double *
    rcond, double *ferr, double *berr, double *work, int *
    iwork, int *info);

int dgetc2_(int *n, double *a, int *lda, int 
    *ipiv, int *jpiv, int *info);

int dgetf2_(int *m, int *n, double *a, int *
    lda, int *ipiv, int *info);




int dggbak_(char *job, char *side, int *n, int *ilo, 
    int *ihi, double *lscale, double *rscale, int *m, 
    double *v, int *ldv, int *info);

int dggbal_(char *job, int *n, double *a, int *
    lda, double *b, int *ldb, int *ilo, int *ihi, 
    double *lscale, double *rscale, double *work, int *
    info);

int dggev_(char *jobvl, char *jobvr, int *n, double *
    a, int *lda, double *b, int *ldb, double *alphar, 
    double *alphai, double *beta, double *vl, int *ldvl, 
    double *vr, int *ldvr, double *work, int *lwork, 
    int *info);

int dggevx_(char *balanc, char *jobvl, char *jobvr, char *
    sense, int *n, double *a, int *lda, double *b, 
    int *ldb, double *alphar, double *alphai, double *
    beta, double *vl, int *ldvl, double *vr, int *ldvr, 
    int *ilo, int *ihi, double *lscale, double *rscale, 
    double *abnrm, double *bbnrm, double *rconde, double *
    rcondv, double *work, int *lwork, int *iwork, int *
    bwork, int *info);

int dggglm_(int *n, int *m, int *p, double *
    a, int *lda, double *b, int *ldb, double *d__, 
    double *x, double *y, double *work, int *lwork, 
    int *info);

int dgghrd_(char *compq, char *compz, int *n, int *
    ilo, int *ihi, double *a, int *lda, double *b, 
    int *ldb, double *q, int *ldq, double *z__, int *
    ldz, int *info);

int dgglse_(int *m, int *n, int *p, double *
    a, int *lda, double *b, int *ldb, double *c__, 
    double *d__, double *x, double *work, int *lwork, 
    int *info);

int dggqrf_(int *n, int *m, int *p, double *
    a, int *lda, double *taua, double *b, int *ldb, 
    double *taub, double *work, int *lwork, int *info);

int dggrqf_(int *m, int *p, int *n, double *
    a, int *lda, double *taua, double *b, int *ldb, 
    double *taub, double *work, int *lwork, int *info);

int dggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
    int *n, int *p, int *k, int *l, double *a, 
    int *lda, double *b, int *ldb, double *alpha, 
    double *beta, double *u, int *ldu, double *v, int 
    *ldv, double *q, int *ldq, double *work, int *iwork, 
    int *info);

int dggsvp_(char *jobu, char *jobv, char *jobq, int *m, 
    int *p, int *n, double *a, int *lda, double *b, 
    int *ldb, double *tola, double *tolb, int *k, int 
    *l, double *u, int *ldu, double *v, int *ldv, 
    double *q, int *ldq, int *iwork, double *tau, 
    double *work, int *info);

int dgtcon_(char *norm, int *n, double *dl, 
    double *d__, double *du, double *du2, int *ipiv, 
    double *anorm, double *rcond, double *work, int *
    iwork, int *info);

int dgtrfs_(char *trans, int *n, int *nrhs, 
    double *dl, double *d__, double *du, double *dlf, 
    double *df, double *duf, double *du2, int *ipiv, 
    double *b, int *ldb, double *x, int *ldx, double *
    ferr, double *berr, double *work, int *iwork, int *
    info);

int dgtsv_(int *n, int *nrhs, double *dl, 
    double *d__, double *du, double *b, int *ldb, int 
    *info);

int dgtsvx_(char *fact, char *trans, int *n, int *
    nrhs, double *dl, double *d__, double *du, double *
    dlf, double *df, double *duf, double *du2, int *ipiv, 
    double *b, int *ldb, double *x, int *ldx, double *
    rcond, double *ferr, double *berr, double *work, int *
    iwork, int *info);



int dgtts2_(int *itrans, int *n, int *nrhs, 
    double *dl, double *d__, double *du, double *du2, 
    int *ipiv, double *b, int *ldb);

int dhgeqz_(char *job, char *compq, char *compz, int *n, 
    int *ilo, int *ihi, double *a, int *lda, double *
    b, int *ldb, double *alphar, double *alphai, double *
    beta, double *q, int *ldq, double *z__, int *ldz, 
    double *work, int *lwork, int *info);

int dhsein_(char *side, char *eigsrc, char *initv, int *
    select, int *n, double *h__, int *ldh, double *wr, 
    double *wi, double *vl, int *ldvl, double *vr, 
    int *ldvr, int *mm, int *m, double *work, int *
    ifaill, int *ifailr, int *info);

int dhseqr_(char *job, char *compz, int *n, int *ilo,
    int *ihi, double *h__, int *ldh, double *wr, 
    double *wi, double *z__, int *ldz, double *work, 
    int *lwork, int *info);

int dlabad_(double *small, double *large);

int dlabrd_(int *m, int *n, int *nb, double *
    a, int *lda, double *d__, double *e, double *tauq, 
    double *taup, double *x, int *ldx, double *y, int 
    *ldy);

int dlacon_(int *n, double *v, double *x, 
    int *isgn, double *est, int *kase);


int dladiv_(double *a, double *b, double *c__, 
    double *d__, double *p, double *q);

int dlae2_(double *a, double *b, double *c__, 
    double *rt1, double *rt2);

int dlaebz_(int *ijob, int *nitmax, int *n, 
    int *mmax, int *minp, int *nbmin, double *abstol, 
    double *reltol, double *pivmin, double *d__, double *
    e, double *e2, int *nval, double *ab, double *c__, 
    int *mout, int *nab, double *work, int *iwork, 
    int *info);

int dlaed0_(int *icompq, int *qsiz, int *n, 
    double *d__, double *e, double *q, int *ldq, 
    double *qstore, int *ldqs, double *work, int *iwork, 
    int *info);

int dlaed1_(int *n, double *d__, double *q, 
    int *ldq, int *indxq, double *rho, int *cutpnt, 
    double *work, int *iwork, int *info);

int dlaed2_(int *k, int *n, int *n1, double *
    d__, double *q, int *ldq, int *indxq, double *rho, 
    double *z__, double *dlamda, double *w, double *q2, 
    int *indx, int *indxc, int *indxp, int *coltyp, 
    int *info);

int dlaed3_(int *k, int *n, int *n1, double *
    d__, double *q, int *ldq, double *rho, double *dlamda,
    double *q2, int *indx, int *ctot, double *w, 
    double *s, int *info);

int dlaed4_(int *n, int *i__, double *d__, 
    double *z__, double *delta, double *rho, double *dlam,
    int *info);

int dlaed5_(int *i__, double *d__, double *z__, 
    double *delta, double *rho, double *dlam);

int dlaed6_(int *kniter, int *orgati, double *
    rho, double *d__, double *z__, double *finit, double *
    tau, int *info);

int dlaed7_(int *icompq, int *n, int *qsiz, 
    int *tlvls, int *curlvl, int *curpbm, double *d__, 
    double *q, int *ldq, int *indxq, double *rho, int 
    *cutpnt, double *qstore, int *qptr, int *prmptr, int *
    perm, int *givptr, int *givcol, double *givnum, 
    double *work, int *iwork, int *info);

int dlaed8_(int *icompq, int *k, int *n, int 
    *qsiz, double *d__, double *q, int *ldq, int *indxq, 
    double *rho, int *cutpnt, double *z__, double *dlamda,
    double *q2, int *ldq2, double *w, int *perm, int 
    *givptr, int *givcol, double *givnum, int *indxp, int 
    *indx, int *info);

int dlaed9_(int *k, int *kstart, int *kstop, 
    int *n, double *d__, double *q, int *ldq, double *
    rho, double *dlamda, double *w, double *s, int *lds, 
    int *info);

int dlaeda_(int *n, int *tlvls, int *curlvl, 
    int *curpbm, int *prmptr, int *perm, int *givptr, 
    int *givcol, double *givnum, double *q, int *qptr, 
    double *z__, double *ztemp, int *info);

int dlaein_(int *rightv, int *noinit, int *n, 
    double *h__, int *ldh, double *wr, double *wi, 
    double *vr, double *vi, double *b, int *ldb, 
    double *work, double *eps3, double *smlnum, double *
    bignum, int *info);

int dlaev2_(double *a, double *b, double *c__, 
    double *rt1, double *rt2, double *cs1, double *sn1);

int dlaexc_(int *wantq, int *n, double *t, 
    int *ldt, double *q, int *ldq, int *j1, int *n1, 
    int *n2, double *work, int *info);

int dlag2_(double *a, int *lda, double *b, 
    int *ldb, double *safmin, double *scale1, double *
    scale2, double *wr1, double *wr2, double *wi);

int dlags2_(int *upper, double *a1, double *a2, 
    double *a3, double *b1, double *b2, double *b3, 
    double *csu, double *snu, double *csv, double *snv, 
    double *csq, double *snq);

int dlagtf_(int *n, double *a, double *lambda, 
    double *b, double *c__, double *tol, double *d__, 
    int *in, int *info);


int dlagts_(int *job, int *n, double *a, 
    double *b, double *c__, double *d__, int *in, 
    double *y, double *tol, int *info);

int dlagv2_(double *a, int *lda, double *b, 
    int *ldb, double *alphar, double *alphai, double *
    beta, double *csl, double *snl, double *csr, double *
    snr);

int dlahqr_(int *wantt, int *wantz, int *n, 
    int *ilo, int *ihi, double *h__, int *ldh, double 
    *wr, double *wi, int *iloz, int *ihiz, double *z__, 
    int *ldz, int *info);

int dlahrd_(int *n, int *k, int *nb, double *
    a, int *lda, double *tau, double *t, int *ldt, 
    double *y, int *ldy);

int dlaic1_(int *job, int *j, double *x, 
    double *sest, double *w, double *gamma, double *
    sestpr, double *s, double *c__);

int dlaln2_(int *ltrans, int *na, int *nw, 
    double *smin, double *ca, double *a, int *lda, 
    double *d1, double *d2, double *b, int *ldb, 
    double *wr, double *wi, double *x, int *ldx, 
    double *scale, double *xnorm, int *info);

int dlals0_(int *icompq, int *nl, int *nr, 
    int *sqre, int *nrhs, double *b, int *ldb, double 
    *bx, int *ldbx, int *perm, int *givptr, int *givcol, 
    int *ldgcol, double *givnum, int *ldgnum, double *
    poles, double *difl, double *difr, double *z__, int *
    k, double *c__, double *s, double *work, int *info);

int dlalsa_(int *icompq, int *smlsiz, int *n, 
    int *nrhs, double *b, int *ldb, double *bx, int *
    ldbx, double *u, int *ldu, double *vt, int *k, 
    double *difl, double *difr, double *z__, double *
    poles, int *givptr, int *givcol, int *ldgcol, int *
    perm, double *givnum, double *c__, double *s, double *
    work, int *iwork, int *info);

int dlalsd_(char *uplo, int *smlsiz, int *n, int 
    *nrhs, double *d__, double *e, double *b, int *ldb, 
    double *rcond, int *rank, double *work, int *iwork, 
    int *info);

int dlamc1_(int *beta, int *t, int *rnd, int 
    *ieee1);

int dlamc2_(int *beta, int *t, int *rnd, 
    double *eps, int *emin, double *rmin, int *emax, 
    double *rmax);

int dlamc4_(int *emin, double *start, int *base);

int dlamc5_(int *beta, int *p, int *emin, 
    int *ieee, int *emax, double *rmax);

int dlamrg_(int *n1, int *n2, double *a, int 
    *dtrd1, int *dtrd2, int *index);

int dlanv2_(double *a, double *b, double *c__, 
    double *d__, double *rt1r, double *rt1i, double *rt2r,
    double *rt2i, double *cs, double *sn);

int dlapll_(int *n, double *x, int *incx, 
    double *y, int *incy, double *ssmin);

int dlapmt_(int *forwrd, int *m, int *n, 
    double *x, int *ldx, int *k);

int dlaqgb_(int *m, int *n, int *kl, int *ku,
    double *ab, int *ldab, double *r__, double *c__, 
    double *rowcnd, double *colcnd, double *amax, char *equed);

int dlaqge_(int *m, int *n, double *a, int *
    lda, double *r__, double *c__, double *rowcnd, double 
    *colcnd, double *amax, char *equed);

int dlaqp2_(int *m, int *n, int *offset, 
    double *a, int *lda, int *jpvt, double *tau, 
    double *vn1, double *vn2, double *work);

int dlaqps_(int *m, int *n, int *offset, int 
    *nb, int *kb, double *a, int *lda, int *jpvt, 
    double *tau, double *vn1, double *vn2, double *auxv, 
    double *f, int *ldf);

int dlaqsb_(char *uplo, int *n, int *kd, double *
    ab, int *ldab, double *s, double *scond, double *amax,
    char *equed);

int dlaqsp_(char *uplo, int *n, double *ap, 
    double *s, double *scond, double *amax, char *equed);

int dlaqsy_(char *uplo, int *n, double *a, int *
    lda, double *s, double *scond, double *amax, char *equed);

int dlaqtr_(int *ltran, int *lfloat, int *n, 
    double *t, int *ldt, double *b, double *w, double 
    *scale, double *x, double *work, int *info);

int dlar1v_(int *n, int *b1, int *bn, double 
    *sigma, double *d__, double *l, double *ld, double *
    lld, double *gersch, double *z__, double *ztz, double 
    *mingma, int *r__, int *isuppz, double *work);

int dlar2v_(int *n, double *x, double *y, 
    double *z__, int *incx, double *c__, double *s, 
    int *incc);

int dlarf_(char *side, int *m, int *n, double *v,
    int *incv, double *tau, double *c__, int *ldc, 
    double *work);

int dlarfb_(char *side, char *trans, char *direct, char *
    storev, int *m, int *n, int *k, double *v, int *
    ldv, double *t, int *ldt, double *c__, int *ldc, 
    double *work, int *ldwork);

int dlarfg_(int *n, double *alpha, double *x, 
    int *incx, double *tau);

int dlarft_(char *direct, char *storev, int *n, int *
    k, double *v, int *ldv, double *tau, double *t, 
    int *ldt);

int dlarfx_(char *side, int *m, int *n, double *
    v, double *tau, double *c__, int *ldc, double *work);

int dlargv_(int *n, double *x, int *incx, 
    double *y, int *incy, double *c__, int *incc);

int dlarnv_(int *idist, int *iseed, int *n, 
    double *x);

int dlarrb_(int *n, double *d__, double *l, 
    double *ld, double *lld, int *ifirst, int *ilast, 
    double *sigma, double *reltol, double *w, double *
    wgap, double *werr, double *work, int *iwork, int *
    info);

int dlarre_(int *n, double *d__, double *e, 
    double *tol, int *nsplit, int *isplit, int *m, 
    double *w, double *woff, double *gersch, double *work,
    int *info);

int dlarrf_(int *n, double *d__, double *l, 
    double *ld, double *lld, int *ifirst, int *ilast, 
    double *w, double *dplus, double *lplus, double *work,
    int *iwork, int *info);

int dlarrv_(int *n, double *d__, double *l, 
    int *isplit, int *m, double *w, int *iblock, 
    double *gersch, double *tol, double *z__, int *ldz, 
    int *isuppz, double *work, int *iwork, int *info);

int dlartg_(double *f, double *g, double *cs, 
    double *sn, double *r__);

int dlartv_(int *n, double *x, int *incx, 
    double *y, int *incy, double *c__, double *s, int 
    *incc);

int dlaruv_(int *iseed, int *n, double *x);

int dlarz_(char *side, int *m, int *n, int *l, 
    double *v, int *incv, double *tau, double *c__, 
    int *ldc, double *work);

int dlarzb_(char *side, char *trans, char *direct, char *
    storev, int *m, int *n, int *k, int *l, double *v,
    int *ldv, double *t, int *ldt, double *c__, int *
    ldc, double *work, int *ldwork);

int dlarzt_(char *direct, char *storev, int *n, int *
    k, double *v, int *ldv, double *tau, double *t, 
    int *ldt);

int dlas2_(double *f, double *g, double *h__, 
    double *ssmin, double *ssmax);

int dlascl_(char *type__, int *kl, int *ku, 
    double *cfrom, double *cto, int *m, int *n, 
    double *a, int *lda, int *info);

int dlasd0_(int *n, int *sqre, double *d__, 
    double *e, double *u, int *ldu, double *vt, int *
    ldvt, int *smlsiz, int *iwork, double *work, int *
    info);

int dlasd1_(int *nl, int *nr, int *sqre, 
    double *d__, double *alpha, double *beta, double *u, 
    int *ldu, double *vt, int *ldvt, int *idxq, int *
    iwork, double *work, int *info);

int dlasd2_(int *nl, int *nr, int *sqre, int 
    *k, double *d__, double *z__, double *alpha, double *
    beta, double *u, int *ldu, double *vt, int *ldvt, 
    double *dsigma, double *u2, int *ldu2, double *vt2, 
    int *ldvt2, int *idxp, int *idx, int *idxc, int *
    idxq, int *coltyp, int *info);

int dlasd3_(int *nl, int *nr, int *sqre, int 
    *k, double *d__, double *q, int *ldq, double *dsigma, 
    double *u, int *ldu, double *u2, int *ldu2, 
    double *vt, int *ldvt, double *vt2, int *ldvt2, 
    int *idxc, int *ctot, double *z__, int *info);

int dlasd4_(int *n, int *i__, double *d__, 
    double *z__, double *delta, double *rho, double *
    sigma, double *work, int *info);

int dlasd5_(int *i__, double *d__, double *z__, 
    double *delta, double *rho, double *dsigma, double *
    work);

int dlasd6_(int *icompq, int *nl, int *nr, 
    int *sqre, double *d__, double *vf, double *vl, 
    double *alpha, double *beta, int *idxq, int *perm, 
    int *givptr, int *givcol, int *ldgcol, double *givnum,
    int *ldgnum, double *poles, double *difl, double *
    difr, double *z__, int *k, double *c__, double *s, 
    double *work, int *iwork, int *info);

int dlasd7_(int *icompq, int *nl, int *nr, 
    int *sqre, int *k, double *d__, double *z__, 
    double *zw, double *vf, double *vfw, double *vl, 
    double *vlw, double *alpha, double *beta, double *
    dsigma, int *idx, int *idxp, int *idxq, int *perm, 
    int *givptr, int *givcol, int *ldgcol, double *givnum,
    int *ldgnum, double *c__, double *s, int *info);

int dlasd8_(int *icompq, int *k, double *d__, 
    double *z__, double *vf, double *vl, double *difl, 
    double *difr, int *lddifr, double *dsigma, double *
    work, int *info);

int dlasd9_(int *icompq, int *ldu, int *k, 
    double *d__, double *z__, double *vf, double *vl, 
    double *difl, double *difr, double *dsigma, double *
    work, int *info);

int dlasda_(int *icompq, int *smlsiz, int *n, 
    int *sqre, double *d__, double *e, double *u, int 
    *ldu, double *vt, int *k, double *difl, double *difr, 
    double *z__, double *poles, int *givptr, int *givcol, 
    int *ldgcol, int *perm, double *givnum, double *c__, 
    double *s, double *work, int *iwork, int *info);

int dlasdq_(char *uplo, int *sqre, int *n, int *
    ncvt, int *nru, int *ncc, double *d__, double *e, 
    double *vt, int *ldvt, double *u, int *ldu, 
    double *c__, int *ldc, double *work, int *info);

int dlasdt_(int *n, int *lvl, int *nd, int *
    inode, int *ndiml, int *ndimr, int *msub);

int dlaset_(char *uplo, int *m, int *n, double *
    alpha, double *beta, double *a, int *lda);

int dlasq1_(int *n, double *d__, double *e, 
    double *work, int *info);

int dlasq2_(int *n, double *z__, int *info);

int dlasq3_(int *i0, int *n0, double *z__, 
    int *pp, double *dmin__, double *sigma, double *desig,
    double *qmax, int *nfail, int *iter, int *ndiv, 
    int *ieee);

int dlasq4_(int *i0, int *n0, double *z__, 
    int *pp, int *n0in, double *dmin__, double *dmin1, 
    double *dmin2, double *dn, double *dn1, double *dn2, 
    double *tau, int *ttype);

int dlasq5_(int *i0, int *n0, double *z__, 
    int *pp, double *tau, double *dmin__, double *dmin1, 
    double *dmin2, double *dn, double *dnm1, double *dnm2,
    int *ieee);

int dlasq6_(int *i0, int *n0, double *z__, 
    int *pp, double *dmin__, double *dmin1, double *dmin2,
    double *dn, double *dnm1, double *dnm2);

int dlasr_(char *side, char *pivot, char *direct, int *m,
    int *n, double *c__, double *s, double *a, int *
    lda);

int dlasrt_(char *id, int *n, double *d__, int *
    info);

int dlassq_(int *n, double *x, int *incx, 
    double *scale, double *sumsq);

int dlasv2_(double *f, double *g, double *h__, 
    double *ssmin, double *ssmax, double *snr, double *
    csr, double *snl, double *csl);

int dlaswp_(int *n, double *a, int *lda, int 
    *k1, int *k2, int *ipiv, int *incx);

int dlasy2_(int *ltranl, int *ltranr, int *isgn, 
    int *n1, int *n2, double *tl, int *ldtl, double *
    tr, int *ldtr, double *b, int *ldb, double *scale, 
    double *x, int *ldx, double *xnorm, int *info);

int dlasyf_(char *uplo, int *n, int *nb, int *kb,
    double *a, int *lda, int *ipiv, double *w, int *
    ldw, int *info);

int dlatbs_(char *uplo, char *trans, char *diag, char *
    normin, int *n, int *kd, double *ab, int *ldab, 
    double *x, double *scale, double *cnorm, int *info);

int dlatdf_(int *ijob, int *n, double *z__, 
    int *ldz, double *rhs, double *rdsum, double *rdscal, 
    int *ipiv, int *jpiv);

int dlatps_(char *uplo, char *trans, char *diag, char *
    normin, int *n, double *ap, double *x, double *scale, 
    double *cnorm, int *info);

int dlatrd_(char *uplo, int *n, int *nb, double *
    a, int *lda, double *e, double *tau, double *w, 
    int *ldw);

int dlatrs_(char *uplo, char *trans, char *diag, char *
    normin, int *n, double *a, int *lda, double *x, 
    double *scale, double *cnorm, int *info);

int dlatrz_(int *m, int *n, int *l, double *
    a, int *lda, double *tau, double *work);

int dlatzm_(char *side, int *m, int *n, double *
    v, int *incv, double *tau, double *c1, double *c2, 
    int *ldc, double *work);

int dlauu2_(char *uplo, int *n, double *a, int *
    lda, int *info);


int dopgtr_(char *uplo, int *n, double *ap, 
    double *tau, double *q, int *ldq, double *work, 
    int *info);

int dopmtr_(char *side, char *uplo, char *trans, int *m, 
    int *n, double *ap, double *tau, double *c__, int 
    *ldc, double *work, int *info);

int dorg2l_(int *m, int *n, int *k, double *
    a, int *lda, double *tau, double *work, int *info);

int dorg2r_(int *m, int *n, int *k, double *
    a, int *lda, double *tau, double *work, int *info);

int dorgbr_(char *vect, int *m, int *n, int *k, 
    double *a, int *lda, double *tau, double *work, 
    int *lwork, int *info);

int dorghr_(int *n, int *ilo, int *ihi, 
    double *a, int *lda, double *tau, double *work, 
    int *lwork, int *info);

int dorgl2_(int *m, int *n, int *k, double *
    a, int *lda, double *tau, double *work, int *info);


int dorgql_(int *m, int *n, int *k, double *
    a, int *lda, double *tau, double *work, int *lwork, 
    int *info);


int dorgr2_(int *m, int *n, int *k, double *
    a, int *lda, double *tau, double *work, int *info);

int dorgrq_(int *m, int *n, int *k, double *
    a, int *lda, double *tau, double *work, int *lwork, 
    int *info);

int dorgtr_(char *uplo, int *n, double *a, int *
    lda, double *tau, double *work, int *lwork, int *info);

int dorm2l_(char *side, char *trans, int *m, int *n, 
    int *k, double *a, int *lda, double *tau, double *
    c__, int *ldc, double *work, int *info);

int dorm2r_(char *side, char *trans, int *m, int *n, 
    int *k, double *a, int *lda, double *tau, double *
    c__, int *ldc, double *work, int *info);

int dormbr_(char *vect, char *side, char *trans, int *m, 
    int *n, int *k, double *a, int *lda, double *tau, 
    double *c__, int *ldc, double *work, int *lwork, 
    int *info);

int dormhr_(char *side, char *trans, int *m, int *n, 
    int *ilo, int *ihi, double *a, int *lda, double *
    tau, double *c__, int *ldc, double *work, int *lwork, 
    int *info);

int dorml2_(char *side, char *trans, int *m, int *n, 
    int *k, double *a, int *lda, double *tau, double *
    c__, int *ldc, double *work, int *info);


int dormql_(char *side, char *trans, int *m, int *n, 
    int *k, double *a, int *lda, double *tau, double *
    c__, int *ldc, double *work, int *lwork, int *info);


int dormr2_(char *side, char *trans, int *m, int *n, 
    int *k, double *a, int *lda, double *tau, double *
    c__, int *ldc, double *work, int *info);

int dormr3_(char *side, char *trans, int *m, int *n, 
    int *k, int *l, double *a, int *lda, double *tau, 
    double *c__, int *ldc, double *work, int *info);

int dormrq_(char *side, char *trans, int *m, int *n, 
    int *k, double *a, int *lda, double *tau, double *
    c__, int *ldc, double *work, int *lwork, int *info);

int dormrz_(char *side, char *trans, int *m, int *n, 
    int *k, int *l, double *a, int *lda, double *tau, 
    double *c__, int *ldc, double *work, int *lwork, 
    int *info);

int dormtr_(char *side, char *uplo, char *trans, int *m, 
    int *n, double *a, int *lda, double *tau, double *
    c__, int *ldc, double *work, int *lwork, int *info);

int dpbcon_(char *uplo, int *n, int *kd, double *
    ab, int *ldab, double *anorm, double *rcond, double *
    work, int *iwork, int *info);

int dpbequ_(char *uplo, int *n, int *kd, double *
    ab, int *ldab, double *s, double *scond, double *amax,
    int *info);

int dpbrfs_(char *uplo, int *n, int *kd, int *
    nrhs, double *ab, int *ldab, double *afb, int *ldafb, 
    double *b, int *ldb, double *x, int *ldx, double *
    ferr, double *berr, double *work, int *iwork, int *
    info);

int dpbstf_(char *uplo, int *n, int *kd, double *
    ab, int *ldab, int *info);

int dpbsv_(char *uplo, int *n, int *kd, int *
    nrhs, double *ab, int *ldab, double *b, int *ldb, 
    int *info);

int dpbsvx_(char *fact, char *uplo, int *n, int *kd, 
    int *nrhs, double *ab, int *ldab, double *afb, 
    int *ldafb, char *equed, double *s, double *b, int *
    ldb, double *x, int *ldx, double *rcond, double *ferr,
    double *berr, double *work, int *iwork, int *info);

int dpbtf2_(char *uplo, int *n, int *kd, double *
    ab, int *ldab, int *info);


int dpbtrs_(char *uplo, int *n, int *kd, int *
    nrhs, double *ab, int *ldab, double *b, int *ldb, 
    int *info);

int dpocon_(char *uplo, int *n, double *a, int *
    lda, double *anorm, double *rcond, double *work, int *
    iwork, int *info);

int dpoequ_(int *n, double *a, int *lda, 
    double *s, double *scond, double *amax, int *info);

int dporfs_(char *uplo, int *n, int *nrhs, 
    double *a, int *lda, double *af, int *ldaf, 
    double *b, int *ldb, double *x, int *ldx, double *
    ferr, double *berr, double *work, int *iwork, int *
    info);

int dposv_(char *uplo, int *n, int *nrhs, double 
    *a, int *lda, double *b, int *ldb, int *info);

int dposvx_(char *fact, char *uplo, int *n, int *
    nrhs, double *a, int *lda, double *af, int *ldaf, 
    char *equed, double *s, double *b, int *ldb, double *
    x, int *ldx, double *rcond, double *ferr, double *
    berr, double *work, int *iwork, int *info);

int dpotf2_(char *uplo, int *n, double *a, int *
    lda, int *info);



int dpotrs_(char *uplo, int *n, int *nrhs, 
    double *a, int *lda, double *b, int *ldb, int *
    info);

int dppcon_(char *uplo, int *n, double *ap, 
    double *anorm, double *rcond, double *work, int *
    iwork, int *info);

int dppequ_(char *uplo, int *n, double *ap, 
    double *s, double *scond, double *amax, int *info);

int dpprfs_(char *uplo, int *n, int *nrhs, 
    double *ap, double *afp, double *b, int *ldb, 
    double *x, int *ldx, double *ferr, double *berr, 
    double *work, int *iwork, int *info);

int dppsv_(char *uplo, int *n, int *nrhs, double 
    *ap, double *b, int *ldb, int *info);

int dppsvx_(char *fact, char *uplo, int *n, int *
    nrhs, double *ap, double *afp, char *equed, double *s, 
    double *b, int *ldb, double *x, int *ldx, double *
    rcond, double *ferr, double *berr, double *work, int *
    iwork, int *info);

int dpptrf_(char *uplo, int *n, double *ap, int *
    info);

int dpptri_(char *uplo, int *n, double *ap, int *
    info);

int dpptrs_(char *uplo, int *n, int *nrhs, 
    double *ap, double *b, int *ldb, int *info);

int dptcon_(int *n, double *d__, double *e, 
    double *anorm, double *rcond, double *work, int *info);

int dpteqr_(char *compz, int *n, double *d__, 
    double *e, double *z__, int *ldz, double *work, 
    int *info);

int dptrfs_(int *n, int *nrhs, double *d__, 
    double *e, double *df, double *ef, double *b, int 
    *ldb, double *x, int *ldx, double *ferr, double *berr,
    double *work, int *info);

int dptsv_(int *n, int *nrhs, double *d__, 
    double *e, double *b, int *ldb, int *info);

int dptsvx_(char *fact, int *n, int *nrhs, 
    double *d__, double *e, double *df, double *ef, 
    double *b, int *ldb, double *x, int *ldx, double *
    rcond, double *ferr, double *berr, double *work, int *
    info);


int dpttrs_(int *n, int *nrhs, double *d__, 
    double *e, double *b, int *ldb, int *info);

int dptts2_(int *n, int *nrhs, double *d__, 
    double *e, double *b, int *ldb);

int drscl_(int *n, double *sa, double *sx, 
    int *incx);

int dsbev_(char *jobz, char *uplo, int *n, int *kd, 
    double *ab, int *ldab, double *w, double *z__, 
    int *ldz, double *work, int *info);

int dsbevd_(char *jobz, char *uplo, int *n, int *kd, 
    double *ab, int *ldab, double *w, double *z__, 
    int *ldz, double *work, int *lwork, int *iwork, 
    int *liwork, int *info);

int dsbevx_(char *jobz, char *range, char *uplo, int *n, 
    int *kd, double *ab, int *ldab, double *q, int *
    ldq, double *vl, double *vu, int *il, int *iu, 
    double *abstol, int *m, double *w, double *z__, 
    int *ldz, double *work, int *iwork, int *ifail, 
    int *info);

int dsbgst_(char *vect, char *uplo, int *n, int *ka, 
    int *kb, double *ab, int *ldab, double *bb, int *
    ldbb, double *x, int *ldx, double *work, int *info);

int dsbgv_(char *jobz, char *uplo, int *n, int *ka, 
    int *kb, double *ab, int *ldab, double *bb, int *
    ldbb, double *w, double *z__, int *ldz, double *work, 
    int *info);

int dsbgvd_(char *jobz, char *uplo, int *n, int *ka, 
    int *kb, double *ab, int *ldab, double *bb, int *
    ldbb, double *w, double *z__, int *ldz, double *work, 
    int *lwork, int *iwork, int *liwork, int *info);

int dsbgvx_(char *jobz, char *range, char *uplo, int *n, 
    int *ka, int *kb, double *ab, int *ldab, double *
    bb, int *ldbb, double *q, int *ldq, double *vl, 
    double *vu, int *il, int *iu, double *abstol, int 
    *m, double *w, double *z__, int *ldz, double *work, 
    int *iwork, int *ifail, int *info);


int dspcon_(char *uplo, int *n, double *ap, int *
    ipiv, double *anorm, double *rcond, double *work, int 
    *iwork, int *info);

int dspev_(char *jobz, char *uplo, int *n, double *
    ap, double *w, double *z__, int *ldz, double *work, 
    int *info);

int dspevd_(char *jobz, char *uplo, int *n, double *
    ap, double *w, double *z__, int *ldz, double *work, 
    int *lwork, int *iwork, int *liwork, int *info);

int dspevx_(char *jobz, char *range, char *uplo, int *n, 
    double *ap, double *vl, double *vu, int *il, int *
    iu, double *abstol, int *m, double *w, double *z__, 
    int *ldz, double *work, int *iwork, int *ifail, 
    int *info);

int dspgst_(int *itype, char *uplo, int *n, 
    double *ap, double *bp, int *info);

int dspgv_(int *itype, char *jobz, char *uplo, int *
    n, double *ap, double *bp, double *w, double *z__, 
    int *ldz, double *work, int *info);

int dspgvd_(int *itype, char *jobz, char *uplo, int *
    n, double *ap, double *bp, double *w, double *z__, 
    int *ldz, double *work, int *lwork, int *iwork, 
    int *liwork, int *info);

int dspgvx_(int *itype, char *jobz, char *range, char *
    uplo, int *n, double *ap, double *bp, double *vl, 
    double *vu, int *il, int *iu, double *abstol, int 
    *m, double *w, double *z__, int *ldz, double *work, 
    int *iwork, int *ifail, int *info);

int dsprfs_(char *uplo, int *n, int *nrhs, 
    double *ap, double *afp, int *ipiv, double *b, 
    int *ldb, double *x, int *ldx, double *ferr, 
    double *berr, double *work, int *iwork, int *info);

int dspsv_(char *uplo, int *n, int *nrhs, double 
    *ap, int *ipiv, double *b, int *ldb, int *info);

int dspsvx_(char *fact, char *uplo, int *n, int *
    nrhs, double *ap, double *afp, int *ipiv, double *b, 
    int *ldb, double *x, int *ldx, double *rcond, 
    double *ferr, double *berr, double *work, int *iwork, 
    int *info);

int dsptrd_(char *uplo, int *n, double *ap, 
    double *d__, double *e, double *tau, int *info);

int dsptrf_(char *uplo, int *n, double *ap, int *
    ipiv, int *info);

int dsptri_(char *uplo, int *n, double *ap, int *
    ipiv, double *work, int *info);

int dsptrs_(char *uplo, int *n, int *nrhs, 
    double *ap, int *ipiv, double *b, int *ldb, int *
    info);

int dstebz_(char *range, char *order, int *n, double 
    *vl, double *vu, int *il, int *iu, double *abstol, 
    double *d__, double *e, int *m, int *nsplit, 
    double *w, int *iblock, int *isplit, double *work, 
    int *iwork, int *info);



int dstein_(int *n, double *d__, double *e, 
    int *m, double *w, int *iblock, int *isplit, 
    double *z__, int *ldz, double *work, int *iwork, 
    int *ifail, int *info);

int dsteqr_(char *compz, int *n, double *d__, 
    double *e, double *z__, int *ldz, double *work, 
    int *info);

int dsterf_(int *n, double *d__, double *e, 
    int *info);

int dstev_(char *jobz, int *n, double *d__, 
    double *e, double *z__, int *ldz, double *work, 
    int *info);

int dstevd_(char *jobz, int *n, double *d__, 
    double *e, double *z__, int *ldz, double *work, 
    int *lwork, int *iwork, int *liwork, int *info);

int dstevr_(char *jobz, char *range, int *n, double *
    d__, double *e, double *vl, double *vu, int *il, 
    int *iu, double *abstol, int *m, double *w, 
    double *z__, int *ldz, int *isuppz, double *work, 
    int *lwork, int *iwork, int *liwork, int *info);

int dstevx_(char *jobz, char *range, int *n, double *
    d__, double *e, double *vl, double *vu, int *il, 
    int *iu, double *abstol, int *m, double *w, 
    double *z__, int *ldz, double *work, int *iwork, 
    int *ifail, int *info);

int dsycon_(char *uplo, int *n, double *a, int *
    lda, int *ipiv, double *anorm, double *rcond, double *
    work, int *iwork, int *info);

int dsyev_(char *jobz, char *uplo, int *n, double *a,
    int *lda, double *w, double *work, int *lwork, 
    int *info);

int dsyevd_(char *jobz, char *uplo, int *n, double *
    a, int *lda, double *w, double *work, int *lwork, 
    int *iwork, int *liwork, int *info);

int dsyevr_(char *jobz, char *range, char *uplo, int *n, 
    double *a, int *lda, double *vl, double *vu, int *
    il, int *iu, double *abstol, int *m, double *w, 
    double *z__, int *ldz, int *isuppz, double *work, 
    int *lwork, int *iwork, int *liwork, int *info);

int dsyevx_(char *jobz, char *range, char *uplo, int *n, 
    double *a, int *lda, double *vl, double *vu, int *
    il, int *iu, double *abstol, int *m, double *w, 
    double *z__, int *ldz, double *work, int *lwork, 
    int *iwork, int *ifail, int *info);

int dsygs2_(int *itype, char *uplo, int *n, 
    double *a, int *lda, double *b, int *ldb, int *
    info);

int dsygst_(int *itype, char *uplo, int *n, 
    double *a, int *lda, double *b, int *ldb, int *
    info);

int dsygv_(int *itype, char *jobz, char *uplo, int *
    n, double *a, int *lda, double *b, int *ldb, 
    double *w, double *work, int *lwork, int *info);

int dsygvd_(int *itype, char *jobz, char *uplo, int *
    n, double *a, int *lda, double *b, int *ldb, 
    double *w, double *work, int *lwork, int *iwork, 
    int *liwork, int *info);

int dsygvx_(int *itype, char *jobz, char *range, char *
    uplo, int *n, double *a, int *lda, double *b, int 
    *ldb, double *vl, double *vu, int *il, int *iu, 
    double *abstol, int *m, double *w, double *z__, 
    int *ldz, double *work, int *lwork, int *iwork, 
    int *ifail, int *info);

int dsyrfs_(char *uplo, int *n, int *nrhs, 
    double *a, int *lda, double *af, int *ldaf, int *
    ipiv, double *b, int *ldb, double *x, int *ldx, 
    double *ferr, double *berr, double *work, int *iwork, 
    int *info);

int dsysv_(char *uplo, int *n, int *nrhs, double 
    *a, int *lda, int *ipiv, double *b, int *ldb, 
    double *work, int *lwork, int *info);

int dsysvx_(char *fact, char *uplo, int *n, int *
    nrhs, double *a, int *lda, double *af, int *ldaf, 
    int *ipiv, double *b, int *ldb, double *x, int *
    ldx, double *rcond, double *ferr, double *berr, 
    double *work, int *lwork, int *iwork, int *info);

int dsytd2_(char *uplo, int *n, double *a, int *
    lda, double *d__, double *e, double *tau, int *info);

int dsytf2_(char *uplo, int *n, double *a, int *
    lda, int *ipiv, int *info);



int dsytri_(char *uplo, int *n, double *a, int *
    lda, int *ipiv, double *work, int *info);

int dsytrs_(char *uplo, int *n, int *nrhs, 
    double *a, int *lda, int *ipiv, double *b, int *
    ldb, int *info);

int dtbcon_(char *norm, char *uplo, char *diag, int *n, 
    int *kd, double *ab, int *ldab, double *rcond, 
    double *work, int *iwork, int *info);

int dtbrfs_(char *uplo, char *trans, char *diag, int *n, 
    int *kd, int *nrhs, double *ab, int *ldab, double 
    *b, int *ldb, double *x, int *ldx, double *ferr, 
    double *berr, double *work, int *iwork, int *info);


int dtgevc_(char *side, char *howmny, int *select, 
    int *n, double *a, int *lda, double *b, int *ldb, 
    double *vl, int *ldvl, double *vr, int *ldvr, int 
    *mm, int *m, double *work, int *info);

int dtgex2_(int *wantq, int *wantz, int *n, 
    double *a, int *lda, double *b, int *ldb, double *
    q, int *ldq, double *z__, int *ldz, int *j1, int *
    n1, int *n2, double *work, int *lwork, int *info);

int dtgexc_(int *wantq, int *wantz, int *n, 
    double *a, int *lda, double *b, int *ldb, double *
    q, int *ldq, double *z__, int *ldz, int *ifst, 
    int *ilst, double *work, int *lwork, int *info);

int dtgsen_(int *ijob, int *wantq, int *wantz, 
    int *select, int *n, double *a, int *lda, double *
    b, int *ldb, double *alphar, double *alphai, double *
    beta, double *q, int *ldq, double *z__, int *ldz, 
    int *m, double *pl, double *pr, double *dif, 
    double *work, int *lwork, int *iwork, int *liwork, 
    int *info);

int dtgsja_(char *jobu, char *jobv, char *jobq, int *m, 
    int *p, int *n, int *k, int *l, double *a, 
    int *lda, double *b, int *ldb, double *tola, 
    double *tolb, double *alpha, double *beta, double *u, 
    int *ldu, double *v, int *ldv, double *q, int *
    ldq, double *work, int *ncycle, int *info);

int dtgsna_(char *job, char *howmny, int *select, 
    int *n, double *a, int *lda, double *b, int *ldb, 
    double *vl, int *ldvl, double *vr, int *ldvr, 
    double *s, double *dif, int *mm, int *m, double *
    work, int *lwork, int *iwork, int *info);

int dtgsy2_(char *trans, int *ijob, int *m, int *
    n, double *a, int *lda, double *b, int *ldb, 
    double *c__, int *ldc, double *d__, int *ldd, 
    double *e, int *lde, double *f, int *ldf, double *
    scale, double *rdsum, double *rdscal, int *iwork, int 
    *pq, int *info);

int dtgsyl_(char *trans, int *ijob, int *m, int *
    n, double *a, int *lda, double *b, int *ldb, 
    double *c__, int *ldc, double *d__, int *ldd, 
    double *e, int *lde, double *f, int *ldf, double *
    scale, double *dif, double *work, int *lwork, int *
    iwork, int *info);

int dtpcon_(char *norm, char *uplo, char *diag, int *n, 
    double *ap, double *rcond, double *work, int *iwork, 
    int *info);

int dtprfs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, double *ap, double *b, int *ldb, 
    double *x, int *ldx, double *ferr, double *berr, 
    double *work, int *iwork, int *info);

int dtptri_(char *uplo, char *diag, int *n, double *
    ap, int *info);

int dtptrs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, double *ap, double *b, int *ldb, int *
    info);

int dtrcon_(char *norm, char *uplo, char *diag, int *n, 
    double *a, int *lda, double *rcond, double *work, 
    int *iwork, int *info);

int dtrevc_(char *side, char *howmny, int *select, 
    int *n, double *t, int *ldt, double *vl, int *
    ldvl, double *vr, int *ldvr, int *mm, int *m, 
    double *work, int *info);

int dtrexc_(char *compq, int *n, double *t, int *
    ldt, double *q, int *ldq, int *ifst, int *ilst, 
    double *work, int *info);

int dtrrfs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, double *a, int *lda, double *b, int *
    ldb, double *x, int *ldx, double *ferr, double *berr, 
    double *work, int *iwork, int *info);

int dtrsen_(char *job, char *compq, int *select, int 
    *n, double *t, int *ldt, double *q, int *ldq, 
    double *wr, double *wi, int *m, double *s, double 
    *sep, double *work, int *lwork, int *iwork, int *
    liwork, int *info);

int dtrsna_(char *job, char *howmny, int *select, 
    int *n, double *t, int *ldt, double *vl, int *
    ldvl, double *vr, int *ldvr, double *s, double *sep, 
    int *mm, int *m, double *work, int *ldwork, int *
    iwork, int *info);

int dtrsyl_(char *trana, char *tranb, int *isgn, int 
    *m, int *n, double *a, int *lda, double *b, int *
    ldb, double *c__, int *ldc, double *scale, int *info);

int dtrti2_(char *uplo, char *diag, int *n, double *
    a, int *lda, int *info);


int dtrtrs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, double *a, int *lda, double *b, int *
    ldb, int *info);

int dtzrqf_(int *m, int *n, double *a, int *
    lda, double *tau, int *info);

int dtzrzf_(int *m, int *n, double *a, int *
    lda, double *tau, double *work, int *lwork, int *info);

int icmax1_(int *n, void *cx, int *incx);

int ieeeck_(int *ispec, float *zero, float *one);

int ilaenv_(int *ispec, char *name__, char *opts, int *n1, 
    int *n2, int *n3, int *n4, int name_len, int opts_len);

int izmax1_(int *n, void *cx, int *incx);


int sbdsqr_(char *uplo, int *n, int *ncvt, int *
    nru, int *ncc, float *d__, float *e, float *vt, int *ldvt, float *
    u, int *ldu, float *c__, int *ldc, float *work, int *info);

int sdisna_(char *job, int *m, int *n, float *d__, 
    float *sep, int *info);


int sgbcon_(char *norm, int *n, int *kl, int *ku,
    float *ab, int *ldab, int *ipiv, float *anorm, float *rcond, 
    float *work, int *iwork, int *info);

int sgbequ_(int *m, int *n, int *kl, int *ku,
    float *ab, int *ldab, float *r__, float *c__, float *rowcnd, float *
    colcnd, float *amax, int *info);

int sgbrfs_(char *trans, int *n, int *kl, int *
    ku, int *nrhs, float *ab, int *ldab, float *afb, int *ldafb,
    int *ipiv, float *b, int *ldb, float *x, int *ldx, float *
    ferr, float *berr, float *work, int *iwork, int *info);

int sgbsv_(int *n, int *kl, int *ku, int *
    nrhs, float *ab, int *ldab, int *ipiv, float *b, int *ldb, 
    int *info);

int sgbsvx_(char *fact, char *trans, int *n, int *kl,
    int *ku, int *nrhs, float *ab, int *ldab, float *afb, 
    int *ldafb, int *ipiv, char *equed, float *r__, float *c__, 
    float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr,
    float *berr, float *work, int *iwork, int *info);

int sgbtf2_(int *m, int *n, int *kl, int *ku,
    float *ab, int *ldab, int *ipiv, int *info);



int sgebak_(char *job, char *side, int *n, int *ilo, 
    int *ihi, float *scale, int *m, float *v, int *ldv, int 
    *info);

int sgebal_(char *job, int *n, float *a, int *lda, 
    int *ilo, int *ihi, float *scale, int *info);

int sgebd2_(int *m, int *n, float *a, int *lda, 
    float *d__, float *e, float *tauq, float *taup, float *work, int *info);


int sgecon_(char *norm, int *n, float *a, int *lda, 
    float *anorm, float *rcond, float *work, int *iwork, int *info);

int sgeequ_(int *m, int *n, float *a, int *lda, 
    float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, int 
    *info);

int sgeev_(char *jobvl, char *jobvr, int *n, float *a, 
    int *lda, float *wr, float *wi, float *vl, int *ldvl, float *vr, 
    int *ldvr, float *work, int *lwork, int *info);

int sgeevx_(char *balanc, char *jobvl, char *jobvr, char *
    sense, int *n, float *a, int *lda, float *wr, float *wi, float *
    vl, int *ldvl, float *vr, int *ldvr, int *ilo, int *
    ihi, float *scale, float *abnrm, float *rconde, float *rcondv, float *work,
    int *lwork, int *iwork, int *info);

int sgegs_(char *jobvsl, char *jobvsr, int *n, float *a, 
    int *lda, float *b, int *ldb, float *alphar, float *alphai, float 
    *beta, float *vsl, int *ldvsl, float *vsr, int *ldvsr, float *
    work, int *lwork, int *info);

int sgegv_(char *jobvl, char *jobvr, int *n, float *a, 
    int *lda, float *b, int *ldb, float *alphar, float *alphai, float 
    *beta, float *vl, int *ldvl, float *vr, int *ldvr, float *work, 
    int *lwork, int *info);

int sgehd2_(int *n, int *ilo, int *ihi, float *a, 
    int *lda, float *tau, float *work, int *info);

int sgehrd_(int *n, int *ilo, int *ihi, float *a, 
    int *lda, float *tau, float *work, int *lwork, int *info);

int sgelq2_(int *m, int *n, float *a, int *lda, 
    float *tau, float *work, int *info);


int sgels_(char *trans, int *m, int *n, int *
    nrhs, float *a, int *lda, float *b, int *ldb, float *work, 
    int *lwork, int *info);

int sgelsd_(int *m, int *n, int *nrhs, float *a, 
    int *lda, float *b, int *ldb, float *s, float *rcond, int *
    rank, float *work, int *lwork, int *iwork, int *info);

int sgelss_(int *m, int *n, int *nrhs, float *a, 
    int *lda, float *b, int *ldb, float *s, float *rcond, int *
    rank, float *work, int *lwork, int *info);

int sgelsx_(int *m, int *n, int *nrhs, float *a, 
    int *lda, float *b, int *ldb, int *jpvt, float *rcond, 
    int *rank, float *work, int *info);

int sgelsy_(int *m, int *n, int *nrhs, float *a, 
    int *lda, float *b, int *ldb, int *jpvt, float *rcond, 
    int *rank, float *work, int *lwork, int *info);

int sgeql2_(int *m, int *n, float *a, int *lda, 
    float *tau, float *work, int *info);

int sgeqlf_(int *m, int *n, float *a, int *lda, 
    float *tau, float *work, int *lwork, int *info);


int sgeqpf_(int *m, int *n, float *a, int *lda, 
    int *jpvt, float *tau, float *work, int *info);

int sgeqr2_(int *m, int *n, float *a, int *lda, 
    float *tau, float *work, int *info);


int sgerfs_(char *trans, int *n, int *nrhs, float *a, 
    int *lda, float *af, int *ldaf, int *ipiv, float *b, 
    int *ldb, float *x, int *ldx, float *ferr, float *berr, float *
    work, int *iwork, int *info);

int sgerq2_(int *m, int *n, float *a, int *lda, 
    float *tau, float *work, int *info);

int sgerqf_(int *m, int *n, float *a, int *lda, 
    float *tau, float *work, int *lwork, int *info);

int sgesc2_(int *n, float *a, int *lda, float *rhs, 
    int *ipiv, int *jpiv, float *scale);

int sgesdd_(char *jobz, int *m, int *n, float *a, 
    int *lda, float *s, float *u, int *ldu, float *vt, int *ldvt,
    float *work, int *lwork, int *iwork, int *info);

int sgesv_(int *n, int *nrhs, float *a, int *lda, 
    int *ipiv, float *b, int *ldb, int *info);

int sgesvd_(char *jobu, char *jobvt, int *m, int *n, 
    float *a, int *lda, float *s, float *u, int *ldu, float *vt, 
    int *ldvt, float *work, int *lwork, int *info);

int sgesvx_(char *fact, char *trans, int *n, int *
    nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, 
    char *equed, float *r__, float *c__, float *b, int *ldb, float *x, 
    int *ldx, float *rcond, float *ferr, float *berr, float *work, 
    int *iwork, int *info);

int sgetc2_(int *n, float *a, int *lda, int *ipiv,
    int *jpiv, int *info);

int sgetf2_(int *m, int *n, float *a, int *lda, 
    int *ipiv, int *info);




int sggbak_(char *job, char *side, int *n, int *ilo, 
    int *ihi, float *lscale, float *rscale, int *m, float *v, 
    int *ldv, int *info);

int sggbal_(char *job, int *n, float *a, int *lda, 
    float *b, int *ldb, int *ilo, int *ihi, float *lscale, float 
    *rscale, float *work, int *info);

int sggev_(char *jobvl, char *jobvr, int *n, float *a, 
    int *lda, float *b, int *ldb, float *alphar, float *alphai, float 
    *beta, float *vl, int *ldvl, float *vr, int *ldvr, float *work, 
    int *lwork, int *info);

int sggevx_(char *balanc, char *jobvl, char *jobvr, char *
    sense, int *n, float *a, int *lda, float *b, int *ldb, float 
    *alphar, float *alphai, float *beta, float *vl, int *ldvl, float *vr, 
    int *ldvr, int *ilo, int *ihi, float *lscale, float *rscale,
    float *abnrm, float *bbnrm, float *rconde, float *rcondv, float *work, 
    int *lwork, int *iwork, int *bwork, int *info);

int sggglm_(int *n, int *m, int *p, float *a, 
    int *lda, float *b, int *ldb, float *d__, float *x, float *y, 
    float *work, int *lwork, int *info);

int sgghrd_(char *compq, char *compz, int *n, int *
    ilo, int *ihi, float *a, int *lda, float *b, int *ldb, float 
    *q, int *ldq, float *z__, int *ldz, int *info);

int sgglse_(int *m, int *n, int *p, float *a, 
    int *lda, float *b, int *ldb, float *c__, float *d__, float *x, 
    float *work, int *lwork, int *info);

int sggqrf_(int *n, int *m, int *p, float *a, 
    int *lda, float *taua, float *b, int *ldb, float *taub, float *
    work, int *lwork, int *info);

int sggrqf_(int *m, int *p, int *n, float *a, 
    int *lda, float *taua, float *b, int *ldb, float *taub, float *
    work, int *lwork, int *info);

int sggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
    int *n, int *p, int *k, int *l, float *a, int *lda,
    float *b, int *ldb, float *alpha, float *beta, float *u, int *
    ldu, float *v, int *ldv, float *q, int *ldq, float *work, 
    int *iwork, int *info);

int sggsvp_(char *jobu, char *jobv, char *jobq, int *m, 
    int *p, int *n, float *a, int *lda, float *b, int *ldb, 
    float *tola, float *tolb, int *k, int *l, float *u, int *ldu,
    float *v, int *ldv, float *q, int *ldq, int *iwork, float *
    tau, float *work, int *info);

int sgtcon_(char *norm, int *n, float *dl, float *d__, 
    float *du, float *du2, int *ipiv, float *anorm, float *rcond, float *
    work, int *iwork, int *info);

int sgtrfs_(char *trans, int *n, int *nrhs, float *dl,
    float *d__, float *du, float *dlf, float *df, float *duf, float *du2, 
    int *ipiv, float *b, int *ldb, float *x, int *ldx, float *
    ferr, float *berr, float *work, int *iwork, int *info);

int sgtsv_(int *n, int *nrhs, float *dl, float *d__, 
    float *du, float *b, int *ldb, int *info);

int sgtsvx_(char *fact, char *trans, int *n, int *
    nrhs, float *dl, float *d__, float *du, float *dlf, float *df, float *duf, 
    float *du2, int *ipiv, float *b, int *ldb, float *x, int *
    ldx, float *rcond, float *ferr, float *berr, float *work, int *iwork, 
    int *info);



int sgtts2_(int *itrans, int *n, int *nrhs, float 
    *dl, float *d__, float *du, float *du2, int *ipiv, float *b, int *
    ldb);

int shgeqz_(char *job, char *compq, char *compz, int *n, 
    int *ilo, int *ihi, float *a, int *lda, float *b, int *
    ldb, float *alphar, float *alphai, float *beta, float *q, int *ldq, 
    float *z__, int *ldz, float *work, int *lwork, int *info);

int shsein_(char *side, char *eigsrc, char *initv, int *
    select, int *n, float *h__, int *ldh, float *wr, float *wi, float 
    *vl, int *ldvl, float *vr, int *ldvr, int *mm, int *m, 
    float *work, int *ifaill, int *ifailr, int *info);

int shseqr_(char *job, char *compz, int *n, int *ilo,
    int *ihi, float *h__, int *ldh, float *wr, float *wi, float *z__,
    int *ldz, float *work, int *lwork, int *info);

int slabad_(float *small, float *large);

int slabrd_(int *m, int *n, int *nb, float *a, 
    int *lda, float *d__, float *e, float *tauq, float *taup, float *x, 
    int *ldx, float *y, int *ldy);

int slacon_(int *n, float *v, float *x, int *isgn, 
    float *est, int *kase);


int sladiv_(float *a, float *b, float *c__, float *d__, float *p, 
    float *q);

int slae2_(float *a, float *b, float *c__, float *rt1, float *rt2);

int slaebz_(int *ijob, int *nitmax, int *n, 
    int *mmax, int *minp, int *nbmin, float *abstol, float *
    reltol, float *pivmin, float *d__, float *e, float *e2, int *nval, 
    float *ab, float *c__, int *mout, int *nab, float *work, int 
    *iwork, int *info);

int slaed0_(int *icompq, int *qsiz, int *n, float 
    *d__, float *e, float *q, int *ldq, float *qstore, int *ldqs, 
    float *work, int *iwork, int *info);

int slaed1_(int *n, float *d__, float *q, int *ldq, 
    int *indxq, float *rho, int *cutpnt, float *work, int *
    iwork, int *info);

int slaed2_(int *k, int *n, int *n1, float *d__, 
    float *q, int *ldq, int *indxq, float *rho, float *z__, float *
    dlamda, float *w, float *q2, int *indx, int *indxc, int *
    indxp, int *coltyp, int *info);

int slaed3_(int *k, int *n, int *n1, float *d__, 
    float *q, int *ldq, float *rho, float *dlamda, float *q2, int *
    indx, int *ctot, float *w, float *s, int *info);

int slaed4_(int *n, int *i__, float *d__, float *z__, 
    float *delta, float *rho, float *dlam, int *info);

int slaed5_(int *i__, float *d__, float *z__, float *delta, 
    float *rho, float *dlam);

int slaed6_(int *kniter, int *orgati, float *rho, 
    float *d__, float *z__, float *finit, float *tau, int *info);

int slaed7_(int *icompq, int *n, int *qsiz, 
    int *tlvls, int *curlvl, int *curpbm, float *d__, float *q, 
    int *ldq, int *indxq, float *rho, int *cutpnt, float *
    qstore, int *qptr, int *prmptr, int *perm, int *
    givptr, int *givcol, float *givnum, float *work, int *iwork, 
    int *info);

int slaed8_(int *icompq, int *k, int *n, int 
    *qsiz, float *d__, float *q, int *ldq, int *indxq, float *rho, 
    int *cutpnt, float *z__, float *dlamda, float *q2, int *ldq2, 
    float *w, int *perm, int *givptr, int *givcol, float *
    givnum, int *indxp, int *indx, int *info);

int slaed9_(int *k, int *kstart, int *kstop, 
    int *n, float *d__, float *q, int *ldq, float *rho, float *dlamda,
    float *w, float *s, int *lds, int *info);

int slaeda_(int *n, int *tlvls, int *curlvl, 
    int *curpbm, int *prmptr, int *perm, int *givptr, 
    int *givcol, float *givnum, float *q, int *qptr, float *z__, 
    float *ztemp, int *info);

int slaein_(int *rightv, int *noinit, int *n, 
    float *h__, int *ldh, float *wr, float *wi, float *vr, float *vi, float 
    *b, int *ldb, float *work, float *eps3, float *smlnum, float *bignum, 
    int *info);

int slaev2_(float *a, float *b, float *c__, float *rt1, float *
    rt2, float *cs1, float *sn1);

int slaexc_(int *wantq, int *n, float *t, int *
    ldt, float *q, int *ldq, int *j1, int *n1, int *n2, 
    float *work, int *info);

int slag2_(float *a, int *lda, float *b, int *ldb, 
    float *safmin, float *scale1, float *scale2, float *wr1, float *wr2, float *
    wi);

int slags2_(int *upper, float *a1, float *a2, float *a3, 
    float *b1, float *b2, float *b3, float *csu, float *snu, float *csv, float *
    snv, float *csq, float *snq);

int slagtf_(int *n, float *a, float *lambda, float *b, float 
    *c__, float *tol, float *d__, int *in, int *info);


int slagts_(int *job, int *n, float *a, float *b, float 
    *c__, float *d__, int *in, float *y, float *tol, int *info);

int slagv2_(float *a, int *lda, float *b, int *ldb, 
    float *alphar, float *alphai, float *beta, float *csl, float *snl, float *
    csr, float *snr);

int slahqr_(int *wantt, int *wantz, int *n, 
    int *ilo, int *ihi, float *h__, int *ldh, float *wr, float *
    wi, int *iloz, int *ihiz, float *z__, int *ldz, int *
    info);

int slahrd_(int *n, int *k, int *nb, float *a, 
    int *lda, float *tau, float *t, int *ldt, float *y, int *ldy);

int slaic1_(int *job, int *j, float *x, float *sest, 
    float *w, float *gamma, float *sestpr, float *s, float *c__);

int slaln2_(int *ltrans, int *na, int *nw, float *
    smin, float *ca, float *a, int *lda, float *d1, float *d2, float *b, 
    int *ldb, float *wr, float *wi, float *x, int *ldx, float *scale, 
    float *xnorm, int *info);

int slals0_(int *icompq, int *nl, int *nr, 
    int *sqre, int *nrhs, float *b, int *ldb, float *bx, 
    int *ldbx, int *perm, int *givptr, int *givcol, 
    int *ldgcol, float *givnum, int *ldgnum, float *poles, float *
    difl, float *difr, float *z__, int *k, float *c__, float *s, float *
    work, int *info);

int slalsa_(int *icompq, int *smlsiz, int *n, 
    int *nrhs, float *b, int *ldb, float *bx, int *ldbx, float *
    u, int *ldu, float *vt, int *k, float *difl, float *difr, float *
    z__, float *poles, int *givptr, int *givcol, int *ldgcol, 
    int *perm, float *givnum, float *c__, float *s, float *work, int *
    iwork, int *info);

int slalsd_(char *uplo, int *smlsiz, int *n, int 
    *nrhs, float *d__, float *e, float *b, int *ldb, float *rcond, 
    int *rank, float *work, int *iwork, int *info);

int slamc1_(int *beta, int *t, int *rnd, int 
    *ieee1);

int slamc2_(int *beta, int *t, int *rnd, float *
    eps, int *emin, float *rmin, int *emax, float *rmax);

int slamc4_(int *emin, float *start, int *base);

int slamc5_(int *beta, int *p, int *emin, 
    int *ieee, int *emax, float *rmax);

int slamrg_(int *n1, int *n2, float *a, int *
    strd1, int *strd2, int *index);

int slanv2_(float *a, float *b, float *c__, float *d__, float *
    rt1r, float *rt1i, float *rt2r, float *rt2i, float *cs, float *sn);

int slapll_(int *n, float *x, int *incx, float *y, 
    int *incy, float *ssmin);

int slapmt_(int *forwrd, int *m, int *n, float *x,
    int *ldx, int *k);

int slaqgb_(int *m, int *n, int *kl, int *ku,
    float *ab, int *ldab, float *r__, float *c__, float *rowcnd, float *
    colcnd, float *amax, char *equed);

int slaqge_(int *m, int *n, float *a, int *lda, 
    float *r__, float *c__, float *rowcnd, float *colcnd, float *amax, char *
    equed);

int slaqp2_(int *m, int *n, int *offset, float *a,
    int *lda, int *jpvt, float *tau, float *vn1, float *vn2, float *
    work);

int slaqps_(int *m, int *n, int *offset, int 
    *nb, int *kb, float *a, int *lda, int *jpvt, float *tau, 
    float *vn1, float *vn2, float *auxv, float *f, int *ldf);

int slaqsb_(char *uplo, int *n, int *kd, float *ab, 
    int *ldab, float *s, float *scond, float *amax, char *equed);

int slaqsp_(char *uplo, int *n, float *ap, float *s, float *
    scond, float *amax, char *equed);

int slaqsy_(char *uplo, int *n, float *a, int *lda, 
    float *s, float *scond, float *amax, char *equed);

int slaqtr_(int *ltran, int *lfloat, int *n, float 
    *t, int *ldt, float *b, float *w, float *scale, float *x, float *work, 
    int *info);

int slar1v_(int *n, int *b1, int *bn, float *
    sigma, float *d__, float *l, float *ld, float *lld, float *gersch, float *
    z__, float *ztz, float *mingma, int *r__, int *isuppz, float *
    work);

int slar2v_(int *n, float *x, float *y, float *z__, int 
    *incx, float *c__, float *s, int *incc);

int slarf_(char *side, int *m, int *n, float *v, 
    int *incv, float *tau, float *c__, int *ldc, float *work);

int slarfb_(char *side, char *trans, char *direct, char *
    storev, int *m, int *n, int *k, float *v, int *ldv, 
    float *t, int *ldt, float *c__, int *ldc, float *work, int *
    ldwork);

int slarfg_(int *n, float *alpha, float *x, int *incx, 
    float *tau);

int slarft_(char *direct, char *storev, int *n, int *
    k, float *v, int *ldv, float *tau, float *t, int *ldt);

int slarfx_(char *side, int *m, int *n, float *v, 
    float *tau, float *c__, int *ldc, float *work);

int slargv_(int *n, float *x, int *incx, float *y, 
    int *incy, float *c__, int *incc);

int slarnv_(int *idist, int *iseed, int *n, float 
    *x);

int slarrb_(int *n, float *d__, float *l, float *ld, float *
    lld, int *ifirst, int *ilast, float *sigma, float *reltol, float 
    *w, float *wgap, float *werr, float *work, int *iwork, int *info);

int slarre_(int *n, float *d__, float *e, float *tol, 
    int *nsplit, int *isplit, int *m, float *w, float *woff, 
    float *gersch, float *work, int *info);

int slarrf_(int *n, float *d__, float *l, float *ld, float *
    lld, int *ifirst, int *ilast, float *w, float *dplus, float *
    lplus, float *work, int *iwork, int *info);

int slarrv_(int *n, float *d__, float *l, int *isplit, 
    int *m, float *w, int *iblock, float *gersch, float *tol, float *
    z__, int *ldz, int *isuppz, float *work, int *iwork, 
    int *info);

int slartg_(float *f, float *g, float *cs, float *sn, float *r__);

int slartv_(int *n, float *x, int *incx, float *y, 
    int *incy, float *c__, float *s, int *incc);

int slaruv_(int *iseed, int *n, float *x);

int slarz_(char *side, int *m, int *n, int *l, 
    float *v, int *incv, float *tau, float *c__, int *ldc, float *
    work);

int slarzb_(char *side, char *trans, char *direct, char *
    storev, int *m, int *n, int *k, int *l, float *v, 
    int *ldv, float *t, int *ldt, float *c__, int *ldc, float *
    work, int *ldwork);

int slarzt_(char *direct, char *storev, int *n, int *
    k, float *v, int *ldv, float *tau, float *t, int *ldt);

int slas2_(float *f, float *g, float *h__, float *ssmin, float *
    ssmax);

int slascl_(char *type__, int *kl, int *ku, float *
    cfrom, float *cto, int *m, int *n, float *a, int *lda, 
    int *info);

int slasd0_(int *n, int *sqre, float *d__, float *e, 
    float *u, int *ldu, float *vt, int *ldvt, int *smlsiz, 
    int *iwork, float *work, int *info);

int slasd1_(int *nl, int *nr, int *sqre, float *
    d__, float *alpha, float *beta, float *u, int *ldu, float *vt, 
    int *ldvt, int *idxq, int *iwork, float *work, int *
    info);

int slasd2_(int *nl, int *nr, int *sqre, int 
    *k, float *d__, float *z__, float *alpha, float *beta, float *u, int *
    ldu, float *vt, int *ldvt, float *dsigma, float *u2, int *ldu2, 
    float *vt2, int *ldvt2, int *idxp, int *idx, int *idxc,
    int *idxq, int *coltyp, int *info);

int slasd3_(int *nl, int *nr, int *sqre, int 
    *k, float *d__, float *q, int *ldq, float *dsigma, float *u, int *
    ldu, float *u2, int *ldu2, float *vt, int *ldvt, float *vt2, 
    int *ldvt2, int *idxc, int *ctot, float *z__, int *
    info);

int slasd4_(int *n, int *i__, float *d__, float *z__, 
    float *delta, float *rho, float *sigma, float *work, int *info);

int slasd5_(int *i__, float *d__, float *z__, float *delta, 
    float *rho, float *dsigma, float *work);

int slasd6_(int *icompq, int *nl, int *nr, 
    int *sqre, float *d__, float *vf, float *vl, float *alpha, float *beta,
    int *idxq, int *perm, int *givptr, int *givcol, 
    int *ldgcol, float *givnum, int *ldgnum, float *poles, float *
    difl, float *difr, float *z__, int *k, float *c__, float *s, float *
    work, int *iwork, int *info);

int slasd7_(int *icompq, int *nl, int *nr, 
    int *sqre, int *k, float *d__, float *z__, float *zw, float *vf, 
    float *vfw, float *vl, float *vlw, float *alpha, float *beta, float *dsigma,
    int *idx, int *idxp, int *idxq, int *perm, int *
    givptr, int *givcol, int *ldgcol, float *givnum, int *
    ldgnum, float *c__, float *s, int *info);

int slasd8_(int *icompq, int *k, float *d__, float *
    z__, float *vf, float *vl, float *difl, float *difr, int *lddifr, 
    float *dsigma, float *work, int *info);

int slasd9_(int *icompq, int *ldu, int *k, float *
    d__, float *z__, float *vf, float *vl, float *difl, float *difr, float *
    dsigma, float *work, int *info);

int slasda_(int *icompq, int *smlsiz, int *n, 
    int *sqre, float *d__, float *e, float *u, int *ldu, float *vt, 
    int *k, float *difl, float *difr, float *z__, float *poles, int *
    givptr, int *givcol, int *ldgcol, int *perm, float *givnum,
    float *c__, float *s, float *work, int *iwork, int *info);

int slasdq_(char *uplo, int *sqre, int *n, int *
    ncvt, int *nru, int *ncc, float *d__, float *e, float *vt, 
    int *ldvt, float *u, int *ldu, float *c__, int *ldc, float *
    work, int *info);

int slasdt_(int *n, int *lvl, int *nd, int *
    inode, int *ndiml, int *ndimr, int *msub);

int slaset_(char *uplo, int *m, int *n, float *alpha, 
    float *beta, float *a, int *lda);

int slasq1_(int *n, float *d__, float *e, float *work, 
    int *info);

int slasq2_(int *n, float *z__, int *info);

int slasq3_(int *i0, int *n0, float *z__, int *pp,
    float *dmin__, float *sigma, float *desig, float *qmax, int *nfail, 
    int *iter, int *ndiv, int *ieee);

int slasq4_(int *i0, int *n0, float *z__, int *pp,
    int *n0in, float *dmin__, float *dmin1, float *dmin2, float *dn, 
    float *dn1, float *dn2, float *tau, int *ttype);

int slasq5_(int *i0, int *n0, float *z__, int *pp,
    float *tau, float *dmin__, float *dmin1, float *dmin2, float *dn, float *
    dnm1, float *dnm2, int *ieee);

int slasq6_(int *i0, int *n0, float *z__, int *pp,
    float *dmin__, float *dmin1, float *dmin2, float *dn, float *dnm1, float *
    dnm2);

int slasr_(char *side, char *pivot, char *direct, int *m,
    int *n, float *c__, float *s, float *a, int *lda);

int slasrt_(char *id, int *n, float *d__, int *info);

int slassq_(int *n, float *x, int *incx, float *scale, 
    float *sumsq);

int slasv2_(float *f, float *g, float *h__, float *ssmin, float *
    ssmax, float *snr, float *csr, float *snl, float *csl);

int slaswp_(int *n, float *a, int *lda, int *k1, 
    int *k2, int *ipiv, int *incx);

int slasy2_(int *ltranl, int *ltranr, int *isgn, 
    int *n1, int *n2, float *tl, int *ldtl, float *tr, int *
    ldtr, float *b, int *ldb, float *scale, float *x, int *ldx, float 
    *xnorm, int *info);

int slasyf_(char *uplo, int *n, int *nb, int *kb,
    float *a, int *lda, int *ipiv, float *w, int *ldw, int 
    *info);

int slatbs_(char *uplo, char *trans, char *diag, char *
    normin, int *n, int *kd, float *ab, int *ldab, float *x, 
    float *scale, float *cnorm, int *info);

int slatdf_(int *ijob, int *n, float *z__, int *
    ldz, float *rhs, float *rdsum, float *rdscal, int *ipiv, int *
    jpiv);

int slatps_(char *uplo, char *trans, char *diag, char *
    normin, int *n, float *ap, float *x, float *scale, float *cnorm, 
    int *info);

int slatrd_(char *uplo, int *n, int *nb, float *a, 
    int *lda, float *e, float *tau, float *w, int *ldw);

int slatrs_(char *uplo, char *trans, char *diag, char *
    normin, int *n, float *a, int *lda, float *x, float *scale, float 
    *cnorm, int *info);

int slatrz_(int *m, int *n, int *l, float *a, 
    int *lda, float *tau, float *work);

int slatzm_(char *side, int *m, int *n, float *v, 
    int *incv, float *tau, float *c1, float *c2, int *ldc, float *
    work);

int slauu2_(char *uplo, int *n, float *a, int *lda, 
    int *info);


int sopgtr_(char *uplo, int *n, float *ap, float *tau, 
    float *q, int *ldq, float *work, int *info);

int sopmtr_(char *side, char *uplo, char *trans, int *m, 
    int *n, float *ap, float *tau, float *c__, int *ldc, float *work, 
    int *info);

int sorg2l_(int *m, int *n, int *k, float *a, 
    int *lda, float *tau, float *work, int *info);

int sorg2r_(int *m, int *n, int *k, float *a, 
    int *lda, float *tau, float *work, int *info);

int sorgbr_(char *vect, int *m, int *n, int *k, 
    float *a, int *lda, float *tau, float *work, int *lwork, int 
    *info);

int sorghr_(int *n, int *ilo, int *ihi, float *a, 
    int *lda, float *tau, float *work, int *lwork, int *info);

int sorgl2_(int *m, int *n, int *k, float *a, 
    int *lda, float *tau, float *work, int *info);


int sorgql_(int *m, int *n, int *k, float *a, 
    int *lda, float *tau, float *work, int *lwork, int *info);


int sorgr2_(int *m, int *n, int *k, float *a, 
    int *lda, float *tau, float *work, int *info);

int sorgrq_(int *m, int *n, int *k, float *a, 
    int *lda, float *tau, float *work, int *lwork, int *info);

int sorgtr_(char *uplo, int *n, float *a, int *lda, 
    float *tau, float *work, int *lwork, int *info);

int sorm2l_(char *side, char *trans, int *m, int *n, 
    int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
    float *work, int *info);

int sorm2r_(char *side, char *trans, int *m, int *n, 
    int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
    float *work, int *info);

int sormbr_(char *vect, char *side, char *trans, int *m, 
    int *n, int *k, float *a, int *lda, float *tau, float *c__, 
    int *ldc, float *work, int *lwork, int *info);

int sormhr_(char *side, char *trans, int *m, int *n, 
    int *ilo, int *ihi, float *a, int *lda, float *tau, float *
    c__, int *ldc, float *work, int *lwork, int *info);

int sorml2_(char *side, char *trans, int *m, int *n, 
    int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
    float *work, int *info);


int sormql_(char *side, char *trans, int *m, int *n, 
    int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
    float *work, int *lwork, int *info);


int sormr2_(char *side, char *trans, int *m, int *n, 
    int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
    float *work, int *info);

int sormr3_(char *side, char *trans, int *m, int *n, 
    int *k, int *l, float *a, int *lda, float *tau, float *c__, 
    int *ldc, float *work, int *info);

int sormrq_(char *side, char *trans, int *m, int *n, 
    int *k, float *a, int *lda, float *tau, float *c__, int *ldc,
    float *work, int *lwork, int *info);

int sormrz_(char *side, char *trans, int *m, int *n, 
    int *k, int *l, float *a, int *lda, float *tau, float *c__, 
    int *ldc, float *work, int *lwork, int *info);

int sormtr_(char *side, char *uplo, char *trans, int *m, 
    int *n, float *a, int *lda, float *tau, float *c__, int *ldc,
    float *work, int *lwork, int *info);

int spbcon_(char *uplo, int *n, int *kd, float *ab, 
    int *ldab, float *anorm, float *rcond, float *work, int *iwork, 
    int *info);

int spbequ_(char *uplo, int *n, int *kd, float *ab, 
    int *ldab, float *s, float *scond, float *amax, int *info);

int spbrfs_(char *uplo, int *n, int *kd, int *
    nrhs, float *ab, int *ldab, float *afb, int *ldafb, float *b, 
    int *ldb, float *x, int *ldx, float *ferr, float *berr, float *
    work, int *iwork, int *info);

int spbstf_(char *uplo, int *n, int *kd, float *ab, 
    int *ldab, int *info);

int spbsv_(char *uplo, int *n, int *kd, int *
    nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);

int spbsvx_(char *fact, char *uplo, int *n, int *kd, 
    int *nrhs, float *ab, int *ldab, float *afb, int *ldafb, 
    char *equed, float *s, float *b, int *ldb, float *x, int *ldx, 
    float *rcond, float *ferr, float *berr, float *work, int *iwork, 
    int *info);

int spbtf2_(char *uplo, int *n, int *kd, float *ab, 
    int *ldab, int *info);


int spbtrs_(char *uplo, int *n, int *kd, int *
    nrhs, float *ab, int *ldab, float *b, int *ldb, int *info);

int spocon_(char *uplo, int *n, float *a, int *lda, 
    float *anorm, float *rcond, float *work, int *iwork, int *info);

int spoequ_(int *n, float *a, int *lda, float *s, float 
    *scond, float *amax, int *info);

int sporfs_(char *uplo, int *n, int *nrhs, float *a, 
    int *lda, float *af, int *ldaf, float *b, int *ldb, float *x,
    int *ldx, float *ferr, float *berr, float *work, int *iwork, 
    int *info);

int sposv_(char *uplo, int *n, int *nrhs, float *a, 
    int *lda, float *b, int *ldb, int *info);

int sposvx_(char *fact, char *uplo, int *n, int *
    nrhs, float *a, int *lda, float *af, int *ldaf, char *equed, 
    float *s, float *b, int *ldb, float *x, int *ldx, float *rcond, 
    float *ferr, float *berr, float *work, int *iwork, int *info);

int spotf2_(char *uplo, int *n, float *a, int *lda, 
    int *info);



int spotrs_(char *uplo, int *n, int *nrhs, float *a, 
    int *lda, float *b, int *ldb, int *info);

int sppcon_(char *uplo, int *n, float *ap, float *anorm, 
    float *rcond, float *work, int *iwork, int *info);

int sppequ_(char *uplo, int *n, float *ap, float *s, float *
    scond, float *amax, int *info);

int spprfs_(char *uplo, int *n, int *nrhs, float *ap, 
    float *afp, float *b, int *ldb, float *x, int *ldx, float *ferr, 
    float *berr, float *work, int *iwork, int *info);

int sppsv_(char *uplo, int *n, int *nrhs, float *ap, 
    float *b, int *ldb, int *info);

int sppsvx_(char *fact, char *uplo, int *n, int *
    nrhs, float *ap, float *afp, char *equed, float *s, float *b, int *
    ldb, float *x, int *ldx, float *rcond, float *ferr, float *berr, float 
    *work, int *iwork, int *info);

int spptrf_(char *uplo, int *n, float *ap, int *info);

int spptri_(char *uplo, int *n, float *ap, int *info);

int spptrs_(char *uplo, int *n, int *nrhs, float *ap, 
    float *b, int *ldb, int *info);

int sptcon_(int *n, float *d__, float *e, float *anorm, 
    float *rcond, float *work, int *info);

int spteqr_(char *compz, int *n, float *d__, float *e, 
    float *z__, int *ldz, float *work, int *info);

int sptrfs_(int *n, int *nrhs, float *d__, float *e, 
    float *df, float *ef, float *b, int *ldb, float *x, int *ldx, 
    float *ferr, float *berr, float *work, int *info);

int sptsv_(int *n, int *nrhs, float *d__, float *e, 
    float *b, int *ldb, int *info);

int sptsvx_(char *fact, int *n, int *nrhs, float *d__,
    float *e, float *df, float *ef, float *b, int *ldb, float *x, int 
    *ldx, float *rcond, float *ferr, float *berr, float *work, int *info);


int spttrs_(int *n, int *nrhs, float *d__, float *e, 
    float *b, int *ldb, int *info);

int sptts2_(int *n, int *nrhs, float *d__, float *e, 
    float *b, int *ldb);

int srscl_(int *n, float *sa, float *sx, int *incx);

int ssbev_(char *jobz, char *uplo, int *n, int *kd, 
    float *ab, int *ldab, float *w, float *z__, int *ldz, float *work,
    int *info);

int ssbevd_(char *jobz, char *uplo, int *n, int *kd, 
    float *ab, int *ldab, float *w, float *z__, int *ldz, float *work,
    int *lwork, int *iwork, int *liwork, int *info);

int ssbevx_(char *jobz, char *range, char *uplo, int *n, 
    int *kd, float *ab, int *ldab, float *q, int *ldq, float *vl,
    float *vu, int *il, int *iu, float *abstol, int *m, float *
    w, float *z__, int *ldz, float *work, int *iwork, int *
    ifail, int *info);

int ssbgst_(char *vect, char *uplo, int *n, int *ka, 
    int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *
    x, int *ldx, float *work, int *info);

int ssbgv_(char *jobz, char *uplo, int *n, int *ka, 
    int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *
    w, float *z__, int *ldz, float *work, int *info);

int ssbgvd_(char *jobz, char *uplo, int *n, int *ka, 
    int *kb, float *ab, int *ldab, float *bb, int *ldbb, float *
    w, float *z__, int *ldz, float *work, int *lwork, int *
    iwork, int *liwork, int *info);

int ssbgvx_(char *jobz, char *range, char *uplo, int *n, 
    int *ka, int *kb, float *ab, int *ldab, float *bb, int *
    ldbb, float *q, int *ldq, float *vl, float *vu, int *il, int 
    *iu, float *abstol, int *m, float *w, float *z__, int *ldz, float 
    *work, int *iwork, int *ifail, int *info);


int sspcon_(char *uplo, int *n, float *ap, int *ipiv, 
    float *anorm, float *rcond, float *work, int *iwork, int *info);

int sspev_(char *jobz, char *uplo, int *n, float *ap, 
    float *w, float *z__, int *ldz, float *work, int *info);

int sspevd_(char *jobz, char *uplo, int *n, float *ap, 
    float *w, float *z__, int *ldz, float *work, int *lwork, int 
    *iwork, int *liwork, int *info);

int sspevx_(char *jobz, char *range, char *uplo, int *n, 
    float *ap, float *vl, float *vu, int *il, int *iu, float *abstol, 
    int *m, float *w, float *z__, int *ldz, float *work, int *
    iwork, int *ifail, int *info);

int sspgst_(int *itype, char *uplo, int *n, float *ap,
    float *bp, int *info);

int sspgv_(int *itype, char *jobz, char *uplo, int *
    n, float *ap, float *bp, float *w, float *z__, int *ldz, float *work, 
    int *info);

int sspgvd_(int *itype, char *jobz, char *uplo, int *
    n, float *ap, float *bp, float *w, float *z__, int *ldz, float *work, 
    int *lwork, int *iwork, int *liwork, int *info);

int sspgvx_(int *itype, char *jobz, char *range, char *
    uplo, int *n, float *ap, float *bp, float *vl, float *vu, int *il,
    int *iu, float *abstol, int *m, float *w, float *z__, int *
    ldz, float *work, int *iwork, int *ifail, int *info);

int ssprfs_(char *uplo, int *n, int *nrhs, float *ap, 
    float *afp, int *ipiv, float *b, int *ldb, float *x, int *
    ldx, float *ferr, float *berr, float *work, int *iwork, int *
    info);

int sspsv_(char *uplo, int *n, int *nrhs, float *ap, 
    int *ipiv, float *b, int *ldb, int *info);

int sspsvx_(char *fact, char *uplo, int *n, int *
    nrhs, float *ap, float *afp, int *ipiv, float *b, int *ldb, float 
    *x, int *ldx, float *rcond, float *ferr, float *berr, float *work, 
    int *iwork, int *info);

int ssptrd_(char *uplo, int *n, float *ap, float *d__, 
    float *e, float *tau, int *info);

int ssptrf_(char *uplo, int *n, float *ap, int *ipiv, 
    int *info);

int ssptri_(char *uplo, int *n, float *ap, int *ipiv, 
    float *work, int *info);

int ssptrs_(char *uplo, int *n, int *nrhs, float *ap, 
    int *ipiv, float *b, int *ldb, int *info);

int sstebz_(char *range, char *order, int *n, float *vl, 
    float *vu, int *il, int *iu, float *abstol, float *d__, float *e, 
    int *m, int *nsplit, float *w, int *iblock, int *
    isplit, float *work, int *iwork, int *info);



int sstein_(int *n, float *d__, float *e, int *m, float 
    *w, int *iblock, int *isplit, float *z__, int *ldz, float *
    work, int *iwork, int *ifail, int *info);

int ssteqr_(char *compz, int *n, float *d__, float *e, 
    float *z__, int *ldz, float *work, int *info);

int ssterf_(int *n, float *d__, float *e, int *info);

int sstev_(char *jobz, int *n, float *d__, float *e, float *
    z__, int *ldz, float *work, int *info);

int sstevd_(char *jobz, int *n, float *d__, float *e, float 
    *z__, int *ldz, float *work, int *lwork, int *iwork, 
    int *liwork, int *info);

int sstevr_(char *jobz, char *range, int *n, float *d__, 
    float *e, float *vl, float *vu, int *il, int *iu, float *abstol, 
    int *m, float *w, float *z__, int *ldz, int *isuppz, float *
    work, int *lwork, int *iwork, int *liwork, int *info);

int sstevx_(char *jobz, char *range, int *n, float *d__, 
    float *e, float *vl, float *vu, int *il, int *iu, float *abstol, 
    int *m, float *w, float *z__, int *ldz, float *work, int *
    iwork, int *ifail, int *info);

int ssycon_(char *uplo, int *n, float *a, int *lda, 
    int *ipiv, float *anorm, float *rcond, float *work, int *iwork, 
    int *info);

int ssyev_(char *jobz, char *uplo, int *n, float *a, 
    int *lda, float *w, float *work, int *lwork, int *info);

int ssyevd_(char *jobz, char *uplo, int *n, float *a, 
    int *lda, float *w, float *work, int *lwork, int *iwork, 
    int *liwork, int *info);

int ssyevr_(char *jobz, char *range, char *uplo, int *n, 
    float *a, int *lda, float *vl, float *vu, int *il, int *iu, 
    float *abstol, int *m, float *w, float *z__, int *ldz, int *
    isuppz, float *work, int *lwork, int *iwork, int *liwork, 
    int *info);

int ssyevx_(char *jobz, char *range, char *uplo, int *n, 
    float *a, int *lda, float *vl, float *vu, int *il, int *iu, 
    float *abstol, int *m, float *w, float *z__, int *ldz, float *
    work, int *lwork, int *iwork, int *ifail, int *info);

int ssygs2_(int *itype, char *uplo, int *n, float *a, 
    int *lda, float *b, int *ldb, int *info);

int ssygst_(int *itype, char *uplo, int *n, float *a, 
    int *lda, float *b, int *ldb, int *info);

int ssygv_(int *itype, char *jobz, char *uplo, int *
    n, float *a, int *lda, float *b, int *ldb, float *w, float *work, 
    int *lwork, int *info);

int ssygvd_(int *itype, char *jobz, char *uplo, int *
    n, float *a, int *lda, float *b, int *ldb, float *w, float *work, 
    int *lwork, int *iwork, int *liwork, int *info);

int ssygvx_(int *itype, char *jobz, char *range, char *
    uplo, int *n, float *a, int *lda, float *b, int *ldb, float *
    vl, float *vu, int *il, int *iu, float *abstol, int *m, 
    float *w, float *z__, int *ldz, float *work, int *lwork, int 
    *iwork, int *ifail, int *info);

int ssyrfs_(char *uplo, int *n, int *nrhs, float *a, 
    int *lda, float *af, int *ldaf, int *ipiv, float *b, 
    int *ldb, float *x, int *ldx, float *ferr, float *berr, float *
    work, int *iwork, int *info);

int ssysv_(char *uplo, int *n, int *nrhs, float *a, 
    int *lda, int *ipiv, float *b, int *ldb, float *work, 
    int *lwork, int *info);

int ssysvx_(char *fact, char *uplo, int *n, int *
    nrhs, float *a, int *lda, float *af, int *ldaf, int *ipiv, 
    float *b, int *ldb, float *x, int *ldx, float *rcond, float *ferr,
    float *berr, float *work, int *lwork, int *iwork, int *
    info);

int ssytd2_(char *uplo, int *n, float *a, int *lda, 
    float *d__, float *e, float *tau, int *info);

int ssytf2_(char *uplo, int *n, float *a, int *lda, 
    int *ipiv, int *info);



int ssytri_(char *uplo, int *n, float *a, int *lda, 
    int *ipiv, float *work, int *info);

int ssytrs_(char *uplo, int *n, int *nrhs, float *a, 
    int *lda, int *ipiv, float *b, int *ldb, int *info);

int stbcon_(char *norm, char *uplo, char *diag, int *n, 
    int *kd, float *ab, int *ldab, float *rcond, float *work, 
    int *iwork, int *info);

int stbrfs_(char *uplo, char *trans, char *diag, int *n, 
    int *kd, int *nrhs, float *ab, int *ldab, float *b, int 
    *ldb, float *x, int *ldx, float *ferr, float *berr, float *work, 
    int *iwork, int *info);


int stgevc_(char *side, char *howmny, int *select, 
    int *n, float *a, int *lda, float *b, int *ldb, float *vl, 
    int *ldvl, float *vr, int *ldvr, int *mm, int *m, float 
    *work, int *info);

int stgex2_(int *wantq, int *wantz, int *n, float 
    *a, int *lda, float *b, int *ldb, float *q, int *ldq, float *
    z__, int *ldz, int *j1, int *n1, int *n2, float *work, 
    int *lwork, int *info);

int stgexc_(int *wantq, int *wantz, int *n, float 
    *a, int *lda, float *b, int *ldb, float *q, int *ldq, float *
    z__, int *ldz, int *ifst, int *ilst, float *work, int *
    lwork, int *info);

int stgsen_(int *ijob, int *wantq, int *wantz, 
    int *select, int *n, float *a, int *lda, float *b, int *
    ldb, float *alphar, float *alphai, float *beta, float *q, int *ldq, 
    float *z__, int *ldz, int *m, float *pl, float *pr, float *dif, 
    float *work, int *lwork, int *iwork, int *liwork, int *
    info);

int stgsja_(char *jobu, char *jobv, char *jobq, int *m, 
    int *p, int *n, int *k, int *l, float *a, int *lda,
    float *b, int *ldb, float *tola, float *tolb, float *alpha, float *
    beta, float *u, int *ldu, float *v, int *ldv, float *q, int *
    ldq, float *work, int *ncycle, int *info);

int stgsna_(char *job, char *howmny, int *select, 
    int *n, float *a, int *lda, float *b, int *ldb, float *vl, 
    int *ldvl, float *vr, int *ldvr, float *s, float *dif, int *
    mm, int *m, float *work, int *lwork, int *iwork, int *
    info);

int stgsy2_(char *trans, int *ijob, int *m, int *
    n, float *a, int *lda, float *b, int *ldb, float *c__, int *
    ldc, float *d__, int *ldd, float *e, int *lde, float *f, int 
    *ldf, float *scale, float *rdsum, float *rdscal, int *iwork, int 
    *pq, int *info);

int stgsyl_(char *trans, int *ijob, int *m, int *
    n, float *a, int *lda, float *b, int *ldb, float *c__, int *
    ldc, float *d__, int *ldd, float *e, int *lde, float *f, int 
    *ldf, float *scale, float *dif, float *work, int *lwork, int *
    iwork, int *info);

int stpcon_(char *norm, char *uplo, char *diag, int *n, 
    float *ap, float *rcond, float *work, int *iwork, int *info);

int stprfs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, float *ap, float *b, int *ldb, float *x, int *ldx,
    float *ferr, float *berr, float *work, int *iwork, int *info);

int stptri_(char *uplo, char *diag, int *n, float *ap, 
    int *info);

int stptrs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, float *ap, float *b, int *ldb, int *info);

int strcon_(char *norm, char *uplo, char *diag, int *n, 
    float *a, int *lda, float *rcond, float *work, int *iwork, 
    int *info);

int strevc_(char *side, char *howmny, int *select, 
    int *n, float *t, int *ldt, float *vl, int *ldvl, float *vr, 
    int *ldvr, int *mm, int *m, float *work, int *info);

int strexc_(char *compq, int *n, float *t, int *ldt, 
    float *q, int *ldq, int *ifst, int *ilst, float *work, 
    int *info);

int strrfs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, float *a, int *lda, float *b, int *ldb, float *x, 
    int *ldx, float *ferr, float *berr, float *work, int *iwork, 
    int *info);

int strsen_(char *job, char *compq, int *select, int 
    *n, float *t, int *ldt, float *q, int *ldq, float *wr, float *wi, 
    int *m, float *s, float *sep, float *work, int *lwork, int *
    iwork, int *liwork, int *info);

int strsna_(char *job, char *howmny, int *select, 
    int *n, float *t, int *ldt, float *vl, int *ldvl, float *vr, 
    int *ldvr, float *s, float *sep, int *mm, int *m, float *
    work, int *ldwork, int *iwork, int *info);

int strsyl_(char *trana, char *tranb, int *isgn, int 
    *m, int *n, float *a, int *lda, float *b, int *ldb, float *
    c__, int *ldc, float *scale, int *info);

int strti2_(char *uplo, char *diag, int *n, float *a, 
    int *lda, int *info);


int strtrs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, float *a, int *lda, float *b, int *ldb, int *
    info);

int stzrqf_(int *m, int *n, float *a, int *lda, 
    float *tau, int *info);

int stzrzf_(int *m, int *n, float *a, int *lda, 
    float *tau, float *work, int *lwork, int *info);

int xerbla_(char *srname, int *info);

int zbdsqr_(char *uplo, int *n, int *ncvt, int *
    nru, int *ncc, double *d__, double *e, void *vt, 
    int *ldvt, void *u, int *ldu, void *c__, 
    int *ldc, double *rwork, int *info);

//int zdrot_(int *n, void *cx, int *incx, 
    //void *cy, int *incy, double *c__, double *s);

int zdrscl_(int *n, double *sa, void *sx, 
    int *incx);


int zgbcon_(char *norm, int *n, int *kl, int *ku,
    void *ab, int *ldab, int *ipiv, double *anorm, 
    double *rcond, void *work, double *rwork, int *
    info);

int zgbequ_(int *m, int *n, int *kl, int *ku,
    void *ab, int *ldab, double *r__, double *c__, 
    double *rowcnd, double *colcnd, double *amax, int *
    info);

int zgbrfs_(char *trans, int *n, int *kl, int *
    ku, int *nrhs, void *ab, int *ldab, void *
    afb, int *ldafb, int *ipiv, void *b, int *ldb, 
    void *x, int *ldx, double *ferr, double *berr, 
    void *work, double *rwork, int *info);

int zgbsv_(int *n, int *kl, int *ku, int *
    nrhs, void *ab, int *ldab, int *ipiv, void *
    b, int *ldb, int *info);

int zgbsvx_(char *fact, char *trans, int *n, int *kl,
    int *ku, int *nrhs, void *ab, int *ldab, 
    void *afb, int *ldafb, int *ipiv, char *equed, 
    double *r__, double *c__, void *b, int *ldb, 
    void *x, int *ldx, double *rcond, double *ferr, 
    double *berr, void *work, double *rwork, int *
    info);

int zgbtf2_(int *m, int *n, int *kl, int *ku,
    void *ab, int *ldab, int *ipiv, int *info);



int zgebak_(char *job, char *side, int *n, int *ilo, 
    int *ihi, double *scale, int *m, void *v, 
    int *ldv, int *info);

int zgebal_(char *job, int *n, void *a, int 
    *lda, int *ilo, int *ihi, double *scale, int *info);

int zgebd2_(int *m, int *n, void *a, 
    int *lda, double *d__, double *e, void *tauq, 
    void *taup, void *work, int *info);


int zgecon_(char *norm, int *n, void *a, 
    int *lda, double *anorm, double *rcond, void *
    work, double *rwork, int *info);

int zgeequ_(int *m, int *n, void *a, 
    int *lda, double *r__, double *c__, double *rowcnd, 
    double *colcnd, double *amax, int *info);

int zgeev_(char *jobvl, char *jobvr, int *n, 
    void *a, int *lda, void *w, void *vl, 
    int *ldvl, void *vr, int *ldvr, void *work, 
    int *lwork, double *rwork, int *info);

int zgeevx_(char *balanc, char *jobvl, char *jobvr, char *
    sense, int *n, void *a, int *lda, void *w, 
    void *vl, int *ldvl, void *vr, int *ldvr, 
    int *ilo, int *ihi, double *scale, double *abnrm, 
    double *rconde, double *rcondv, void *work, int *
    lwork, double *rwork, int *info);

int zgegs_(char *jobvsl, char *jobvsr, int *n, 
    void *a, int *lda, void *b, int *ldb, 
    void *alpha, void *beta, void *vsl, 
    int *ldvsl, void *vsr, int *ldvsr, void *
    work, int *lwork, double *rwork, int *info);

int zgegv_(char *jobvl, char *jobvr, int *n, 
    void *a, int *lda, void *b, int *ldb, 
    void *alpha, void *beta, void *vl, int 
    *ldvl, void *vr, int *ldvr, void *work, int 
    *lwork, double *rwork, int *info);

int zgehd2_(int *n, int *ilo, int *ihi, 
    void *a, int *lda, void *tau, void *
    work, int *info);

int zgehrd_(int *n, int *ilo, int *ihi, 
    void *a, int *lda, void *tau, void *
    work, int *lwork, int *info);

int zgelq2_(int *m, int *n, void *a, 
    int *lda, void *tau, void *work, int *info);


int zgels_(char *trans, int *m, int *n, int *
    nrhs, void *a, int *lda, void *b, int *ldb, 
    void *work, int *lwork, int *info);

int zgelsx_(int *m, int *n, int *nrhs, 
    void *a, int *lda, void *b, int *ldb, 
    int *jpvt, double *rcond, int *rank, void *work, 
    double *rwork, int *info);

int zgelsy_(int *m, int *n, int *nrhs, 
    void *a, int *lda, void *b, int *ldb, 
    int *jpvt, double *rcond, int *rank, void *work, 
    int *lwork, double *rwork, int *info);

int zgeql2_(int *m, int *n, void *a, 
    int *lda, void *tau, void *work, int *info);

int zgeqlf_(int *m, int *n, void *a, 
    int *lda, void *tau, void *work, int *lwork,
    int *info);


int zgeqpf_(int *m, int *n, void *a, 
    int *lda, int *jpvt, void *tau, void *work, 
    double *rwork, int *info);

int zgeqr2_(int *m, int *n, void *a, 
    int *lda, void *tau, void *work, int *info);


int zgerfs_(char *trans, int *n, int *nrhs, 
    void *a, int *lda, void *af, int *ldaf, 
    int *ipiv, void *b, int *ldb, void *x, 
    int *ldx, double *ferr, double *berr, void *work,
    double *rwork, int *info);

int zgerq2_(int *m, int *n, void *a, 
    int *lda, void *tau, void *work, int *info);

int zgerqf_(int *m, int *n, void *a, 
    int *lda, void *tau, void *work, int *lwork,
    int *info);

int zgesc2_(int *n, void *a, int *lda, 
    void *rhs, int *ipiv, int *jpiv, double *scale);

int zgesv_(int *n, int *nrhs, void *a, 
    int *lda, int *ipiv, void *b, int *ldb, int *
    info);

int zgesvx_(char *fact, char *trans, int *n, int *
    nrhs, void *a, int *lda, void *af, int *
    ldaf, int *ipiv, char *equed, double *r__, double *c__, 
    void *b, int *ldb, void *x, int *ldx, 
    double *rcond, double *ferr, double *berr, void *
    work, double *rwork, int *info);

int zgetc2_(int *n, void *a, int *lda, 
    int *ipiv, int *jpiv, int *info);

int zgetf2_(int *m, int *n, void *a, 
    int *lda, int *ipiv, int *info);




int zggbak_(char *job, char *side, int *n, int *ilo, 
    int *ihi, double *lscale, double *rscale, int *m, 
    void *v, int *ldv, int *info);

int zggbal_(char *job, int *n, void *a, int 
    *lda, void *b, int *ldb, int *ilo, int *ihi, 
    double *lscale, double *rscale, double *work, int *
    info);

int zggev_(char *jobvl, char *jobvr, int *n, 
    void *a, int *lda, void *b, int *ldb, 
    void *alpha, void *beta, void *vl, int 
    *ldvl, void *vr, int *ldvr, void *work, int 
    *lwork, double *rwork, int *info);

int zggevx_(char *balanc, char *jobvl, char *jobvr, char *
    sense, int *n, void *a, int *lda, void *b, 
    int *ldb, void *alpha, void *beta, 
    void *vl, int *ldvl, void *vr, int *ldvr, 
    int *ilo, int *ihi, double *lscale, double *rscale, 
    double *abnrm, double *bbnrm, double *rconde, double *
    rcondv, void *work, int *lwork, double *rwork, 
    int *iwork, int *bwork, int *info);

int zggglm_(int *n, int *m, int *p, 
    void *a, int *lda, void *b, int *ldb, 
    void *d__, void *x, void *y, void 
    *work, int *lwork, int *info);

int zgghrd_(char *compq, char *compz, int *n, int *
    ilo, int *ihi, void *a, int *lda, void *b, 
    int *ldb, void *q, int *ldq, void *z__, 
    int *ldz, int *info);

int zgglse_(int *m, int *n, int *p, 
    void *a, int *lda, void *b, int *ldb, 
    void *c__, void *d__, void *x, 
    void *work, int *lwork, int *info);

int zggqrf_(int *n, int *m, int *p, 
    void *a, int *lda, void *taua, void *b,
    int *ldb, void *taub, void *work, int *
    lwork, int *info);

int zggrqf_(int *m, int *p, int *n, 
    void *a, int *lda, void *taua, void *b,
    int *ldb, void *taub, void *work, int *
    lwork, int *info);

int zggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
    int *n, int *p, int *k, int *l, void *a, 
    int *lda, void *b, int *ldb, double *alpha, 
    double *beta, void *u, int *ldu, void *v, 
    int *ldv, void *q, int *ldq, void *work, 
    double *rwork, int *iwork, int *info);

int zggsvp_(char *jobu, char *jobv, char *jobq, int *m, 
    int *p, int *n, void *a, int *lda, void 
    *b, int *ldb, double *tola, double *tolb, int *k, 
    int *l, void *u, int *ldu, void *v, int 
    *ldv, void *q, int *ldq, int *iwork, double *
    rwork, void *tau, void *work, int *info);

int zgtcon_(char *norm, int *n, void *dl, 
    void *d__, void *du, void *du2, int *
    ipiv, double *anorm, double *rcond, void *work, 
    int *info);

int zgtrfs_(char *trans, int *n, int *nrhs, 
    void *dl, void *d__, void *du, 
    void *dlf, void *df, void *duf, 
    void *du2, int *ipiv, void *b, int *ldb, 
    void *x, int *ldx, double *ferr, double *berr, 
    void *work, double *rwork, int *info);

int zgtsv_(int *n, int *nrhs, void *dl, 
    void *d__, void *du, void *b, int *ldb,
    int *info);

int zgtsvx_(char *fact, char *trans, int *n, int *
    nrhs, void *dl, void *d__, void *du, 
    void *dlf, void *df, void *duf, 
    void *du2, int *ipiv, void *b, int *ldb, 
    void *x, int *ldx, double *rcond, double *ferr, 
    double *berr, void *work, double *rwork, int *
    info);



int zgtts2_(int *itrans, int *n, int *nrhs, 
    void *dl, void *d__, void *du, 
    void *du2, int *ipiv, void *b, int *ldb);

int zhbev_(char *jobz, char *uplo, int *n, int *kd, 
    void *ab, int *ldab, double *w, void *z__, 
    int *ldz, void *work, double *rwork, int *info);

int zhbevd_(char *jobz, char *uplo, int *n, int *kd, 
    void *ab, int *ldab, double *w, void *z__, 
    int *ldz, void *work, int *lwork, double *rwork, 
    int *lrwork, int *iwork, int *liwork, int *info);

int zhbevx_(char *jobz, char *range, char *uplo, int *n, 
    int *kd, void *ab, int *ldab, void *q, 
    int *ldq, double *vl, double *vu, int *il, int *
    iu, double *abstol, int *m, double *w, void *z__,
    int *ldz, void *work, double *rwork, int *iwork,
    int *ifail, int *info);

int zhbgst_(char *vect, char *uplo, int *n, int *ka, 
    int *kb, void *ab, int *ldab, void *bb, 
    int *ldbb, void *x, int *ldx, void *work, 
    double *rwork, int *info);

int zhbgv_(char *jobz, char *uplo, int *n, int *ka, 
    int *kb, void *ab, int *ldab, void *bb, 
    int *ldbb, double *w, void *z__, int *ldz, 
    void *work, double *rwork, int *info);

int zhbgvx_(char *jobz, char *range, char *uplo, int *n, 
    int *ka, int *kb, void *ab, int *ldab, 
    void *bb, int *ldbb, void *q, int *ldq, 
    double *vl, double *vu, int *il, int *iu, double *
    abstol, int *m, double *w, void *z__, int *ldz, 
    void *work, double *rwork, int *iwork, int *
    ifail, int *info);


int zhecon_(char *uplo, int *n, void *a, 
    int *lda, int *ipiv, double *anorm, double *rcond, 
    void *work, int *info);

int zheev_(char *jobz, char *uplo, int *n, void 
    *a, int *lda, double *w, void *work, int *lwork, 
    double *rwork, int *info);

int zheevd_(char *jobz, char *uplo, int *n, 
    void *a, int *lda, double *w, void *work, 
    int *lwork, double *rwork, int *lrwork, int *iwork, 
    int *liwork, int *info);

int zheevr_(char *jobz, char *range, char *uplo, int *n, 
    void *a, int *lda, double *vl, double *vu, 
    int *il, int *iu, double *abstol, int *m, double *
    w, void *z__, int *ldz, int *isuppz, void *
    work, int *lwork, double *rwork, int *lrwork, int *
    iwork, int *liwork, int *info);

int zheevx_(char *jobz, char *range, char *uplo, int *n, 
    void *a, int *lda, double *vl, double *vu, 
    int *il, int *iu, double *abstol, int *m, double *
    w, void *z__, int *ldz, void *work, int *
    lwork, double *rwork, int *iwork, int *ifail, int *
    info);

int zhegs2_(int *itype, char *uplo, int *n, 
    void *a, int *lda, void *b, int *ldb, 
    int *info);

int zhegst_(int *itype, char *uplo, int *n, 
    void *a, int *lda, void *b, int *ldb, 
    int *info);

int zhegv_(int *itype, char *jobz, char *uplo, int *
    n, void *a, int *lda, void *b, int *ldb, 
    double *w, void *work, int *lwork, double *rwork,
    int *info);

int zhegvd_(int *itype, char *jobz, char *uplo, int *
    n, void *a, int *lda, void *b, int *ldb, 
    double *w, void *work, int *lwork, double *rwork,
    int *lrwork, int *iwork, int *liwork, int *info);

int zhegvx_(int *itype, char *jobz, char *range, char *
    uplo, int *n, void *a, int *lda, void *b, 
    int *ldb, double *vl, double *vu, int *il, int *
    iu, double *abstol, int *m, double *w, void *z__,
    int *ldz, void *work, int *lwork, double *rwork,
    int *iwork, int *ifail, int *info);

int zherfs_(char *uplo, int *n, int *nrhs, 
    void *a, int *lda, void *af, int *ldaf, 
    int *ipiv, void *b, int *ldb, void *x, 
    int *ldx, double *ferr, double *berr, void *work,
    double *rwork, int *info);

int zhesv_(char *uplo, int *n, int *nrhs, 
    void *a, int *lda, int *ipiv, void *b, 
    int *ldb, void *work, int *lwork, int *info);

int zhesvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *a, int *lda, void *af, int *
    ldaf, int *ipiv, void *b, int *ldb, void *x,
    int *ldx, double *rcond, double *ferr, double *berr, 
    void *work, int *lwork, double *rwork, int *info);

int zhetf2_(char *uplo, int *n, void *a, 
    int *lda, int *ipiv, int *info);



int zhetri_(char *uplo, int *n, void *a, 
    int *lda, int *ipiv, void *work, int *info);

int zhetrs_(char *uplo, int *n, int *nrhs, 
    void *a, int *lda, int *ipiv, void *b, 
    int *ldb, int *info);

int zhgeqz_(char *job, char *compq, char *compz, int *n, 
    int *ilo, int *ihi, void *a, int *lda, 
    void *b, int *ldb, void *alpha, void *
    beta, void *q, int *ldq, void *z__, int *
    ldz, void *work, int *lwork, double *rwork, int *
    info);

int zhpcon_(char *uplo, int *n, void *ap, 
    int *ipiv, double *anorm, double *rcond, void *
    work, int *info);

int zhpev_(char *jobz, char *uplo, int *n, void 
    *ap, double *w, void *z__, int *ldz, void *
    work, double *rwork, int *info);

int zhpevd_(char *jobz, char *uplo, int *n, 
    void *ap, double *w, void *z__, int *ldz, 
    void *work, int *lwork, double *rwork, int *
    lrwork, int *iwork, int *liwork, int *info);

int zhpevx_(char *jobz, char *range, char *uplo, int *n, 
    void *ap, double *vl, double *vu, int *il, 
    int *iu, double *abstol, int *m, double *w, 
    void *z__, int *ldz, void *work, double *
    rwork, int *iwork, int *ifail, int *info);

int zhpgst_(int *itype, char *uplo, int *n, 
    void *ap, void *bp, int *info);

int zhpgv_(int *itype, char *jobz, char *uplo, int *
    n, void *ap, void *bp, double *w, void 
    *z__, int *ldz, void *work, double *rwork, int *
    info);

int zhpgvd_(int *itype, char *jobz, char *uplo, int *
    n, void *ap, void *bp, double *w, void 
    *z__, int *ldz, void *work, int *lwork, double *
    rwork, int *lrwork, int *iwork, int *liwork, int *
    info);

int zhpgvx_(int *itype, char *jobz, char *range, char *
    uplo, int *n, void *ap, void *bp, double *
    vl, double *vu, int *il, int *iu, double *abstol, 
    int *m, double *w, void *z__, int *ldz, 
    void *work, double *rwork, int *iwork, int *
    ifail, int *info);

int zhprfs_(char *uplo, int *n, int *nrhs, 
    void *ap, void *afp, int *ipiv, void *
    b, int *ldb, void *x, int *ldx, double *ferr, 
    double *berr, void *work, double *rwork, int *
    info);

int zhpsv_(char *uplo, int *n, int *nrhs, 
    void *ap, int *ipiv, void *b, int *ldb, 
    int *info);

int zhpsvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *ap, void *afp, int *ipiv, 
    void *b, int *ldb, void *x, int *ldx, 
    double *rcond, double *ferr, double *berr, void *
    work, double *rwork, int *info);

int zhptrd_(char *uplo, int *n, void *ap, 
    double *d__, double *e, void *tau, int *info);

int zhptrf_(char *uplo, int *n, void *ap, 
    int *ipiv, int *info);

int zhptri_(char *uplo, int *n, void *ap, 
    int *ipiv, void *work, int *info);

int zhptrs_(char *uplo, int *n, int *nrhs, 
    void *ap, int *ipiv, void *b, int *ldb, 
    int *info);

int zhsein_(char *side, char *eigsrc, char *initv, int *
    select, int *n, void *h__, int *ldh, void *
    w, void *vl, int *ldvl, void *vr, int *ldvr,
    int *mm, int *m, void *work, double *rwork, 
    int *ifaill, int *ifailr, int *info);

int zhseqr_(char *job, char *compz, int *n, int *ilo,
    int *ihi, void *h__, int *ldh, void *w, 
    void *z__, int *ldz, void *work, int *lwork,
    int *info);

int zlabrd_(int *m, int *n, int *nb, 
    void *a, int *lda, double *d__, double *e, 
    void *tauq, void *taup, void *x, int *
    ldx, void *y, int *ldy);


int zlacon_(int *n, void *v, void *x, 
    double *est, int *kase);

int zlacp2_(char *uplo, int *m, int *n, double *
    a, int *lda, void *b, int *ldb);


int zlacrm_(int *m, int *n, void *a, 
    int *lda, double *b, int *ldb, void *c__, 
    int *ldc, double *rwork);

int zlacrt_(int *n, void *cx, int *incx, 
    void *cy, int *incy, void *c__, void *
    s);

int zlaed0_(int *qsiz, int *n, double *d__, 
    double *e, void *q, int *ldq, void *qstore, 
    int *ldqs, double *rwork, int *iwork, int *info);

int zlaed7_(int *n, int *cutpnt, int *qsiz, 
    int *tlvls, int *curlvl, int *curpbm, double *d__, 
    void *q, int *ldq, double *rho, int *indxq, 
    double *qstore, int *qptr, int *prmptr, int *perm, 
    int *givptr, int *givcol, double *givnum, void *
    work, double *rwork, int *iwork, int *info);

int zlaed8_(int *k, int *n, int *qsiz, 
    void *q, int *ldq, double *d__, double *rho, 
    int *cutpnt, double *z__, double *dlamda, void *
    q2, int *ldq2, double *w, int *indxp, int *indx, 
    int *indxq, int *perm, int *givptr, int *givcol, 
    double *givnum, int *info);

int zlaein_(int *rightv, int *noinit, int *n, 
    void *h__, int *ldh, void *w, void *v, 
    void *b, int *ldb, double *rwork, double *eps3, 
    double *smlnum, int *info);

int zlaesy_(void *a, void *b, 
    void *c__, void *rt1, void *rt2, 
    void *evscal, void *cs1, void *sn1);

int zlaev2_(void *a, void *b, 
    void *c__, double *rt1, double *rt2, double *cs1,
    void *sn1);

int zlags2_(int *upper, double *a1, void *
    a2, double *a3, double *b1, void *b2, double *b3,
    double *csu, void *snu, double *csv, void *
    snv, double *csq, void *snq);


int zlahef_(char *uplo, int *n, int *nb, int *kb,
    void *a, int *lda, int *ipiv, void *w, 
    int *ldw, int *info);

int zlahqr_(int *wantt, int *wantz, int *n, 
    int *ilo, int *ihi, void *h__, int *ldh, 
    void *w, int *iloz, int *ihiz, void *z__, 
    int *ldz, int *info);

int zlahrd_(int *n, int *k, int *nb, 
    void *a, int *lda, void *tau, void *t, 
    int *ldt, void *y, int *ldy);

int zlaic1_(int *job, int *j, void *x, 
    double *sest, void *w, void *gamma, double *
    sestpr, void *s, void *c__);

int zlals0_(int *icompq, int *nl, int *nr, 
    int *sqre, int *nrhs, void *b, int *ldb, 
    void *bx, int *ldbx, int *perm, int *givptr, 
    int *givcol, int *ldgcol, double *givnum, int *ldgnum,
    double *poles, double *difl, double *difr, double *
    z__, int *k, double *c__, double *s, double *rwork, 
    int *info);

int zlalsa_(int *icompq, int *smlsiz, int *n, 
    int *nrhs, void *b, int *ldb, void *bx, 
    int *ldbx, double *u, int *ldu, double *vt, int *
    k, double *difl, double *difr, double *z__, double *
    poles, int *givptr, int *givcol, int *ldgcol, int *
    perm, double *givnum, double *c__, double *s, double *
    rwork, int *iwork, int *info);

int zlapll_(int *n, void *x, int *incx, 
    void *y, int *incy, double *ssmin);

int zlapmt_(int *forwrd, int *m, int *n, 
    void *x, int *ldx, int *k);

int zlaqgb_(int *m, int *n, int *kl, int *ku,
    void *ab, int *ldab, double *r__, double *c__, 
    double *rowcnd, double *colcnd, double *amax, char *equed);

int zlaqge_(int *m, int *n, void *a, 
    int *lda, double *r__, double *c__, double *rowcnd, 
    double *colcnd, double *amax, char *equed);

int zlaqhb_(char *uplo, int *n, int *kd, 
    void *ab, int *ldab, double *s, double *scond, 
    double *amax, char *equed);

int zlaqhe_(char *uplo, int *n, void *a, 
    int *lda, double *s, double *scond, double *amax, 
    char *equed);

int zlaqhp_(char *uplo, int *n, void *ap, 
    double *s, double *scond, double *amax, char *equed);

int zlaqp2_(int *m, int *n, int *offset, 
    void *a, int *lda, int *jpvt, void *tau, 
    double *vn1, double *vn2, void *work);

int zlaqps_(int *m, int *n, int *offset, int 
    *nb, int *kb, void *a, int *lda, int *jpvt, 
    void *tau, double *vn1, double *vn2, void *
    auxv, void *f, int *ldf);

int zlaqsb_(char *uplo, int *n, int *kd, 
    void *ab, int *ldab, double *s, double *scond, 
    double *amax, char *equed);

int zlaqsp_(char *uplo, int *n, void *ap, 
    double *s, double *scond, double *amax, char *equed);

int zlaqsy_(char *uplo, int *n, void *a, 
    int *lda, double *s, double *scond, double *amax, 
    char *equed);

int zlar1v_(int *n, int *b1, int *bn, double 
    *sigma, double *d__, double *l, double *ld, double *
    lld, double *gersch, void *z__, double *ztz, 
    double *mingma, int *r__, int *isuppz, double *work);

int zlar2v_(int *n, void *x, void *y, 
    void *z__, int *incx, double *c__, void *s, 
    int *incc);

int zlarcm_(int *m, int *n, double *a, int *
    lda, void *b, int *ldb, void *c__, int *ldc,
    double *rwork);

int zlarf_(char *side, int *m, int *n, void 
    *v, int *incv, void *tau, void *c__, int *
    ldc, void *work);

int zlarfb_(char *side, char *trans, char *direct, char *
    storev, int *m, int *n, int *k, void *v, int 
    *ldv, void *t, int *ldt, void *c__, int *
    ldc, void *work, int *ldwork);

int zlarfg_(int *n, void *alpha, void *
    x, int *incx, void *tau);

int zlarft_(char *direct, char *storev, int *n, int *
    k, void *v, int *ldv, void *tau, void *
    t, int *ldt);

int zlarfx_(char *side, int *m, int *n, 
    void *v, void *tau, void *c__, int *
    ldc, void *work);

int zlargv_(int *n, void *x, int *incx, 
    void *y, int *incy, double *c__, int *incc);

int zlarnv_(int *idist, int *iseed, int *n, 
    void *x);

int zlarrv_(int *n, double *d__, double *l, 
    int *isplit, int *m, double *w, int *iblock, 
    double *gersch, double *tol, void *z__, int *ldz,
    int *isuppz, double *work, int *iwork, int *info);

int zlartg_(void *f, void *g, double *
    cs, void *sn, void *r__);

int zlartv_(int *n, void *x, int *incx, 
    void *y, int *incy, double *c__, void *s, 
    int *incc);

int zlarz_(char *side, int *m, int *n, int *l, 
    void *v, int *incv, void *tau, void *
    c__, int *ldc, void *work);

int zlarzb_(char *side, char *trans, char *direct, char *
    storev, int *m, int *n, int *k, int *l, void 
    *v, int *ldv, void *t, int *ldt, void *c__, 
    int *ldc, void *work, int *ldwork);

int zlarzt_(char *direct, char *storev, int *n, int *
    k, void *v, int *ldv, void *tau, void *
    t, int *ldt);

int zlascl_(char *type__, int *kl, int *ku, 
    double *cfrom, double *cto, int *m, int *n, 
    void *a, int *lda, int *info);

int zlaset_(char *uplo, int *m, int *n, 
    void *alpha, void *beta, void *a, int *
    lda);

int zlasr_(char *side, char *pivot, char *direct, int *m,
    int *n, double *c__, double *s, void *a, 
    int *lda);

int zlassq_(int *n, void *x, int *incx, 
    double *scale, double *sumsq);

int zlaswp_(int *n, void *a, int *lda, 
    int *k1, int *k2, int *ipiv, int *incx);

int zlasyf_(char *uplo, int *n, int *nb, int *kb,
    void *a, int *lda, int *ipiv, void *w, 
    int *ldw, int *info);

int zlatbs_(char *uplo, char *trans, char *diag, char *
    normin, int *n, int *kd, void *ab, int *ldab, 
    void *x, double *scale, double *cnorm, int *info);

int zlatdf_(int *ijob, int *n, void *z__, 
    int *ldz, void *rhs, double *rdsum, double *
    rdscal, int *ipiv, int *jpiv);

int zlatps_(char *uplo, char *trans, char *diag, char *
    normin, int *n, void *ap, void *x, double *
    scale, double *cnorm, int *info);

int zlatrd_(char *uplo, int *n, int *nb, 
    void *a, int *lda, double *e, void *tau, 
    void *w, int *ldw);

int zlatrs_(char *uplo, char *trans, char *diag, char *
    normin, int *n, void *a, int *lda, void *x, 
    double *scale, double *cnorm, int *info);

int zlatrz_(int *m, int *n, int *l, 
    void *a, int *lda, void *tau, void *
    work);

int zlatzm_(char *side, int *m, int *n, 
    void *v, int *incv, void *tau, void *
    c1, void *c2, int *ldc, void *work);

int zlauu2_(char *uplo, int *n, void *a, 
    int *lda, int *info);


int zpbcon_(char *uplo, int *n, int *kd, 
    void *ab, int *ldab, double *anorm, double *
    rcond, void *work, double *rwork, int *info);

int zpbequ_(char *uplo, int *n, int *kd, 
    void *ab, int *ldab, double *s, double *scond, 
    double *amax, int *info);

int zpbrfs_(char *uplo, int *n, int *kd, int *
    nrhs, void *ab, int *ldab, void *afb, int *
    ldafb, void *b, int *ldb, void *x, int *ldx,
    double *ferr, double *berr, void *work, double *
    rwork, int *info);

int zpbstf_(char *uplo, int *n, int *kd, 
    void *ab, int *ldab, int *info);

int zpbsv_(char *uplo, int *n, int *kd, int *
    nrhs, void *ab, int *ldab, void *b, int *
    ldb, int *info);

int zpbsvx_(char *fact, char *uplo, int *n, int *kd, 
    int *nrhs, void *ab, int *ldab, void *afb, 
    int *ldafb, char *equed, double *s, void *b, int 
    *ldb, void *x, int *ldx, double *rcond, double *
    ferr, double *berr, void *work, double *rwork, 
    int *info);

int zpbtf2_(char *uplo, int *n, int *kd, 
    void *ab, int *ldab, int *info);


int zpbtrs_(char *uplo, int *n, int *kd, int *
    nrhs, void *ab, int *ldab, void *b, int *
    ldb, int *info);

int zpocon_(char *uplo, int *n, void *a, 
    int *lda, double *anorm, double *rcond, void *
    work, double *rwork, int *info);

int zpoequ_(int *n, void *a, int *lda, 
    double *s, double *scond, double *amax, int *info);

int zporfs_(char *uplo, int *n, int *nrhs, 
    void *a, int *lda, void *af, int *ldaf, 
    void *b, int *ldb, void *x, int *ldx, 
    double *ferr, double *berr, void *work, double *
    rwork, int *info);

int zposv_(char *uplo, int *n, int *nrhs, 
    void *a, int *lda, void *b, int *ldb, 
    int *info);

int zposvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *a, int *lda, void *af, int *
    ldaf, char *equed, double *s, void *b, int *ldb, 
    void *x, int *ldx, double *rcond, double *ferr, 
    double *berr, void *work, double *rwork, int *
    info);

int zpotf2_(char *uplo, int *n, void *a, 
    int *lda, int *info);



int zpotrs_(char *uplo, int *n, int *nrhs, 
    void *a, int *lda, void *b, int *ldb, 
    int *info);

int zppcon_(char *uplo, int *n, void *ap, 
    double *anorm, double *rcond, void *work, double 
    *rwork, int *info);

int zppequ_(char *uplo, int *n, void *ap, 
    double *s, double *scond, double *amax, int *info);

int zpprfs_(char *uplo, int *n, int *nrhs, 
    void *ap, void *afp, void *b, int *ldb,
    void *x, int *ldx, double *ferr, double *berr, 
    void *work, double *rwork, int *info);

int zppsv_(char *uplo, int *n, int *nrhs, 
    void *ap, void *b, int *ldb, int *info);

int zppsvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *ap, void *afp, char *equed, double *
    s, void *b, int *ldb, void *x, int *ldx, 
    double *rcond, double *ferr, double *berr, void *
    work, double *rwork, int *info);

int zpptrf_(char *uplo, int *n, void *ap, 
    int *info);

int zpptri_(char *uplo, int *n, void *ap, 
    int *info);

int zpptrs_(char *uplo, int *n, int *nrhs, 
    void *ap, void *b, int *ldb, int *info);

int zptcon_(int *n, double *d__, void *e, 
    double *anorm, double *rcond, double *rwork, int *
    info);

int zptrfs_(char *uplo, int *n, int *nrhs, 
    double *d__, void *e, double *df, void *ef, 
    void *b, int *ldb, void *x, int *ldx, 
    double *ferr, double *berr, void *work, double *
    rwork, int *info);

int zptsv_(int *n, int *nrhs, double *d__, 
    void *e, void *b, int *ldb, int *info);

int zptsvx_(char *fact, int *n, int *nrhs, 
    double *d__, void *e, double *df, void *ef, 
    void *b, int *ldb, void *x, int *ldx, 
    double *rcond, double *ferr, double *berr, void *
    work, double *rwork, int *info);


int zpttrs_(char *uplo, int *n, int *nrhs, 
    double *d__, void *e, void *b, int *ldb, 
    int *info);

int zptts2_(int *iuplo, int *n, int *nrhs, 
    double *d__, void *e, void *b, int *ldb);


int zspcon_(char *uplo, int *n, void *ap, 
    int *ipiv, double *anorm, double *rcond, void *
    work, int *info);

int zspmv_(char *uplo, int *n, void *alpha, 
    void *ap, void *x, int *incx, void *
    beta, void *y, int *incy);

int zspr_(char *uplo, int *n, void *alpha, 
    void *x, int *incx, void *ap);

int zsprfs_(char *uplo, int *n, int *nrhs, 
    void *ap, void *afp, int *ipiv, void *
    b, int *ldb, void *x, int *ldx, double *ferr, 
    double *berr, void *work, double *rwork, int *
    info);

int zspsv_(char *uplo, int *n, int *nrhs, 
    void *ap, int *ipiv, void *b, int *ldb, 
    int *info);

int zspsvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *ap, void *afp, int *ipiv, 
    void *b, int *ldb, void *x, int *ldx, 
    double *rcond, double *ferr, double *berr, void *
    work, double *rwork, int *info);

int zsptrf_(char *uplo, int *n, void *ap, 
    int *ipiv, int *info);

int zsptri_(char *uplo, int *n, void *ap, 
    int *ipiv, void *work, int *info);

int zsptrs_(char *uplo, int *n, int *nrhs, 
    void *ap, int *ipiv, void *b, int *ldb, 
    int *info);


int zstein_(int *n, double *d__, double *e, 
    int *m, double *w, int *iblock, int *isplit, 
    void *z__, int *ldz, double *work, int *iwork, 
    int *ifail, int *info);

int zsteqr_(char *compz, int *n, double *d__, 
    double *e, void *z__, int *ldz, double *work, 
    int *info);

int zsycon_(char *uplo, int *n, void *a, 
    int *lda, int *ipiv, double *anorm, double *rcond, 
    void *work, int *info);



int zsyrfs_(char *uplo, int *n, int *nrhs, 
    void *a, int *lda, void *af, int *ldaf, 
    int *ipiv, void *b, int *ldb, void *x, 
    int *ldx, double *ferr, double *berr, void *work,
    double *rwork, int *info);

int zsysv_(char *uplo, int *n, int *nrhs, 
    void *a, int *lda, int *ipiv, void *b, 
    int *ldb, void *work, int *lwork, int *info);

int zsysvx_(char *fact, char *uplo, int *n, int *
    nrhs, void *a, int *lda, void *af, int *
    ldaf, int *ipiv, void *b, int *ldb, void *x,
    int *ldx, double *rcond, double *ferr, double *berr, 
    void *work, int *lwork, double *rwork, int *info);

int zsytf2_(char *uplo, int *n, void *a, 
    int *lda, int *ipiv, int *info);


int zsytri_(char *uplo, int *n, void *a, 
    int *lda, int *ipiv, void *work, int *info);

int zsytrs_(char *uplo, int *n, int *nrhs, 
    void *a, int *lda, int *ipiv, void *b, 
    int *ldb, int *info);

int ztbcon_(char *norm, char *uplo, char *diag, int *n, 
    int *kd, void *ab, int *ldab, double *rcond, 
    void *work, double *rwork, int *info);

int ztbrfs_(char *uplo, char *trans, char *diag, int *n, 
    int *kd, int *nrhs, void *ab, int *ldab, 
    void *b, int *ldb, void *x, int *ldx, 
    double *ferr, double *berr, void *work, double *
    rwork, int *info);


int ztgevc_(char *side, char *howmny, int *select, 
    int *n, void *a, int *lda, void *b, int 
    *ldb, void *vl, int *ldvl, void *vr, int *
    ldvr, int *mm, int *m, void *work, double *rwork,
    int *info);

int ztgex2_(int *wantq, int *wantz, int *n, 
    void *a, int *lda, void *b, int *ldb, 
    void *q, int *ldq, void *z__, int *ldz, 
    int *j1, int *info);

int ztgexc_(int *wantq, int *wantz, int *n, 
    void *a, int *lda, void *b, int *ldb, 
    void *q, int *ldq, void *z__, int *ldz, 
    int *ifst, int *ilst, int *info);

int ztgsen_(int *ijob, int *wantq, int *wantz, 
    int *select, int *n, void *a, int *lda, 
    void *b, int *ldb, void *alpha, void *
    beta, void *q, int *ldq, void *z__, int *
    ldz, int *m, double *pl, double *pr, double *dif, 
    void *work, int *lwork, int *iwork, int *liwork, 
    int *info);

int ztgsja_(char *jobu, char *jobv, char *jobq, int *m, 
    int *p, int *n, int *k, int *l, void *a, 
    int *lda, void *b, int *ldb, double *tola, 
    double *tolb, double *alpha, double *beta, void *
    u, int *ldu, void *v, int *ldv, void *q, 
    int *ldq, void *work, int *ncycle, int *info);

int ztgsna_(char *job, char *howmny, int *select, 
    int *n, void *a, int *lda, void *b, int 
    *ldb, void *vl, int *ldvl, void *vr, int *
    ldvr, double *s, double *dif, int *mm, int *m, 
    void *work, int *lwork, int *iwork, int *info);

int ztgsy2_(char *trans, int *ijob, int *m, int *
    n, void *a, int *lda, void *b, int *ldb, 
    void *c__, int *ldc, void *d__, int *ldd, 
    void *e, int *lde, void *f, int *ldf, 
    double *scale, double *rdsum, double *rdscal, int *
    info);

int ztgsyl_(char *trans, int *ijob, int *m, int *
    n, void *a, int *lda, void *b, int *ldb, 
    void *c__, int *ldc, void *d__, int *ldd, 
    void *e, int *lde, void *f, int *ldf, 
    double *scale, double *dif, void *work, int *
    lwork, int *iwork, int *info);

int ztpcon_(char *norm, char *uplo, char *diag, int *n, 
    void *ap, double *rcond, void *work, double 
    *rwork, int *info);

int ztprfs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, void *ap, void *b, int *ldb, 
    void *x, int *ldx, double *ferr, double *berr, 
    void *work, double *rwork, int *info);

int ztptri_(char *uplo, char *diag, int *n, 
    void *ap, int *info);

int ztptrs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, void *ap, void *b, int *ldb, 
    int *info);

int ztrcon_(char *norm, char *uplo, char *diag, int *n, 
    void *a, int *lda, double *rcond, void *
    work, double *rwork, int *info);

int ztrevc_(char *side, char *howmny, int *select, 
    int *n, void *t, int *ldt, void *vl, 
    int *ldvl, void *vr, int *ldvr, int *mm, int 
    *m, void *work, double *rwork, int *info);

int ztrexc_(char *compq, int *n, void *t, 
    int *ldt, void *q, int *ldq, int *ifst, int *
    ilst, int *info);

int ztrrfs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, void *a, int *lda, void *b, 
    int *ldb, void *x, int *ldx, double *ferr, 
    double *berr, void *work, double *rwork, int *
    info);

int ztrsen_(char *job, char *compq, int *select, int 
    *n, void *t, int *ldt, void *q, int *ldq, 
    void *w, int *m, double *s, double *sep, 
    void *work, int *lwork, int *info);

int ztrsna_(char *job, char *howmny, int *select, 
    int *n, void *t, int *ldt, void *vl, 
    int *ldvl, void *vr, int *ldvr, double *s, 
    double *sep, int *mm, int *m, void *work, 
    int *ldwork, double *rwork, int *info);

int ztrsyl_(char *trana, char *tranb, int *isgn, int 
    *m, int *n, void *a, int *lda, void *b, 
    int *ldb, void *c__, int *ldc, double *scale, 
    int *info);

int ztrti2_(char *uplo, char *diag, int *n, 
    void *a, int *lda, int *info);


int ztrtrs_(char *uplo, char *trans, char *diag, int *n, 
    int *nrhs, void *a, int *lda, void *b, 
    int *ldb, int *info);

int ztzrqf_(int *m, int *n, void *a, 
    int *lda, void *tau, int *info);

int ztzrzf_(int *m, int *n, void *a, 
    int *lda, void *tau, void *work, int *lwork,
    int *info);

int zung2l_(int *m, int *n, int *k, 
    void *a, int *lda, void *tau, void *
    work, int *info);

int zung2r_(int *m, int *n, int *k, 
    void *a, int *lda, void *tau, void *
    work, int *info);

int zungbr_(char *vect, int *m, int *n, int *k, 
    void *a, int *lda, void *tau, void *
    work, int *lwork, int *info);

int zunghr_(int *n, int *ilo, int *ihi, 
    void *a, int *lda, void *tau, void *
    work, int *lwork, int *info);

int zungl2_(int *m, int *n, int *k, 
    void *a, int *lda, void *tau, void *
    work, int *info);


int zungql_(int *m, int *n, int *k, 
    void *a, int *lda, void *tau, void *
    work, int *lwork, int *info);


int zungr2_(int *m, int *n, int *k, 
    void *a, int *lda, void *tau, void *
    work, int *info);

int zungrq_(int *m, int *n, int *k, 
    void *a, int *lda, void *tau, void *
    work, int *lwork, int *info);

int zungtr_(char *uplo, int *n, void *a, 
    int *lda, void *tau, void *work, int *lwork,
    int *info);

int zunm2l_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *info);

int zunm2r_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *info);

int zunmbr_(char *vect, char *side, char *trans, int *m, 
    int *n, int *k, void *a, int *lda, void 
    *tau, void *c__, int *ldc, void *work, int *
    lwork, int *info);

int zunmhr_(char *side, char *trans, int *m, int *n, 
    int *ilo, int *ihi, void *a, int *lda, 
    void *tau, void *c__, int *ldc, void *
    work, int *lwork, int *info);

int zunml2_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *info);


int zunmql_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *lwork,
    int *info);


int zunmr2_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *info);

int zunmr3_(char *side, char *trans, int *m, int *n, 
    int *k, int *l, void *a, int *lda, void 
    *tau, void *c__, int *ldc, void *work, int *
    info);

int zunmrq_(char *side, char *trans, int *m, int *n, 
    int *k, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *lwork,
    int *info);

int zunmrz_(char *side, char *trans, int *m, int *n, 
    int *k, int *l, void *a, int *lda, void 
    *tau, void *c__, int *ldc, void *work, int *
    lwork, int *info);

int zunmtr_(char *side, char *uplo, char *trans, int *m, 
    int *n, void *a, int *lda, void *tau, 
    void *c__, int *ldc, void *work, int *lwork,
    int *info);

int zupgtr_(char *uplo, int *n, void *ap, 
    void *tau, void *q, int *ldq, void *
    work, int *info);

int zupmtr_(char *side, char *uplo, char *trans, int *m, 
    int *n, void *ap, void *tau, void *c__,
    int *ldc, void *work, int *info);

#endif /* __CLAPACK_H */
