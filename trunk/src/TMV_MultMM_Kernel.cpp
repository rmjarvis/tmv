
#define TMV_CompilingLibrary

#include "TMV_Blas.h"
#include "tmv/TMV_MultMM_Kernel.h"

namespace tmv {

#if !defined(BLAS) || !defined(TMV_INST_MIX)
#ifdef TMV_INST_FLOAT
#ifdef __SSE__
    template void multmm_16_16_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_16_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_M_16_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_M_16_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_N_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_N_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);



    template void multmm_16_16_64(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_16_64(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_M_16_64(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_M_16_64(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_N_64(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<1,float>& x,
        const float* A, const float* B, float* C0);
    template void multmm_16_N_64(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<-1,float>& x,
        const float* A, const float* B, float* C0);
#endif
#endif

#ifdef TMV_INST_DOUBLE
#ifdef __SSE2__
    template void multmm_16_16_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<1,double>& x,
        const double* A, const double* B, double* C0);
    template void multmm_16_16_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<-1,double>& x,
        const double* A, const double* B, double* C0);
    template void multmm_M_16_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<1,double>& x,
        const double* A, const double* B, double* C0);
    template void multmm_M_16_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<-1,double>& x,
        const double* A, const double* B, double* C0);
    template void multmm_16_N_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<1,double>& x,
        const double* A, const double* B, double* C0);
    template void multmm_16_N_32(
        const ptrdiff_t M, const ptrdiff_t N, const ptrdiff_t K,
        const Scaling<-1,double>& x,
        const double* A, const double* B, double* C0);

#endif
#endif
#endif

#define InstFile "TMV_MultMM_Kernel.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


