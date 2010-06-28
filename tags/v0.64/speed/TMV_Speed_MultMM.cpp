
#include "TMV.h"

#ifdef NOBLAS
#undef NOBLAS
#endif

#define CBLAS
//#define FBLAS

#include "../src/TMV_Blas.h"

#include <iostream>
#include <sys/time.h>
#include <fstream>
#include <assert.h>


const int NLOOPS = 3;

static void ClearCache()
{
    tmv::Matrix<double> M(2000,2000,8.);
    for(int i=0;i<2000;i++) for(int j=0;j<2000;j++) {
        M(i,j) = 8.*i-6.*j-748.;
    }
}

template <class T>
void CallGemm(
    tmv::StorageType sa, tmv::StorageType sb, 
    int xm, int xn, int xk, T xa, 
    const T* a, int xlda, const T* b, int xldb, 
    T xb, T* c, int xldc)
{ assert(false); }

template <>
void CallGemm<double>(
    tmv::StorageType sa, tmv::StorageType sb, 
    int xm, int xn, int xk, double xa, 
    const double* a, int xlda, const double* b, int xldb, 
    double xb, double* c, int xldc)
{
    BLASNAME(dgemm) (
        BLASCM 
        sa==tmv::RowMajor ? BLASCH_T : BLASCH_NT,
        sb==tmv::RowMajor ? BLASCH_T : BLASCH_NT,
        BLASV(xm),BLASV(xn),BLASV(xk),BLASV(xa),
        BLASP(a),BLASV(xlda),BLASP(b),BLASV(xldb),
        BLASV(xb),BLASP(c),BLASV(xldc) BLAS1 BLAS1);
}

template <>
void CallGemm<float>(
    tmv::StorageType sa, tmv::StorageType sb, 
    int xm, int xn, int xk, float xa, 
    const float* a, int xlda, const float* b, int xldb, 
    float xb, float* c, int xldc)
{
    BLASNAME(sgemm) (
        BLASCM 
        sa==tmv::RowMajor ? BLASCH_T : BLASCH_NT,
        sb==tmv::RowMajor ? BLASCH_T : BLASCH_NT,
        BLASV(xm),BLASV(xn),BLASV(xk),BLASV(xa),
        BLASP(a),BLASV(xlda),BLASP(b),BLASV(xldb),
        BLASV(xb),BLASP(c),BLASV(xldc) BLAS1 BLAS1);
}

template <>
void CallGemm<std::complex<double> >(
    tmv::StorageType sa, tmv::StorageType sb, 
    int xm, int xn, int xk, std::complex<double> xa, 
    const std::complex<double>* a, int xlda, 
    const std::complex<double>* b, int xldb, 
    std::complex<double> xb, std::complex<double>* c, int xldc)
{
    BLASNAME(zgemm) (
        BLASCM 
        sa==tmv::RowMajor ? BLASCH_T : BLASCH_NT,
        sb==tmv::RowMajor ? BLASCH_T : BLASCH_NT,
        BLASV(xm),BLASV(xn),BLASV(xk),BLASP(&xa),
        BLASP(a),BLASV(xlda),BLASP(b),BLASV(xldb),
        BLASP(&xb),BLASP(c),BLASV(xldc) BLAS1 BLAS1);
}

template <>
void CallGemm<std::complex<float> >(
    tmv::StorageType sa, tmv::StorageType sb, 
    int xm, int xn, int xk, std::complex<float> xa, 
    const std::complex<float>* a, int xlda,
    const std::complex<float>* b, int xldb, 
    std::complex<float> xb, std::complex<float>* c, int xldc)
{
    BLASNAME(cgemm) (
        BLASCM 
        sa==tmv::RowMajor ? BLASCH_T : BLASCH_NT,
        sb==tmv::RowMajor ? BLASCH_T : BLASCH_NT,
        BLASV(xm),BLASV(xn),BLASV(xk),BLASP(&xa),
        BLASP(a),BLASV(xlda),BLASP(b),BLASV(xldb),
        BLASP(&xb),BLASP(c),BLASV(xldc) BLAS1 BLAS1);
}

template <class T, tmv::StorageType S> 
static void Complexify(tmv::Matrix<T,S>& M)
{}


template <class T, tmv::StorageType S> 
static void Complexify(tmv::Matrix<std::complex<T>,S>& M)
{ 
    for(size_t i=0;i<M.colsize();i++) for(size_t j=0;j<M.rowsize();j++) {
        M(i,j) *= std::complex<T>(0.3,0.1);
    }
}


template <class T, tmv::StorageType SA, tmv::StorageType SB>
static void Speed_MultMM(int Mmult, int Nmult, int Kmult, const char* file)
{
    std::cout<<tmv::TMV_Text(T())<<"  "<<
        TMV_Text(SA)<<"  "<<TMV_Text(SB)<<std::endl;
    std::cout<<file<<"  "<<Mmult<<" "<<Nmult<<" "<<Kmult<<std::endl;
    std::ofstream os(file);

    os<<"# N  TMV  BLAS\n";

    for(int I=2;I<=2048;I*=2) {
        std::cout<<I<<std::endl;

        const int M = Mmult*I+3;
        const int N = Nmult*I+3;
        const int K = Kmult*I+3;

        tmv::Matrix<T,SA> A1(M,K);
        tmv::Matrix<T,SB> B1(K,N);
        tmv::Matrix<T,tmv::ColMajor> C1(M,N);
        tmv::Matrix<T,SA> A2(M,K);
        tmv::Matrix<T,SB> B2(K,N);
        tmv::Matrix<T,tmv::ColMajor> C2(M,N);

        timeval tp;

        for(int i=0;i<M;i++) for(int k=0;k<K;k++) {
            A1(i,k) = (i+3.*k-43.)/(i+k+11.);
        }
        for(int k=0;k<K;k++) for(int j=0;j<N;j++) {
            B1(k,j) = (-2.*k+2.*k*j-133.)/(2.*j+N*k+121.);
        }
        Complexify(A1);
        Complexify(B1);
        A2 = A1;
        B2 = B1;

        double tmvtime=1.e100;
        double blastime=1.e100;

        for(int i=0;i<NLOOPS;i++) {
            ClearCache();

            gettimeofday(&tp,0);
            double t1 = tp.tv_sec + tp.tv_usec/1.e6;

            C1 = A1*B1;

            gettimeofday(&tp,0);
            double t2 = tp.tv_sec + tp.tv_usec/1.e6;

            double time = t2-t1;
            if (time < tmvtime) tmvtime = time;
            std::cout<<time<<"  ";
        }
        std::cout<<tmvtime<<"  (TMV)\n";

        for(int i=0;i<NLOOPS;i++) {
            ClearCache();

            gettimeofday(&tp,0);
            double t3 = tp.tv_sec + tp.tv_usec/1.e6;

            int xm = M;
            int xn = N;
            int xk = K;
            T xa = 1.;
            T xb = 0.;
            int xlda = SA==tmv::RowMajor ? K : M;
            int xldb = SB==tmv::RowMajor ? N : K;
            int xldc = M;
            CallGemm<T>(SA,SB,xm,xn,xk,xa,A2.cptr(),xlda,B2.cptr(),
                        xldb,xb,C2.ptr(),xldc);

            gettimeofday(&tp,0);
            double t4 = tp.tv_sec + tp.tv_usec/1.e6;

            double time = t4-t3;
            if (time < blastime) blastime = time;
            std::cout<<time<<"  ";
        }
        std::cout<<blastime<<"  (BLAS)\n";

        std::cout<<"Norm(C1-C2) = "<<Norm(C1-C2)<<
            "  Norm(C) = "<<Norm(C1)<<std::endl;
        assert(Norm(C1-C2) < 1.e-6*Norm(C1));

        double ops = double(M)*double(N)*double(K);
        if (tmv::isComplex(T())) ops *= 4.;
        os<<N<<"  "<<ops/tmvtime<<"  "<<ops/blastime<<std::endl;
    }
}

int main() try 
{
    Speed_MultMM<double,tmv::ColMajor,tmv::ColMajor>(
        1,1,1, "speed_multmm_double_sq_ccc.data");
    Speed_MultMM<float,tmv::ColMajor,tmv::ColMajor>(
        1,1,1, "speed_multmm_float_sq_ccc.data");
    Speed_MultMM<std::complex<double>,tmv::ColMajor,tmv::ColMajor>(
        1,1,1, "speed_multmm_complexdouble_sq_ccc.data");
    Speed_MultMM<std::complex<float>,tmv::ColMajor,tmv::ColMajor>(
        1,1,1, "speed_multmm_complexfloat_sq_ccc.data");

    Speed_MultMM<double,tmv::RowMajor,tmv::ColMajor>(
        1,1,1, "speed_multmm_double_sq_rcc.data");
    Speed_MultMM<float,tmv::RowMajor,tmv::ColMajor>(
        1,1,1, "speed_multmm_float_sq_rcc.data");
    Speed_MultMM<std::complex<double>,tmv::RowMajor,tmv::ColMajor>(
        1,1,1, "speed_multmm_complexdouble_sq_rcc.data");
    Speed_MultMM<std::complex<float>,tmv::RowMajor,tmv::ColMajor>(
        1,1,1, "speed_multmm_complexfloat_sq_rcc.data");

    Speed_MultMM<double,tmv::ColMajor,tmv::RowMajor>(
        1,1,1, "speed_multmm_double_sq_crc.data");
    Speed_MultMM<float,tmv::ColMajor,tmv::RowMajor>(
        1,1,1, "speed_multmm_float_sq_crc.data");
    Speed_MultMM<std::complex<double>,tmv::ColMajor,tmv::RowMajor>(
        1,1,1, "speed_multmm_complexdouble_sq_crc.data");
    Speed_MultMM<std::complex<float>,tmv::ColMajor,tmv::RowMajor>(
        1,1,1, "speed_multmm_complexfloat_sq_crc.data");

    Speed_MultMM<double,tmv::RowMajor,tmv::RowMajor>(
        1,1,1, "speed_multmm_double_sq_rrc.data");
    Speed_MultMM<float,tmv::RowMajor,tmv::RowMajor>(
        1,1,1, "speed_multmm_float_sq_rrc.data");
    Speed_MultMM<std::complex<double>,tmv::RowMajor,tmv::RowMajor>(
        1,1,1, "speed_multmm_complexdouble_sq_rrc.data");
    Speed_MultMM<std::complex<float>,tmv::RowMajor,tmv::RowMajor>(
        1,1,1, "speed_multmm_complexfloat_sq_rrc.data");

    return 0;

} catch (tmv::Error& e) {
    std::cerr<<e<<std::endl;
    exit(1);
}
