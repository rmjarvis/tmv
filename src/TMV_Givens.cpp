
#include "TMV_Blas.h"
#include "tmv/TMV_Givens.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"

namespace tmv {

    template <class T, class RT>
    void InstGivensRotate(T& x, T& y, RT& c, T& s)
    { InlineGivensRotate(x,y,c,s); }

    template <class T, class RT, class T2>
    static void DoGivensMult(
        RT c, T s, VectorView<T2>& v1, VectorView<T2>& v2)
    {
        if (v1.step() == 1 && v2.step() == 1) {
            VectorView<T2,Unit> v1u = v1.unitView();
            VectorView<T2,Unit> v2u = v2.unitView();
            InlineGivensMultV(c,s,v1u,v2u);
        } else {
            InlineGivensMultV(c,s,v1,v2);
        }
    }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
    static void DoGivensMult(
        double c, double s,
        VectorView<double>& v1, VectorView<double>& v2)
    { 
        int n=v1.size();
        int s1=v1.step();
        int s2=v2.step();
        double* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        double* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        BLASNAME(drot) (
            BLASV(n),BLASP(v1p),BLASV(s1),
            BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
    }
#ifdef BLASZDROT
    static void DoGivensMult(
        double c, double s,
        VectorView<std::complex<double> >& v1,
        VectorView<std::complex<double> >& v2)
    { 
        int n=v1.size();
        int s1=v1.step();
        int s2=v2.step();
        std::complex<double>* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        std::complex<double>* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        BLASNAME(zdrot) (
            BLASV(n),BLASP(v1p),BLASV(s1),
            BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
    }
#endif
#ifdef ELAP
    static void DoGivensMult(
        double c, std::complex<double> s,
        VectorView<std::complex<double> >& v1,
        VectorView<std::complex<double> >& v2)
    {
        int n = v1.size();
        int s1 = v1.step();
        int s2 = v2.step();
        std::complex<double>* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        std::complex<double>* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        LAPNAME(zrot) (
            LAPV(n),LAPP(v1p),LAPV(s1),LAPP(v2p),LAPV(s2),
            LAPV(c),LAPP(&s)); 
    }
#endif // ELAP
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
    static void DoGivensMult(
        float c, float s,
        VectorView<float>& v1, VectorView<float>& v2)
    { 
        int n=v1.size();
        int s1=v1.step();
        int s2=v2.step();
        float* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        float* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        BLASNAME(srot) (
            BLASV(n),BLASP(v1p),BLASV(s1),
            BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
    }
#ifdef BLASZDROT
    static void DoGivensMult(
        float c, float s,
        VectorView<std::complex<float> >& v1,
        VectorView<std::complex<float> >& v2)
    { 
        int n=v1.size();
        int s1=v1.step();
        int s2=v2.step();
        std::complex<float>* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        std::complex<float>* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        BLASNAME(csrot) (
            BLASV(n),BLASP(v1p),BLASV(s1),
            BLASP(v2p),BLASV(s2),BLASV(c),BLASV(s)); 
    }
#endif
#ifdef ELAP
    static void DoGivensMult(
        float c, std::complex<float> s,
        VectorView<std::complex<float> >& v1,
        VectorView<std::complex<float> >& v2)
    {
        int n = v1.size();
        int s1 = v1.step();
        int s2 = v2.step();
        std::complex<float>* v1p = v1.ptr();
        if (s1<0) v1p += (n-1)*s1;
        std::complex<float>* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        LAPNAME(crot) (
            LAPV(n),LAPP(v1p),LAPV(s1),LAPP(v2p),LAPV(s2),
            LAPV(c),LAPP(&s)); 
    }
#endif // ELAP
#endif // DOUBLE
#endif // BLAS

    template <class T, class RT, class T2>
    void InstGivensMultV(RT c, T s, VectorView<T2> v1, VectorView<T2> v2)
    {
        if (v1.size() > 0 && s != T(0)) {
            DoGivensMult(c,s,v1,v2);
        }
    }

#define InstFile "TMV_Givens.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


