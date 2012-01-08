
//#define PRINTALGO_XV

#include "TMV_Blas.h"
#include "tmv/TMV_MultXV.h"
#include "tmv/TMV_ScaleV.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"

namespace tmv {

    template <bool add, class T, class V1, class V2>
    static void DoMultXV3(const T x, const V1& v1, V2& v2)
    { 
        if (x == T(1))
            InlineMultXV<add>(Scaling<1,T>(x),v1,v2); 
        else if (x == T(-1))
            InlineMultXV<add>(Scaling<-1,T>(x),v1,v2); 
        else if (x == T(0))
            Maybe<!add>::zero(v2);
        else
            InlineMultXV<add>(Scaling<0,T>(x),v1,v2); 
    }

    template <bool add, class T, class V1, class V2>
    static void DoMultXV3(const std::complex<T> x, const V1& v1, V2& v2)
    {
        if (imag(x) == T(0)) {
            if (real(x) == T(1))
                InlineMultXV<add>(Scaling<1,T>(real(x)),v1,v2);
            else if (real(x) == T(-1))
                InlineMultXV<add>(Scaling<-1,T>(real(x)),v1,v2);
            else if (real(x) == T(0))
                Maybe<!add>::zero(v2);
            else
                InlineMultXV<add>(Scaling<0,T>(real(x)),v1,v2);
        } else 
            InlineMultXV<add>(Scaling<0,std::complex<T> >(x),v1,v2); 
    }

    template <bool add, class T, class V1>
    static void DoMultXV2(const T x, const V1& v1, VectorView<T> v2)
    {
        if (v1.step() == 1 && v2.step() == 1) {
            VectorView<T,Unit> v2u = v2.unitView();
            DoMultXV3<add>(x,v1.unitView(),v2u);
        } else 
            DoMultXV3<add>(x,v1,v2); 
    }

    template <class T1, int C1, class T2>
    static void DoMultXV(
        const T2 x, const ConstVectorView<T1,C1>& v1,
        VectorView<T2> v2)
    { DoMultXV2<false>(x,v1,v2); }
    template <class T1, int C1, class T2>
    static void DoAddMultXV(
        const T2 x, const ConstVectorView<T1,C1>& v1,
        VectorView<T2> v2)
    { DoMultXV2<true>(x,v1,v2); }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
    static void DoMultXV(
        const double x,
        const ConstVectorView<double>& v1, VectorView<double> v2)
    { InstScale(x,v2=v1); }
    static void DoMultXV(
        const std::complex<double> x, 
        const ConstVectorView<std::complex<double> >& v1, 
        VectorView<std::complex<double> > v2)
    { InstScale(x,v2=v1); }
#ifdef TMV_INST_MIX
    static void DoMultXV(
        const std::complex<double> x, 
        const ConstVectorView<double>& v1, 
        VectorView<std::complex<double> > v2)
    { InstScale(x,v2=v1); }
#endif
    static void DoAddMultXV(
        const double x,
        const ConstVectorView<double>& v1, VectorView<double> v2)
    { 
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        const double* v1p = v1.cptr();
        if (s1<0) v1p += (n-1)*s1;
        double* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        BLASNAME(daxpy) (
            BLASV(n),BLASV(x),BLASP(v1p),BLASV(s1),
            BLASP(v2p),BLASV(s2));
    }
    static void DoAddMultXV(
        const std::complex<double> x, 
        const ConstVectorView<std::complex<double> >& v1, 
        VectorView<std::complex<double> > v2)
    {
        if (imag(x) == 0.) {
            if (v1.step() == 1 && v2.step() == 1) {
                DoAddMultXV(
                    real(x),v1.flatten().xView(),v2.flatten().xView());
            } else {
                DoAddMultXV(real(x),v1.realPart(),v2.realPart());
                DoAddMultXV(real(x),v1.imagPart(),v2.imagPart());
            }
        } else {
            int n=v2.size();
            int s1=v1.step();
            int s2=v2.step();
            const std::complex<double>* v1p = v1.cptr();
            if (s1<0) v1p += (n-1)*s1;
            std::complex<double>* v2p = v2.ptr();
            if (s2<0) v2p += (n-1)*s2;
            BLASNAME(zaxpy) (
                BLASV(n),BLASP(&x),BLASP(v1p),BLASV(s1),
                BLASP(v2p),BLASV(s2));
        }
    }
#ifdef TMV_INST_MIX
    static void DoAddMultXV(
        const std::complex<double> x, 
        const ConstVectorView<double>& v1, 
        VectorView<std::complex<double> > v2)
    {
        double xr = real(x);
        if (xr != 0.) DoAddMultXV(xr,v1,v2.realPart());
        double xi = imag(x);
        if (xi != 0.) DoAddMultXV(xi,v1,v2.imagPart());
    }
#endif
#endif
#ifdef TMV_INST_FLOAT
    static void DoMultXV(
        const float x,
        const ConstVectorView<float>& v1, VectorView<float> v2)
    { InstScale(x,v2=v1); }
    static void DoMultXV(
        const std::complex<float> x, 
        const ConstVectorView<std::complex<float> >& v1, 
        VectorView<std::complex<float> > v2)
    { InstScale(x,v2=v1); }
#ifdef TMV_INST_MIX
    static void DoMultXV(
        const std::complex<float> x, 
        const ConstVectorView<float>& v1, 
        VectorView<std::complex<float> > v2)
    { InstScale(x,v2=v1); }
#endif
    static void DoAddMultXV(
        const float x,
        const ConstVectorView<float>& v1, VectorView<float> v2)
    {
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        const float* v1p = v1.cptr();
        if (s1<0) v1p += (n-1)*s1;
        float* v2p = v2.ptr();
        if (s2<0) v2p += (n-1)*s2;
        BLASNAME(saxpy) (
            BLASV(n),BLASV(x),BLASP(v1p),BLASV(s1),
            BLASP(v2p),BLASV(s2));
    }
    static void DoAddMultXV(
        const std::complex<float> x, 
        const ConstVectorView<std::complex<float> >& v1, 
        VectorView<std::complex<float> > v2)
    {
        if (imag(x) == 0.F) {
            if (v1.step() == 1 && v2.step() == 1) 
                DoAddMultXV(
                    real(x),v1.flatten().xView(),v2.flatten().xView());
            else {
                DoAddMultXV(real(x),v1.realPart(),v2.realPart());
                DoAddMultXV(real(x),v1.imagPart(),v2.imagPart());
            }
        } else {
            int n=v2.size();
            int s1=v1.step();
            int s2=v2.step();
            const std::complex<float>* v1p = v1.cptr();
            if (s1<0) v1p += (n-1)*s1;
            std::complex<float>* v2p = v2.ptr();
            if (s2<0) v2p += (n-1)*s2;
            BLASNAME(caxpy) (
                BLASV(n),BLASP(&x),BLASP(v1p),BLASV(s1),
                BLASP(v2p),BLASV(s2));
        }
    }
#ifdef TMV_INST_MIX
    static void DoAddMultXV(
        const std::complex<float> x, 
        const ConstVectorView<float>& v1, 
        VectorView<std::complex<float> > v2)
    {
        float xr = real(x);
        if (xr != 0.F) DoAddMultXV(xr,v1,v2.realPart());
        float xi = imag(x);
        if (xi != 0.F) DoAddMultXV(xi,v1,v2.imagPart());
    }
#endif
#endif
#endif // BLAS

    template <class T1, int C1, class T2>
    void InstMultXV(
        const T2 x, const ConstVectorView<T1,C1>& v1, VectorView<T2> v2)
    { DoMultXV(x,v1,v2); }
    template <class T1, int C1, class T2>
    void InstAddMultXV(
        const T2 x, const ConstVectorView<T1,C1>& v1, VectorView<T2> v2)
    { DoAddMultXV(x,v1,v2); }

    template <class T1, int C1, class T2>
    void InstAliasMultXV(
        const T2 x, const ConstVectorView<T1,C1>& v1, VectorView<T2> v2)
    { InlineAliasMultXV<false>(Scaling<0,T2>(x),v1,v2); }
    template <class T1, int C1, class T2>
    void InstAliasAddMultXV(
        const T2 x, const ConstVectorView<T1,C1>& v1, VectorView<T2> v2)
    { InlineAliasMultXV<true>(Scaling<0,T2>(x),v1,v2); }

#define InstFile "TMV_MultXV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


