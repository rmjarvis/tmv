

#include "TMV_Blas.h"
#include "tmv/TMV_ScaleV.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_CopyV.h"
//#include "tmv/TMV_ProdXV.h"

namespace tmv {

    template <class T, class V> 
    static void DoScale2(const T x, V& v)
    { 
        if (x == T(-1))
            InlineScale(Scaling<-1,T>(x),v); 
        else if (x == T(0))
            v.setZero();
        else if (x != T(1))
            InlineScale(Scaling<0,T>(x),v); 
    }

    template <class T, class V> 
    static void DoScale2(const std::complex<T> x, V& v)
    {
        if (imag(x) == T(0)) {
            if (real(x) == T(-1))
                InlineScale(Scaling<-1,T>(real(x)),v);
            else if (real(x) == T(0))
                v.setZero();
            else if (real(x) != T(1))
                InlineScale(Scaling<0,T>(real(x)),v);
        } else 
            InlineScale(Scaling<0,std::complex<T> >(x),v); 
    }

    template <class T> 
    static void DoScale(const T x, VectorView<T> v)
    { 
        if (v.step() == 1) {
            VectorView<T,Unit> vu = v.unitView();
            DoScale2(x,vu);
        } else 
            DoScale2(x,v); 
    }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
    static void DoScale(const double x, VectorView<double> v)
    {
        int n=v.size();
        if (n==0) return;
        TMVAssert(v.size()>0);
        int s=v.step();
        if (s == 0) { v[0] *= x; return; }
        double* vp = v.ptr();
        if (s < 0) vp += (n-1)*s;
        BLASNAME(dscal) (BLASV(n),BLASV(x),BLASP(vp),BLASV(s));
    }
    static void DoScale(
        const std::complex<double> x, VectorView<std::complex<double> > v)
    {
        if (imag(x) == double(0)) {
            int n=v.size();
            if (n==0) return;
            TMVAssert(v.size()>0);
            int s=v.step();
            double xr = real(x);
            if (s == 0) { v[0] *= x; return; }
            std::complex<double>* vp = v.ptr();
            if (s < 0) vp += (n-1)*s;
            BLASNAME(zdscal) (BLASV(n),BLASV(xr),BLASP(vp),BLASV(s));
        } else {
            int n=v.size();
            if (n==0) return;
            TMVAssert(v.size()>0);
            int s=v.step();
            if (s == 0) { v[0] *= x; return; }
            std::complex<double>* vp = v.ptr();
            if (s < 0) vp += (n-1)*s;
            BLASNAME(zscal) (BLASV(n),BLASP(&x),BLASP(vp),BLASV(s));
        }
    }
#endif
#ifdef TMV_INST_FLOAT
    static void DoScale(const float x, VectorView<float> v)
    {
        int n=v.size();
        if (n==0) return;
        TMVAssert(v.size()>0);
        int s=v.step();
        if (s == 0) { v[0] *= x; return; }
        float* vp = v.ptr();
        if (s < 0) vp += (n-1)*s;
        BLASNAME(sscal) (BLASV(n),BLASV(x),BLASP(vp),BLASV(s));
    }
    static void DoScale(
        const std::complex<float> x, VectorView<std::complex<float> > v)
    {
        if (imag(x) == float(0)) {
            int n=v.size();
            if (n==0) return;
            TMVAssert(v.size()>0);
            int s=v.step();
            float xr = real(x);
            if (s == 0) { v[0] *= x; return; }
            std::complex<float>* vp = v.ptr();
            if (s < 0) vp += (n-1)*s;
            BLASNAME(csscal) (BLASV(n),BLASV(xr),BLASP(vp),BLASV(s));
        } else {
            int n=v.size();
            if (n==0) return;
            TMVAssert(v.size()>0);
            int s=v.step();
            if (s == 0) { v[0] *= x; return; }
            std::complex<float>* vp = v.ptr();
            if (s < 0) vp += (n-1)*s;
            BLASNAME(cscal) (BLASV(n),BLASP(&x),BLASP(vp),BLASV(s));
        }
    }
#endif
#endif

    template <class T> 
    void InstScale(const T x, VectorView<T> v)
    { DoScale(x,v); }

#define InstFile "TMV_ScaleV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


