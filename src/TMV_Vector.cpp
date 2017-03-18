
//#define PRINTALGO_MinMax

#include "TMV_Blas.h"
#include <iostream>
#include <cstring> // for memmove
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SmallVector.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_NormV.h"
#include "tmv/TMV_Norm.h"
#include "tmv/TMV_SortV.h"
#include "tmv/TMV_SwapV.h"
#include "tmv/TMV_CopyV.h"
#include "tmv/TMV_ReverseV.h"
#include "tmv/TMV_ConjugateV.h"
#include "tmv/TMV_VectorIO.h"
#include "tmv/TMV_ScaleV.h"

namespace tmv {

    // First some things from TMV_Blas.h
    void LAP_Results(int Lap_info, const char* fn)
    {
        if (Lap_info < 0) {
#ifdef NOTHROW
            std::cerr<<"info < 0 returned by LAPACK function "<<fn<<std::endl;
            exit(1);
#else
            throw Error("info < 0 returned by LAPACK function ",fn);
#endif
        }
    }

    void LAP_Results(
        int Lap_info, int lwork_opt, int m, int n, int lwork, const char* fn)
    {
        LAP_Results(Lap_info,fn);
        if (lwork_opt > lwork) {
            std::ostringstream s;
            s << "LAPACK function " << fn << 
                " requested more workspace than provided";
            TMV_Warning(s.str());
            s.str(std::string());
            s<<"for matrix with m,n = "<<m<<','<<n<<std::endl;
            TMV_Warning(s.str());
            s.str(std::string());
            s<<"Given: "<<lwork<<", requested "<<lwork_opt<<std::endl;
            TMV_Warning(s.str());
        }
    }




    //
    // Copy Vectors
    //

    template <class T1, int C1, class T2> 
    static void DoInstCopy(
        const ConstVectorView<T1,C1>& v1, VectorView<T2> v2)
    {
        if (v2.step() == 1) {
            VectorView<T2,Unit> v2u = v2.unitView();
            if (v1.step() == 1) InlineCopy(v1.unitView(),v2u);
            else InlineCopy(v1,v2u);
        } else 
            if (v1.step() == 1) InlineCopy(v1.unitView(),v2);
            else InlineCopy(v1,v2);
    }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
    static void DoInstCopy(
        const ConstVectorView<double>& v1, VectorView<double> v2)
    {
        int n=v2.size();
        if (n == 0) return;
        int s1=v1.step();
        int s2=v2.step();
        const double* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        double* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(dcopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    static void DoInstCopy(
        const ConstVectorView<std::complex<double> >& v1,
        VectorView<std::complex<double> > v2)
    {
        int n=v2.size();
        if (n == 0) return;
        int s1=v1.step();
        int s2=v2.step();
        const std::complex<double>* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        std::complex<double>* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(zcopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    static void DoInstCopy(
        const ConstVectorView<std::complex<double>,Conj>& v1, 
        VectorView<std::complex<double> > v2)
    {
        DoInstCopy(v1.conjugate(),v2);
        v2.conjugateSelf();
    }
#endif
#ifdef TMV_INST_FLOAT
    static void DoInstCopy(
        const ConstVectorView<float>& v1, VectorView<float> v2)
    {
        int n=v2.size();
        if (n == 0) return;
        int s1=v1.step();
        int s2=v2.step();
        const float* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        float* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(scopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    static void DoInstCopy(
        const ConstVectorView<std::complex<float> >& v1,
        VectorView<std::complex<float> > v2)
    {
        int n=v2.size();
        if (n == 0) return;
        int s1=v1.step();
        int s2=v2.step();
        const std::complex<float>* v1p = v1.cptr();
        if (s1 < 0) v1p += (n-1)*s1;
        std::complex<float>* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(ccopy) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    static void DoInstCopy(
        const ConstVectorView<std::complex<float>,Conj>& v1, 
        VectorView<std::complex<float> > v2)
    {
        DoInstCopy(v1.conjugate(),v2);
        v2.conjugateSelf();
    }
#endif // FLOAT
#endif // BLAS

    static inline bool shouldReverse(const ptrdiff_t step1, const ptrdiff_t step2)
    {
        return ( (step2 < 0 && (step1 != 1 || step2 == -1)) ||
                 (step1 == -1 && step2 != 1) );
    }

    template <class T1, int C1, class T2> 
    void InstCopy(const ConstVectorView<T1,C1>& v1, VectorView<T2> v2)
    { 
        if (shouldReverse(v1.step(),v2.step())) 
            DoInstCopy(v1.reverse(),v2.reverse());
        else 
            DoInstCopy(v1,v2); 
    }

    template <class T1, int C1, class T2> 
    void InstAliasCopy(
        const ConstVectorView<T1,C1>& v1, VectorView<T2> v2)
    { InlineAliasCopy(v1,v2); }


    //
    // Swap Vectors
    //

    template <class T, int C> 
    static void DoInstSwap(VectorView<T,C> v1, VectorView<T> v2)
    {
        if (v1.step() == 1) {
            VectorView<T,C|Unit> v1u = v1.unitView();
            if (v2.step() == 1) {
                VectorView<T,Unit> v2u = v2.unitView();
                InlineSwap(v1u,v2u);
            } else 
                InlineSwap(v1u,v2);
        } else {
            if (v2.step() == 1) {
                VectorView<T,Unit> v2u = v2.unitView();
                InlineSwap(v1,v2u);
            } else 
                InlineSwap(v1,v2);
        }
    }

#ifdef BLAS
#ifdef TMV_INST_DOUBLE
    static void DoInstSwap(VectorView<double> v1, VectorView<double> v2)
    {
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        double* v1p = v1.ptr();
        if (s1 < 0) v1p += (n-1)*s1;
        double* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(dswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    static void DoInstSwap(
        VectorView<std::complex<double> > v1, 
        VectorView<std::complex<double> > v2)
    {
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        std::complex<double>* v1p = v1.ptr();
        if (s1 < 0) v1p += (n-1)*s1;
        std::complex<double>* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(zswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    static void DoInstSwap(
        VectorView<std::complex<double>,Conj> v1, 
        VectorView<std::complex<double> > v2)
    {
        DoInstSwap(v1.conjugate(),v2);
        v1.conjugateSelf();
        v2.conjugateSelf();
    }
#endif
#ifdef TMV_INST_FLOAT
    static void DoInstSwap(VectorView<float> v1, VectorView<float> v2)
    {
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        float* v1p = v1.ptr();
        if (s1 < 0) v1p += (n-1)*s1;
        float* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(sswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    static void DoInstSwap(
        VectorView<std::complex<float> > v1, 
        VectorView<std::complex<float> > v2)
    {
        int n=v2.size();
        int s1=v1.step();
        int s2=v2.step();
        std::complex<float>* v1p = v1.ptr();
        if (s1 < 0) v1p += (n-1)*s1;
        std::complex<float>* v2p = v2.ptr();
        if (s2 < 0) v2p += (n-1)*s2;
        BLASNAME(cswap) (BLASV(n),BLASP(v1p),BLASV(s1),BLASP(v2p),BLASV(s2));
    }
    static void DoInstSwap(
        VectorView<std::complex<float>,Conj> v1, 
        VectorView<std::complex<float> > v2)
    {
        DoInstSwap(v1.conjugate(),v2);
        v1.conjugateSelf();
        v2.conjugateSelf();
    }
#endif
#endif // BLAS

    template <class T, int C> 
    void InstSwap(VectorView<T,C> v1, VectorView<T> v2)
    { 
        if (shouldReverse(v1.step(),v2.step())) 
            DoInstSwap(v1.reverse(),v2.reverse());
        else 
            DoInstSwap(v1,v2); 
    }

    template <class T, int C> 
    void InstAliasSwap(VectorView<T,C> v1, VectorView<T> v2)
    { InlineAliasSwap(v1,v2); }


    //
    // ReverseSelf
    //

    template <class T> 
    void InstReverseSelf(VectorView<T> v)
    {
        if (v.step() == 1) {
            VectorView<T,Unit> vu = v.unitView();
            InlineReverseSelf(vu);
        } else if (v.step() < 0) {
            InstReverseSelf(v.reverse());
        } else 
            InlineReverseSelf(v);
    }


    //
    // conjugateSelf
    //

    template <class T> 
    static void DoInstConjugateSelf(VectorView<T> v)
    {
        if (v.step() == 1) {
            VectorView<T,Unit> vu = v.unitView();
            InlineConjugateSelf(vu);
        } else 
            InlineConjugateSelf(v);
    }

#ifdef ELAP
#ifdef TMV_INST_DOUBLE
    static void DoInstConjugateSelf(VectorView<std::complex<double> > v)
    {
        int n = v.size();
        int s = v.step();
        std::complex<double>* vp = v.ptr();
        //if (s < 0) vp += (n-1)*s;
        LAPNAME(zlacgv) (LAPV(n),LAPP(vp),LAPV(s)); 
    }
#endif
#ifdef TMV_INST_FLOAT
    static void DoInstConjugateSelf(VectorView<std::complex<float> > v)
    {
        int n = v.size();
        int s = v.step();
        std::complex<float>* vp = v.ptr();
        //if (s < 0) vp += (n-1)*s;
        LAPNAME(clacgv) (LAPV(n),LAPP(vp),LAPV(s)); 
    }
#endif
#endif // ELAP

    template <class T> 
    void InstConjugateSelf(VectorView<T> v)
    {
        if (v.step() >= 0) DoInstConjugateSelf(v); 
        else DoInstConjugateSelf(v.reverse()); 
    }

    //
    // SumElements
    //

    template <class T> 
    T InstSumElements(const ConstVectorView<T>& v)
    {
        if (v.step() == 1) return InlineSumElements(v.unitView());
        else if (v.step() < 0) return InstSumElements(v.reverse());
        else return InlineSumElements(v);
    }


    //
    // SumAbsElements
    //

    template <class T> 
    static typename ConstVectorView<T>::float_type DoInstSumAbsElements(
        const ConstVectorView<T>& v)
    {
        if (v.step() == 1) return InlineSumAbsElements(v.unitView());
        else return InlineSumAbsElements(v);
    }

    template <class T> 
    static typename Traits<T>::real_type DoInstSumAbs2Elements(
        const ConstVectorView<T>& v)
    {
        if (v.step() == 1) return InlineSumAbs2Elements(v.unitView());
        else return InlineSumAbs2Elements(v);
    }

#ifdef BLAS
#ifndef BLASNORETURN
#ifdef TMV_INST_DOUBLE
    static double DoInstSumAbs2Elements(const ConstVectorView<double>& v)
    {
        int n = v.size();
        if (n == 0) return double(0);
        int s = v.step();
        // If s == 0, the BLAS standard is to return 0 for the sum,
        // rather than the correct value.  Weird.
        // The non-blas version works correctly, but just do it here anyway.
        if (s == 0) return n*TMV_ABS(v[0]);
        const double* vp = v.cptr();
        //if (s < 0) vp += (n-1)*s;
        return BLASNAME(dasum) (BLASV(n),BLASP(vp),BLASV(s));
    }
    static double DoInstSumAbsElements(const ConstVectorView<double>& v)
    { return DoInstSumAbs2Elements(v); }
    static double DoInstSumAbs2Elements(
        const ConstVectorView<std::complex<double> >& v)
    {
        int n = v.size();
        if (n == 0) return double(0);
        int s = v.step();
        if (s == 0) return n*TMV_ABS(v[0]);
        const std::complex<double>* vp = v.cptr();
        //if (s < 0) vp += (n-1)*s;
        return BLASNAME(dzasum) (BLASV(n),BLASP(vp),BLASV(s));
    }
#endif
#ifdef TMV_INST_FLOAT
    static float DoInstSumAbs2Elements(const ConstVectorView<float>& v)
    {
        int n = v.size();
        if (n == 0) return float(0);
        int s = v.step();
        if (s == 0) return n*TMV_ABS(v[0]);
        const float* vp = v.cptr();
        //if (s < 0) vp += (n-1)*s;
        return BLASNAME(sasum) (BLASV(n),BLASP(vp),BLASV(s));
    }
    static float DoInstSumAbsElements(const ConstVectorView<float>& v)
    { return DoInstSumAbs2Elements(v); }
    static float DoInstSumAbs2Elements(
        const ConstVectorView<std::complex<float> >& v)
    {
        int n = v.size();
        if (n == 0) return float(0);
        int s = v.step();
        if (s == 0) return n*TMV_ABS(v[0]);
        const std::complex<float>* vp = v.cptr();
        //if (s < 0) vp += (n-1)*s;
        return BLASNAME(scasum) (BLASV(n),BLASP(vp),BLASV(s));
    }
#endif
#ifdef XLAP
#ifdef TMV_INST_DOUBLE
    static double DoInstSumAbsElements(
        const ConstVectorView<std::complex<double> >& v)
    {
        int n = v.size();
        if (n == 0) return double(0);
        int s = v.step();
        if (s == 0) return n*TMV_ABS(v[0]);
        const std::complex<double>* vp = v.cptr();
        //if (s < 0) vp += (n-1)*s;
        return LAPNAME(dzsum1) (LAPV(n),LAPP(vp),LAPV(s)); 
    }
#endif
#ifdef TMV_INST_FLOAT
    static float DoInstSumAbsElements(
        const ConstVectorView<std::complex<float> >& v)
    {
        int n = v.size();
        if (n == 0) return float(0);
        int s = v.step();
        if (s == 0) return n*TMV_ABS(v[0]);
        const std::complex<float>* vp = v.cptr();
        //if (s < 0) vp += (n-1)*s;
        return LAPNAME(scsum1) (LAPV(n),LAPP(vp),LAPV(s)); 
    }
#endif
#endif // XLAP
#endif // BLASNORETURN
#endif // BLAS

    template <class T> 
    typename ConstVectorView<T>::float_type InstSumAbsElements(
        const ConstVectorView<T>& v)
    { 
        if (v.step() >= 0) return DoInstSumAbsElements(v); 
        else return DoInstSumAbsElements(v.reverse()); 
    }

    template <class T> 
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstVectorView<T>& v)
    {
        if (v.step() >= 0) return DoInstSumAbs2Elements(v); 
        else return DoInstSumAbs2Elements(v.reverse()); 
    }

    //
    // MaxElement
    //

    template <class T> 
    T InstMaxElement(const ConstVectorView<T>& v, ptrdiff_t* imax)
    {
        if (v.step() == 1) return InlineMaxElement(v.unitView(),imax);
        else return InlineMaxElement(v,imax);
    }

    //
    // MaxAbsElement
    //

    template <class T> 
    static typename ConstVectorView<T>::float_type DoInstMaxAbsElement(
        const ConstVectorView<T>& v, ptrdiff_t* imax)
    {
        if (v.step() == 1) return InlineMaxAbsElement(v.unitView(),imax);
        else return InlineMaxAbsElement(v,imax);
    }

    template <class T> 
    static typename Traits<T>::real_type DoInstMaxAbs2Element(
        const ConstVectorView<T>& v, ptrdiff_t* imax)
    {
        if (v.step() == 1) return InlineMaxAbs2Element(v.unitView(),imax);
        else return InlineMaxAbs2Element(v,imax);
    }

#ifdef BLAS
    // These return values seem to work, so I don't guard this segment 
    // with BLASNORETURN
#ifdef TMV_INST_DOUBLE
    static double DoInstMaxAbs2Element(
        const ConstVectorView<double>& v, ptrdiff_t* imax)
    {
        int n=v.size();
        if (n == 0) return double(0);
        int s=v.step();
        if (s == 0) {
            if (imax) *imax = 0;
            return TMV_ABS2(v[0]);
        } else if (s < 0) {
            double max = DoInstMaxAbs2Element(v.reverse(),imax);
            if (imax) *imax = n-1-(*imax);
            return max;
        } else {
            int i = BLASNAME(idamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
            --i;
#endif
            if (imax) *imax = i;
            return TMV_ABS2(v[i]);
        }
    }
    static double DoInstMaxAbsElement(
        const ConstVectorView<double>& v, ptrdiff_t* imax)
    { return DoInstMaxAbs2Element(v,imax); }
    static double DoInstMaxAbs2Element(
        const ConstVectorView<std::complex<double> >& v, ptrdiff_t* imax)
    {
        int n=v.size();
        if (n == 0) return double(0);
        int s=v.step();
        if (s == 0) {
            if (imax) *imax = 0;
            return TMV_ABS2(v[0]);
        } else if (s < 0) {
            double max = DoInstMaxAbs2Element(v.reverse(),imax);
            if (imax) *imax = n-1-(*imax);
            return max;
        } else {
            int i = BLASNAME(izamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
            --i;
#endif
            if (imax) *imax = i;
            return TMV_ABS2(v[i]);
        }
    }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
    static float DoInstMaxAbs2Element(
        const ConstVectorView<float>& v, ptrdiff_t* imax)
    {
        int n=v.size();
        if (n == 0) return float(0);
        int s=v.step();
        if (s == 0) {
            if (imax) *imax = 0;
            return TMV_ABS2(v[0]);
        } else if (s < 0) {
            float max = DoInstMaxAbs2Element(v.reverse(),imax);
            if (imax) *imax = n-1-(*imax);
            return max;
        } else {
            int i = BLASNAME(isamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
            --i;
#endif
            if (imax) *imax = i;
            return TMV_ABS2(v[i]);
        }
    }
    static float DoInstMaxAbsElement(
        const ConstVectorView<float>& v, ptrdiff_t* imax)
    { return DoInstMaxAbs2Element(v,imax); }
    static float DoInstMaxAbs2Element(
        const ConstVectorView<std::complex<float> >& v, ptrdiff_t* imax)
    {
        int n=v.size();
        if (n == 0) return float(0);
        int s=v.step();
        if (s == 0) {
            if (imax) *imax = 0;
            return TMV_ABS2(v[0]);
        } else if (s < 0) {
            float max = DoInstMaxAbs2Element(v.reverse(),imax);
            if (imax) *imax = n-1-(*imax);
            return max;
        } else {
            int i = BLASNAME(icamax) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
            --i;
#endif
            if (imax) *imax = i;
            return TMV_ABS2(v[i]);
        }
    }
#endif // FLOAT
#endif // BLAS

    template <class T> 
    typename ConstVectorView<T>::float_type InstMaxAbsElement(
        const ConstVectorView<T>& v, ptrdiff_t* imax)
    { return DoInstMaxAbsElement(v,imax); }

    template <class T> 
    typename Traits<T>::real_type InstMaxAbs2Element(
        const ConstVectorView<T>& v, ptrdiff_t* imax)
    { return DoInstMaxAbs2Element(v,imax); }


    //
    // MinElement
    //

    template <class T> 
    T InstMinElement(const ConstVectorView<T>& v, ptrdiff_t* imin)
    {
        if (v.step() == 1) return InlineMinElement(v.unitView(),imin);
        else return InlineMinElement(v,imin);
    }


    //
    // MinAbsElement
    //

    template <class T> 
    static typename ConstVectorView<T>::float_type DoInstMinAbsElement(
        const ConstVectorView<T>& v, ptrdiff_t* imin)
    {
        if (v.step() == 1) return InlineMinAbsElement(v.unitView(),imin);
        else return InlineMinAbsElement(v,imin);
    }

    template <class T> 
    typename Traits<T>::real_type DoInstMinAbs2Element(
        const ConstVectorView<T>& v, ptrdiff_t* imin)
    {
        if (v.step() == 1) return InlineMinAbs2Element(v.unitView(),imin);
        else return InlineMinAbs2Element(v,imin);
    }

#ifdef BLAS
#ifdef BLASIDAMIN
#ifdef TMV_INST_DOUBLE
    double DoInstMinAbs2Element(
        const ConstVectorView<double>& v, ptrdiff_t* imin)
    {
        int n=v.size();
        if (n == 0) return double(0); 
        int s=v.step();
        if (s == 0) {
            if (imin) *imin = 0;
            return TMV_ABS2(v[0]);
        } else if (s < 0) {
            double min = v.reverse().minAbsElement(imin);
            if (imin) *imin = n-1-(*imin);
            return min;
        } else {
            int i = BLASNAME(idamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
            --i;
#endif
            if (imin) *imin = i;
            return TMV_ABS2(v[i]);
        }
    }
    double DoInstMinAbsElement(
        const ConstVectorView<double>& v, ptrdiff_t* imin)
    { DoInstMinAbs2Element(v,imin); }
    double DoInstMinAbs2Element(
        const ConstVectorView<std::complex<double> >& v, ptrdiff_t* imin)
    {
        int n=v.size();
        if (n == 0) return double(0); 
        int s=v.step();
        if (s == 0) {
            if (imin) *imin = 0;
            return TMV_ABS2(v[0]);
        } else if (s < 0) {
            double min = v.reverse().minAbsElement(imin);
            if (imin) *imin = n-1-(*imin);
            return min;
        } else {
            int i = BLASNAME(izamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
            --i;
#endif
            if (imin) *imin = i;
            return TMV_ABS2(v[i]);
        }
    }
#endif // DOUBLE
#ifdef TMV_INST_FLOAT
    float DoInstMinAbs2Element(
        const ConstVectorView<float>& v, ptrdiff_t* imin)
    {
        int n=v.size();
        if (n == 0) return float(0); 
        int s=v.step();
        if (s == 0) {
            if (imin) *imin = 0;
            return TMV_ABS2(v[0]);
        } else if (s < 0) {
            float min = v.reverse().minAbsElement(imin);
            if (imin) *imin = n-1-(*imin);
            return min;
        } else {
            int i = BLASNAME(isamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
            --i;
#endif
            if (imin) *imin = i;
            return TMV_ABS2(v[i]);
        }
    }
    float DoInstMinAbsElement(
        const ConstVectorView<float>& v, ptrdiff_t* imin)
    { DoInstMinAbs2Element(v,imin); }
    float DoInstMinAbs2Element(
        const ConstVectorView<std::complex<float> >& v, ptrdiff_t* imin)
    {
        int n=v.size();
        if (n == 0) return float(0); 
        int s=v.step();
        if (s == 0) {
            if (imin) *imin = 0;
            return TMV_ABS2(v[0]);
        } else if (s < 0) {
            float min = v.reverse().minAbsElement(imin);
            if (imin) *imin = n-1-(*imin);
            return min;
        } else {
            int i = BLASNAME(icamin) (BLASV(n),BLASP(v.cptr()),BLASV(s));
#ifndef CBLAS
            --i;
#endif
            if (imin) *imin = i;
            return TMV_ABS2(v[i]);
        }
    }
#endif // FLOAT
#endif // BLASIDAMIN
#endif // BLAS

    template <class T> 
    typename ConstVectorView<T>::float_type InstMinAbsElement(
        const ConstVectorView<T>& v, ptrdiff_t* imin)
    { return DoInstMinAbsElement(v,imin); }

    template <class T> 
    typename Traits<T>::real_type InstMinAbs2Element(
        const ConstVectorView<T>& v, ptrdiff_t* imin)
    { return DoInstMinAbs2Element(v,imin); }


    //
    // NormSq
    //

    template <class T> 
    typename Traits<T>::real_type InstNormSq(const ConstVectorView<T>& v)
    {
        if (v.step() == 1) return InlineNormSq(v.unitView());
        else if (v.step() < 0) return InstNormSq(v.reverse());
        else return InlineNormSq(v);
    }

    template <class T> 
    typename ConstVectorView<T>::float_type InstNormSq(
        const ConstVectorView<T>& v,
        const typename ConstVectorView<T>::float_type scale)
    {
        if (v.step() == 1) return InlineNormSq(v.unitView(),scale);
        else if (v.step() < 0) return InstNormSq(v.reverse(),scale);
        else return InlineNormSq(v,scale);
    }


    //
    // Norm2
    //

    template <class T> 
    static typename ConstVectorView<T>::float_type DoInstNorm2(
        const ConstVectorView<T>& v)
    {
        if (v.step() == 1) return InlineNorm2(v.unitView());
        else return InlineNorm2(v);
    }

#ifdef BLAS
#ifndef BLASNORETURN
#ifdef TMV_INST_DOUBLE
    static double DoInstNorm2(const ConstVectorView<double>& v)
    {
        int n=v.size();
        int s=v.step();
        const double* vp = v.cptr();
        //if (s < 0) vp += (n-1)*s;
        return BLASNAME(dnrm2) (BLASV(n),BLASP(vp),BLASV(s));
    }
    static double DoInstNorm2(const ConstVectorView<std::complex<double> >& v)
    {
        int n=v.size();
        int s=v.step();
        const std::complex<double>* vp = v.cptr();
        //if (s < 0) vp += (n-1)*s;
        return BLASNAME(dznrm2) (BLASV(n),BLASP(vp),BLASV(s));
    }
#endif
#ifdef TMV_INST_FLOAT
    static float DoInstNorm2(const ConstVectorView<float>& v)
    {
        int n=v.size();
        int s=v.step();
        const float* vp = v.cptr();
        //if (s < 0) vp += (n-1)*s;
        return BLASNAME(snrm2) (BLASV(n),BLASP(vp),BLASV(s));
    }
    static float DoInstNorm2(const ConstVectorView<std::complex<float> >& v)
    {
        int n=v.size();
        int s=v.step();
        const std::complex<float>* vp = v.cptr();
        //if (s < 0) vp += (n-1)*s;
        return BLASNAME(scnrm2) (BLASV(n),BLASP(vp),BLASV(s));
    }
#endif
#endif // BLASNORETURN
#endif // BLAS

    template <class T> 
    typename ConstVectorView<T>::float_type InstNorm2(
        const ConstVectorView<T>& v)
    { 
        if (v.step() >= 0) return DoInstNorm2(v); 
        else return DoInstNorm2(v.reverse()); 
    }

    //
    // Sort
    //

    template <class T> 
    void InstSort(VectorView<T> v, ADType ad, CompType comp)
    {
        if (v.step() == 1) {
            VectorView<T,Unit> vu = v.unitView();
            InlineSort(vu,ad,comp);
        } else {
            InlineSort(v,ad,comp); 
        }
    }
    template <class T> 
    void InstSort(VectorView<T> v, ptrdiff_t* P, ADType ad, CompType comp)
    {
        if (v.step() == 1) {
            VectorView<T,Unit> vu = v.unitView();
            InlineSort(vu,P,ad,comp);
        } else 
            InlineSort(v,P,ad,comp); 
    }


    //
    // I/O
    //

    template <class T, int C> 
    void InstWrite(const TMV_Writer& writer, const ConstVectorView<T,C>& v)
    {
        if (v.step() == 1) InlineWrite(writer,v.unitView()); 
        else InlineWrite(writer,v); 
    }

    template <class T>
    void InstRead(const TMV_Reader& reader, VectorView<T> v)
    {
        if (v.step() == 1) {
            VectorView<T,Unit> vu = v.unitView();
            InlineRead(reader,vu); 
        } else 
            InlineRead(reader,v); 
    }

#define InstFile "TMV_Vector.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


