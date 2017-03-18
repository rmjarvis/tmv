
//#define PRINTALGO_NormB
//#define PRINTALGO_NormM

#include "TMV_Blas.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_SmallBandMatrix.h"
#include "tmv/TMV_BandMatrixIO.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_CopyB.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_ScaleB.h"
#include "tmv/TMV_SwapB.h"
#include "tmv/TMV_TransposeB.h"
#include "tmv/TMV_NormB.h"
#include "tmv/TMV_Norm.h"

//#include "tmv/TMV_BandLUD.h"
//#include "tmv/TMV_BandQRD.h"
//#include "tmv/TMV_BandSVD.h"

namespace tmv {

    // 
    // Copy BandMatrix
    //

    template <class T1, int C1, class T2>
    void InstCopy(
        const ConstBandMatrixView<T1,C1>& m1, BandMatrixView<T2> m2)
    {
        if (m2.iscm()) {
            BandMatrixView<T2,ColMajor> m2cm = m2;
            if (m1.iscm()) {
                InlineCopy(m1.cmView(),m2cm);
            } else if (m1.isrm()) {
                InlineCopy(m1.rmView(),m2cm);
            } else if (m1.isdm()) {
                InlineCopy(m1.dmView(),m2cm);
            } else
                InlineCopy(m1,m2cm);
        } else if (m2.isrm()) {
            InstCopy(m1.transpose(),m2.transpose());
        } else if (m2.isdm()) {
            BandMatrixView<T2,DiagMajor> m2dm = m2;
            if (m1.iscm())
                InlineCopy(m1.cmView(),m2dm);
            else if (m1.isrm())
                InlineCopy(m1.rmView(),m2dm);
            else if (m1.isdm())
                InlineCopy(m1.dmView(),m2dm);
            else
                InlineCopy(m1,m2dm);
        } else {
            if (m1.iscm())
                InlineCopy(m1.cmView(),m2);
            else if (m1.isrm())
                InlineCopy(m1.rmView(),m2);
            else if (m1.isdm())
                InlineCopy(m1.dmView(),m2);
            else
                InlineCopy(m1,m2);
        }
    }

    template <class T1, int C1, class T2>
    void InstAliasCopy(
        const ConstBandMatrixView<T1,C1>& m1, BandMatrixView<T2> m2)
    { InlineAliasCopy(m1,m2); }


    // 
    // Swap BandMatrix
    //

    template <class T, int C>
    void InstSwap(BandMatrixView<T,C> m1, BandMatrixView<T> m2)
    {
        if (m2.iscm()) {
            BandMatrixView<T,ColMajor> m2cm = m2;
            if (m1.iscm()) {
                BandMatrixView<T,C|ColMajor> m1cm = m1;
                InlineSwap(m1cm,m2cm);
            } else if (m1.isrm()) {
                BandMatrixView<T,C|RowMajor> m1rm = m1;
                InlineSwap(m1rm,m2cm);
            } else if (m1.isdm()) {
                BandMatrixView<T,C|DiagMajor> m1dm = m1;
                InlineSwap(m1dm,m2cm);
            } else
                InlineSwap(m1,m2cm);
        } else if (m2.isrm()) {
            InstSwap(m1.transpose(),m2.transpose());
        } else if (m2.isdm()) {
            BandMatrixView<T,DiagMajor> m2dm = m2;
            if (m1.iscm()) {
                BandMatrixView<T,C|ColMajor> m1cm = m1;
                InlineSwap(m1cm,m2dm);
            } else if (m1.isrm()) {
                BandMatrixView<T,C|RowMajor> m1rm = m1;
                InlineSwap(m1rm,m2dm);
            } else if (m1.isdm()) {
                BandMatrixView<T,C|DiagMajor> m1dm = m1;
                InlineSwap(m1dm,m2dm);
            } else
                InlineSwap(m1,m2dm);
        } else {
            if (m1.iscm()) {
                BandMatrixView<T,C|ColMajor> m1cm = m1;
                InlineSwap(m1cm,m2);
            } else if (m1.isrm()) {
                BandMatrixView<T,C|RowMajor> m1rm = m1;
                InlineSwap(m1rm,m2);
            } else if (m1.isdm()) {
                BandMatrixView<T,C|DiagMajor> m1dm = m1;
                InlineSwap(m1dm,m2);
            } else
                InlineSwap(m1,m2);
        }
    }

    template <class T, int C>
    void InstAliasSwap(BandMatrixView<T,C> m1, BandMatrixView<T> m2)
    { InlineAliasSwap(m1,m2); }


    // 
    // TransposeSelf
    //

    template <class T>
    void InstTransposeSelf(BandMatrixView<T> m)
    {
        if (m.iscm()) {
            BandMatrixView<T,ColMajor> mcm = m;
            InlineTransposeSelf(mcm);
        } else if (m.isrm()) {
            BandMatrixView<T,ColMajor> mt = m.transpose();
            InlineTransposeSelf(mt);
        } else if (m.isdm()) {
            BandMatrixView<T,DiagMajor> mdm = m;
            InlineTransposeSelf(mdm);
        } else {
            InlineTransposeSelf(m);
        }
    }

    //
    // Norms
    //

    template <class T>
    static T DoInstSumElements(const ConstBandMatrixView<T>& m)
    {
        if (m.iscm()) return InlineSumElements(m.cmView());
        else if (m.isrm()) return InlineSumElements(m.transpose().cmView());
        else if (m.isdm()) return InlineSumElements(m.dmView());
        else return InlineSumElements(m);
    }

    template <class T>
    static typename ConstBandMatrixView<T>::float_type DoInstSumAbsElements(
        const ConstBandMatrixView<T>& m)
    {
        if (m.iscm()) return InlineSumAbsElements(m.cmView());
        else if (m.isrm()) return InlineSumAbsElements(m.transpose().cmView());
        else if (m.isdm()) return InlineSumAbsElements(m.dmView());
        else return InlineSumAbsElements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstSumAbs2Elements(
        const ConstBandMatrixView<T>& m)
    {
        if (m.iscm()) return InlineSumAbs2Elements(m.cmView());
        else if (m.isrm()) return InlineSumAbs2Elements(
            m.transpose().cmView());
        else if (m.isdm()) return InlineSumAbs2Elements(m.dmView());
        else return InlineSumAbs2Elements(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstNormSq(
        const ConstBandMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNormSq(m.cmView());
        else if (m.isrm()) return InlineNormSq(m.transpose().cmView());
        else if (m.isdm()) return InlineNormSq(m.dmView());
        else return InlineNormSq(m);
    }

    template <class T>
    static typename ConstBandMatrixView<T>::float_type DoInstNormSq(
        const ConstBandMatrixView<T>& m,
        const typename ConstBandMatrixView<T>::float_type scale)
    {
        if (m.iscm()) return InlineNormSq(m.cmView(),scale);
        else if (m.isrm()) return InlineNormSq(m.transpose().cmView(),scale);
        else if (m.isdm()) return InlineNormSq(m.dmView(),scale);
        else return InlineNormSq(m,scale);
    }

    template <class T>
    static typename ConstBandMatrixView<T>::float_type DoInstNormF(
        const ConstBandMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNormF(m.cmView());
        else if (m.isrm()) return InlineNormF(m.transpose().cmView());
        else if (m.isdm()) return InlineNormF(m.dmView());
        else return InlineNormF(m);
    }

    template <class T>
    static typename ConstBandMatrixView<T>::float_type DoInstMaxAbsElement(
        const ConstBandMatrixView<T>& m)
    {
        if (m.iscm()) return InlineMaxAbsElement(m.cmView());
        else if (m.isrm()) return InlineMaxAbsElement(m.transpose().cmView());
        else if (m.isdm()) return InlineMaxAbsElement(m.dmView());
        else return InlineMaxAbsElement(m);
    }

    template <class T>
    static typename Traits<T>::real_type DoInstMaxAbs2Element(
        const ConstBandMatrixView<T>& m)
    {
        if (m.iscm()) return InlineMaxAbs2Element(m.cmView());
        else if (m.isrm()) return InlineMaxAbs2Element(m.transpose().cmView());
        else if (m.isdm()) return InlineMaxAbs2Element(m.dmView());
        else return InlineMaxAbs2Element(m);
    }

    template <class T>
    static typename ConstBandMatrixView<T>::float_type DoInstNorm1(
        const ConstBandMatrixView<T>& m)
    {
        if (m.iscm()) return InlineNorm1(m.cmView());
        else if (m.isrm()) return InlineNorm1(m.rmView());
        else if (m.isdm()) return InlineNorm1(m.dmView());
        else return InlineNorm1(m);
    }


#ifdef XLAP
#ifdef TMV_INST_DOUBLE
    static double LapNorm(
        char c, const ConstBandMatrixView<double>& m)
    {
        if (m.colsize() > m.rowsize()) {
            BandMatrix<double,ColMajor> mc(
                m.colsize(),m.colsize(),m.nlo(),m.nhi(),0.);
            mc.colRange(0,m.rowsize()) = m;
            return LapNorm(c,mc.view());
        }
        if (m.rowsize() > m.colsize()) {
            BandMatrix<double,ColMajor> mc(
                m.rowsize(),m.rowsize(),m.nlo(),m.nhi(),0.);
            mc.rowRange(0,m.colsize()) = m;
            return LapNorm(c,mc.view());
        }
        if (!m.iscm() && !(m.isdm() && m.nlo()==1 && m.nhi()==1))
            return LapNorm(c,m.copy().view());

        if (m1.iscm() || m1.isrm()) {
            int n = m.colsize();
            int kl = m.nlo();
            int ku = m.nhi();
            int lda = m.stepj() + 1;
#ifndef LAPNOWORK
            int lwork = c=='I' : N : 0;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            double norm = LAPNAME(dlangb) (
                LAPCM LAPV(c),LAPV(n),LAPV(kl),LAPV(ku),
                LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()) LAP1);
            return norm;
        } else {
            char c = '1';
            int n = m.colsize();
            double norm = LAPNAME(dlangt) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),
                LAPP(m.cptr()+m.stepj()));
            return norm;
        }
    }
    static double LapNorm(
        char c, const ConstBandMatrixView<std::complex<double> >& m)
    {
        if (m.colsize() > m.rowsize()) {
            BandMatrix<std::complex<double>,ColMajor> mc(
                m.colsize(),m.colsize(),m.nlo(),m.nhi(),0.);
            mc.colRange(0,m.rowsize()) = m;
            return LapNorm(c,mc.view());
        }
        if (m.rowsize() > m.colsize()) {
            BandMatrix<std::complex<double>,ColMajor> mc(
                m.rowsize(),m.rowsize(),m.nlo(),m.nhi(),0.);
            mc.rowRange(0,m.colsize()) = m;
            return LapNorm(c,mc.view());
        }
        if (!m.iscm() && m.isrm()) {
            char c2 = c == 'I' : '1' : c=='1' : 'I' : c;
            return LapNorm(c,m.transpose());
        }
        if (!m.iscm() && !(m.isdm() && m.nlo()==1 && m.nhi()==1)) {
            return LapNorm(c,m.copy().view());
        }

        if (m1.iscm()) {
            int n = m.colsize();
            int kl = m.nlo();
            int ku = m.nhi();
            int lda = m.stepj() + 1;
#ifndef LAPNOWORK
            int lwork = c=='I' : N : 0;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            double norm = LAPNAME(zlangb) (
                LAPCM LAPV(c),LAPV(n),LAPV(kl),LAPV(ku),
                LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()) LAP1);
            return norm;
        } else {
            char c = '1';
            int n = m.colsize();
            double norm = LAPNAME(zlangt) (
                LAPCM LAPV(c),LAPV(n),
                LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),
                LAPP(m.cptr()+m.stepj()));
            return norm;
        }
    }
    static double DoInstNormF(const ConstBandMatrixView<double>& m)
    { return LapNorm('F',m); }
    static double DoInstNormF(
        const ConstBandMatrixView<std::complex<double> >& m)
    { return LapNorm('F',m); }
    static double DoInstMaxAbsElement(const ConstBandMatrixView<double>& m)
    { return LapNorm('M',m); }
    static double DoInstMaxAbsElement(
        const ConstBandMatrixView<std::complex<double> >& m)
    { return LapNorm('M',m); }
    static double DoInstMaxAbs2Element(const ConstBandMatrixView<double>& m)
    { return DoInstMaxAbsElement(m); }
    static double DoInstMaxAbs2Element(
        const ConstBandMatrixView<std::complex<double> >& m)
    { return InlineMaxAbs2Element(m); }
    static double DoInstNorm1(const ConstBandMatrixView<double>& m)
    { return LapNorm('1',m); }
    static double DoInstNorm1(
        const ConstBandMatrixView<std::complex<double> >& m)
    { return LapNorm('1',m); }
#endif // DOUBLE
#endif // XLAP

    template <class T>
    T InstSumElements(const ConstBandMatrixView<T>& m)
    { return DoInstSumElements(m); }

    template <class T>
    typename ConstBandMatrixView<T>::float_type InstSumAbsElements(
        const ConstBandMatrixView<T>& m)
    { return DoInstSumAbsElements(m); }

    template <class T>
    typename Traits<T>::real_type InstSumAbs2Elements(
        const ConstBandMatrixView<T>& m)
    { return DoInstSumAbs2Elements(m); }

    template <class T>
    typename Traits<T>::real_type InstNormSq(const ConstBandMatrixView<T>& m)
    { return DoInstNormSq(m); }

    template <class T>
    typename ConstBandMatrixView<T>::float_type InstNormSq(
        const ConstBandMatrixView<T>& m,
        const typename ConstBandMatrixView<T>::float_type scale)
    { return DoInstNormSq(m,scale); }

    template <class T>
    typename ConstBandMatrixView<T>::float_type InstNormF(
        const ConstBandMatrixView<T>& m)
    { return DoInstNormF(m); }

    template <class T>
    typename ConstBandMatrixView<T>::float_type InstMaxAbsElement(
        const ConstBandMatrixView<T>& m)
    { return DoInstMaxAbsElement(m); }

    template <class T>
    typename Traits<T>::real_type InstMaxAbs2Element(
        const ConstBandMatrixView<T>& m)
    { return DoInstMaxAbs2Element(m); }

    template <class T>
    typename ConstBandMatrixView<T>::float_type InstNorm1(
        const ConstBandMatrixView<T>& m)
    { return DoInstNorm1(m); }


    //
    // I/O
    //

    template <class T, int C>
    void InstWrite(
        const TMV_Writer& writer, const ConstBandMatrixView<T,C>& m)
    {
        if (m.isrm()) InlineWrite(writer,m.rmView());
        else InlineWrite(writer,m);
    }

    template <class T>
    void InstRead(const TMV_Reader& reader, BandMatrixView<T> m)
    {
        if (m.isrm()) {
            BandMatrixView<T,RowMajor> mrm = m.rmView();
            InlineRead(reader,mrm);
        } else
            InlineRead(reader,m);
    }


    //
    // BandMatrixDivHelper
    //

    template <class T>
    BandMatrixDivHelper2<T>::BandMatrixDivHelper2() {}

    template <class T>
    BandMatrixDivHelper2<T>::~BandMatrixDivHelper2() {}

    template <class T>
    bool BandMatrixDivHelper2<T>::mIsSquare() const 
    { return getConstView().isSquare(); }

    template <class T>
    Matrix<T> BandMatrixDivHelper2<T>::getMatrix() const
    { return Matrix<T>(getConstView()); }

    template <class T> 
    void BandMatrixDivHelper2<T>::setDiv() const
    {
        TMVStaticAssert(!Traits<T>::isinteger);
        if (!this->divIsSet()) {
            DivType dt = this->getDivType();
            TMVAssert(dt == tmv::LU || dt == tmv::QR ||
                      dt == tmv::SV);
            switch (dt) {
#if 0
              case tmv::LU :
                   this->divider.reset(new InstBandLUD<T>(
                           this->getConstView(),this->divIsInPlace()));
                   break;
              case tmv::QR :
                   this->divider.reset(new InstBandQRD<T>(
                           this->getConstView(),this->divIsInPlace()));
                   break;
              case tmv::SV :
                   this->divider.reset(new InstBandSVD<T>(
                           this->getConstView(),this->divIsInPlace()));
                   break;
#endif
              default :
                   // The above assert should have already failed.
                   // So go ahead and fall through.
                   break;
            }
        }
    }

#if 0
    template <class T>
    const InstBandLUD<T>& BandMatrixDivHelper2<T>::lud() const
    {
        this->divideUsing(LU);
        setDiv();
        TMVAssert(dynamic_cast<const InstBandLUD<T>*>(this->getDiv()));
        return static_cast<const InstBandLUD<T>&>(*this->getDiv());
    }

    template <class T>
    const InstBandQRD<T>& BandMatrixDivHelper2<T>::qrd() const
    {
        this->divideUsing(QR);
        setDiv();
        TMVAssert(dynamic_cast<const InstBandQRD<T>*>(this->getDiv()));
        return static_cast<const InstBandQRD<T>&>(*this->getDiv());
    }

    template <class T>
    const InstBandSVD<T>& BandMatrixDivHelper2<T>::svd() const
    {
        this->divideUsing(SV);
        setDiv();
        TMVAssert(dynamic_cast<const InstBandSVD<T>*>(this->getDiv()));
        return static_cast<const InstBandSVD<T>&>(*this->getDiv());
    }
#endif

#define InstFile "TMV_BandMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


