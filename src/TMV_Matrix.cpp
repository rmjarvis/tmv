///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////



#include "TMV_Blas.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Divider.h"
#include "tmv/TMV_LUD.h"
#include "tmv/TMV_QRD.h"
#include "tmv/TMV_QRPD.h"
#include "tmv/TMV_SVD.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VIt.h"
#include "TMV_IntegerDet.h"
#include <iostream>

namespace tmv {

#define RT TMV_RealType(T)

#ifdef TMV_BLOCKSIZE
#define PERM_BLOCKSIZE TMV_BLOCKSIZE/2
#else
#define PERM_BLOCKSIZE 32
#endif

    //
    // Access
    //

    template <class T>
    T GenMatrix<T>::cref(ptrdiff_t i, ptrdiff_t j) const
    {
        const T* mi = cptr() + i*stepi() + j*stepj();
        return isconj() ? TMV_CONJ(*mi) : *mi;
    }

    template <class T, int A>
    typename MatrixView<T,A>::reference MatrixView<T,A>::ref(ptrdiff_t i, ptrdiff_t j) 
    {
        T* mi = ptr() + i*itssi + j*stepj();
        return RefHelper<T>::makeRef(mi,ct());
    }

    template <class T>
    void GenMatrix<T>::setDiv() const
    {
        if (!this->divIsSet()) {
            DivType dt = this->getDivType();
            TMVAssert(dt == tmv::LU || dt == tmv::QR ||
                      dt == tmv::QRP || dt == tmv::SV);
            switch (dt) {
              case LU : 
                   this->divider.reset(
                       new LUDiv<T>(*this,this->divIsInPlace())); 
                   break;
              case QR : 
                   this->divider.reset(
                       new QRDiv<T>(*this,this->divIsInPlace())); 
                   break;
              case QRP : 
                   this->divider.reset(
                       new QRPDiv<T>(*this,this->divIsInPlace())); 
                   break;
              case SV : 
                   this->divider.reset(
                       new SVDiv<T>(*this,this->divIsInPlace())); 
                   break;
              default : 
                   // The above assert should have already failed
                   // so go ahead and fall through.
                   break;
            }
        }
    }

#ifdef INST_INT
    template <>
    void GenMatrix<int>::setDiv() const
    { TMVAssert(TMV_FALSE); }
    template <>
    void GenMatrix<std::complex<int> >::setDiv() const
    { TMVAssert(TMV_FALSE); }
#endif

    // Note: These need to be in the .cpp file, not the .h file for
    // dynamic libraries.  Apparently the typeinfo used for dynamic_cast
    // doesn't get shared correctly by different modules, so the 
    // dynamic_cast fails when called in one module for an object that 
    // was created in a different module.
    //
    // So putting these functions here puts the dynamic cast in the shared
    // library, which is also where it is created (by setDiv() above).
    template <class T>
    bool GenMatrix<T>::divIsLUDiv() const
    { return dynamic_cast<const LUDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenMatrix<T>::divIsQRDiv() const
    { return dynamic_cast<const QRDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenMatrix<T>::divIsQRPDiv() const
    { return dynamic_cast<const QRPDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenMatrix<T>::divIsSVDiv() const
    { return dynamic_cast<const SVDiv<T>*>(this->getDiv()); }


#ifdef INST_INT
    template <>
    bool GenMatrix<int>::divIsLUDiv() const
    { return false; }
    template <>
    bool GenMatrix<int>::divIsQRDiv() const
    { return false; }
    template <>
    bool GenMatrix<int>::divIsQRPDiv() const
    { return false; }
    template <>
    bool GenMatrix<int>::divIsSVDiv() const
    { return false; }

    template <>
    bool GenMatrix<std::complex<int> >::divIsLUDiv() const
    { return false; }
    template <>
    bool GenMatrix<std::complex<int> >::divIsQRDiv() const
    { return false; }
    template <>
    bool GenMatrix<std::complex<int> >::divIsQRPDiv() const
    { return false; }
    template <>
    bool GenMatrix<std::complex<int> >::divIsSVDiv() const
    { return false; }
#endif

    //
    // OK? (SubMatrix, SubVector)
    //

    template <class T>
    bool GenMatrix<T>::hasSubMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
    {
        if (i1==i2 || j1==j2) return true; // no elements, so whatever...
        bool ok = true;
        if (istep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1 < 0 || i1 >= colsize()) {
            ok = false;
            std::cerr<<"first col element ("<<i1<<") must be in 0 -- ";
            std::cerr<<colsize()-1<<std::endl;
        }
        if (i2-istep < 0 || i2-istep >= colsize()) {
            ok = false;
            std::cerr<<"last col element ("<<i2-istep<<") must be in 0 -- ";
            std::cerr<<colsize()-1<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cerr<<"col range ("<<i2-i1<<") must be multiple of istep (";
            std::cerr<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cerr<<"n col elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
        }
        if (jstep == 0) {
            ok = false;
            std::cerr<<"jstep ("<<jstep<<") can not be 0\n";
        }
        if (j1 < 0 || j1 >= rowsize()) {
            ok = false;
            std::cerr<<"first row element ("<<j1<<") must be in 0 -- ";
            std::cerr<<rowsize()-1<<std::endl;
        }
        if (j2-jstep < 0 || j2-jstep >= rowsize()) {
            ok = false;
            std::cerr<<"last row element ("<<j2-jstep<<") must be in 0 -- ";
            std::cerr<<rowsize()-1<<std::endl;
        }
        if ((j2-j1)%jstep != 0) {
            ok = false;
            std::cerr<<"row range ("<<j2-j1<<") must be multiple of istep (";
            std::cerr<<jstep<<")\n";
        }
        if ((j2-j1)/jstep < 0) {
            ok = false;
            std::cerr<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
        }
        return ok;
    }

    template <class T>
    bool GenMatrix<T>::hasSubVector(
        ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const 
    {
        if (size==0) return true;
        bool ok = true;
        if (istep == 0 && jstep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") and jstep ("<<jstep;
            std::cerr<<") can not both be 0\n";
        }
        if (i < 0 || i >= colsize()) {
            ok = false;
            std::cerr<<"i ("<<i<<") must be in 0 -- "<<colsize()-1<<std::endl;
        }
        if (j < 0 || j >= rowsize()) {
            ok = false;
            std::cerr<<"j ("<<j<<") must be in 0 -- "<<rowsize()-1<<std::endl;
        }
        ptrdiff_t i2 = i+istep*(size-1);
        ptrdiff_t j2 = j+jstep*(size-1);
        if (i2 < 0 || i2 >= colsize()) {
            ok = false;
            std::cerr<<"last element's i ("<<i2<<") must be in 0 -- ";
            std::cerr<<colsize()-1<<std::endl;
        }
        if (j2 < 0 || j2 >= rowsize()) {
            ok = false;
            std::cerr<<"last element's j ("<<j2<<") must be in 0 -- ";
            std::cerr<<rowsize()-1<<std::endl;
        }
        return ok;
    }

    template <class T>
    bool ConstMatrixView<T,FortranStyle>::hasSubMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
    {
        if (i1==i2 || j1==j2) return true; // no elements, so whatever...
        bool ok = true;
        if (istep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1 < 1 || i1 > this->colsize()) {
            ok = false;
            std::cerr<<"first col element ("<<i1<<") must be in 1 -- ";
            std::cerr<<this->colsize()<<std::endl;
        }
        if (i2 < 1 || i2 > this->colsize()) {
            ok = false;
            std::cerr<<"last col element ("<<i2<<") must be in 1 -- ";
            std::cerr<<this->colsize()<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cerr<<"col range ("<<i2-i1<<") must be multiple of istep (";
            std::cerr<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cerr<<"n col elements ("<<(i2-i1)/istep+1<<") must be positive\n";
        }
        if (jstep == 0) {
            ok = false;
            std::cerr<<"jstep ("<<jstep<<") can not be 0\n";
        }
        if (j1 < 1 || j1 > this->rowsize()) {
            ok = false;
            std::cerr<<"first row element ("<<j1<<") must be in 1 -- ";
            std::cerr<<this->rowsize()<<std::endl;
        }
        if (j2 < 1 || j2 > this->rowsize()) {
            ok = false;
            std::cerr<<"last row element ("<<j2<<") must be in 1 -- ";
            std::cerr<<this->rowsize()<<std::endl;
        }
        if ((j2-j1)%jstep != 0) {
            ok = false;
            std::cerr<<"row range ("<<j2-j1<<") must be multiple of istep (";
            std::cerr<<jstep<<")\n";
        }
        if ((j2-j1)/jstep < 0) {
            ok = false;
            std::cerr<<"n row elements ("<<(j2-j1)/jstep+1<<") must be positive\n";
        }
        return ok;
    }

    template <class T>
    bool ConstMatrixView<T,FortranStyle>::hasSubVector(
        ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const 
    {
        if (size==0) return true;
        bool ok = true;
        if (istep == 0 && jstep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") and jstep ("<<jstep;
            std::cerr<<") can not both be 0\n";
        }
        if (i < 1 || i > this->colsize()) {
            ok = false;
            std::cerr<<"i ("<<i<<") must be in 1 -- "<<this->colsize()<<std::endl;
        }
        if (j < 1 || j > this->rowsize()) {
            ok = false;
            std::cerr<<"j ("<<j<<") must be in 1 -- "<<this->rowsize()<<std::endl;
        }
        ptrdiff_t i2 = i+istep*(size-1);
        ptrdiff_t j2 = j+jstep*(size-1);
        if (i2 < 1 || i2 > this->colsize()) {
            ok = false;
            std::cerr<<"last element's i ("<<i2<<") must be in 1 -- ";
            std::cerr<<this->colsize()<<std::endl;
        }
        if (j2 < 1 || j2 > this->rowsize()) {
            ok = false;
            std::cerr<<"last element's j ("<<j2<<") must be in 1 -- ";
            std::cerr<<this->rowsize()<<std::endl;
        }
        return ok;
    }


    //
    // Norms, det, etc.
    //

    template <class T>
    T GenMatrix<T>::det() const
    { return DivHelper<T>::det(); }

    template <class T>
    RT GenMatrix<T>::logDet(T* sign) const
    { return DivHelper<T>::logDet(sign); }

    template <class T>
    bool GenMatrix<T>::isSingular() const
    { return DivHelper<T>::isSingular(); }

#ifdef INST_INT
    template <>
    int GenMatrix<int>::det() const
    { return IntegerDet(*this); }

    template <>
    std::complex<int> GenMatrix<std::complex<int> >::det() const
    { return IntegerDet(*this); }

    template <>
    int GenMatrix<int>::logDet(int* ) const
    { TMVAssert(TMV_FALSE); return 0; }

    template <>
    int GenMatrix<std::complex<int> >::logDet(std::complex<int>* ) const
    { TMVAssert(TMV_FALSE); return 0; }

    template <>
    bool GenMatrix<int>::isSingular() const
    { return det() == 0; }

    template <>
    bool GenMatrix<std::complex<int> >::isSingular() const
    { return det() == 0; }
#endif

    template <class T>
    T GenMatrix<T>::sumElements() const
    {
        if (canLinearize()) return constLinearView().sumElements();
        else {
            T sum(0);
            if (iscm()) {
                const ptrdiff_t N = rowsize();
                for(ptrdiff_t j=0;j<N;++j) sum += col(j).sumElements();
            } else {
                const ptrdiff_t M = colsize();
                for(ptrdiff_t i=0;i<M;++i) sum += row(i).sumElements();
            }
            return sum;
        }
    }

    template <class T>
    static RT DoSumAbsElements(const GenMatrix<T>& m)
    {
        if (m.canLinearize()) return m.constLinearView().sumAbsElements();
        else {
            RT sum(0);
            if (m.iscm()) {
                const ptrdiff_t N = m.rowsize();
                for(ptrdiff_t j=0;j<N;++j) sum += m.col(j).sumAbsElements();
            } else {
                const ptrdiff_t M = m.colsize();
                for(ptrdiff_t i=0;i<M;++i) sum += m.row(i).sumAbsElements();
            }
            return sum;
        }
    }

    template <class T>
    static RT DoSumAbs2Elements(const GenMatrix<T>& m)
    {
        if (m.canLinearize()) return m.constLinearView().sumAbs2Elements();
        else {
            RT sum(0);
            if (m.iscm()) {
                const ptrdiff_t N = m.rowsize();
                for(ptrdiff_t j=0;j<N;++j) sum += m.col(j).sumAbs2Elements();
            } else {
                const ptrdiff_t M = m.colsize();
                for(ptrdiff_t i=0;i<M;++i) sum += m.row(i).sumAbs2Elements();
            }
            return sum;
        }
    }

#ifdef INST_INT
    static int DoSumAbsElements(const GenMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T>
    RT GenMatrix<T>::sumAbsElements() const
    { return DoSumAbsElements(*this); }

    template <class T>
    RT GenMatrix<T>::sumAbs2Elements() const
    { return DoSumAbs2Elements(*this); }

    template <class T>
    RT GenMatrix<T>::normSq(const RT scale) const
    {
        if (canLinearize()) return constLinearView().normSq(scale);
        else {
            RT sum(0);
            if (isrm()) {
                const ptrdiff_t M = colsize();
                for(ptrdiff_t i=0;i<M;++i) sum += row(i).normSq(scale);
            } else {
                const ptrdiff_t N = rowsize();
                for(ptrdiff_t j=0;j<N;++j) sum += col(j).normSq(scale);
            }
            return sum;
        }
    }

    template <class T>
    static RT NonLapMaxAbsElement(const GenMatrix<T>& m)
    {
        if (m.canLinearize()) return m.constLinearView().maxAbsElement();
        else {
            RT max(0);
            if (m.iscm()) {
                const ptrdiff_t N = m.rowsize();
                for(ptrdiff_t j=0;j<N;++j) {
                    RT temp = m.col(j).normInf();
                    if (temp > max) max = temp;
                }
            } else {
                const ptrdiff_t M = m.colsize();
                for(ptrdiff_t i=0;i<M;++i) {
                    RT temp = m.row(i).normInf();
                    if (temp > max) max = temp;
                }
            }
            return max;
        }
    }

    template <class T>
    static RT NonLapMaxAbs2Element(const GenMatrix<T>& m)
    {
        if (m.canLinearize()) return m.constLinearView().maxAbs2Element();
        else {
            RT max(0);
            if (m.iscm()) {
                const ptrdiff_t N = m.rowsize();
                for(ptrdiff_t j=0;j<N;++j) {
                    RT temp = m.col(j).maxAbs2Element();
                    if (temp > max) max = temp;
                }
            } else {
                const ptrdiff_t M = m.colsize();
                for(ptrdiff_t i=0;i<M;++i) {
                    RT temp = m.row(i).maxAbs2Element();
                    if (temp > max) max = temp;
                }
            }
            return max;
        }
    }

    template <class T>
    static RT NonLapNorm1(const GenMatrix<T>& m)
    {
        RT max(0);
        const ptrdiff_t N = m.rowsize();
        for(ptrdiff_t j=0;j<N;++j) {
            RT temp = m.col(j).norm1();
            if (temp > max) max = temp;
        }
        return max;
    }

    template <class T>
    static RT NonLapNormF(const GenMatrix<T>& m)
    {
        const RT eps = TMV_Epsilon<T>();

        RT mmax = m.maxAbs2Element();
        RT norm;
        if (mmax == RT(0)) {
            // Then norm is also 0
            norm = RT(0);
        } else if (TMV_Underflow(mmax * mmax)) {
            // Then we need to rescale, since underflow has caused 
            // rounding errors.
            // Epsilon is a pure power of 2, so no rounding errors from 
            // rescaling.
            const RT inveps = RT(1)/eps;
            RT scale = inveps;
            mmax *= scale;
            const RT eps2 = eps*eps;
            while (mmax < eps2) { scale *= inveps; mmax *= inveps; }
            norm = TMV_SQRT(m.normSq(scale))/scale;
        } else if (RT(1) / mmax == RT(0)) {
            // Then mmax is already inf, so no hope of making it more accurate.
            norm = mmax;
        } else if (RT(1) / (mmax*mmax) == RT(0)) {
            // Then we have overflow, so we need to rescale:
            const RT inveps = RT(1)/eps;
            RT scale = eps;
            mmax *= scale;
            while (mmax > inveps) { scale *= eps; mmax *= eps; }
            norm = TMV_SQRT(m.normSq(scale))/scale;
        } else {
            norm = TMV_SQRT(m.normSq());
        }
        return norm;
    }

    template <class T>
    static inline RT NonLapNormInf(const GenMatrix<T>& m)
    { return NonLapNorm1(m.transpose()); }

#ifdef XLAP
    template <class T>
    static RT LapNorm(const char c, const GenMatrix<T>& m)
    {
        switch(c) {
          case 'M' : return NonLapMaxAbsElement(m);
          case '1' : return NonLapNorm1(m);
          case 'F' : return NonLapNormF(m);
          case 'I' : return NonLapNormInf(m);
          default : TMVAssert(TMV_FALSE); 
        }
        return RT(0);
    }
#ifdef INST_DOUBLE
    template <>
    double LapNorm(const char c, const GenMatrix<double>& m)
    {
        TMVAssert(m.iscm());
        char cc = c;
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
#ifndef LAPNOWORK
        int lwork = c=='I' ? M : 0;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        double norm = LAPNAME(dlange) (
            LAPCM LAPV(cc),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
        return norm;
    }
    template <>
    double LapNorm(const char c, const GenMatrix<std::complex<double> >& m)
    {
        TMVAssert(m.iscm());
        char cc = c;
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
#ifndef LAPNOWORK
        int lwork = c=='I' ? M : 0;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        double norm = LAPNAME(zlange) (
            LAPCM LAPV(cc),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
        return norm;
    }
#endif
#ifdef INST_FLOAT
    template <>
    float LapNorm(const char c, const GenMatrix<float>& m)
    {
        TMVAssert(m.iscm());
        char cc = c;
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
#ifndef LAPNOWORK
        int lwork = c=='I' ? M : 0;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        double norm = LAPNAME(slange) (
            LAPCM LAPV(cc),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
        return norm;
    }
    template <>
    float LapNorm(const char c, const GenMatrix<std::complex<float> >& m)
    {
        TMVAssert(m.iscm());
        char cc = c;
        int M = m.colsize();
        int N = m.rowsize();
        int lda = m.stepj();
#ifndef LAPNOWORK
        int lwork = c=='I' ? M : 0;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        double norm = LAPNAME(clange) (
            LAPCM LAPV(cc),LAPV(M),LAPV(N),
            LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
        return norm;
    }
#endif
#endif // XLAP

#ifdef INST_INT
    static int NonLapNormF(const GenMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNormF(const GenMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNorm1(const GenMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapMaxAbsElement(const GenMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T>
    RT GenMatrix<T>::maxAbsElement() const
    {
#ifdef XLAP
        if (isrm() && stepi() > 0) return LapNorm('M',transpose());
        else if (iscm() && stepj() > 0) return LapNorm('M',*this);
        else
#endif
            return NonLapMaxAbsElement(*this);
    }
    template <class T>
    RT GenMatrix<T>::maxAbs2Element() const
    {
#ifdef XLAP
        if (Traits<T>::iscomplex) return NonLapMaxAbs2Element(*this);
        else if (isrm() && stepi() > 0) return LapNorm('M',transpose());
        else if (iscm() && stepj() > 0) return LapNorm('M',*this);
        else
#endif
            return NonLapMaxAbs2Element(*this);
    }
    template <class T>
    RT GenMatrix<T>::norm1() const
    {
#ifdef XLAP
        if (isrm() && stepi() > 0) return LapNorm('I',transpose());
        else if (iscm() && stepj() > 0) return LapNorm('1',*this);
        else
#endif
            return NonLapNorm1(*this);
    }
    template <class T>
    RT GenMatrix<T>::normF() const
    {
#ifdef XLAP
        if (isrm() && stepi() > 0) return LapNorm('F',transpose());
        else if (iscm() && stepj() > 0) return LapNorm('F',*this);
        else
#endif
            return NonLapNormF(*this);
    }

    template <class T>
    static RT DoNorm2(const GenMatrix<T>& m)
    {
        if (m.colsize() < m.rowsize()) return DoNorm2(m.transpose());
        if (m.rowsize() == 0) return RT(0);
        Matrix<T> m2(m);
        DiagMatrix<RT> S(m.rowsize());
        SV_Decompose(m2.view(),S.view(),false);
        return S(0);
    }

    template <class T>
    static RT DoCondition(const GenMatrix<T>& m) 
    {
        if (m.colsize() < m.rowsize()) return DoCondition(m.transpose());
        if (m.rowsize() == 0) return RT(1);
        Matrix<T> m2(m);
        DiagMatrix<RT> S(m.rowsize());
        SV_Decompose(m2.view(),S.view(),false);
        return S(0)/S(S.size()-1);
    }

#ifdef INST_INT
    static int DoNorm2(const GenMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoCondition(const GenMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoNorm2(const GenMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoCondition(const GenMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T>
    RT GenMatrix<T>::doNorm2() const
    { return tmv::DoNorm2(*this); }

    template <class T>
    RT GenMatrix<T>::doCondition() const
    { return tmv::DoCondition(*this); }

    template <class T>
    QuotXM<T,T> GenMatrix<T>::QInverse() const
    { return QuotXM<T,T>(T(1),*this); }

    //
    // Modifying Functions
    //

    template <class T, int A>
    MatrixView<T,A>& MatrixView<T,A>::clip(RT thresh) 
    {
        TMVAssert(A==CStyle);
        if (this->canLinearize()) linearView().clip(thresh);
        else {
            if (this->isrm()) {
                const ptrdiff_t M = colsize();
                for(ptrdiff_t i=0;i<M;++i) row(i).clip(thresh);
            } else {
                const ptrdiff_t N = rowsize();
                for(ptrdiff_t j=0;j<N;++j) col(j).clip(thresh);
            }
        }
        return *this; 
    }

    template <class T, int A>
    MatrixView<T,A>& MatrixView<T,A>::setZero() 
    {
        TMVAssert(A==CStyle);
        if (this->canLinearize()) linearView().setZero();
        else {
            if (this->isrm()) {
                const ptrdiff_t M = colsize();
                for(ptrdiff_t i=0;i<M;++i) row(i).setZero();
            } else {
                const ptrdiff_t N = rowsize();
                for(ptrdiff_t j=0;j<N;++j) col(j).setZero();
            }
        }
        return *this; 
    }

    template <class T, int A>
    MatrixView<T,A>& MatrixView<T,A>::setAllTo(const T& x) 
    {
        TMVAssert(A==CStyle);
        if (this->canLinearize()) linearView().setAllTo(x);
        else {
            if (this->isrm()) {
                const ptrdiff_t M = colsize();
                for(ptrdiff_t i=0;i<M;++i) row(i).setAllTo(x); 
            } else  {
                const ptrdiff_t N = rowsize();
                for(ptrdiff_t j=0;j<N;++j) col(j).setAllTo(x); 
            }
        }
        return *this; 
    }

    template <class T, int A>
    MatrixView<T,A>& MatrixView<T,A>::addToAll(const T& x) 
    {
        TMVAssert(A==CStyle);
        if (this->canLinearize()) linearView().addToAll(x);
        else {
            if (this->isrm()) {
                const ptrdiff_t M = colsize();
                for(ptrdiff_t i=0;i<M;++i) row(i).addToAll(x); 
            } else  {
                const ptrdiff_t N = rowsize();
                for(ptrdiff_t j=0;j<N;++j) col(j).addToAll(x); 
            }
        }
        return *this; 
    }

    template <class T, int A>
    MatrixView<T,A>& MatrixView<T,A>::transposeSelf() 
    {
        TMVAssert(A==CStyle);
        TMVAssert(colsize() == rowsize());
        const ptrdiff_t M = colsize();
        for(ptrdiff_t i=1;i<M;++i) Swap(row(i,0,i),col(i,0,i));
        return *this; 
    }

    template <class T, int A>
    MatrixView<T,A>& MatrixView<T,A>::conjugateSelf() 
    {
        TMVAssert(A==CStyle);
        if (isComplex(T())) {
            if (this->canLinearize()) linearView().conjugateSelf();
            else {
                if (this->isrm()) {
                    const ptrdiff_t M = colsize();
                    for(ptrdiff_t i=0;i<M;++i) row(i).conjugateSelf();
                } else  {
                    const ptrdiff_t N = rowsize();
                    for(ptrdiff_t j=0;j<N;++j) col(j).conjugateSelf();
                }
            }
        }
        return *this; 
    }

    template <class T, int A>
    MatrixView<T,A>& MatrixView<T,A>::setToIdentity(const T& x) 
    {
        TMVAssert(A==CStyle);
        TMVAssert(colsize() == rowsize());
        setZero();
        diag().setAllTo(x);
        return *this;
    }

    template <class T, int A>
    MatrixView<T,A>& MatrixView<T,A>::permuteRows(
        const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2) 
    {
        TMVAssert(A==CStyle);
        TMVAssert(i2<=colsize());
        TMVAssert(i1<=i2);
        // This idea of doing the permutations a block at a time
        // is cribbed from the LAPack code.  It does speed things up 
        // quite a bit for large matrices.  On my machine where BLOCKSIZE=64
        // is optimal for most routines, blocks of 32 were optimal here,
        // so I use BLOCKSIZE/2 in general.
        const ptrdiff_t N = rowsize();
        const ptrdiff_t Nx = N/PERM_BLOCKSIZE*PERM_BLOCKSIZE;
        if (Nx != 0) {
            for(ptrdiff_t j=0;j<Nx;) {
                ptrdiff_t j2 = j+PERM_BLOCKSIZE;
                const ptrdiff_t* pi = p+i1;
                for(ptrdiff_t i=i1;i<i2;++i,++pi) {
                    TMVAssert(*pi < colsize());
                    colRange(j,j2).swapRows(i,*pi);
                }
                j = j2;
            }
        }
        if (Nx != N) {
            const ptrdiff_t* pi = p+i1;
            for(ptrdiff_t i=i1;i<i2;++i,++pi) {
                TMVAssert(*pi < colsize());
                colRange(Nx,N).swapRows(i,*pi);
            }
        }
        return *this;
    }

    template <class T, int A>
    MatrixView<T,A>& MatrixView<T,A>::reversePermuteRows(
        const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2) 
    {
        TMVAssert(A==CStyle);
        TMVAssert(i2<=colsize());
        TMVAssert(i1<=i2);
        const ptrdiff_t N = rowsize();
        const ptrdiff_t Nx = N/PERM_BLOCKSIZE*PERM_BLOCKSIZE;
        if (Nx != 0) {
            for(ptrdiff_t j=0;j<Nx;) {
                ptrdiff_t j2 = j+PERM_BLOCKSIZE;
                const ptrdiff_t* pi = p+i2;
                for(ptrdiff_t i=i2;i>i1;) {
                    --i; --pi;
                    TMVAssert(*pi < colsize());
                    colRange(j,j2).swapRows(i,*pi);
                }
                j = j2;
            }
        }
        if (Nx != N) {
            const ptrdiff_t* pi = p+i2;
            for(ptrdiff_t i=i2;i>i1;) {
                --i; --pi;
                TMVAssert(*pi < colsize());
                colRange(Nx,N).swapRows(i,*pi);
            }
        }
        return *this;
    }

    //
    // Copy Matrices
    //

    template <class T>
    static void NonLapCopy(const GenMatrix<T>& m1, MatrixView<T> m2)
    {
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m1.ct()==NonConj);
        TMVAssert(m2.ct()==NonConj);
        const ptrdiff_t M = m2.colsize();
        const ptrdiff_t N = m2.rowsize();

        if (m1.iscm() && m2.iscm()) {
            const T* p1 = m1.cptr();
            T* p2 = m2.ptr();
            const ptrdiff_t s1 = m1.stepj();
            const ptrdiff_t s2 = m2.stepj();
            for(ptrdiff_t j=0;j<N;++j,p1+=s1,p2+=s2) {
                std::copy(p1,p1+M,p2);
            }
        } else if (M > N) {
            if (shouldReverse(m1.stepi(),m2.stepi()))
                for(ptrdiff_t j=0;j<N;++j) 
                    DoCopySameType(m1.col(j).reverse(),m2.col(j).reverse());
            else
                for(ptrdiff_t j=0;j<N;++j) 
                    DoCopySameType(m1.col(j),m2.col(j));
        } else {
            if (shouldReverse(m1.stepj(),m2.stepj()))
                for(ptrdiff_t i=0;i<M;++i) 
                    DoCopySameType(m1.row(i).reverse(),m2.row(i).reverse());
            else
                for(ptrdiff_t i=0;i<M;++i) 
                    DoCopySameType(m1.row(i),m2.row(i));
        }
    }
#ifdef ELAP
    template <class T>
    static inline void LapCopy(const GenMatrix<T>& m1, MatrixView<T> m2)
    { NonLapCopy(m1,m2); }
#ifdef INST_DOUBLE
    template <>
    void LapCopy(const GenMatrix<double>& m1, MatrixView<double> m2)
    {
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m1.iscm());
        TMVAssert(m2.iscm());
        char c = 'A';
        int m = m1.colsize();
        int n = m1.rowsize();
        int ld1 = m1.stepj();
        int ld2 = m2.stepj();
        LAPNAME(dlacpy) (
            LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
            LAPP(m2.ptr()),LAPV(ld2));
    }
    template <>
    void LapCopy(
        const GenMatrix<std::complex<double> >& m1,
        MatrixView<std::complex<double> > m2)
    {
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m1.iscm());
        TMVAssert(m2.iscm());
        TMVAssert(m1.ct() == NonConj);
        TMVAssert(m2.ct() == NonConj);
        char c = 'A';
        int m = m1.colsize();
        int n = m1.rowsize();
        int ld1 = m1.stepj();
        int ld2 = m2.stepj();
        LAPNAME(zlacpy) (
            LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
            LAPP(m2.ptr()),LAPV(ld2));
    }
#endif
#ifdef INST_FLOAT
    template <>
    void LapCopy(const GenMatrix<float>& m1, MatrixView<float> m2)
    {
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m1.iscm());
        TMVAssert(m2.iscm());
        char c = 'A';
        int m = m1.colsize();
        int n = m1.rowsize();
        int ld1 = m1.stepj();
        int ld2 = m2.stepj();
        LAPNAME(slacpy) (
            LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
            LAPP(m2.ptr()),LAPV(ld2));
    }
    template <>
    void LapCopy(
        const GenMatrix<std::complex<float> >& m1,
        MatrixView<std::complex<float> > m2)
    {
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m1.iscm());
        TMVAssert(m2.iscm());
        TMVAssert(m1.ct() == NonConj);
        TMVAssert(m2.ct() == NonConj);
        char c = 'A';
        int m = m1.colsize();
        int n = m1.rowsize();
        int ld1 = m1.stepj();
        int ld2 = m2.stepj();
        LAPNAME(clacpy) (
            LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
            LAPP(m2.ptr()),LAPV(ld2));
    }
#endif
#endif
    template <class T>
    void DoCopySameType(const GenMatrix<T>& m1, MatrixView<T> m2)
    {
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m2.rowsize() > 0);
        TMVAssert(m2.colsize() > 0);
        TMVAssert(m1.ct() == NonConj);
        TMVAssert(m2.ct() == NonConj);
        TMVAssert(!m2.isSameAs(m1));

#ifdef ELAP
        if (m1.iscm() && m2.iscm() && m1.stepj()>0 && m2.stepj()>0) 
            LapCopy(m1,m2);
        else
#endif
            NonLapCopy(m1,m2);
    }

    // 
    // Swap
    //

    template <class T>
    void Swap(MatrixView<T> m1, MatrixView<T> m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        if (m1.canLinearize() && m2.canLinearize() &&
            m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) {
            TMVAssert(m1.linearView().size() == m2.linearView().size());
            TMVAssert(m1.stepi()==m2.stepi() && m1.stepj()==m2.stepj());
            Swap(m1.linearView(),m2.linearView());
        } else {
            if (m1.isrm() && m2.isrm()) {
                const ptrdiff_t M = m1.colsize();
                for(ptrdiff_t i=0;i<M;++i) Swap(m1.row(i),m2.row(i)); 
            } else {
                const ptrdiff_t N = m1.rowsize();
                for(ptrdiff_t j=0;j<N;++j) Swap(m1.col(j),m2.col(j)); 
            }
        }
    }

    //
    // m1 == m2
    //

    template <class T1, class T2>
    bool operator==(const GenMatrix<T1>& m1, const GenMatrix<T2>& m2)
    {
        if (m1.colsize() != m2.colsize()) return false;
        else if (m1.rowsize() != m2.rowsize()) return false;
        else if (m1.isSameAs(m2)) return true;
        else if (m1.stepi()==m2.stepi() && m1.stepj()==m2.stepj() &&
                 m1.canLinearize() && m2.canLinearize())
            return m1.constLinearView() == m2.constLinearView();
        else {
            const ptrdiff_t M = m1.colsize();
            for(ptrdiff_t i=0;i<M;++i) 
                if (m1.row(i) != m2.row(i)) return false;
            return true;  
        }
    }

    //
    // I/O
    //

    template <class T>
    void GenMatrix<T>::write(const TMV_Writer& writer) const
    {
        const ptrdiff_t M = colsize();
        const ptrdiff_t N = rowsize();
        writer.begin();
        writer.writeCode("M");
        writer.writeSize(M);
        writer.writeSize(N);
        writer.writeStart();
        for(ptrdiff_t i=0;i<M;++i) {
            writer.writeLParen();
            for(ptrdiff_t j=0;j<N;++j) {
                if (j > 0) writer.writeSpace();
                writer.writeValue(cref(i,j));
            }
            writer.writeRParen();
            if (i < M-1) writer.writeRowEnd();
        }
        writer.writeFinal();
        writer.end();
    }

#ifndef NOTHROW
    template <class T>
    class MatrixReadError : public ReadError
    {
    public :
        Matrix<T> m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t cs,rs;
        bool is,iseof,isbad;

        MatrixReadError(std::istream& _is) throw() :
            ReadError("Matrix."),
            i(0), j(0), cs(0), rs(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        MatrixReadError(
            std::istream& _is, 
            const std::string& _e, const std::string& _g) throw() :
            ReadError("Matrix."),
            i(0), j(0), exp(_e), got(_g), cs(0), rs(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        MatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenMatrix<T>& _m, std::istream& _is) throw() :
            ReadError("Matrix."),
            m(_m), i(_i), j(_j),
            cs(m.colsize()), rs(m.rowsize()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        MatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenMatrix<T>& _m, std::istream& _is, 
            const std::string& _e, const std::string& _g) throw() :
            ReadError("Matrix."),
            m(_m), i(_i), j(_j), exp(_e), got(_g),
            cs(m.colsize()), rs(m.rowsize()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        MatrixReadError(
            const GenMatrix<T>& _m,
            std::istream& _is, ptrdiff_t _cs, ptrdiff_t _rs) throw() :
            ReadError("Matrix."),
            m(_m), i(0), j(0), cs(_cs), rs(_rs),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        MatrixReadError(const MatrixReadError<T>& rhs) :
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got),
            cs(rhs.cs), rs(rhs.rs),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~MatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for Matrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (cs != m.colsize()) {
                os<<"Wrong column size: expected "<<m.colsize()<<
                    ", got "<<cs<<".\n";
            }
            if (rs != m.rowsize()) {
                os<<"Wrong row size: expected "<<m.rowsize()<<
                    ", got "<<rs<<".\n";
            }
            if (!is) {
                if (iseof) {
                    os<<"Input stream reached end-of-file prematurely.\n";
                } else if (isbad) {
                    os<<"Input stream is corrupted.\n";
                } else {
                    os<<"Input stream cannot read next character.\n";
                }
            }
            if (m.colsize() > 0 || m.rowsize() > 0) {
                const ptrdiff_t N = m.rowsize();
                os<<"The portion of the Matrix which was successfully "
                    "read is: \n";
                for(ptrdiff_t ii=0;ii<i;++ii) {
                    os<<"( ";
                    for(ptrdiff_t jj=0;jj<N;++jj) os<<' '<<m.cref(ii,jj)<<' ';
                    os<<" )\n";
                }
                os<<"( ";
                for(ptrdiff_t jj=0;jj<j;++jj) os<<' '<<m.cref(i,jj)<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    template <class T>
    static void FinishRead(const TMV_Reader& reader, MatrixView<T> m) 
    {
        const ptrdiff_t M = m.colsize();
        const ptrdiff_t N = m.rowsize();
        std::string exp, got;
        T temp;
        if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw MatrixReadError<T>(0,0,m,reader.getis(),exp,got);
#endif
        }
        for(ptrdiff_t i=0;i<M;++i) {
            if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw MatrixReadError<T>(i,0,m,reader.getis(),exp,got);
#endif
            }        
            for(ptrdiff_t j=0;j<N;++j) {
                if (j>0) {
                    if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw MatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                }
                if (!reader.readValue(temp)) {
#ifdef NOTHROW
                    std::cerr<<"Matrix Read Error: reading value\n";
                    exit(1);
#else
                    throw MatrixReadError<T>(i,j,m,reader.getis());
#endif
                }
                m.ref(i,j) = temp;
            }
            if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw MatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
            if (i < M-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw MatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
        }
        if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw MatrixReadError<T>(M,0,m,reader.getis(),exp,got);
#endif
        }
    }

    template <class T, int A>
    void Matrix<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        if (!reader.readCode("M",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw MatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t cs=colsize(), rs=rowsize();
        if (!reader.readSize(cs,exp,got) || 
            !reader.readSize(rs,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: reading size\n";
            exit(1);
#else
            throw MatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (cs != colsize() || rs != rowsize()) resize(cs,rs);
        MatrixView<T> v = view();
        FinishRead(reader,v);
    }

    template <class T, int A>
    void MatrixView<T,A>::read(const TMV_Reader& reader) 
    {
        std::string exp,got;
        if (!reader.readCode("M",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw MatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t cs=colsize(), rs=rowsize();
        if (!reader.readSize(cs,exp,got) || 
            !reader.readSize(rs,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: reading size\n";
            exit(1);
#else
            throw MatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (cs != colsize() || rs != rowsize()) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: wrong size\n";
            exit(1);
#else
            throw MatrixReadError<T>(*this,reader.getis(),cs,rs);
#endif
        }
        FinishRead(reader,*this);
    }

#undef RT

#define InstFile "TMV_Matrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


