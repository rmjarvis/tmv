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


//#define XDEBUG


#include "TMV_Blas.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_SymLDLD.h"
#include "tmv/TMV_SymCHD.h"
#include "tmv/TMV_SymSVD.h"
#include "tmv/TMV_VIt.h"
#include "tmv/TMV_SymMatrixArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "TMV_IntegerDet.h"
#include <iostream>

#ifdef XDEBUG
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT TMV_RealType(T)

    //
    // Access
    //

    template <class T>
    T GenSymMatrix<T>::cref(ptrdiff_t i, ptrdiff_t j) const
    {
        if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
            const T* mi = cptr() + i*stepi() + j*stepj();
            return isconj() ? TMV_CONJ(*mi) : *mi;
        } else {
            const T* mi = cptr() + j*stepi() + i*stepj();
            return issym() != isconj() ? *mi : TMV_CONJ(*mi);
        }
    }

    template <class T, int A>
    typename SymMatrixView<T,A>::reference SymMatrixView<T,A>::ref(
        ptrdiff_t i, ptrdiff_t j)
    {
        if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
            T* mi = ptr() + i*stepi() + j*stepj();
            return RefHelper<T>::makeRef(mi,ct());
        } else {
            T* mi = ptr() + j*stepi() + i*stepj();
            return RefHelper<T>::makeRef(
                mi, this->issym() != this->isconj() ? NonConj : Conj);
        }
    }

    template <class T>
    void GenSymMatrix<T>::setDiv() const
    {
        if (!this->divIsSet()) {
            DivType dt = this->getDivType();
            TMVAssert(dt == tmv::LU || dt == tmv::CH || dt == tmv::SV);
            TMVAssert(isherm() || dt != tmv::CH);
            switch (dt) {
              case LU : 
                   this->divider.reset(
                       new SymLDLDiv<T>(*this,this->divIsInPlace())); 
                   break;
              case CH : 
                   this->divider.reset(
                       new HermCHDiv<T>(*this,this->divIsInPlace())); 
                   break;
              case SV : 
                   if (isherm()) 
                       this->divider.reset(
                           new HermSVDiv<T>(*this,this->divIsInPlace()));
                   else
                       this->divider.reset(
                           new SymSVDiv<T>(*this,this->divIsInPlace()));
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
    void GenSymMatrix<int>::setDiv() const
    { TMVAssert(TMV_FALSE); }
    template <>
    void GenSymMatrix<std::complex<int> >::setDiv() const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T> 
    QuotXS<T,T> GenSymMatrix<T>::QInverse() const
    { return QuotXS<T,T>(T(1),*this); }

    template <class T> template <class T1> 
    void GenSymMatrix<T>::doMakeInverse(SymMatrixView<T1> sinv) const
    {
        TMVAssert(issym() == sinv.issym());
        TMVAssert(isherm() == sinv.isherm());

        this->setDiv();
        const SymDivider<T>* sdiv = dynamic_cast<const SymDivider<T>*>(
            this->getDiv());
        TMVAssert(sdiv);
        sdiv->makeInverse(sinv);
    }

    template <class T>
    bool GenSymMatrix<T>::divIsLUDiv() const
    { return dynamic_cast<const SymLDLDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenSymMatrix<T>::divIsCHDiv() const
    { return dynamic_cast<const HermCHDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenSymMatrix<T>::divIsHermSVDiv() const
    { return dynamic_cast<const HermSVDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenSymMatrix<T>::divIsSymSVDiv() const
    { return dynamic_cast<const SymSVDiv<T>*>(this->getDiv()); }

#ifdef INST_INT
    template <>
    bool GenSymMatrix<int>::divIsLUDiv() const
    { return false; }
    template <>
    bool GenSymMatrix<int>::divIsCHDiv() const
    { return false; }
    template <>
    bool GenSymMatrix<int>::divIsHermSVDiv() const
    { return false; }
    template <>
    bool GenSymMatrix<int>::divIsSymSVDiv() const
    { return false; }

    template <>
    bool GenSymMatrix<std::complex<int> >::divIsLUDiv() const
    { return false; }
    template <>
    bool GenSymMatrix<std::complex<int> >::divIsCHDiv() const
    { return false; }
    template <>
    bool GenSymMatrix<std::complex<int> >::divIsHermSVDiv() const
    { return false; }
    template <>
    bool GenSymMatrix<std::complex<int> >::divIsSymSVDiv() const
    { return false; }
#endif

    //
    // OK? (SubMatrix, etc.)
    //

    template <class T> 
    bool GenSymMatrix<T>::hasSubMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
    {
        if (i1==i2 || j1==j2) return true; // no elements, so whatever...
        bool ok = true;
        ptrdiff_t i2x = i2-istep;
        ptrdiff_t j2x = j2-jstep;
        if (istep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1 < 0 || i1 >= size()) {
            ok = false;
            std::cerr<<"first col element ("<<i1<<") must be in 0 -- ";
            std::cerr<<size()-1<<std::endl;
        }
        if (i2x < 0 || i2x >= size()) {
            ok = false;
            std::cerr<<"last col element ("<<i2x<<") must be in 0 -- ";
            std::cerr<<size()-1<<std::endl;
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
        if (j1 < 0 || j1 >= size()) {
            ok = false;
            std::cerr<<"first row element ("<<j1<<") must be in 0 -- ";
            std::cerr<<size()-1<<std::endl;
        }
        if (j2-jstep < 0 || j2-jstep >= size()) {
            ok = false;
            std::cerr<<"last row element ("<<j2-jstep<<") must be in 0 -- ";
            std::cerr<<size()-1<<std::endl;
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
        if ((i1<j1 && i2x>j2x) || (i1>j1 && i2x<j2x)) {
            ok = false;
            std::cerr<<"Upper left ("<<i1<<','<<j1<<") and lower right (";
            std::cerr<<i2x<<','<<j2x<<") corners must be in same triangle\n";
        }
        if ((i2x<j1 && i1>j2x) || (i2x>j1 && i1<j2x)) {
            ok = false;
            std::cerr<<"Upper right ("<<i1<<','<<j2x<<") and lower left (";
            std::cerr<<i2x<<','<<j1<<") corners must be in same triangle\n";
        }
        return ok;
    }

    template <class T> 
    bool GenSymMatrix<T>::hasSubVector(
        ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const 
    {
        if (n==0) return true;
        bool ok = true;
        if (istep == 0 && jstep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") and jstep ("<<jstep;
            std::cerr<<") can not both be 0\n";
        }
        if (i<0 || i >= size()) {
            ok = false;
            std::cerr<<"i ("<<i<<") must be in 0 -- "<<size()-1<<std::endl;
        }
        if (j<0 || j >= size()) {
            ok = false;
            std::cerr<<"j ("<<j<<") must be in 0 -- "<<size()-1<<std::endl;
        }
        ptrdiff_t i2 = i+istep*(n-1);
        ptrdiff_t j2 = j+jstep*(n-1);
        if (i2 < 0 || i2 >= size()) {
            ok = false;
            std::cerr<<"last element's i ("<<i2<<") must be in 0 -- ";
            std::cerr<<size()-1<<std::endl;
        }
        if (j2 < 0 || j2 >= size()) {
            ok = false;
            std::cerr<<"last element's j ("<<j2<<") must be in 0 -- ";
            std::cerr<<size()-1<<std::endl;
        }
        if ((i<j && i2>j2) || (i>j && i2<j2)) {
            ok = false;
            std::cerr<<"First ("<<i<<','<<j<<") and last ("<<i2<<','<<j2;
            std::cerr<<") elements must be in same triangle\n";
        }
        return ok;
    }

    template <class T> 
    bool GenSymMatrix<T>::hasSubSymMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const 
    {
        if (i1==i2) return true;
        bool ok=true;
        if (istep == 0) {
            ok = false; 
            std::cerr<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1<0 || i1 >= size()) {
            ok = false;
            std::cerr<<"first diag element ("<<i1<<") must be in 0 -- ";
            std::cerr<<size()-1<<std::endl;
        }
        if (i2-istep<0 || i2-istep >= size()) {
            ok = false;
            std::cerr<<"last diag element ("<<i2-istep<<") must be in 0 -- ";
            std::cerr<<size()-1<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cerr<<"range ("<<i2-i1<<") must be multiple of istep (";
            std::cerr<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cerr<<"n diag elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymMatrixView<T,FortranStyle>::hasSubMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
    {
        if (i1==i2 || j1==j2) return true; // no elements, so whatever...
        bool ok = true;
        if (istep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1 < 1 || i1 > this->size()) {
            ok = false;
            std::cerr<<"first col element ("<<i1<<") must be in 1 -- ";
            std::cerr<<this->size()<<std::endl;
        }
        if (i2 < 1 || i2 > this->size()) {
            ok = false;
            std::cerr<<"last col element ("<<i2-istep<<") must be in 1 -- ";
            std::cerr<<this->size()<<std::endl;
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
        if (j1 < 1 || j1 > this->size()) {
            ok = false;
            std::cerr<<"first row element ("<<j1<<") must be in 1 -- ";
            std::cerr<<this->size()<<std::endl;
        }
        if (j2 < 1 || j2 > this->size()) {
            ok = false;
            std::cerr<<"last row element ("<<j2-jstep<<") must be in 1 -- ";
            std::cerr<<this->size()<<std::endl;
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
        if ((i1<j1 && i2>j2) || (i1>j1 && i2<j2)) {
            ok = false;
            std::cerr<<"Upper left ("<<i1<<','<<j1<<") and lower right (";
            std::cerr<<i2<<','<<j2<<") corners must be in same triangle\n";
        }
        if ((i2<j1 && i1>j2) || (i2>j1 && i1<j2)) {
            ok = false;
            std::cerr<<"Upper right ("<<i1<<','<<j2<<") and lower left (";
            std::cerr<<i2<<','<<j1<<") corners must be in same triangle\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymMatrixView<T,FortranStyle>::hasSubVector(
        ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t n) const 
    {
        if (n==0) return true;
        bool ok = true;
        if (istep == 0 && jstep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") and jstep ("<<jstep;
            std::cerr<<") can not both be 0\n";
        }
        if (i < 1 || i > this->size()) {
            ok = false;
            std::cerr<<"i ("<<i<<") must be in 1 -- "<<this->size()<<std::endl;
        }
        if (i < 1 || j > this->size()) {
            ok = false;
            std::cerr<<"j ("<<j<<") must be in 1 -- "<<this->size()<<std::endl;
        }
        ptrdiff_t i2 = i+istep*(n-1);
        ptrdiff_t j2 = j+jstep*(n-1);
        if (i2 < 1 || i2 > this->size()) {
            ok = false;
            std::cerr<<"last element's i ("<<i2<<") must be in 1 -- ";
            std::cerr<<this->size()<<std::endl;
        }
        if (j2 < 1 || j2 > this->size()) {
            ok = false;
            std::cerr<<"last element's j ("<<j2<<") must be in 1 -- ";
            std::cerr<<this->size()<<std::endl;
        }
        if ((i<j && i2>j2) || (i>j && i2<j2)) {
            ok = false;
            std::cerr<<"First ("<<i<<','<<j<<") and last ("<<i2<<','<<j2;
            std::cerr<<") elements must be in same triangle\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymMatrixView<T,FortranStyle>::hasSubSymMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const 
    {
        if (i1==i2) return true;
        bool ok=true;
        if (istep == 0) {
            ok = false; 
            std::cerr<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1<1 || i1 > this->size()) {
            ok = false;
            std::cerr<<"first diag element ("<<i1<<") must be in 1 -- ";
            std::cerr<<this->size()<<std::endl;
        }
        if (i2-istep<1 || i2-istep > this->size()) {
            ok = false;
            std::cerr<<"last diag element ("<<i2-istep<<") must be in 1 -- ";
            std::cerr<<this->size()<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cerr<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep;
            std::cerr<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cerr<<"n diag elements ("<<(i2-i1)/istep+1<<") must be positive\n";
        }
        return ok;
    }

    // 
    // SwapRowsCols
    //

    template <class T, int A>
    SymMatrixView<T,A>& SymMatrixView<T,A>::swapRowsCols(ptrdiff_t i1, ptrdiff_t i2)
    {
        TMVAssert(i1<size());
        TMVAssert(i2<size());
        if (i1 == i2) return *this;
        else {
#ifdef XDEBUG
            Matrix<T> M(*this);
            Matrix<T> M0(*this);
            M.swapRows(i1,i2);
            M.swapCols(i1,i2);
#endif
            if (i1 > i2) TMV_SWAP(i1,i2);
            // Now i1 < i2
            if (uplo() == Upper) transpose().swapRowsCols(i1,i2);
            else {
                Swap(row(i1,0,i1),row(i2,0,i1));
                Swap(row(i2,i1+1,i2),col(i1,i1+1,i2));
                if (!this->issym()) {
                    row(i2,i1,i2).conjugateSelf(); // Conj m(i2,i1) as well
                    col(i1,i1+1,i2).conjugateSelf();
                }
                Swap(col(i1,i2+1,size()),col(i2,i2+1,size()));
                diag().swap(i1,i2);
            }
#ifdef XDEBUG
            if (Norm(M-*this) > 1.e-5*TMV_MAX(RT(1),Norm(M))) {
                cerr<<"swapRowsCols: i1,i2 = "<<i1<<','<<i2<<endl;
                cerr<<"M0 = "<<TMV_Text(*this)<<"  "<<M0<<endl;
                cerr<<"M = "<<M<<endl;
                cerr<<"S = "<<*this<<endl;
                abort();
            }
#endif
            return *this;
        }
    }

    template <class T, int A>
    SymMatrixView<T,A>& SymMatrixView<T,A>::permuteRowsCols(
        const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2) 
    {
        const ptrdiff_t* pi = p+i1;
        for(ptrdiff_t i=i1;i<i2;++i,++pi) {
            TMVAssert(*pi < size());
            swapRowsCols(i,*pi);
        }
        return *this;
    }

    template <class T, int A>
    SymMatrixView<T,A>& SymMatrixView<T,A>::reversePermuteRowsCols(
        const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
    {
        const ptrdiff_t* pi = p+i2;
        for(ptrdiff_t i=i2;i>i1;) {
            --i; --pi;
            TMVAssert(*pi < size());
            swapRowsCols(i,*pi);
        }
        return *this;
    }


    //
    // Norms
    //

    template <class T>
    T GenSymMatrix<T>::det() const
    { return DivHelper<T>::det(); }

    template <class T>
    RT GenSymMatrix<T>::logDet(T* sign) const
    { return DivHelper<T>::logDet(sign); }

    template <class T>
    bool GenSymMatrix<T>::isSingular() const
    { return DivHelper<T>::isSingular(); }

#ifdef INST_INT
    template <>
    int GenSymMatrix<int>::det() const
    { return IntegerDet(*this); }

    template <>
    std::complex<int> GenSymMatrix<std::complex<int> >::det() const
    { return IntegerDet(*this); }

    template <>
    int GenSymMatrix<int>::logDet(int* ) const
    { TMVAssert(TMV_FALSE); return 0; }

    template <>
    int GenSymMatrix<std::complex<int> >::logDet(std::complex<int>* ) const
    { TMVAssert(TMV_FALSE); return 0; }

    template <>
    bool GenSymMatrix<int>::isSingular() const
    { return det() == 0; }

    template <>
    bool GenSymMatrix<std::complex<int> >::isSingular() const
    { return det() == 0; }
#endif

    template <class T> 
    T GenSymMatrix<T>::sumElements() const
    {
        T sum = diag().sumElements();
        if (size() > 1) {
            T temp = upperTri().offDiag().sumElements();
            if (issym()) {
                sum += RT(2) * temp;
            } else {
                // temp + conj(temp) = 2*real(temp)
                sum += RT(2) * TMV_REAL(temp);
            }
        }
        return sum;
    }

    template <class T> 
    static RT DoSumAbsElements(const GenSymMatrix<T>& m)
    {
        RT sum = m.diag().sumAbsElements();
        if (m.size() > 1) 
            sum += RT(2) * m.upperTri().offDiag().sumAbsElements();
        return sum;
    }

    template <class T> 
    static RT DoSumAbs2Elements(const GenSymMatrix<T>& m)
    {
        RT sum = m.diag().sumAbs2Elements();
        if (m.size() > 1) 
            sum += RT(2) * m.upperTri().offDiag().sumAbs2Elements();
        return sum;
    }

#ifdef INST_INT
    static int DoSumAbsElements(const GenSymMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T> 
    RT GenSymMatrix<T>::sumAbsElements() const
    { return DoSumAbsElements(*this); }

    template <class T> 
    RT GenSymMatrix<T>::sumAbs2Elements() const
    { return DoSumAbs2Elements(*this); }

    template <class T> 
    RT GenSymMatrix<T>::normSq(const RT scale) const
    {
        RT sum = diag().normSq(scale);
        if (size() > 1) sum += RT(2) * upperTri().offDiag().normSq(scale);
        return sum;
    }

    template <class T> 
    static RT NonLapNorm1(const GenSymMatrix<T>& m) 
    {
        RT max(0);
        for(ptrdiff_t j=0;j<m.size();++j) {
            RT temp = m.col(j,0,j).norm1();
            temp += m.col(j,j,m.size()).norm1();
            if (temp > max) max = temp;
        }
        return max;
    } 

    template <class T> 
    static RT NonLapNormF(const GenSymMatrix<T>& m)
    {
        const RT eps = TMV_Epsilon<T>();

        RT mmax = m.maxAbs2Element();
        if (mmax == RT(0)) return RT(0);
        else if (TMV_Underflow(mmax * mmax)) {
            // Then we need to rescale, since underflow has caused 
            // rounding errors.
            // Epsilon is a pure power of 2, so no rounding errors from 
            // rescaling.
            const RT inveps = RT(1)/eps;
            RT scale = inveps;
            mmax *= scale;
            const RT eps2 = eps*eps;
            while (mmax < eps2) { scale *= inveps; mmax *= inveps; }
            return TMV_SQRT(m.normSq(scale))/scale;
        } else if (RT(1) / mmax == RT(0)) {
            // Then mmax is already inf, so no hope of making it more accurate.
            return mmax;
        } else if (RT(1) / (mmax*mmax) == RT(0)) {
            // Then we have overflow, so we need to rescale:
            const RT inveps = RT(1)/eps;
            RT scale = eps;
            mmax *= scale;
            while (mmax > inveps) { scale *= eps; mmax *= eps; }
            return TMV_SQRT(m.normSq(scale))/scale;
        } else {
            return TMV_SQRT(m.normSq());
        }
    }

    template <class T> 
    static inline RT NonLapMaxAbsElement(const GenSymMatrix<T>& m)
    { return m.upperTri().maxAbsElement(); }

    template <class T> 
    static inline RT NonLapMaxAbs2Element(const GenSymMatrix<T>& m)
    { return m.upperTri().maxAbs2Element(); }

#ifdef XLAP
    template <class T> 
    static RT LapNorm(const char c, const GenSymMatrix<T>& m)
    {
        switch(c) {
          case 'M' : return NonLapMaxAbsElement(m);
          case '1' : return NonLapNorm1(m);
          case 'F' : return NonLapNormF(m);
          default : TMVAssert(TMV_FALSE);
        }
        return RT(0);
    }
#ifdef INST_DOUBLE
    template <> 
    double LapNorm(const char c, const GenSymMatrix<double>& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        char cc = c;
        int N = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
        int lwork = c=='1' ? N : 0;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        double norm = LAPNAME(dlansy) (
            LAPCM LAPV(cc),
            (m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
            LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1);
        return norm;
    }
    template <> 
    double LapNorm(const char c, const GenSymMatrix<std::complex<double> >& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        char cc = c;
        int N = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
        int lwork = c=='1' ? N : 0;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        double norm = m.isherm() ?
            LAPNAME(zlanhe) (
                LAPCM LAPV(cc),
                (m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
                LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1) :
            LAPNAME(zlansy) (
                LAPCM LAPV(cc),
                (m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
                LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1);
        return norm;
    }
#endif
#ifdef INST_FLOAT
    template <> 
    float LapNorm(const char c, const GenSymMatrix<float>& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        char cc = c;
        int N = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
        int lwork = c=='1' ? N : 0;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        float norm = LAPNAME(slansy) (
            LAPCM LAPV(cc),
            (m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
            LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1);
        return norm;
    }
    template <> 
    float LapNorm(const char c, const GenSymMatrix<std::complex<float> >& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        char cc = c;
        int N = m.size();
        int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
        int lwork = c=='1' ? N : 0;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        float norm = m.isherm() ?
            LAPNAME(clanhe) (
                LAPCM LAPV(cc),
                (m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
                LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1) :
            LAPNAME(clansy) (
                LAPCM LAPV(cc),
                (m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
                LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1);
        return norm;
    }
#endif
#endif // XLAP

#ifdef INST_INT
    static int NonLapNormF(const GenSymMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNormF(const GenSymMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNorm1(const GenSymMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapMaxAbsElement(const GenSymMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T> 
    RT GenSymMatrix<T>::maxAbsElement() const
    {
#ifdef XLAP
        if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
            return LapNorm('M',*this);
        else
#endif
            return NonLapMaxAbsElement(*this);
    }
    template <class T> 
    RT GenSymMatrix<T>::maxAbs2Element() const
    {
#ifdef XLAP
        if (Traits<T>::iscomplex) return NonLapMaxAbs2Element(*this);
        else if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
            return LapNorm('M',*this);
        else
#endif
            return NonLapMaxAbs2Element(*this);
    }
    template <class T> 
    RT GenSymMatrix<T>::norm1() const
    {
#ifdef XLAP
        if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
            return LapNorm('1',*this);
        else
#endif
            return NonLapNorm1(*this);
    }
    template <class T> 
    RT GenSymMatrix<T>::normF() const
    {
#ifdef XLAP
        if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
            return LapNorm('F',*this);
        else
#endif
            return NonLapNormF(*this);
    }

    template <class T> 
    static RT DoNorm2(const GenSymMatrix<T>& m)
    {
        if (m.size() == 0) return RT(0);
        DiagMatrix<RT> S(m.size());
        if (m.isherm()) {
            HermMatrix<T> m2(m);
            SV_Decompose(m2.view(),S.view());
        } else {
            SymMatrix<T> m2(m);
            SV_Decompose(m2.view(),S.view());
        }
        return S(0);
    }

    template <class T> 
    static RT DoCondition(const GenSymMatrix<T>& m)
    {
        if (m.size() == 0) return RT(1);
        DiagMatrix<RT> S(m.size());
        if (m.isherm()) {
            HermMatrix<T> m2(m);
            SV_Decompose(m2.view(),S.view());
        } else {
            SymMatrix<T> m2(m);
            SV_Decompose(m2.view(),S.view());
        }
        return TMV_ABS(S(0)/S(S.size()-1));
    }

#ifdef INST_INT
    static int DoNorm2(const GenSymMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoCondition(const GenSymMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoNorm2(const GenSymMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoCondition(const GenSymMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T> 
    RT GenSymMatrix<T>::doNorm2() const
    { return tmv::DoNorm2(*this); }
    template <class T> 
    RT GenSymMatrix<T>::doCondition() const
    { return tmv::DoCondition(*this); }


    //
    // I/O
    //

    template <class T> 
    void GenSymMatrix<T>::write(const TMV_Writer& writer) const
    {
        const ptrdiff_t N = size();
        writer.begin();
        writer.writeCode(issym() ? "S" : "H");
        writer.writeSize(N);
        writer.writeSimpleSize(N);
        writer.writeStart();
        for(ptrdiff_t i=0;i<N;++i) {
            writer.writeLParen();
            for(ptrdiff_t j=0;j<i+1;++j) {
                if (j > 0) writer.writeSpace();
                writer.writeValue(cref(i,j));
            }
            if (!writer.isCompact()) {
                for(ptrdiff_t j=i+1;j<N;++j) {
                    writer.writeSpace();
                    writer.writeValue(cref(i,j));
                }
            }
            writer.writeRParen();
            if (i < N-1) writer.writeRowEnd();
        }
        writer.writeFinal();
        writer.end();
    }

#ifndef NOTHROW
    template <class T> 
    class SymMatrixReadError : public ReadError
    {
    public :
        SymMatrix<T> m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t s;
        T v1, v2;
        bool is, iseof, isbad;

        SymMatrixReadError(std::istream& _is) throw() :
            ReadError("SymMatrix."),
            i(0), j(0), exp(0), got(0), s(0), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        SymMatrixReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("SymMatrix."),
            i(0), j(0), exp(_e), got(_g), s(0), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        SymMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenSymMatrix<T>& _m,
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("SymMatrix."),
            m(_m), i(_i), j(_j), exp(_e), got(_g), s(m.size()), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        SymMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenSymMatrix<T>& _m,
            std::istream& _is, T _v1=0, T _v2=0) throw() :
            ReadError("SymMatrix."),
            m(_m), i(_i), j(_j), s(m.size()), v1(_v1), v2(_v2),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        SymMatrixReadError(
            const GenSymMatrix<T>& _m, std::istream& _is, ptrdiff_t _s) throw() :
            ReadError("SymMatrix."),
            m(_m), i(0), j(0), exp(0), got(0), s(_s), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        SymMatrixReadError(const SymMatrixReadError<T>& rhs) :
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got), 
            s(rhs.s), v1(rhs.v1), v2(rhs.v2),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        virtual ~SymMatrixReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for SymMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"'";
                if (isReal(T()) && exp == "S") os<<" (or 'H')";
                os<<", got '"<<got<<"'.\n";
            }
            if (s != m.size()) {
                os<<"Wrong size: expected "<<m.size()<<", got "<<s<<".\n";
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
            if (v1 != v2) {
                os<<"Input matrix is not symmetric.\n";
                os<<"Lower triangle has the value "<<v1<<" at ("<<i<<","<<j<<")\n";
                os<<"Upper triangle has the value "<<v2<<" at ("<<j<<","<<i<<")\n";
            }
            if (m.size() > 0) {
                os<<"The portion of the SymMatrix which was successfully "
                    "read is: \n";
                const ptrdiff_t N = m.size();
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

    template <class T> 
    class HermMatrixReadError : public ReadError
    {
    public :
        HermMatrix<T> m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t s;
        T v1, v2;
        bool is, iseof, isbad;

        HermMatrixReadError(std::istream& _is) throw() :
            ReadError("HermMatrix."),
            i(0), j(0), exp(0), got(0), s(0), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermMatrixReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("HermMatrix."),
            i(0), j(0), exp(_e), got(_g), s(0), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        HermMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenSymMatrix<T>& _m,
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("HermMatrix."),
            m(_m), i(_i), j(_j), exp(_e), got(_g), 
            s(m.size()), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenSymMatrix<T>& _m,
            std::istream& _is, T _v1=0, T _v2=0) throw() :
            ReadError("HermMatrix."),
            m(_m), i(_i), j(_j), exp(0), got(0), 
            s(m.size()), v1(_v1), v2(_v2),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermMatrixReadError(
            const GenSymMatrix<T>& _m, std::istream& _is, ptrdiff_t _s) throw() :
            ReadError("HermMatrix."),
            m(_m), i(0), j(0), exp(0), got(0), s(_s), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        HermMatrixReadError(const HermMatrixReadError<T>& rhs) :
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got), 
            s(rhs.s), v1(rhs.v1), v2(rhs.v2),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        virtual ~HermMatrixReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for HermMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"'";
                if (isReal(T()) && exp == "H") os<<" (or 'S')";
                os<<", got '"<<got<<"'.\n";
            }
            if (s != m.size()) {
                os<<"Wrong size: expected "<<m.size()<<", got "<<s<<".\n";
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
            if (i==j && v1 != T(0)) {
                os<<"Non-real value found on diagonal: "<<v1<<std::endl;
            }
            if (i!=j && v1 != v2) {
                os<<"Input matrix is not Hermitian.\n";
                os<<"Lower triangle has the value "<<v1<<" at ("<<i<<","<<j<<")\n";
                os<<"Upper triangle has the value "<<v2<<" at ("<<j<<","<<i<<")\n";
            }
            if (m.size() > 0) {
                os<<"The portion of the HermMatrix which was successfully "
                    "read is: \n";
                const ptrdiff_t N = m.size();
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
    static void FinishRead(const TMV_Reader& reader, SymMatrixView<T> m) 
    {
        const ptrdiff_t N = m.size();
        std::string exp, got;
        T temp;
        if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            if (m.issym())
                throw SymMatrixReadError<T>(0,0,m,reader.getis(),exp,got);
            else
                throw HermMatrixReadError<T>(0,0,m,reader.getis(),exp,got);
#endif
        }
        for(ptrdiff_t i=0;i<N;++i) {
            if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"SymMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                if (m.issym())
                    throw SymMatrixReadError<T>(i,0,m,reader.getis(),exp,got);
                else
                    throw HermMatrixReadError<T>(i,0,m,reader.getis(),exp,got);
#endif
            }
            for(ptrdiff_t j=0;j<i+1;++j) {
                if (j>0 && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"SymMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    if (m.issym())
                        throw SymMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
                    else
                        throw HermMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                }
                if (!reader.readValue(temp)) {
#ifdef NOTHROW
                    std::cerr<<"SymMatrix Read Error: reading value\n";
                    exit(1);
#else
                    if (m.issym())
                        throw SymMatrixReadError<T>(i,j,m,reader.getis());
                    else
                        throw HermMatrixReadError<T>(i,j,m,reader.getis());
#endif
                }
                if (j == i && !m.issym() && TMV_IMAG(temp) != RT(0)) {
#ifdef NOTHROW
                    std::cerr<<"SymMatrix Read Error: "<<temp<<" not real\n";
                    exit(1);
#else
                    throw HermMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                }
                if (j < i && !reader.isCompact()) {
                    T temp2 = m.issym() ? m.cref(j,i) : TMV_CONJ(m.cref(j,i));
                    if (temp != temp2) {
#ifdef NOTHROW
                        std::cerr<<"SymMatrix Read Error: "<<temp<<" != "<<temp2<<"\n";
                        exit(1);
#else
                        if (m.issym())
                            throw SymMatrixReadError<T>(i,j,m,reader.getis(),temp,m.cref(j,i));
                        else
                            throw HermMatrixReadError<T>(i,j,m,reader.getis(),temp,m.cref(j,i));
#endif
                    }
                } else {
                    m.ref(i,j) = temp;
                }
            }
            if (!reader.isCompact()) {
                for(ptrdiff_t j=i+1;j<N;++j) {
                    if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"SymMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        if (m.issym())
                            throw SymMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
                        else
                            throw HermMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"SymMatrix Read Error: reading value\n";
                        exit(1);
#else
                        if (m.issym())
                            throw SymMatrixReadError<T>(i,j,m,reader.getis());
                        else
                            throw HermMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    m.ref(i,j) = temp;
                }
            }
            if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"SymMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                if (m.issym())
                    throw SymMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
                else
                    throw HermMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
            if (i < N-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"SymMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                if (m.issym())
                    throw SymMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
                else
                    throw HermMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
        }
        if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            if (m.issym())
                throw SymMatrixReadError<T>(N,0,m,reader.getis(),exp,got);
            else
                throw HermMatrixReadError<T>(N,0,m,reader.getis(),exp,got);
#endif
        }
    }

    template <class T, int A>
    void SymMatrix<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        bool ok = isReal(T()) ? 
            reader.readCode("S","H",exp,got) : reader.readCode("S",exp,got);
        if (!ok) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw SymMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: reading size\n";
            exit(1);
#else
            throw SymMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) resize(s);
        s=size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: reading size\n";
            exit(1);
#else
            throw SymMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw SymMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        SymMatrixView<T> v = view();
        FinishRead(reader,v);
    }

    template <class T, int A>
    void HermMatrix<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        bool ok = isReal(T()) ? 
            reader.readCode("S","H",exp,got) : reader.readCode("H",exp,got);
        if (!ok) {
#ifdef NOTHROW
            std::cerr<<"HermMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw HermMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"HermMatrix Read Error: reading size\n";
            exit(1);
#else
            throw HermMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) resize(s);
        s=size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"HermMatrix Read Error: reading size\n";
            exit(1);
#else
            throw HermMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"HermMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw HermMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        SymMatrixView<T> v = view();
        FinishRead(reader,v);
    }

    template <class T, int A>
    void SymMatrixView<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        bool ok = isReal(T()) ? 
            reader.readCode("S","H",exp,got) : 
            reader.readCode(issym()?"S":"H",exp,got);
        if (!ok) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            if (issym())
                throw SymMatrixReadError<T>(reader.getis(),exp,got);
            else
                throw HermMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: reading size\n";
            exit(1);
#else
            if (issym())
                throw SymMatrixReadError<T>(reader.getis(),exp,got);
            else
                throw HermMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: Wrong size\n";
            exit(1);
#else
            if (issym())
                throw SymMatrixReadError<T>(*this,reader.getis(),s);
            else
                throw HermMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        s=size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: reading size\n";
            exit(1);
#else
            if (issym())
                throw SymMatrixReadError<T>(reader.getis(),exp,got);
            else
                throw HermMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"SymMatrix Read Error: Wrong size\n";
            exit(1);
#else
            if (issym())
                throw SymMatrixReadError<T>(*this,reader.getis(),s);
            else
                throw HermMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        FinishRead(reader,*this);
    }



#define InstFile "TMV_SymMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


