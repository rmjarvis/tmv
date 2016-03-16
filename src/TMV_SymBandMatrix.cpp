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
#include "tmv/TMV_SymBandMatrix.h"
#include "tmv/TMV_SymBandCHD.h"
#include "tmv/TMV_BandLUD.h"
#include "tmv/TMV_SymBandSVD.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VIt.h"
#include "tmv/TMV_SymBandMatrixArith.h"
#include "tmv/TMV_SymMatrix.h"
#include "tmv/TMV_BandMatrix.h"
#include "TMV_IntegerDet.h"
#include <iostream>
#include <string>

namespace tmv {

#define RT TMV_RealType(T)

    //
    // Access
    //

    template <class T> 
    T GenSymBandMatrix<T>::cref(ptrdiff_t i, ptrdiff_t j) const
    {
        if (okij(i,j)) {
            if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
                const T* mi = cptr() + i*stepi() + j*stepj();
                return isconj() ? TMV_CONJ(*mi) : *mi;
            } else {
                const T* mi = cptr() + j*stepi() + i*stepj();
                return issym() != isconj() ? *mi : TMV_CONJ(*mi);
            }
        } else {
            return T(0);
        }
    }

    template <class T, int A>
    typename SymBandMatrixView<T,A>::reference SymBandMatrixView<T,A>::ref(
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
    void GenSymBandMatrix<T>::setDiv() const
    {
        if (!this->divIsSet()) {
            DivType dt = this->getDivType();
            TMVAssert(dt == tmv::LU || dt == tmv::CH || dt == tmv::SV);
            TMVAssert(isherm() || dt != tmv::CH);
            switch (dt) {
              case LU : 
                   this->divider.reset(new BandLUDiv<T>(*this)); 
                   break;
              case CH : 
                   this->divider.reset(
                       new HermBandCHDiv<T>(*this,this->divIsInPlace())); 
                   break;
              case SV : 
                   if (isherm()) 
                       this->divider.reset(new HermBandSVDiv<T>(*this));
                   else
                       this->divider.reset(new SymBandSVDiv<T>(*this));
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
    void GenSymBandMatrix<int>::setDiv() const
    { TMVAssert(TMV_FALSE); }
    template <>
    void GenSymBandMatrix<std::complex<int> >::setDiv() const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T>
    bool GenSymBandMatrix<T>::divIsLUDiv() const
    { return dynamic_cast<const BandLUDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenSymBandMatrix<T>::divIsCHDiv() const
    { return dynamic_cast<const HermBandCHDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenSymBandMatrix<T>::divIsHermSVDiv() const
    { return dynamic_cast<const HermBandSVDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenSymBandMatrix<T>::divIsSymSVDiv() const
    { return dynamic_cast<const SymBandSVDiv<T>*>(this->getDiv()); }

#ifdef INST_INT
    template <>
    bool GenSymBandMatrix<int>::divIsLUDiv() const
    { return false; }
    template <>
    bool GenSymBandMatrix<int>::divIsCHDiv() const
    { return false; }
    template <>
    bool GenSymBandMatrix<int>::divIsHermSVDiv() const
    { return false; }
    template <>
    bool GenSymBandMatrix<int>::divIsSymSVDiv() const
    { return false; }

    template <>
    bool GenSymBandMatrix<std::complex<int> >::divIsLUDiv() const
    { return false; }
    template <>
    bool GenSymBandMatrix<std::complex<int> >::divIsCHDiv() const
    { return false; }
    template <>
    bool GenSymBandMatrix<std::complex<int> >::divIsHermSVDiv() const
    { return false; }
    template <>
    bool GenSymBandMatrix<std::complex<int> >::divIsSymSVDiv() const
    { return false; }
#endif

    template <class T> 
    QuotXsB<T,T> GenSymBandMatrix<T>::QInverse() const
    { return QuotXsB<T,T>(T(1),*this); }

    template <class T> template <class T1> 
    void GenSymBandMatrix<T>::doMakeInverse(SymMatrixView<T1> sinv) const
    {
        TMVAssert(issym() == sinv.issym());
        TMVAssert(isherm() == sinv.isherm());

        this->setDiv();
        const SymDivider<T>* sdiv = dynamic_cast<const SymDivider<T>*>(
            this->getDiv());
        TMVAssert(sdiv);
        sdiv->makeInverse(sinv);
    }

    //
    // OK? (SubBandMatrix, etc.)
    //

    template <class T> 
    bool GenSymBandMatrix<T>::hasSubMatrix(
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
            std::cerr<<"row range ("<<j2-j1<<") must be multiple of jstep (";
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
        if (!okij(i1,j2x)) {
            ok = false;
            std::cerr<<"Upper right ("<<i1<<','<<j2x<<") corner must be in band.\n";
        }
        if (!okij(i2x,j1)) {
            ok = false;
            std::cerr<<"Lower left ("<<i2x<<','<<j1<<") corner must be in band.\n";
        }
        return ok;
    }

    template <class T> 
    bool GenSymBandMatrix<T>::hasSubBandMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
        ptrdiff_t istep, ptrdiff_t jstep) const
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
        if (!okij(i1,j1)) {
            ok = false;
            std::cerr<<"Upper left corner ("<<i1<<','<<j1<<") must be in band\n";
        }
        if (!okij(i1,j1+newnhi)) {
            ok = false;
            std::cerr<<"Start of top diagonal ("<<i1<<','<<j1+newnhi;
            std::cerr<<") must be in band\n";
        }
        if (!okij(i1+newnlo,j1)) {
            ok = false;
            std::cerr<<"Start of bottom diagonal ("<<i1+newnlo<<','<<j1;
            std::cerr<<") must be in band\n";
        }
        if (newnhi >= j2-j1) {
            ok = false;
            std::cerr<<"new nhi ("<<newnhi<<") must be less than the new rowsize (";
            std::cerr<<j2-j1<<")\n";
        }
        if (newnlo >= i2-i1) {
            ok = false;
            std::cerr<<"new nlo ("<<newnlo<<") must be less than the new colsize (";
            std::cerr<<i2-i1<<")\n";
        }
        if ((i1<j1+newnhi && i1+newnlo>j1) || (i1>j1+newnhi && i1+newnlo<j1)) {
            ok = false;
            std::cerr<<"Start of top ("<<i1<<','<<j1+newnhi<<") and bottom (";
            std::cerr<<i1+newnlo<<','<<j1<<") diagonals must be in same triangle\n";
        }
        return ok;
    }

    template <class T> 
    bool GenSymBandMatrix<T>::hasSubVector(
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
        if (!okij(i,j)) {
            ok = false;
            std::cerr<<"First ("<<i<<','<<j<<") element must be in band\n";
        }
        if (!okij(i2,j2)) {
            ok = false;
            std::cerr<<"Last ("<<i2<<','<<j2<<") element must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool GenSymBandMatrix<T>::hasSubSymMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t istep) const 
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
        if (!okij(i1,i2-istep)) {
            ok = false;
            std::cerr<<"Upper right ("<<i1<<','<<i2-istep;
            std::cerr<<") corner must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool GenSymBandMatrix<T>::hasSubSymBandMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t newnlo, ptrdiff_t istep) const 
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
        if (newnlo > nlo()) {
            ok = false;
            std::cerr<<"new number of off-diagonals ("<<newnlo<<") must be less ";
            std::cerr<<"than or equal to the current value ("<<nlo()<<")\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymBandMatrixView<T,FortranStyle>::hasSubMatrix(
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
        if (!this->okij(i1-1,j2-1)) {
            ok = false;
            std::cerr<<"Upper right ("<<i1<<','<<j2<<") corner must be in band.\n";
        }
        if (!this->okij(i2-1,j1-1)) {
            ok = false;
            std::cerr<<"Lower left ("<<i2<<','<<j1<<") corner must be in band.\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymBandMatrixView<T,FortranStyle>::hasSubBandMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t newnlo, ptrdiff_t newnhi,
        ptrdiff_t istep, ptrdiff_t jstep) const
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
        if (!this->okij(i1-1,j1-1)) {
            ok = false;
            std::cerr<<"Upper left corner ("<<i1<<','<<j1<<") must be in band\n";
        }
        if (!this->okij(i1-1,j1-1+newnhi)) {
            ok = false;
            std::cerr<<"Start of top diagonal ("<<i1<<','<<j1+newnhi;
            std::cerr<<") must be in band\n";
        }
        if (!this->okij(i1-1+newnlo,j1-1)) {
            ok = false;
            std::cerr<<"Start of bottom diagonal ("<<i1+newnlo<<','<<j1;
            std::cerr<<") must be in band\n";
        }
        if (newnhi >= j2-j1+1) {
            ok = false;
            std::cerr<<"new nhi ("<<newnhi<<") must be less than the new rowsize (";
            std::cerr<<j2-j1+1<<")\n";
        }
        if (newnlo >= i2-i1+1) {
            ok = false;
            std::cerr<<"new nlo ("<<newnlo<<") must be less than the new colsize (";
            std::cerr<<i2-i1+1<<")\n";
        }
        if ((i1<j1+newnhi && i1+newnlo>j1) || (i1>j1+newnhi && i1+newnlo<j1)) {
            ok = false;
            std::cerr<<"Start of top ("<<i1<<','<<j1+newnhi<<") and bottom (";
            std::cerr<<i1+newnlo<<','<<j1<<") diagonals must be in same triangle\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymBandMatrixView<T,FortranStyle>::hasSubVector(
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
        if (!this->okij(i-1,j-1)) {
            ok = false;
            std::cerr<<"First ("<<i<<','<<j<<") element must be in band\n";
        }
        if (!this->okij(i2-1,j2-1)) {
            ok = false;
            std::cerr<<"Last ("<<i2<<','<<j2<<") element must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymBandMatrixView<T,FortranStyle>::hasSubSymMatrix(
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
        if (!this->okij(i1-1,i2-1)) {
            ok = false;
            std::cerr<<"Upper right ("<<i1<<','<<i2<<") corner must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymBandMatrixView<T,FortranStyle>::hasSubSymBandMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t newnlo, ptrdiff_t istep) const 
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
        if (newnlo > this->nlo()) {
            ok = false;
            std::cerr<<"new number of off-diagonals ("<<newnlo<<") must be less ";
            std::cerr<<"than or equal to the current value ("<<this->nlo()<<")\n";
        }
        return ok;
    }


    //
    // Norms
    //

    template <class T>
    T GenSymBandMatrix<T>::det() const
    { return DivHelper<T>::det(); }

    template <class T>
    RT GenSymBandMatrix<T>::logDet(T* sign) const
    { return DivHelper<T>::logDet(sign); }

    template <class T>
    bool GenSymBandMatrix<T>::isSingular() const
    { return DivHelper<T>::isSingular(); }

#ifdef INST_INT
    template <>
    int GenSymBandMatrix<int>::det() const
    { return IntegerDet(*this); }

    template <>
    std::complex<int> GenSymBandMatrix<std::complex<int> >::det() const
    { return IntegerDet(*this); }

    template <>
    int GenSymBandMatrix<int>::logDet(int* ) const
    { TMVAssert(TMV_FALSE); return 0; }

    template <>
    int GenSymBandMatrix<std::complex<int> >::logDet(std::complex<int>* ) const
    { TMVAssert(TMV_FALSE); return 0; }

    template <>
    bool GenSymBandMatrix<int>::isSingular() const
    { return det() == 0; }

    template <>
    bool GenSymBandMatrix<std::complex<int> >::isSingular() const
    { return det() == 0; }
#endif

    template <class T> 
    T GenSymBandMatrix<T>::sumElements() const
    {
        T sum = diag().sumElements();
        if (size() > 1 && nlo() > 0) {
            T temp = upperBandOff().sumElements();
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
    static RT DoSumAbsElements(const GenSymBandMatrix<T>& m)
    {
        RT sum = m.diag().sumAbsElements();
        if (m.size() > 1 && m.nlo() > 0) 
            sum += RT(2) * m.upperBandOff().sumAbsElements();
        return sum;
    }

    template <class T> 
    static RT DoSumAbs2Elements(const GenSymBandMatrix<T>& m)
    {
        RT sum = m.diag().sumAbs2Elements();
        if (m.size() > 1 && m.nlo() > 0) 
            sum += RT(2) * m.upperBandOff().sumAbs2Elements();
        return sum;
    }

#ifdef INST_INT
    static int DoSumAbsElements(const GenSymBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T> 
    RT GenSymBandMatrix<T>::sumAbsElements() const
    { return DoSumAbsElements(*this); }

    template <class T> 
    RT GenSymBandMatrix<T>::sumAbs2Elements() const
    { return DoSumAbs2Elements(*this); }

    template <class T> 
    RT GenSymBandMatrix<T>::normSq(const RT scale) const
    {
        RT sum = diag().normSq(scale);
        if (size() > 1 && nlo() > 0) 
            sum += RT(2) * upperBandOff().normSq(scale);
        return sum;
    }

    template <class T> 
    static RT NonLapNorm1(const GenSymBandMatrix<T>& m) 
    {
        if (m.nlo() > 0) {
            RT max(0);
            const ptrdiff_t N = m.size();
            if (N > 0) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m.nlo()+1;
                ptrdiff_t k=m.nlo();
                for(ptrdiff_t j=0;j<N;++j) {
                    RT temp = m.col(j,i1,j).norm1();
                    temp += m.col(j,j,i2).norm1();
                    if (temp > max) max = temp;
                    if (k>0) --k; else ++i1;
                    if (i2 < N) ++i2;
                }
            }
            return max;
        } else {
            return m.diag().normInf();
        }
    } 

    template <class T> 
    static RT NonLapNormF(const GenSymBandMatrix<T>& m)
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
        }  else {
            return TMV_SQRT(m.normSq());
        }
    }

    template <class T> 
    static inline RT NonLapMaxAbsElement(const GenSymBandMatrix<T>& m)
    { return m.upperBand().maxAbsElement(); }

    template <class T> 
    static inline RT NonLapMaxAbs2Element(const GenSymBandMatrix<T>& m)
    { return m.upperBand().maxAbs2Element(); }

#ifdef XLAP
    template <class T> 
    static RT LapNorm(const char c, const GenSymBandMatrix<T>& m)
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
    double LapNorm(const char c, const GenSymBandMatrix<double>& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        if (m.iscm()) {
            char cc = c;
            int N = m.size();
            int noff = m.nlo();
            int lda = m.diagstep();
#ifndef LAPNOWORK
            int lwork = c=='1' ? N : 0;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            const double* mp = m.cptr();
            if (m.uplo()==Upper) mp -= m.nhi();
            return LAPNAME(dlansb) (
                LAPCM LAPV(cc), (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
                LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
                LAP1 LAP1);
        } else {
            return LapNorm(c,m.transpose());
        }
    }
    template <> 
    double LapNorm(
        const char c, const GenSymBandMatrix<std::complex<double> >& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        if (m.iscm()) {
            char cc = c;
            int N = m.size();
            int noff = m.nlo();
            int lda = m.diagstep();
#ifndef LAPNOWORK
            int lwork = c=='1' ? N : 0;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            const std::complex<double>* mp = m.cptr();
            if (m.uplo()==Upper) mp -= m.nhi();
            if (m.isherm()) 
                return LAPNAME(zlanhb) (
                    LAPCM LAPV(cc), (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
                    LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
                    LAP1 LAP1);
            else 
                return LAPNAME(zlansb) (
                    LAPCM LAPV(cc), (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
                    LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
                    LAP1 LAP1);
        } else {
            return LapNorm(c,m.transpose());
        }
    }
#endif
#ifdef INST_FLOAT
    template <>
    float LapNorm(const char c, const GenSymBandMatrix<float>& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        if (m.iscm()) {
            char cc = c;
            int N = m.size();
            int noff = m.nlo();
            int lda = m.diagstep();
#ifndef LAPNOWORK
            int lwork = c=='1' ? N : 0;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            const float* mp = m.cptr();
            if (m.uplo()==Upper) mp -= m.nhi();
            return LAPNAME(slansb) (
                LAPCM LAPV(cc), (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
                LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
                LAP1 LAP1);
        } else {
            return LapNorm(c,m.transpose());
        }
    }
    template <>
    float LapNorm(
        const char c, const GenSymBandMatrix<std::complex<float> >& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        if (m.iscm()) {
            char cc = c;
            int N = m.size();
            int noff = m.nlo();
            int lda = m.diagstep();
#ifndef LAPNOWORK
            int lwork = c=='1' ? N : 0;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            const std::complex<float>* mp = m.cptr();
            if (m.uplo()==Upper) mp -= m.nhi();
            if (m.isherm()) 
                return LAPNAME(clanhb) (
                    LAPCM LAPV(cc), (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
                    LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
                    LAP1 LAP1);
            else 
                return LAPNAME(clansb) (
                    LAPCM LAPV(cc), (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO,
                    LAPV(N),LAPV(noff),LAPP(mp),LAPV(lda) LAPWK(work.get()) 
                    LAP1 LAP1);
        } else {
            return LapNorm(c,m.transpose());
        }
    }
#endif
#endif // XLAP

#ifdef INST_INT
    static int NonLapNormF(const GenSymBandMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNormF(const GenSymBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNorm1(const GenSymBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapMaxAbsElement(const GenSymBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T> 
    RT GenSymBandMatrix<T>::maxAbsElement() const
    {
#ifdef XLAP
        if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
            return LapNorm('M',*this);
        else
#endif
            return NonLapMaxAbsElement(*this);
    }
    template <class T> 
    RT GenSymBandMatrix<T>::maxAbs2Element() const
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
    RT GenSymBandMatrix<T>::norm1() const
    {
#ifdef XLAP
        if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
            return LapNorm('1',*this);
        else
#endif
            return NonLapNorm1(*this);
    }
    template <class T> 
    RT GenSymBandMatrix<T>::normF() const
    {
#ifdef XLAP
        if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
            return LapNorm('F',*this);
        else
#endif
            return NonLapNormF(*this);
    }

    template <class T> 
    static RT DoNorm2(const GenSymBandMatrix<T>& m)
    {
        if (m.size() == 0) return RT(0);
        DiagMatrix<RT> S(m.size());
        SV_Decompose(m,S.view());
        return TMV_ABS(S(0));
    }

    template <class T> 
    static RT DoCondition(const GenSymBandMatrix<T>& m)
    {
        if (m.size() == 0) return RT(1);
        DiagMatrix<RT> S(m.size());
        SV_Decompose(m,S.view());
        return TMV_ABS(S(0)/S(S.size()-1));
    }

#ifdef INST_INT
    static int DoNorm2(const GenSymBandMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoCondition(const GenSymBandMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoNorm2(const GenSymBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoCondition(const GenSymBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T> 
    RT GenSymBandMatrix<T>::doNorm2() const
    { return tmv::DoNorm2(*this); }
    template <class T> 
    RT GenSymBandMatrix<T>::doCondition() const
    { return tmv::DoCondition(*this); }

    template <class T> 
    SymBandMatrix<T,Upper|DiagMajor> SymTriDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2)
    {
        TMVAssert2(v2.size() == v1.size()-1);
        SymBandMatrix<T,Upper|DiagMajor> temp(v1.size(),1);
        temp.diag() = v1;
        temp.diag(1) = v2;
        return temp;
    }

    template <class T> 
    HermBandMatrix<T,Upper|DiagMajor> HermTriDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2, UpLoType uplo)
    {
        TMVAssert2(v2.size() == v1.size()-1);
        HermBandMatrix<T,Upper|DiagMajor> temp(v1.size(),1);
        temp.diag() = v1;
        if (uplo == Upper)
            temp.diag(1) = v2;
        else
            temp.diag(-1) = v2;
        return temp;
    }

    template <class T> 
    HermBandMatrix<std::complex<T>,Upper|DiagMajor> HermTriDiagMatrix(
        const GenVector<T>& v1, const GenVector<std::complex<T> >& v2,
        UpLoType uplo)
    {
        TMVAssert2(v2.size() == v1.size()-1);
        HermBandMatrix<std::complex<T>,Upper|DiagMajor> temp(v1.size(),1);
        temp.diag() = v1;
        if (uplo == Upper)
            temp.diag(1) = v2;
        else
            temp.diag(-1) = v2;
        return temp;
    }


    //
    // I/O
    //

    template <class T> 
    void GenSymBandMatrix<T>::write(const TMV_Writer& writer) const
    {
        const ptrdiff_t N = rowsize();
        ptrdiff_t j1=0;
        ptrdiff_t j2=nlo()+1;

        writer.begin();
        writer.writeCode(issym()?"sB":"hB");
        writer.writeSize(N);
        writer.writeSimpleSize(N);
        writer.writeFullSize(nlo());
        writer.writeStart();

        for(ptrdiff_t i=0;i<N;++i) {
            writer.writeLParen();
            if (!writer.isCompact()) {
                for(ptrdiff_t j=0;j<j1;++j) {
                    writer.writeValue(T(0));
                    writer.writeSpace();
                }
            }
            for(ptrdiff_t j=j1;j<i+1;++j) {
                if (j > j1) writer.writeSpace();
                writer.writeValue(cref(i,j));
            }
            if (!writer.isCompact()) {
                for(ptrdiff_t j=i+1;j<j2;++j) {
                    writer.writeSpace();
                    writer.writeValue(cref(i,j));
                }
                for(ptrdiff_t j=j2;j<N;++j) {
                    writer.writeSpace();
                    writer.writeValue(T(0));
                }
                if (j2 < N) ++j2;
            }
            writer.writeRParen();
            if (i < N-1) writer.writeRowEnd();
            if (i >= nlo()) ++j1;
        }
        writer.writeFinal();
        writer.end();
    }

#ifndef NOTHROW
    template <class T> 
    class SymBandMatrixReadError : public ReadError
    {
    public :
        SymBandMatrix<T> m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t s;
        ptrdiff_t lo;
        T v1, v2;
        bool is, iseof, isbad;

        SymBandMatrixReadError(std::istream& _is) throw() :
            ReadError("SymBandMatrix."),
            i(0), j(0), s(0), lo(0), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        SymBandMatrixReadError(
            std::istream& _is, std::string _e, std::string _g) throw() :
            ReadError("SymBandMatrix."),
            i(0), j(0), exp(_e), got(_g), s(0), lo(0), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        SymBandMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenSymBandMatrix<T>& _m,
            std::istream& _is, std::string _e, std::string _g) throw() :
            ReadError("SymBandMatrix."),
            m(_m), i(_i), j(_j), exp(_e), got(_g),
            s(m.size()), lo(m.nlo()), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        SymBandMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenSymBandMatrix<T>& _m,
            std::istream& _is, T _v1=0, T _v2=0) throw() :
            ReadError("SymBandMatrix."),
            m(_m), i(_i), j(_j), s(m.size()), lo(m.nlo()), v1(_v1), v2(_v2),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        SymBandMatrixReadError(
            const GenSymBandMatrix<T>& _m, std::istream& _is, 
            ptrdiff_t _s, ptrdiff_t _lo) throw() :
            ReadError("SymBandMatrix."),
            m(_m), i(0), j(0), s(_s), lo(_lo), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        SymBandMatrixReadError(const SymBandMatrixReadError<T>& rhs) :
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got), 
            s(rhs.s), lo(rhs.lo), v1(rhs.v1), v2(rhs.v2),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}

        virtual ~SymBandMatrixReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for SymBandMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"'";
                if (isReal(T()) && exp == "sB") os<<" (or 'hB')";
                os<<", got '"<<got<<"'.\n";
            }
            if (s != m.size()) {
                os<<"Wrong size: expected "<<m.size()<<", got "<<s<<".\n";
            }
            if (lo != m.nlo()) {
                os<<"Write nlo: expected "<<m.nlo()<<", got "<<lo<<".\n";
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
            if (std::abs(i-j) > m.nlo() && v1 != T(0)) {
                os<<"Invalid input.  Expected 0, got "<<v1<<".\n";
            }
            if (std::abs(i-j) <= m.nlo() && v1 != v2) {
                os<<"Input matrix is not symmetric.\n";
                os<<"Lower triangle has the value "<<v1<<" at ("<<i<<","<<j<<")\n";
                os<<"Upper triangle has the value "<<v2<<" at ("<<j<<","<<i<<")\n";
            }
            if (m.size() > 0) {
                os<<"The portion of the SymBandMatrix which was successfully read is: \n";
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
    class HermBandMatrixReadError : public ReadError
    {
    public :
        HermBandMatrix<T> m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t s;
        ptrdiff_t lo;
        T v1, v2;
        bool is, iseof, isbad;

        HermBandMatrixReadError(std::istream& _is) throw() :
            ReadError("HermBandMatrix."),
            i(0), j(0), s(0), lo(0), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermBandMatrixReadError(
            std::istream& _is, std::string _e, std::string _g) throw() :
            ReadError("HermBandMatrix."),
            i(0), j(0), exp(_e), got(_g), s(0), lo(0), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        HermBandMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenSymBandMatrix<T>& _m,
            std::istream& _is, std::string _e, std::string _g) throw() :
            ReadError("HermBandMatrix."),
            m(_m), i(_i), j(_j), exp(_e), got(_g),
            s(m.size()), lo(m.nlo()), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermBandMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenSymBandMatrix<T>& _m,
            std::istream& _is, T _v1=0, T _v2=0) throw() :
            ReadError("HermBandMatrix."),
            m(_m), i(_i), j(_j), s(m.size()), lo(m.nlo()), v1(_v1), v2(_v2),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermBandMatrixReadError(
            const GenSymBandMatrix<T>& _m, std::istream& _is, 
            ptrdiff_t _s, ptrdiff_t _lo) throw() :
            ReadError("HermBandMatrix."),
            m(_m), i(0), j(0), s(_s), lo(_lo), v1(0), v2(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        HermBandMatrixReadError(const HermBandMatrixReadError<T>& rhs) :
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got),
            s(rhs.s), lo(rhs.lo), v1(rhs.v1), v2(rhs.v2),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}

        virtual ~HermBandMatrixReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for HermBandMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"'";
                if (isReal(T()) && exp == "hB") os<<" (or 'sB')";
                os<<", got '"<<got<<"'.\n";
            }
            if (s != m.size()) {
                os<<"Wrong size: expected "<<m.size()<<", got "<<s<<".\n";
            }
            if (lo != m.nlo()) {
                os<<"Wrong nlo: expected "<<m.nlo()<<", got "<<lo<<".\n";
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
            if (std::abs(i-j) > m.nlo() && v1 != T(0)) {
                os<<"Invalid input.  Expected 0, got "<<v1<<".\n";
            }
            if (i==j && v1 != T(0)) {
                os<<"Non-real value found on diagonal: "<<v1<<std::endl;
            }
            if (std::abs(i-j) <= m.nlo() && i!=j && v1 != v1) {
                os<<"Input matrix is not symmetric.\n";
                os<<"Lower triangle has the value "<<v1<<" at ("<<i<<","<<j<<")\n";
                os<<"Upper triangle has the value "<<v2<<" at ("<<j<<","<<i<<")\n";
            }
            if (m.size() > 0) {
                os<<"The portion of the HermBandMatrix which was successfully read is: \n";
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
    static void FinishRead(const TMV_Reader& reader, SymBandMatrixView<T> m) 
    {
        const ptrdiff_t N = m.size();
        ptrdiff_t j1=0;
        ptrdiff_t j2=m.nlo()+1;
        std::string exp, got;
        T temp;
        if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            if (m.issym())
                throw SymBandMatrixReadError<T>(0,0,m,reader.getis(),exp,got);
            else
                throw HermBandMatrixReadError<T>(0,0,m,reader.getis(),exp,got);
#endif
        }
        for(ptrdiff_t i=0;i<N;++i) {
            if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                if (m.issym())
                    throw SymBandMatrixReadError<T>(i,0,m,reader.getis(),exp,got);
                else
                    throw HermBandMatrixReadError<T>(i,0,m,reader.getis(),exp,got);
#endif
            }
            if (!reader.isCompact()) {
                for(ptrdiff_t j=0;j<j1;++j) {
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix Read Error: reading value\n";
                        exit(1);
#else
                        if (m.issym())
                            throw SymBandMatrixReadError<T>(i,j,m,reader.getis());
                        else
                            throw HermBandMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    if (temp != T(0)) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix Read Error: "<<temp<<" != 0\n";
                        exit(1);
#else
                        if (m.issym())
                            throw SymBandMatrixReadError<T>(i,j,m,reader.getis(),temp);
                        else
                            throw HermBandMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                    }
                    if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        if (m.issym())
                            throw SymBandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
                        else
                            throw HermBandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                }
            }
            for(ptrdiff_t j=j1;j<i+1;++j) {
                if (j>j1 && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    if (m.issym())
                        throw SymBandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
                    else
                        throw HermBandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                }
                if (!reader.readValue(temp)) {
#ifdef NOTHROW
                    std::cerr<<"SymBandMatrix Read Error: reading value\n";
                    exit(1);
#else
                    if (m.issym())
                        throw SymBandMatrixReadError<T>(i,j,m,reader.getis());
                    else
                        throw HermBandMatrixReadError<T>(i,j,m,reader.getis());
#endif
                }
                if (j == i && !m.issym() && TMV_IMAG(temp) != RT(0)) {
#ifdef NOTHROW
                    std::cerr<<"SymBandMatrix Read Error: "<<temp<<" not real\n";
                    exit(1);
#else
                    throw HermBandMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                }
                if (j < i && !reader.isCompact()) {
                    T temp2 = m.issym() ? m.cref(j,i) : TMV_CONJ(m.cref(j,i));
                    if (temp != temp2) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix Read Error: "<<temp<<" != "<<temp2<<"\n";
                        exit(1);
#else
                        if (m.issym())
                            throw SymBandMatrixReadError<T>(i,j,m,reader.getis(),temp,m.cref(j,i));
                        else
                            throw HermBandMatrixReadError<T>(i,j,m,reader.getis(),temp,m.cref(j,i));
#endif
                    }
                } else {
                    m.ref(i,j) = temp;
                }
            }
            if (!reader.isCompact()) {
                for(ptrdiff_t j=i+1;j<j2;++j) {
                    if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        if (m.issym())
                            throw SymBandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
                        else
                            throw HermBandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix Read Error: reading value\n";
                        exit(1);
#else
                        if (m.issym())
                            throw SymBandMatrixReadError<T>(i,j,m,reader.getis());
                        else
                            throw HermBandMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    m.ref(i,j) = temp;
                }
                for(ptrdiff_t j=j2;j<N;++j) {
                    if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        if (m.issym())
                            throw SymBandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
                        else
                            throw HermBandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix Read Error: reading value\n";
                        exit(1);
#else
                        if (m.issym())
                            throw SymBandMatrixReadError<T>(i,j,m,reader.getis());
                        else
                            throw HermBandMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    if (temp != T(0)) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix Read Error: "<<temp<<" != 0\n";
                        exit(1);
#else
                        if (m.issym())
                            throw SymBandMatrixReadError<T>(i,j,m,reader.getis(),temp);
                        else
                            throw HermBandMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                    }
                }
                if (j2 < N) ++j2;
            }
            if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                if (m.issym())
                    throw SymBandMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
                else
                    throw HermBandMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
            if (i < N-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                if (m.issym())
                    throw SymBandMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
                else
                    throw HermBandMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
            if (i >= m.nlo()) ++j1;
        }
        if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            if (m.issym())
                throw SymBandMatrixReadError<T>(N,0,m,reader.getis(),exp,got);
            else
                throw HermBandMatrixReadError<T>(N,0,m,reader.getis(),exp,got);
#endif
        }
    }

    template <class T, int A>
    void SymBandMatrix<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        bool ok = isReal(T()) ? 
            reader.readCode("sB","hB",exp,got) : reader.readCode("sB",exp,got);
        if (!ok) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw SymBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size(), lo=nlo();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: reading size\n";
            exit(1);
#else
            throw SymBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s2=s;
        if (!reader.readSimpleSize(s2,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: reading size\n";
            exit(1);
#else
            throw SymBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s2 != s) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw SymBandMatrixReadError<T>(*this,reader.getis(),s,lo);
#endif
        }
        if (!reader.readFullSize(lo,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: reading size\n";
            exit(1);
#else
            throw SymBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size() || lo != nlo()) resize(s,lo);
        SymBandMatrixView<T> v = view();
        FinishRead(reader,v);
    }

    template <class T, int A>
    void HermBandMatrix<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        bool ok = isReal(T()) ? 
            reader.readCode("sB","hB",exp,got) : reader.readCode("hB",exp,got);
        if (!ok) {
#ifdef NOTHROW
            std::cerr<<"HermBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw HermBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size(), lo=nlo();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"HermBandMatrix Read Error: reading size\n";
            exit(1);
#else
            throw HermBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s2=s;
        if (!reader.readSimpleSize(s2,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"HermBandMatrix Read Error: reading size\n";
            exit(1);
#else
            throw HermBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s2 != s) {
#ifdef NOTHROW
            std::cerr<<"HermBandMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw HermBandMatrixReadError<T>(*this,reader.getis(),s,lo);
#endif
        }
        if (!reader.readFullSize(lo,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"HermBandMatrix Read Error: reading size\n";
            exit(1);
#else
            throw HermBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size() || lo != nlo()) resize(s,lo);
        SymBandMatrixView<T> v = view();
        FinishRead(reader,v);
    }

    template <class T, int A>
    void SymBandMatrixView<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        bool ok = isReal(T()) ? 
            reader.readCode("sB","hB",exp,got) : 
            reader.readCode(issym()?"sB":"hB",exp,got);
        if (!ok) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            if (issym())
                throw SymBandMatrixReadError<T>(reader.getis(),exp,got);
            else
                throw HermBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size(), lo=nlo();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: reading size\n";
            exit(1);
#else
            if (issym())
                throw SymBandMatrixReadError<T>(reader.getis(),exp,got);
            else
                throw HermBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: Wrong size\n";
            exit(1);
#else
            if (issym())
                throw SymBandMatrixReadError<T>(*this,reader.getis(),s,lo);
            else
                throw HermBandMatrixReadError<T>(*this,reader.getis(),s,lo);
#endif
        }
        s=size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: reading size\n";
            exit(1);
#else
            if (issym())
                throw SymBandMatrixReadError<T>(reader.getis(),exp,got);
            else
                throw HermBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: Wrong size\n";
            exit(1);
#else
            if (issym())
                throw SymBandMatrixReadError<T>(*this,reader.getis(),s,lo);
            else
                throw HermBandMatrixReadError<T>(*this,reader.getis(),s,lo);
#endif
        }
        if (!reader.readFullSize(lo,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: reading size\n";
            exit(1);
#else
            if (issym())
                throw SymBandMatrixReadError<T>(reader.getis(),exp,got);
            else
                throw HermBandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (lo != nlo()) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix Read Error: Wrong nlo\n";
            exit(1);
#else
            if (issym())
                throw SymBandMatrixReadError<T>(*this,reader.getis(),s,lo);
            else
                throw HermBandMatrixReadError<T>(*this,reader.getis(),s,lo);
#endif
        }
        FinishRead(reader,*this);
    }

#undef RT

#define InstFile "TMV_SymBandMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


