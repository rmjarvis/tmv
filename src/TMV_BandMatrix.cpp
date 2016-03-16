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
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_BandLUD.h"
#include "tmv/TMV_BandQRD.h"
#include "tmv/TMV_BandSVD.h"
#include "tmv/TMV_VIt.h"
#include "tmv/TMV_BandMatrixArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include "TMV_IntegerDet.h"
#include <iostream>

namespace tmv {

#define RT TMV_RealType(T)

    // 
    // Access
    //

    template <class T> 
    T GenBandMatrix<T>::cref(ptrdiff_t i, ptrdiff_t j) const
    {
        if (okij(i,j)) {
            const T* mi = cptr()+i*stepi()+j*stepj();
            return isconj() ? TMV_CONJ(*mi) : *mi;
        } else {
            return T(0);
        }
    }

    template <class T, int A>
    typename BandMatrixView<T,A>::reference BandMatrixView<T,A>::ref(
        ptrdiff_t i, ptrdiff_t j) 
    {
        T* mi = ptr()+i*stepi()+j*stepj();
        return RefHelper<T>::makeRef(mi,ct());
    }

    template <class T> 
    void GenBandMatrix<T>::setDiv() const
    {
        if (!this->divIsSet()) {
            DivType dt = this->getDivType();
            TMVAssert(dt == tmv::LU || dt == tmv::QR || dt == tmv::SV);
            switch (dt) {
              case LU : 
                   this->divider.reset(
                       new BandLUDiv<T>(*this,this->divIsInPlace())); break;
              case QR : 
                   this->divider.reset(
                       new BandQRDiv<T>(*this,this->divIsInPlace())); break;
              case SV : 
                   this->divider.reset(new BandSVDiv<T>(*this)); break;
              default : 
                   // The above assert should have already failed
                   // so go ahead and fall through.
                   break;
            }
        }
    }

#ifdef INST_INT
    template <>
    void GenBandMatrix<int>::setDiv() const
    { TMVAssert(TMV_FALSE); }
    template <>
    void GenBandMatrix<std::complex<int> >::setDiv() const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T>
    bool GenBandMatrix<T>::divIsLUDiv() const
    { return dynamic_cast<const BandLUDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenBandMatrix<T>::divIsQRDiv() const
    { return dynamic_cast<const BandQRDiv<T>*>(this->getDiv()); }

    template <class T>
    bool GenBandMatrix<T>::divIsSVDiv() const
    { return dynamic_cast<const BandSVDiv<T>*>(this->getDiv()); }

#ifdef INST_INT
    template <>
    bool GenBandMatrix<int>::divIsLUDiv() const
    { return false; }
    template <>
    bool GenBandMatrix<int>::divIsQRDiv() const
    { return false; }
    template <>
    bool GenBandMatrix<int>::divIsSVDiv() const
    { return false; }

    template <>
    bool GenBandMatrix<std::complex<int> >::divIsLUDiv() const
    { return false; }
    template <>
    bool GenBandMatrix<std::complex<int> >::divIsQRDiv() const
    { return false; }
    template <>
    bool GenBandMatrix<std::complex<int> >::divIsSVDiv() const
    { return false; }
#endif

    ptrdiff_t BandStorageLength(StorageType s, ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi)
    {
        TMVAssert(s == RowMajor || s == ColMajor || s == DiagMajor);
        if (cs == 0 || rs == 0) return 0;
        else if (cs == rs) return (cs-1)*(lo+hi)+cs;
        else {
            // correct cs, rs to be actual end of data
            if (cs > rs+lo) cs = rs+lo;
            if (rs > cs+hi) rs = cs+hi;

            if (s == RowMajor) 
                // si = lo+hi, sj = 1, size = (cs-1)*si + (rs-1)*sj + 1
                return (cs-1)*(lo+hi) + rs;
            else if (s == ColMajor) 
                // si = 1, sj = lo+hi, size = (cs-1)*si + (rs-1)*sj + 1
                return (rs-1)*(lo+hi) + cs;
            else if (cs > rs) 
                // si = -rs, sj = 1+rs, size = (rs-lo-hi-1)*si + (rs-1)*sj + 1
                // size = (rs-lo-hi-1)(-rs) + (rs-1)(1+rs) + 1
                //      = rs(lo+hi+1)
                return rs*(lo+hi+1);
            else
                // si = 1-cs, sj = cs, size = (rs-lo-hi-1)*si + (rs-1)*sj + 1
                // size = (cs-lo-hi-1+rs-cs)(1-cs) + (rs-1)(cs) + 1
                //      = cs(lo+hi+1-rs+rs-1) + rs-lo-hi-1 + 1
                //      = (cs-1)(lo+hi) + rs
                return (cs-1)*(lo+hi) + rs;
        }
    }

    ptrdiff_t BandNumElements(ptrdiff_t cs, ptrdiff_t rs, ptrdiff_t lo, ptrdiff_t hi)
    {
        if (cs == 0 || rs == 0) return 0;
        else if (cs == rs) {
            return cs*(lo+hi+1) - (lo*(lo+1)/2) - (hi*(hi+1)/2);
        } else if (cs > rs) {
            // lox = number of subdiagonals that are clipped.
            ptrdiff_t lox = TMV_MAX(rs+lo-cs,ptrdiff_t(0));
            return rs*(lo+hi+1) - (lox*(lox+1)/2) - (hi*(hi+1)/2);
        } else {
            // hix = number of superdiagonals that are clipped.
            ptrdiff_t hix = TMV_MAX(cs+hi-rs,ptrdiff_t(0));
            return cs*(lo+hi+1) - (lo*(lo+1)/2) - (hix*(hix+1)/2);
        }
    }

    //
    // OK? (SubMatrix, etc.)
    //

    template <class T> 
    bool GenBandMatrix<T>::hasSubMatrix(
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
            std::cerr<<"row range ("<<j2-j1<<") must be multiple of jstep (";
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
        if (!okij(i1,j2-jstep)) {
            ok = false;
            std::cerr<<"Upper right corner ("<<i1<<','<<j2-jstep;
            std::cerr<<") must be in band\n";
        }
        if (!okij(i2-istep,j1)) {
            ok = false;
            std::cerr<<"Lower left corner ("<<i2-istep<<','<<j1;
            std::cerr<<") must be in band\n";
        }
        if (!okij(i2-istep,j2-jstep)) {
            ok = false;
            std::cerr<<"Lower right corner ("<<i2-istep<<','<<j2-jstep;
            std::cerr<<") must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool GenBandMatrix<T>::hasSubVector(
        ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const
    {
        if (size==0) return true;
        bool ok = true;
        if (istep == 0 && jstep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") and jstep ("<<jstep;
            std::cerr<<") can not both be 0\n";
        }
        if (i<0 || i >= colsize()) {
            ok = false;
            std::cerr<<"i ("<<i<<") must be in 0 -- "<<colsize()-1<<std::endl;
        }
        if (j<0 || j >= rowsize()) {
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
        if (!okij(i,j)) {
            ok = false;
            std::cerr<<"First element ("<<i<<','<<j<<") must be in band\n";
        }
        if (!okij(i2,j2)) {
            ok = false;
            std::cerr<<"Last element ("<<i2<<','<<j2<<") must be in band\n";
        }
        return ok;
    } 
    template <class T> 
    bool GenBandMatrix<T>::hasSubBandMatrix(
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
        return ok;
    } 

    template <class T> 
    bool ConstBandMatrixView<T,FortranStyle>::hasSubMatrix(
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
        if (j1 < 0 || j1 >= this->rowsize()) {
            ok = false;
            std::cerr<<"first row element ("<<j1<<") must be in 1 -- ";
            std::cerr<<this->rowsize()<<std::endl;
        }
        if (j2 < 0 || j2 >= this->rowsize()) {
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
        if (!this->okij(i1-1,j2-1)) {
            ok = false;
            std::cerr<<"Upper right corner ("<<i1<<','<<j2<<") must be in band\n";
        }
        if (!this->okij(i2-1,j1-1)) {
            ok = false;
            std::cerr<<"Lower left corner ("<<i2<<','<<j1<<") must be in band\n";
        }
        if (!this->okij(i2-1,j2-1)) {
            ok = false;
            std::cerr<<"Lower right corner ("<<i2<<','<<j2<<") must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstBandMatrixView<T,FortranStyle>::hasSubVector(
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
        if (!this->okij(i-1,j-1)) {
            ok = false;
            std::cerr<<"First element ("<<i<<','<<j<<") must be in band\n";
        }
        if (!this->okij(i2-1,j2-1)) {
            ok = false;
            std::cerr<<"Last element ("<<i2<<','<<j2<<") must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstBandMatrixView<T,FortranStyle>::hasSubBandMatrix(
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
        return ok;
    } 

    template <class T, int A>
    bool ConstBandMatrixView<T,A>::canLinearize() const
    {
        if (linsize == -1) {
            ptrdiff_t rs = this->rowsize();
            ptrdiff_t cs = this->colsize();
            ptrdiff_t lo = this->nlo();
            ptrdiff_t hi = this->nhi();
            if (rs > cs+this->nhi()) rs = cs+this->nhi();
            if (cs > rs+this->nlo()) cs = rs+this->nlo();

            if (rs == 0 || cs == 0) linsize = 0;
            else if (stepi() == 1 && stepj() == (lo+hi))  
                linsize = cs + (rs-1)*(lo+hi);
            else if (stepj() == 1 && stepi() == (lo+hi))
                linsize = rs + (cs-1)*(lo+hi);
        }
        return linsize > 0;
    }

    template <class T, int A>
    bool BandMatrixView<T,A>::canLinearize() const
    {
        if (linsize == -1) {
            ptrdiff_t rs = this->rowsize();
            ptrdiff_t cs = this->colsize();
            ptrdiff_t lo = this->nlo();
            ptrdiff_t hi = this->nhi();
            if (rs > cs+this->nhi()) rs = cs+this->nhi();
            if (cs > rs+this->nlo()) cs = rs+this->nlo();

            if (rs == 0 || cs == 0) linsize = 0;
            else if (stepi() == 1 && stepj() == (lo+hi))  
                linsize = cs + (rs-1)*(lo+hi);
            else if (stepj() == 1 && stepi() == (lo+hi))
                linsize = rs + (cs-1)*(lo+hi);
        }
        return linsize > 0;
    }

    template <class T> QuotXB<T,T> GenBandMatrix<T>::QInverse() const
    { return QuotXB<T,T>(T(1),*this); }

    //
    // Norms
    //
    
    template <class T>
    T GenBandMatrix<T>::det() const
    { return DivHelper<T>::det(); }

    template <class T>
    RT GenBandMatrix<T>::logDet(T* sign) const
    { return DivHelper<T>::logDet(sign); }

    template <class T>
    bool GenBandMatrix<T>::isSingular() const
    { return DivHelper<T>::isSingular(); }

#ifdef INST_INT
    template <>
    int GenBandMatrix<int>::det() const
    { return IntegerDet(*this); }

    template <>
    std::complex<int> GenBandMatrix<std::complex<int> >::det() const
    { return IntegerDet(*this); }

    template <>
    int GenBandMatrix<int>::logDet(int* ) const
    { TMVAssert(TMV_FALSE); return 0; }

    template <>
    int GenBandMatrix<std::complex<int> >::logDet(std::complex<int>* ) const
    { TMVAssert(TMV_FALSE); return 0; }

    template <>
    bool GenBandMatrix<int>::isSingular() const
    { return det() == 0; }

    template <>
    bool GenBandMatrix<std::complex<int> >::isSingular() const
    { return det() == 0; }
#endif

    template <class T>
    T GenBandMatrix<T>::sumElements() const
    {
        const ptrdiff_t M = colsize();
        const ptrdiff_t N = rowsize();
        T sum = 0;
        if (M > 0 && N > 0) {
            if (isrm()) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=nhi()+1;
                ptrdiff_t k=nlo();
                for(ptrdiff_t i=0;i<M;++i) {
                    sum += row(i,j1,j2).sumElements();
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (iscm()) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=nlo()+1;
                ptrdiff_t k=nhi();
                for(ptrdiff_t j=0;j<N;++j) {
                    sum += col(j,i1,i2).sumElements();
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(ptrdiff_t i=-nlo();i<=nhi();++i) sum += diag(i).sumElements();
            }
        }
        return sum;
    }

    template <class T>
    static RT DoSumAbsElements(const GenBandMatrix<T>& m)
    {
        const ptrdiff_t M = m.colsize();
        const ptrdiff_t N = m.rowsize();
        RT sum = 0;
        if (M > 0 && N > 0) {
            if (m.isrm()) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=m.nhi()+1;
                ptrdiff_t k=m.nlo();
                for(ptrdiff_t i=0;i<M;++i) {
                    sum += m.row(i,j1,j2).sumAbsElements();
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (m.iscm()) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m.nlo()+1;
                ptrdiff_t k=m.nhi();
                for(ptrdiff_t j=0;j<N;++j) {
                    sum += m.col(j,i1,i2).sumAbsElements();
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(ptrdiff_t i=-m.nlo();i<=m.nhi();++i) 
                    sum += m.diag(i).sumAbsElements();
            }
        }
        return sum;
    }

    template <class T>
    static RT DoSumAbs2Elements(const GenBandMatrix<T>& m)
    {
        const ptrdiff_t M = m.colsize();
        const ptrdiff_t N = m.rowsize();
        RT sum = 0;
        if (M > 0 && N > 0) {
            if (m.isrm()) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=m.nhi()+1;
                ptrdiff_t k=m.nlo();
                for(ptrdiff_t i=0;i<M;++i) {
                    sum += m.row(i,j1,j2).sumAbs2Elements();
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (m.iscm()) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m.nlo()+1;
                ptrdiff_t k=m.nhi();
                for(ptrdiff_t j=0;j<N;++j) {
                    sum += m.col(j,i1,i2).sumAbs2Elements();
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(ptrdiff_t i=-m.nlo();i<=m.nhi();++i) 
                    sum += m.diag(i).sumAbs2Elements();
            }
        }
        return sum;
    }

#ifdef INST_INT
    static int DoSumAbsElements(const GenBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T>
    RT GenBandMatrix<T>::sumAbsElements() const
    { return DoSumAbsElements(*this); }
 
    template <class T>
    RT GenBandMatrix<T>::sumAbs2Elements() const
    { return DoSumAbs2Elements(*this); }
 
    template <class T> 
    RT GenBandMatrix<T>::normSq(RT scale) const
    {
        const ptrdiff_t M = colsize();
        const ptrdiff_t N = rowsize();
        RT sum = 0;
        if (M > 0 && N > 0) {
            if (isrm()) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=nhi()+1;
                ptrdiff_t k=nlo();
                for(ptrdiff_t i=0;i<M;++i) {
                    sum += row(i,j1,j2).normSq(scale);
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (iscm()) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=nlo()+1;
                ptrdiff_t k=nhi();
                for(ptrdiff_t j=0;j<N;++j) {
                    sum += col(j,i1,i2).normSq(scale);
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(ptrdiff_t i=-nlo();i<=nhi();++i) sum += diag(i).normSq(scale);
            }
        }
        return sum;
    }

    template <class T> 
    static RT NonLapMaxAbsElement(const GenBandMatrix<T>& m)
    {
        RT max = 0;
        const ptrdiff_t M = m.colsize();
        const ptrdiff_t N = m.rowsize();
        if (M > 0 && N > 0) {
            if (m.isrm()) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=m.nhi()+1;
                ptrdiff_t k=m.nlo();
                for(ptrdiff_t i=0;i<M;++i) {
                    RT temp = m.row(i,j1,j2).maxAbsElement();
                    if (temp > max) max = temp;
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (m.iscm()) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m.nlo()+1;
                ptrdiff_t k=m.nhi();
                for(ptrdiff_t j=0;j<N;++j) {
                    RT temp = m.col(j,i1,i2).maxAbsElement();
                    if (temp > max) max = temp;
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(ptrdiff_t i=-m.nlo();i<=m.nhi();++i) {
                    RT temp = m.diag(i).maxAbsElement();
                    if (temp > max) max = temp;
                }
            }
        }
        return max;
    }

    template <class T> 
    static RT NonLapMaxAbs2Element(const GenBandMatrix<T>& m)
    {
        RT max = 0;
        const ptrdiff_t M = m.colsize();
        const ptrdiff_t N = m.rowsize();
        if (M > 0 && N > 0) {
            if (m.isrm()) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=m.nhi()+1;
                ptrdiff_t k=m.nlo();
                for(ptrdiff_t i=0;i<M;++i) {
                    RT temp = m.row(i,j1,j2).maxAbs2Element();
                    if (temp > max) max = temp;
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (m.iscm()) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m.nlo()+1;
                ptrdiff_t k=m.nhi();
                for(ptrdiff_t j=0;j<N;++j) {
                    RT temp = m.col(j,i1,i2).maxAbs2Element();
                    if (temp > max) max = temp;
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(ptrdiff_t i=-m.nlo();i<=m.nhi();++i) {
                    RT temp = m.diag(i).maxAbs2Element();
                    if (temp > max) max = temp;
                }
            }
        }
        return max;
    }

    // 1-norm = max_j (sum_i |a_ij|)
    template <class T> 
    static RT NonLapNorm1(const GenBandMatrix<T>& m)
    {
        const ptrdiff_t M = m.colsize();
        const ptrdiff_t N = m.rowsize();
        RT max = 0;
        if (M > 0 && N > 0) {
            ptrdiff_t i1=0;
            ptrdiff_t i2=m.nlo()+1;
            ptrdiff_t k=m.nhi();
            for(ptrdiff_t j=0;j<N;++j) {
                RT temp = m.col(j,i1,i2).norm1();
                if (temp > max) max = temp;
                if (k>0) --k; else ++i1;
                if (i2<M) ++i2;
                else if (i1==M) break;
            }
        }
        return max;
    }

    template <class T> 
    static RT NonLapNormF(const GenBandMatrix<T>& m)
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
    static inline RT NonLapNormInf(const GenBandMatrix<T>& m)
    { return NonLapNorm1(m.transpose()); }

#ifdef XLAP
    template <class T> 
    static RT LapNorm(const char c, const GenBandMatrix<T>& m)
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
    double LapNorm(const char c, const GenBandMatrix<double>& m)
    {
        TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
        TMVAssert(m.isSquare());
        char cc = c;
        int n = m.colsize();
        double norm;
        if (BlasIsCM(m)) {
            int kl = m.nlo();
            int ku = m.nhi();
            int lda = m.stepj()+1;
#ifndef LAPNOWORK
            int lwork = c=='I' ? n : 0;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            norm = LAPNAME(dlangb) (
                LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
                LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
        } else {
            norm = LAPNAME(dlangt) (
                LAPCM LAPV(cc),LAPV(n),
                LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),
                LAPP(m.cptr()+m.stepj()));
        }
        return norm;
    }
    template <> 
    double LapNorm(const char c, const GenBandMatrix<std::complex<double> >& m)
    {
        TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
        TMVAssert(m.isSquare());
        char cc = c;
        int n = m.colsize();
        double norm;
        if (BlasIsCM(m)) {
            int kl = m.nlo();
            int ku = m.nhi();
            int lda = m.stepj()+1;
#ifndef LAPNOWORK
            int lwork = c=='I' ? n : 0;
            AlignedArray<double> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            norm = LAPNAME(zlangb) (
                LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
                LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
        } else {
            norm = LAPNAME(zlangt) (
                LAPCM LAPV(cc),LAPV(n),
                LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),
                LAPP(m.cptr()+m.stepj()));
        }
        return norm;
    }
#endif
#ifdef INST_FLOAT
    template <> 
    float LapNorm(const char c, const GenBandMatrix<float>& m)
    {
        TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
        TMVAssert(m.isSquare());
        char cc = c;
        int n = m.colsize();
        float norm;
        if (BlasIsCM(m)) {
            int kl = m.nlo();
            int ku = m.nhi();
            int lda = m.stepj()+1;
#ifndef LAPNOWORK
            int lwork = c=='I' ? n : 0;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            norm = LAPNAME(slangb) (
                LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
                LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
        } else {
            norm = LAPNAME(slangt) (
                LAPCM LAPV(cc),LAPV(n),
                LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),
                LAPP(m.cptr()+m.stepj()));
        }
        return norm;
    }
    template <> 
    float LapNorm(const char c, const GenBandMatrix<std::complex<float> >& m)
    {
        TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
        TMVAssert(m.isSquare());
        char cc = c;
        int n = m.colsize();
        float norm;
        if (BlasIsCM(m)) {
            int kl = m.nlo();
            int ku = m.nhi();
            int lda = m.stepj()+1;
#ifndef LAPNOWORK
            int lwork = c=='I' ? n : 0;
            AlignedArray<float> work(lwork);
            VectorViewOf(work.get(),lwork).setZero();
#endif
            norm = LAPNAME(clangb) (
                LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
                LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
        } else {
            norm = LAPNAME(clangt) (
                LAPCM LAPV(cc),LAPV(n),
                LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),
                LAPP(m.cptr()+m.stepj()));
        }
        return norm;
    }
#endif
#endif // XLAP

#ifdef INST_INT
    static int NonLapNormF(const GenBandMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNormF(const GenBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNorm1(const GenBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapMaxAbsElement(const GenBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T> 
    RT GenBandMatrix<T>::maxAbsElement() const
    {
#ifdef XLAP
        if (BlasIsRM(*this) && this->isSquare()) 
            return LapNorm('M',transpose());
        else if (BlasIsCM(*this) && this->isSquare()) return LapNorm('M',*this);
        else if (isdm() && this->isSquare() && nlo()==1 && nhi()==1)
            return LapNorm('M',*this);
        else
#endif
            return NonLapMaxAbsElement(*this);
    }
    template <class T> 
    RT GenBandMatrix<T>::maxAbs2Element() const
    {
#ifdef XLAP
        if (Traits<T>::iscomplex) return NonLapMaxAbs2Element(*this);
        else if (BlasIsRM(*this) && this->isSquare()) 
            return LapNorm('M',transpose());
        else if (BlasIsCM(*this) && this->isSquare()) return LapNorm('M',*this);
        else if (isdm() && this->isSquare() && nlo()==1 && nhi()==1)
            return LapNorm('M',*this);
        else
#endif
            return NonLapMaxAbs2Element(*this);
    }
    template <class T> 
    RT GenBandMatrix<T>::norm1() const
    {
#ifdef XLAP
        if (BlasIsRM(*this) && this->isSquare()) return LapNorm('I',transpose());
        else if (BlasIsCM(*this) && this->isSquare()) return LapNorm('1',*this);
        else if (isdm() && this->isSquare() && nlo()==1 && nhi()==1)
            return LapNorm('1',*this);
        else
#endif
            return NonLapNorm1(*this);
    }
    template <class T> 
    RT GenBandMatrix<T>::normF() const
    {
#ifdef XLAP
        if (BlasIsRM(*this) && this->isSquare()) return LapNorm('F',transpose());
        else if (BlasIsCM(*this) && this->isSquare()) return LapNorm('F',*this);
        else if (isdm() && this->isSquare() && nlo()==1 && nhi()==1)
            return LapNorm('F',*this);
        else
#endif
            return NonLapNormF(*this);
    }

    template <class T> 
    static RT DoNorm2(const GenBandMatrix<T>& m)
    {
        if (m.colsize() < m.rowsize()) return DoNorm2(m.transpose());
        if (m.rowsize() == 0) return RT(0);
        DiagMatrix<RT> S(m.rowsize());
        SV_Decompose(m,S.view());
        return S(0);
    }

    template <class T> 
    static RT DoCondition(const GenBandMatrix<T>& m)
    {
        if (m.colsize() < m.rowsize()) return DoCondition(m.transpose());
        if (m.rowsize() == 0) return RT(1);
        DiagMatrix<RT> S(m.rowsize());
        SV_Decompose(m,S.view());
        return S(0)/S(S.size()-1);
    }

#ifdef INST_INT
    static int DoNorm2(const GenBandMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoCondition(const GenBandMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoNorm2(const GenBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoCondition(const GenBandMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T> 
    RT GenBandMatrix<T>::doNorm2() const
    { return tmv::DoNorm2(*this); }

    template <class T> 
    RT GenBandMatrix<T>::doCondition() const
    { return tmv::DoCondition(*this); }

    //
    // Modifying Functions
    //

    template <class T, int A>
    BandMatrixView<T,A>& BandMatrixView<T,A>::clip(RT thresh) 
    {
        if (this->canLinearize()) linearView().clip(thresh);
        else {
            const ptrdiff_t M = colsize();
            const ptrdiff_t N = rowsize();
            if (M > 0 && N > 0) {
                if (isrm()) {
                    ptrdiff_t j1=0;
                    ptrdiff_t j2=nhi()+1;
                    ptrdiff_t k=nlo();
                    for(ptrdiff_t i=0;i<M;++i) {
                        row(i,j1,j2).clip(thresh);
                        if (k>0) --k; else ++j1;
                        if (j2<N) ++j2;
                        else if (j1==N) break;
                    }
                } else if (iscm()) {
                    ptrdiff_t i1=0;
                    ptrdiff_t i2=nlo()+1;
                    ptrdiff_t k=nhi();
                    for(ptrdiff_t j=0;j<N;++j) {
                        col(j,i1,i2).clip(thresh);
                        if (k>0) --k; else ++i1;
                        if (i2<M) ++i2;
                        else if (i1==M) break;
                    }
                } else {
                    for(ptrdiff_t i=-nlo();i<=nhi();++i) diag(i).clip(thresh);
                }
            }
        }
        return *this;
    }

    template <class T, int A>
    BandMatrixView<T,A>& BandMatrixView<T,A>::setZero() 
    {
        if (this->canLinearize()) linearView().setZero();
        else {
            const ptrdiff_t M = colsize();
            const ptrdiff_t N = rowsize();
            if (M > 0 && N > 0) {
                if (isrm()) {
                    ptrdiff_t j1=0;
                    ptrdiff_t j2=nhi()+1;
                    ptrdiff_t k=nlo();
                    for(ptrdiff_t i=0;i<M;++i) {
                        row(i,j1,j2).setZero();
                        if (k>0) --k; else ++j1;
                        if (j2<N) ++j2;
                        else if (j1==N) break;
                    }
                } else if (iscm()) {
                    ptrdiff_t i1=0;
                    ptrdiff_t i2=nlo()+1;
                    ptrdiff_t k=nhi();
                    for(ptrdiff_t j=0;j<N;++j) {
                        col(j,i1,i2).setZero();
                        if (k>0) --k; else ++i1;
                        if (i2<M) ++i2;
                        else if (i1==M) break;
                    }
                } else {
                    for(ptrdiff_t i=-nlo();i<=nhi();++i) diag(i).setZero();
                }
            }
        }
        return *this;
    }

    template <class T, int A>
    BandMatrixView<T,A>& BandMatrixView<T,A>::setAllTo(const T& x) 
    {
        if (this->canLinearize()) linearView().setAllTo(x);
        else {
            const ptrdiff_t M = colsize();
            const ptrdiff_t N = rowsize();
            if (M > 0 && N > 0) {
                if (isrm()) {
                    ptrdiff_t j1=0;
                    ptrdiff_t j2=nhi()+1;
                    ptrdiff_t k=nlo();
                    for(ptrdiff_t i=0;i<M;++i) {
                        row(i,j1,j2).setAllTo(x);
                        if (k>0) --k; else ++j1;
                        if (j2<N) ++j2;
                        else if (j1==N) break;
                    }
                } else if (iscm()) {
                    ptrdiff_t i1=0;
                    ptrdiff_t i2=nlo()+1;
                    ptrdiff_t k=nhi();
                    for(ptrdiff_t j=0;j<N;++j) {
                        col(j,i1,i2).setAllTo(x);
                        if (k>0) --k; else ++i1;
                        if (i2<M) ++i2;
                        else if (i1==M) break;
                    }
                } else {
                    for(ptrdiff_t i=-nlo();i<=nhi();++i) diag(i).setAllTo(x);
                }
            }
        }
        return *this;
    }

    template <class T, int A>
    BandMatrixView<T,A>& BandMatrixView<T,A>::addToAll(const T& x) 
    {
        if (this->canLinearize()) linearView().addToAll(x);
        else {
            const ptrdiff_t M = colsize();
            const ptrdiff_t N = rowsize();
            if (M > 0 && N > 0) {
                if (isrm()) {
                    ptrdiff_t j1=0;
                    ptrdiff_t j2=nhi()+1;
                    ptrdiff_t k=nlo();
                    for(ptrdiff_t i=0;i<M;++i) {
                        row(i,j1,j2).addToAll(x);
                        if (k>0) --k; else ++j1;
                        if (j2<N) ++j2;
                        else if (j1==N) break;
                    }
                } else if (iscm()) {
                    ptrdiff_t i1=0;
                    ptrdiff_t i2=nlo()+1;
                    ptrdiff_t k=nhi();
                    for(ptrdiff_t j=0;j<N;++j) {
                        col(j,i1,i2).addToAll(x);
                        if (k>0) --k; else ++i1;
                        if (i2<M) ++i2;
                        else if (i1==M) break;
                    }
                } else {
                    for(ptrdiff_t i=-nlo();i<=nhi();++i) diag(i).addToAll(x);
                }
            }
        }
        return *this;
    }

    template <class T, int A>
    void BandMatrixView<T,A>::doTransposeSelf() 
    {
        TMVAssert(colsize() == rowsize());
        TMVAssert(nlo() == nhi());
        for(ptrdiff_t i=1;i<=nhi();++i) Swap(diag(-i),diag(i));
    }

    template <class T, int A>
    BandMatrixView<T,A>& BandMatrixView<T,A>::conjugateSelf() 
    {
        if (this->canLinearize()) linearView().conjugateSelf();
        else {
            const ptrdiff_t M = colsize();
            const ptrdiff_t N = rowsize();
            if (isComplex(T()) && M > 0 && N > 0) {
                if (isrm()) {
                    ptrdiff_t j1=0;
                    ptrdiff_t j2=nhi()+1;
                    ptrdiff_t k=nlo();
                    for(ptrdiff_t i=0;i<M;++i) {
                        row(i,j1,j2).conjugateSelf();
                        if (k>0) --k; else ++j1;
                        if (j2<N) ++j2;
                        else if (j1==N) break;
                    }
                } else if (iscm()) {
                    ptrdiff_t i1=0;
                    ptrdiff_t i2=nlo()+1;
                    ptrdiff_t k=nhi();
                    for(ptrdiff_t j=0;j<N;++j) {
                        col(j,i1,i2).conjugateSelf();
                        if (k>0) --k; else ++i1;
                        if (i2<M) ++i2;
                        else if (i1==M) break;
                    }
                } else {
                    for(ptrdiff_t i=-nlo();i<=nhi();++i) diag(i).conjugateSelf();
                }
            }
        }
        return *this;
    }


    //
    // Special Constructors
    //

    template <class T> 
    BandMatrix<T,DiagMajor> UpperBiDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2)
    {
        if (v1.size() == v2.size()) {
            BandMatrix<T,DiagMajor> temp(v1.size(),v1.size()+1,0,1);
            temp.diag(0) = v1;
            temp.diag(1) = v2;
            return temp;
        } else {
            TMVAssert2(v2.size() == v1.size()-1);
            BandMatrix<T,DiagMajor> temp(v1.size(),v1.size(),0,1);
            temp.diag(0) = v1;
            temp.diag(1) = v2;
            return temp;
        }
    }

    template <class T> 
    BandMatrix<T,DiagMajor> LowerBiDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2)
    {
        if (v1.size() == v2.size()) {
            BandMatrix<T,DiagMajor> temp(v2.size()+1,v2.size(),1,0);
            temp.diag(-1) = v1;
            temp.diag(0) = v2;
            return temp;
        } else {
            TMVAssert2(v1.size() == v2.size()-1);
            BandMatrix<T,DiagMajor> temp(v2.size(),v2.size(),1,0);
            temp.diag(-1) = v1;
            temp.diag(0) = v2;
            return temp;
        }
    }

    template <class T> 
    BandMatrix<T,DiagMajor> TriDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2,
        const GenVector<T>& v3)
    {
        if (v1.size() == v2.size()) {
            TMVAssert2(v3.size() == v2.size()-1);
            BandMatrix<T,DiagMajor> temp(v2.size()+1,v2.size(),1,1);
            temp.diag(-1) = v1;
            temp.diag(0) = v2;
            temp.diag(1) = v3;
            return temp;
        } else if (v2.size() == v3.size()) {
            TMVAssert2(v1.size() == v2.size()-1);
            BandMatrix<T,DiagMajor> temp(v2.size(),v2.size()+1,1,1);
            temp.diag(-1) = v1;
            temp.diag(0) = v2;
            temp.diag(1) = v3;
            return temp;
        } else {
            TMVAssert2(v1.size() == v2.size()-1);
            TMVAssert2(v3.size() == v2.size()-1);
            BandMatrix<T,DiagMajor> temp(v2.size(),v2.size(),1,1);
            temp.diag(-1) = v1;
            temp.diag(0) = v2;
            temp.diag(1) = v3;
            return temp;
        }
    }


    //
    // Swap
    //

    template <class T> 
    void Swap(BandMatrixView<T> m1, BandMatrixView<T> m2)
    {
        TMVAssert2(m1.colsize() == m2.colsize());
        TMVAssert2(m1.rowsize() == m2.rowsize());
        TMVAssert2(m1.nlo() == m2.nlo());
        TMVAssert2(m1.nhi() == m2.nhi());
        const ptrdiff_t M = m1.colsize();
        const ptrdiff_t N = m1.rowsize();
        if (M > 0 && N > 0) {
            if (m1.isrm() && m2.isrm()) {
                ptrdiff_t j1=0;
                ptrdiff_t j2=m1.nhi()+1;
                ptrdiff_t k=m1.nlo();
                for(ptrdiff_t i=0;i<M;++i) {
                    Swap(m1.row(i,j1,j2),m2.row(i,j1,j2));
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (m1.iscm() && m2.iscm()) {
                ptrdiff_t i1=0;
                ptrdiff_t i2=m1.nlo()+1;
                ptrdiff_t k=m1.nhi();
                for(ptrdiff_t j=0;j<N;++j) {
                    Swap(m1.col(j,i1,i2),m2.col(j,i1,i2));
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(ptrdiff_t i=-m1.nlo();i<=m1.nhi();++i) 
                    Swap(m1.diag(i),m2.diag(i));
            }
        }
    }

    //
    // m1 == m2
    //

    template <class T1, class T2> 
    bool operator==(const GenBandMatrix<T1>& m1, const GenBandMatrix<T2>& m2)
    {
        if (m1.colsize() != m2.colsize()) return false;
        else if (m1.rowsize() != m2.rowsize()) return false;
        else if (m1.isSameAs(m2)) return true;
        else {
            ptrdiff_t lo = TMV_MIN(m1.nlo(),m2.nlo());
            ptrdiff_t hi = TMV_MIN(m1.nhi(),m2.nhi());
            for(ptrdiff_t i=-lo;i<=hi;++i) 
                if (m1.diag(i) != m2.diag(i)) return false;
            for(ptrdiff_t i=-m1.nlo();i<-lo;++i)
                if (m1.diag(i).maxAbs2Element() != T1(0)) return false;
            for(ptrdiff_t i=-m2.nlo();i<-lo;++i)
                if (m2.diag(i).maxAbs2Element() != T2(0)) return false;
            for(ptrdiff_t i=hi+1;i<=m1.nhi();++i)
                if (m1.diag(i).maxAbs2Element() != T1(0)) return false;
            for(ptrdiff_t i=hi+1;i<=m2.nhi();++i)
                if (m2.diag(i).maxAbs2Element() != T2(0)) return false;
            return true;
        }
    }

    template <class T1, class T2> 
    bool operator==(const GenBandMatrix<T1>& m1, const GenMatrix<T2>& m2)
    {
        if (m1.colsize() != m2.colsize()) return false;
        else if (m1.rowsize() != m2.rowsize()) return false;
        else {
            ConstBandMatrixView<T2> m2b = 
                BandMatrixViewOf(m2,m2.colsize()-1,m2.rowsize()-1);
            if ( m1.diagRange(-m1.nlo(),m1.nhi()+1) !=
                 m2b.diagRange(-m1.nlo(),m1.nhi()+1)) 
                return false;
            if ( m1.nhi()+1 < m1.rowsize() &&
                 m2b.diagRange(m1.nhi()+1,m1.rowsize()).maxAbs2Element() != T2(0)) 
                return false;
            if ( m1.nlo()+1 < m1.colsize() &&
                 m2b.diagRange(-m1.colsize()+1,-m1.nlo()).maxAbs2Element() != T2(0)) 
                return false;
            return true;
        }
    }


    //
    // I/O
    //

    template <class T>
    void GenBandMatrix<T>::write(const TMV_Writer& writer) const
    {
        const ptrdiff_t M = colsize();
        const ptrdiff_t N = rowsize();
        ptrdiff_t j1=0;
        ptrdiff_t j2=nhi()+1;

        writer.begin();
        writer.writeCode("B");
        writer.writeSize(M);
        writer.writeSize(N);
        writer.writeFullSize(nlo());
        writer.writeFullSize(nhi());
        writer.writeStart();

        for(ptrdiff_t i=0;i<M;++i) {
            writer.writeLParen();
            if (!writer.isCompact()) {
                for(ptrdiff_t j=0;j<j1;++j) {
                    writer.writeValue(T(0));
                    if (j<N-1) writer.writeSpace();
                }
            }
            for(ptrdiff_t j=j1;j<j2;++j) {
                if (j > j1) writer.writeSpace();
                writer.writeValue(cref(i,j));
            }
            if (!writer.isCompact()) {
                for(ptrdiff_t j=j2;j<N;++j) {
                    writer.writeSpace();
                    writer.writeValue(T(0));
                }
            }
            writer.writeRParen();
            if (i < M-1) writer.writeRowEnd();
            if (j2 < N) ++j2;
            if (i >= nlo() && j1 < N) ++j1;
        }
        writer.writeFinal();
        writer.end();
    }

#ifndef NOTHROW
    template <class T> 
    class BandMatrixReadError : public ReadError
    {
    public :
        BandMatrix<T> m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t cs,rs;
        ptrdiff_t lo,hi;
        T v1;
        bool is, iseof, isbad;

        BandMatrixReadError(std::istream& _is) throw() :
            ReadError("BandMatrix."),
            i(0), j(0), cs(0), rs(0), lo(0), hi(0), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        BandMatrixReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("BandMatrix."),
            i(0), j(0), exp(_e), got(_g), cs(0), rs(0), lo(0), hi(0), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        BandMatrixReadError(
            const GenBandMatrix<T>& _m, std::istream& _is,
            ptrdiff_t _cs, ptrdiff_t _rs, ptrdiff_t _lo, ptrdiff_t _hi) throw() :
            ReadError("BandMatrix."),
            m(_m), i(0), j(0), cs(_cs), rs(_rs), lo(_lo), hi(_hi), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        BandMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenBandMatrix<T>& _m,
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("BandMatrix."),
            m(_m), i(_i), j(_j), exp(_e), got(_g), 
            cs(m.colsize()), rs(m.rowsize()), lo(m.nlo()), hi(m.nhi()), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        BandMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenBandMatrix<T>& _m,
            std::istream& _is, T _v1=0) throw() :
            ReadError("BandMatrix."),
            m(_m), i(_i), j(_j),
            cs(m.colsize()), rs(m.rowsize()), lo(m.nlo()), hi(m.nhi()), v1(_v1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        BandMatrixReadError(const BandMatrixReadError<T>& rhs) :
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got), 
            cs(rhs.cs), rs(rhs.rs), lo(rhs.lo), hi(rhs.hi), v1(rhs.v1),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~BandMatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for BandMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (cs != m.colsize()) {
                os<<"Wrong colsize: expected "<<m.colsize()<<
                    ", got "<<cs<<".\n";
            }
            if (rs != m.rowsize()) {
                os<<"Wrong rowsize: expected "<<m.rowsize()<<
                    ", got "<<rs<<".\n";
            }
            if (lo != m.nlo()) {
                os<<"Wrong nlo: expected "<<m.nlo()<<", got "<<lo<<".\n";
            }
            if (hi != m.nhi()) {
                os<<"Wrong nhi: expected "<<m.nhi()<<", got "<<hi<<".\n";
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
            if (v1 != T(0)) {
                os<<"Invalid input.  Expected 0, got "<<v1<<".\n";
            }
            if (m.colsize() > 0 || m.rowsize() > 0) {
                os<<"The portion of the BandMatrix which was successfully "
                    "read is: \n";
                const ptrdiff_t N = m.rowsize();
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
    static void FinishRead(const TMV_Reader& reader, BandMatrixView<T> m)
    {
        const ptrdiff_t M = m.colsize();
        const ptrdiff_t N = m.rowsize();
        std::string exp,got;
        T temp;
        ptrdiff_t j1=0;
        ptrdiff_t j2=m.nhi()+1;

        if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw BandMatrixReadError<T>(0,0,m,reader.getis(),exp,got);
#endif
        }
        for(ptrdiff_t i=0;i<M;++i) {
            if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw BandMatrixReadError<T>(i,0,m,reader.getis(),exp,got);
#endif
            }
            if (!reader.isCompact()) {
                for(ptrdiff_t j=0;j<j1;++j) {
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"BandMatrix Read Error: reading value\n";
                        exit(1);
#else
                        throw BandMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    if (temp != T(0)) {
#ifdef NOTHROW
                        std::cerr<<"BandMatrix Read Error: "<<temp<<" != 0\n";
                        exit(1);
#else
                        throw BandMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                    }
                    if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw BandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                }
            }
            for(ptrdiff_t j=j1;j<j2;++j) {
                if (j>j1 && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw BandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                }
                if (!reader.readValue(temp)) {
#ifdef NOTHROW
                    std::cerr<<"BandMatrix Read Error: reading value\n";
                    exit(1);
#else
                    throw BandMatrixReadError<T>(i,j,m,reader.getis());
#endif
                }
                m.ref(i,j) = temp;
            }
            if (!reader.isCompact()) {
                for(ptrdiff_t j=j2;j<N;++j) {
                    if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw BandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"BandMatrix Read Error: reading value\n";
                        exit(1);
#else
                        throw BandMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    if (temp != T(0)) {
#ifdef NOTHROW
                        std::cerr<<"BandMatrix Read Error: "<<temp<<" != 0\n";
                        exit(1);
#else
                        throw BandMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                    }
                }
            }
            if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw BandMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
            if (i < M-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw BandMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
            if (j2 < N) ++j2;
            if (i >= m.nlo() && j1 < N) ++j1;
        }
        if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw BandMatrixReadError<T>(N,0,m,reader.getis(),exp,got);
#endif
        }
    }

    template <class T, int A>
    void BandMatrix<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        if (!reader.readCode("B",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw BandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t cs=colsize(), rs=rowsize(), lo=nlo(), hi=nhi();
        if (!reader.readSize(cs,exp,got) || 
            !reader.readSize(rs,exp,got) || 
            !reader.readFullSize(lo,exp,got) || 
            !reader.readFullSize(hi,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: reading size\n";
            exit(1);
#else
            throw BandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (cs != colsize() || rs != rowsize() || lo != nlo() || hi != nhi()) {
            resize(cs,rs,lo,hi);
        }
        BandMatrixView<T> v = view();
        FinishRead(reader,v);
    }

    template <class T, int A>
    void BandMatrixView<T,A>::read(const TMV_Reader& reader) 
    {
        std::string exp,got;
        if (!reader.readCode("B",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw BandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t cs=colsize(), rs=rowsize(), lo=nlo(), hi=nhi();
        if (!reader.readSize(cs,exp,got) || 
            !reader.readSize(rs,exp,got) || 
            !reader.readFullSize(lo,exp,got) || 
            !reader.readFullSize(hi,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: reading size\n";
            exit(1);
#else
            throw BandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (cs != colsize() || rs != rowsize() || lo != nlo() || hi != nhi()) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw BandMatrixReadError<T>(*this,reader.getis(),cs,rs,lo,hi);
#endif
        }
        FinishRead(reader,*this);
    }

#undef RT

#define InstFile "TMV_BandMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


