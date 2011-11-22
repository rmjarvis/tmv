///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
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
#include "portable_platform.h"

namespace tmv {

#define RT TMV_RealType(T)
#define CT TMV_ComplexType(T)

    //
    // Access
    //

    template <class T> 
    T GenSymBandMatrix<T>::cref(int i, int j) const
    {
        if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
            const T* mi = cptr() + int(i)*stepi() + int(j)*stepj();
            return isconj() ? TMV_CONJ(*mi) : *mi;
        } else {
            const T* mi = cptr() + int(j)*stepi() + int(i)*stepj();
            return issym() != isconj() ? *mi : TMV_CONJ(*mi);
        }
    }

    template <class T, IndexStyle I> 
    typename SymBandMatrixView<T,I>::reference SymBandMatrixView<T,I>::ref(
        int i, int j) const
    {
        if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
            T* mi = ptr() + int(i)*stepi() + int(j)*stepj();
            return TMV_REF(mi,ct());
        } else {
            T* mi = ptr() + int(j)*stepi() + int(i)*stepj();
            return this->issym() != this->isconj() ? TMV_REF(mi,NonConj) : TMV_REF(mi,Conj);
        }
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::newCopy() const
    {
        auto_ptr<BaseMatrix<T> > a;
        if (issym()) {
            if (uplo() == Upper) {
                if (isrm()) a.reset(
                    new SymBandMatrix<T,Upper,RowMajor>(*this));
                else if (iscm()) a.reset(
                    new SymBandMatrix<T,Upper,ColMajor>(*this));
                else a.reset(
                    new SymBandMatrix<T,Upper,DiagMajor>(*this));
            } else {
                if (isrm()) a.reset(
                    new SymBandMatrix<T,Lower,RowMajor>(*this));
                else if (iscm()) a.reset(
                    new SymBandMatrix<T,Lower,ColMajor>(*this));
                else a.reset(
                    new SymBandMatrix<T,Lower,DiagMajor>(*this));
            }
        } else {
            if (uplo() == Upper) {
                if (isrm()) a.reset(
                    new HermBandMatrix<T,Upper,RowMajor>(*this));
                else if (iscm()) a.reset(
                    new HermBandMatrix<T,Upper,ColMajor>(*this));
                else a.reset(
                    new HermBandMatrix<T,Upper,DiagMajor>(*this));
            } else {
                if (isrm()) a.reset(
                    new HermBandMatrix<T,Lower,RowMajor>(*this));
                else if (iscm()) a.reset(
                    new HermBandMatrix<T,Lower,ColMajor>(*this));
                else a.reset(
                    new HermBandMatrix<T,Lower,DiagMajor>(*this));
            }
        }
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::newView() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstSymBandMatrixView<T>(view()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::newTranspose() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstSymBandMatrixView<T>(transpose()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::newConjugate() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstSymBandMatrixView<T>(conjugate()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::newAdjoint() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstSymBandMatrixView<T>(adjoint()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenSymBandMatrix<T>::newInverse() const
    {
        if (issym()) {
            auto_ptr<SymMatrix<T,Upper,ColMajor> > a(
                new SymMatrix<T,Upper,ColMajor>(size()));
            makeInverse(a->view());
            BaseMatrix<T>* ret1 = a.release();
            auto_ptr<BaseMatrix<T> > ret(ret1);
            return ret;
        } else {
            auto_ptr<HermMatrix<T,Upper,ColMajor> > a(
                new HermMatrix<T,Upper,ColMajor>(size()));
            makeInverse(a->view());
            BaseMatrix<T>* ret1 = a.release();
            auto_ptr<BaseMatrix<T> > ret(ret1);
            return ret;
        }
    }
#ifdef INST_INT
    template <>
    auto_ptr<BaseMatrix<int> > GenSymBandMatrix<int>::newInverse() const
    { TMVAssert(TMV_FALSE); return auto_ptr<BaseMatrix<int> >(); }
    template <>
    auto_ptr<BaseMatrix<std::complex<int> > > 
    GenSymBandMatrix<std::complex<int> >::newInverse() const
    { 
        TMVAssert(TMV_FALSE); 
        return auto_ptr<BaseMatrix<std::complex<int> > >(); 
    }
#endif

    template <class T> 
    void GenSymBandMatrix<T>::newDivider() const
    {
        switch(this->getDivType()) {
          case LU : {
              this->setDiv(new BandLUDiv<T>(*this));
              break; 
          }
          case SV : {
              if (isherm()) this->setDiv(new HermBandSVDiv<T>(*this)); 
              else this->setDiv(new SymBandSVDiv<T>(*this)); 
              break; 
          }
          case CH : {
              this->setDiv(
                  new HermBandCHDiv<T>(*this,this->isDivInPlace()));
              break;
          }
          default : TMVAssert(TMV_FALSE); 
        }
    }

#ifdef INST_INT
    template <>
    void GenSymBandMatrix<int>::newDivider() const
    { TMVAssert(TMV_FALSE); }
    template <>
    void GenSymBandMatrix<std::complex<int> >::newDivider() const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T>
    bool GenSymBandMatrix<T>::divIsLUDiv() const
    { return static_cast<bool>(dynamic_cast<const BandLUDiv<T>*>(getDiv())); }

    template <class T>
    bool GenSymBandMatrix<T>::divIsCHDiv() const
    {
        return static_cast<bool>(
            dynamic_cast<const HermBandCHDiv<T>*>(getDiv())); 
    }

    template <class T>
    bool GenSymBandMatrix<T>::divIsHermSVDiv() const
    {
        return static_cast<bool>(
            dynamic_cast<const HermBandSVDiv<T>*>(getDiv())); 
    }

    template <class T>
    bool GenSymBandMatrix<T>::divIsSymSVDiv() const
    { 
        return static_cast<bool>(
            dynamic_cast<const SymBandSVDiv<T>*>(getDiv())); 
    }

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
    void GenSymBandMatrix<T>::doMakeInverse(const SymMatrixView<T1>& sinv) const
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
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
        if (i1==i2 || j1==j2) return true; // no elements, so whatever...
        bool ok = true;
        int i2x = i2-istep;
        int j2x = j2-jstep;
        if (istep == 0) {
            ok = false;
            std::cout<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1 < 0 || i1 >= int(size())) {
            ok = false;
            std::cout<<"first col element ("<<i1<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if (i2x < 0 || i2x >= int(size())) {
            ok = false;
            std::cout<<"last col element ("<<i2x<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cout<<"col range ("<<i2-i1<<") must be multiple of istep (";
            std::cout<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cout<<"n col elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
        }
        if (jstep == 0) {
            ok = false;
            std::cout<<"jstep ("<<jstep<<") can not be 0\n";
        }
        if (j1 < 0 || j1 >= int(size())) {
            ok = false;
            std::cout<<"first row element ("<<j1<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if (j2-jstep < 0 || j2-jstep >= int(size())) {
            ok = false;
            std::cout<<"last row element ("<<j2-jstep<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if ((j2-j1)%jstep != 0) {
            ok = false;
            std::cout<<"row range ("<<j2-j1<<") must be multiple of jstep (";
            std::cout<<jstep<<")\n";
        }
        if ((j2-j1)/jstep < 0) {
            ok = false;
            std::cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
        }
        if ((i1<j1 && i2x>j2x) || (i1>j1 && i2x<j2x)) {
            ok = false;
            std::cout<<"Upper left ("<<i1<<','<<j1<<") and lower right (";
            std::cout<<i2x<<','<<j2x<<") corners must be in same triangle\n";
        }
        if ((i2x<j1 && i1>j2x) || (i2x>j1 && i1<j2x)) {
            ok = false;
            std::cout<<"Upper right ("<<i1<<','<<j2x<<") and lower left (";
            std::cout<<i2x<<','<<j1<<") corners must be in same triangle\n";
        }
        if (!okij(i1,j2x)) {
            ok = false;
            std::cout<<"Upper right ("<<i1<<','<<j2x<<") corner must be in band.\n";
        }
        if (!okij(i2x,j1)) {
            ok = false;
            std::cout<<"Lower left ("<<i2x<<','<<j1<<") corner must be in band.\n";
        }
        return ok;
    }

    template <class T> 
    bool GenSymBandMatrix<T>::hasSubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi,
        int istep, int jstep) const
    {
        if (i1==i2 || j1==j2) return true; // no elements, so whatever...
        bool ok = true;
        if (istep == 0) {
            ok = false;
            std::cout<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1 < 0 || i1 >= int(colsize())) {
            ok = false;
            std::cout<<"first col element ("<<i1<<") must be in 0 -- ";
            std::cout<<colsize()-1<<std::endl;
        }
        if (i2-istep < 0 || i2-istep >= int(colsize())) {
            ok = false;
            std::cout<<"last col element ("<<i2-istep<<") must be in 0 -- ";
            std::cout<<colsize()-1<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cout<<"col range ("<<i2-i1<<") must be multiple of istep (";
            std::cout<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cout<<"n col elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
        }
        if (jstep == 0) {
            ok = false;
            std::cout<<"jstep ("<<jstep<<") can not be 0\n";
        }
        if (j1 < 0 || j1 >= int(rowsize())) {
            ok = false;
            std::cout<<"first row element ("<<j1<<") must be in 0 -- ";
            std::cout<<rowsize()-1<<std::endl;
        }
        if (j2-jstep < 0 || j2-jstep >= int(rowsize())) {
            ok = false;
            std::cout<<"last row element ("<<j2-jstep<<") must be in 0 -- ";
            std::cout<<rowsize()-1<<std::endl;
        }
        if ((j2-j1)%jstep != 0) {
            ok = false;
            std::cout<<"row range ("<<j2-j1<<") must be multiple of istep (";
            std::cout<<jstep<<")\n";
        }
        if ((j2-j1)/jstep < 0) {
            ok = false;
            std::cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
        }
        if (!okij(i1,j1)) {
            ok = false;
            std::cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in band\n";
        }
        if (!okij(i1,j1+newnhi)) {
            ok = false;
            std::cout<<"Start of top diagonal ("<<i1<<','<<j1+newnhi;
            std::cout<<") must be in band\n";
        }
        if (!okij(i1+newnlo,j1)) {
            ok = false;
            std::cout<<"Start of bottom diagonal ("<<i1+newnlo<<','<<j1;
            std::cout<<") must be in band\n";
        }
        if (newnhi >= j2-j1) {
            ok = false;
            std::cout<<"new nhi ("<<newnhi<<") must be less than the new rowsize (";
            std::cout<<j2-j1<<")\n";
        }
        if (newnlo >= i2-i1) {
            ok = false;
            std::cout<<"new nlo ("<<newnlo<<") must be less than the new colsize (";
            std::cout<<i2-i1<<")\n";
        }
        if ((i1<j1+newnhi && i1+newnlo>j1) || (i1>j1+newnhi && i1+newnlo<j1)) {
            ok = false;
            std::cout<<"Start of top ("<<i1<<','<<j1+newnhi<<") and bottom (";
            std::cout<<i1+newnlo<<','<<j1<<") diagonals must be in same triangle\n";
        }
        return ok;
    }

    template <class T> 
    bool GenSymBandMatrix<T>::hasSubVector(
        int i, int j, int istep, int jstep, int n) const 
    {
        if (n==0) return true;
        bool ok = true;
        if (istep == 0 && jstep == 0) {
            ok = false;
            std::cout<<"istep ("<<istep<<") and jstep ("<<jstep;
            std::cout<<") can not both be 0\n";
        }
        if (i<0 || i >= int(size())) {
            ok = false;
            std::cout<<"i ("<<i<<") must be in 0 -- "<<size()-1<<std::endl;
        }
        if (j<0 || j >= int(size())) {
            ok = false;
            std::cout<<"j ("<<j<<") must be in 0 -- "<<size()-1<<std::endl;
        }
        int i2 = int(i)+istep*int(n-1);
        int j2 = int(j)+jstep*int(n-1);
        if (i2 < 0 || i2 >= int(size())) {
            ok = false;
            std::cout<<"last element's i ("<<i2<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if (j2 < 0 || j2 >= int(size())) {
            ok = false;
            std::cout<<"last element's j ("<<j2<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if ((i<j && i2>j2) || (i>j && i2<j2)) {
            ok = false;
            std::cout<<"First ("<<i<<','<<j<<") and last ("<<i2<<','<<j2;
            std::cout<<") elements must be in same triangle\n";
        }
        if (!okij(i,j)) {
            ok = false;
            std::cout<<"First ("<<i<<','<<j<<") element must be in band\n";
        }
        if (!okij(i2,j2)) {
            ok = false;
            std::cout<<"Last ("<<i2<<','<<j2<<") element must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool GenSymBandMatrix<T>::hasSubSymMatrix(
        int i1, int i2, int istep) const 
    {
        if (i1==i2) return true;
        bool ok=true;
        if (istep == 0) {
            ok = false; 
            std::cout<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1<0 || i1 >= int(size())) {
            ok = false;
            std::cout<<"first diag element ("<<i1<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if (i2-istep<0 || i2-istep >= int(size())) {
            ok = false;
            std::cout<<"last diag element ("<<i2-istep<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cout<<"range ("<<i2-i1<<") must be multiple of istep (";
            std::cout<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cout<<"n diag elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
        }
        if (!okij(i1,i2-istep)) {
            ok = false;
            std::cout<<"Upper right ("<<i1<<','<<i2-istep;
            std::cout<<") corner must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool GenSymBandMatrix<T>::hasSubSymBandMatrix(
        int i1, int i2, int newnlo, int istep) const 
    {
        if (i1==i2) return true;
        bool ok=true;
        if (istep == 0) {
            ok = false; 
            std::cout<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1<0 || i1 >= int(size())) {
            ok = false;
            std::cout<<"first diag element ("<<i1<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if (i2-istep<0 || i2-istep >= int(size())) {
            ok = false;
            std::cout<<"last diag element ("<<i2-istep<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cout<<"range ("<<i2-i1<<") must be multiple of istep (";
            std::cout<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cout<<"n diag elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
        }
        if (newnlo > nlo()) {
            ok = false;
            std::cout<<"new number of off-diagonals ("<<newnlo<<") must be less ";
            std::cout<<"than or equal to the current value ("<<nlo()<<")\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymBandMatrixView<T,FortranStyle>::hasSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
        if (i1==i2 || j1==j2) return true; // no elements, so whatever...
        bool ok = true;
        if (istep == 0) {
            ok = false;
            std::cout<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1 < 1 || i1 > int(this->size())) {
            ok = false;
            std::cout<<"first col element ("<<i1<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if (i2 < 1 || i2 > int(this->size())) {
            ok = false;
            std::cout<<"last col element ("<<i2-istep<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cout<<"col range ("<<i2-i1<<") must be multiple of istep (";
            std::cout<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cout<<"n col elements ("<<(i2-i1)/istep+1<<") must be positive\n";
        }
        if (jstep == 0) {
            ok = false;
            std::cout<<"jstep ("<<jstep<<") can not be 0\n";
        }
        if (j1 < 1 || j1 > int(this->size())) {
            ok = false;
            std::cout<<"first row element ("<<j1<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if (j2 < 1 || j2 > int(this->size())) {
            ok = false;
            std::cout<<"last row element ("<<j2-jstep<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if ((j2-j1)%jstep != 0) {
            ok = false;
            std::cout<<"row range ("<<j2-j1<<") must be multiple of istep (";
            std::cout<<jstep<<")\n";
        }
        if ((j2-j1)/jstep < 0) {
            ok = false;
            std::cout<<"n row elements ("<<(j2-j1)/jstep+1<<") must be positive\n";
        }
        if ((i1<j1 && i2>j2) || (i1>j1 && i2<j2)) {
            ok = false;
            std::cout<<"Upper left ("<<i1<<','<<j1<<") and lower right (";
            std::cout<<i2<<','<<j2<<") corners must be in same triangle\n";
        }
        if ((i2<j1 && i1>j2) || (i2>j1 && i1<j2)) {
            ok = false;
            std::cout<<"Upper right ("<<i1<<','<<j2<<") and lower left (";
            std::cout<<i2<<','<<j1<<") corners must be in same triangle\n";
        }
        if (!this->okij(i1-1,j2-1)) {
            ok = false;
            std::cout<<"Upper right ("<<i1<<','<<j2<<") corner must be in band.\n";
        }
        if (!this->okij(i2-1,j1-1)) {
            ok = false;
            std::cout<<"Lower left ("<<i2<<','<<j1<<") corner must be in band.\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymBandMatrixView<T,FortranStyle>::hasSubBandMatrix(
        int i1, int i2, int j1, int j2, int newnlo, int newnhi,
        int istep, int jstep) const
    {
        if (i1==i2 || j1==j2) return true; // no elements, so whatever...
        bool ok = true;
        if (istep == 0) {
            ok = false;
            std::cout<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1 < 1 || i1 > int(this->colsize())) {
            ok = false;
            std::cout<<"first col element ("<<i1<<") must be in 1 -- ";
            std::cout<<this->colsize()<<std::endl;
        }
        if (i2 < 1 || i2 > int(this->colsize())) {
            ok = false;
            std::cout<<"last col element ("<<i2<<") must be in 1 -- ";
            std::cout<<this->colsize()<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cout<<"col range ("<<i2-i1<<") must be multiple of istep (";
            std::cout<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cout<<"n col elements ("<<(i2-i1)/istep+1<<") must be positive\n";
        }
        if (jstep == 0) {
            ok = false;
            std::cout<<"jstep ("<<jstep<<") can not be 0\n";
        }
        if (j1 < 1 || j1 > int(this->rowsize())) {
            ok = false;
            std::cout<<"first row element ("<<j1<<") must be in 1 -- ";
            std::cout<<this->rowsize()<<std::endl;
        }
        if (j2 < 1 || j2 > int(this->rowsize())) {
            ok = false;
            std::cout<<"last row element ("<<j2<<") must be in 1 -- ";
            std::cout<<this->rowsize()<<std::endl;
        }
        if ((j2-j1)%jstep != 0) {
            ok = false;
            std::cout<<"row range ("<<j2-j1<<") must be multiple of istep (";
            std::cout<<jstep<<")\n";
        }
        if ((j2-j1)/jstep < 0) {
            ok = false;
            std::cout<<"n row elements ("<<(j2-j1)/jstep+1<<") must be positive\n";
        }
        if (!this->okij(i1-1,j1-1)) {
            ok = false;
            std::cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in band\n";
        }
        if (!this->okij(i1-1,j1-1+newnhi)) {
            ok = false;
            std::cout<<"Start of top diagonal ("<<i1<<','<<j1+newnhi;
            std::cout<<") must be in band\n";
        }
        if (!this->okij(i1-1+newnlo,j1-1)) {
            ok = false;
            std::cout<<"Start of bottom diagonal ("<<i1+newnlo<<','<<j1;
            std::cout<<") must be in band\n";
        }
        if (newnhi >= j2-j1+1) {
            ok = false;
            std::cout<<"new nhi ("<<newnhi<<") must be less than the new rowsize (";
            std::cout<<j2-j1+1<<")\n";
        }
        if (newnlo >= i2-i1+1) {
            ok = false;
            std::cout<<"new nlo ("<<newnlo<<") must be less than the new colsize (";
            std::cout<<i2-i1+1<<")\n";
        }
        if ((i1<j1+newnhi && i1+newnlo>j1) || (i1>j1+newnhi && i1+newnlo<j1)) {
            ok = false;
            std::cout<<"Start of top ("<<i1<<','<<j1+newnhi<<") and bottom (";
            std::cout<<i1+newnlo<<','<<j1<<") diagonals must be in same triangle\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymBandMatrixView<T,FortranStyle>::hasSubVector(
        int i, int j, int istep, int jstep, int n) const 
    {
        if (n==0) return true;
        bool ok = true;
        if (istep == 0 && jstep == 0) {
            ok = false;
            std::cout<<"istep ("<<istep<<") and jstep ("<<jstep;
            std::cout<<") can not both be 0\n";
        }
        if (i < 1 || i > int(this->size())) {
            ok = false;
            std::cout<<"i ("<<i<<") must be in 1 -- "<<this->size()<<std::endl;
        }
        if (i < 1 || j > int(this->size())) {
            ok = false;
            std::cout<<"j ("<<j<<") must be in 1 -- "<<this->size()<<std::endl;
        }
        int i2 = int(i)+istep*int(n-1);
        int j2 = int(j)+jstep*int(n-1);
        if (i2 < 1 || i2 > int(this->size())) {
            ok = false;
            std::cout<<"last element's i ("<<i2<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if (j2 < 1 || j2 > int(this->size())) {
            ok = false;
            std::cout<<"last element's j ("<<j2<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if ((i<j && i2>j2) || (i>j && i2<j2)) {
            ok = false;
            std::cout<<"First ("<<i<<','<<j<<") and last ("<<i2<<','<<j2;
            std::cout<<") elements must be in same triangle\n";
        }
        if (!this->okij(i-1,j-1)) {
            ok = false;
            std::cout<<"First ("<<i<<','<<j<<") element must be in band\n";
        }
        if (!this->okij(i2-1,j2-1)) {
            ok = false;
            std::cout<<"Last ("<<i2<<','<<j2<<") element must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymBandMatrixView<T,FortranStyle>::hasSubSymMatrix(
        int i1, int i2, int istep) const 
    {
        if (i1==i2) return true;
        bool ok=true;
        if (istep == 0) {
            ok = false; 
            std::cout<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1<1 || i1 > int(this->size())) {
            ok = false;
            std::cout<<"first diag element ("<<i1<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if (i2-istep<1 || i2-istep > int(this->size())) {
            ok = false;
            std::cout<<"last diag element ("<<i2-istep<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cout<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep;
            std::cout<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cout<<"n diag elements ("<<(i2-i1)/istep+1<<") must be positive\n";
        }
        if (!this->okij(i1-1,i2-1)) {
            ok = false;
            std::cout<<"Upper right ("<<i1<<','<<i2<<") corner must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstSymBandMatrixView<T,FortranStyle>::hasSubSymBandMatrix(
        int i1, int i2, int newnlo, int istep) const 
    {
        if (i1==i2) return true;
        bool ok=true;
        if (istep == 0) {
            ok = false; 
            std::cout<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1<1 || i1 > int(this->size())) {
            ok = false;
            std::cout<<"first diag element ("<<i1<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if (i2-istep<1 || i2-istep > int(this->size())) {
            ok = false;
            std::cout<<"last diag element ("<<i2-istep<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cout<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep;
            std::cout<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cout<<"n diag elements ("<<(i2-i1)/istep+1<<") must be positive\n";
        }
        if (newnlo > this->nlo()) {
            ok = false;
            std::cout<<"new number of off-diagonals ("<<newnlo<<") must be less ";
            std::cout<<"than or equal to the current value ("<<this->nlo()<<")\n";
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
            const int N = m.size();
            if (N > 0) {
                int i1=0;
                int i2=m.nlo()+1;
                int k=m.nlo();
                for(int j=0;j<N;++j) {
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
        if (m.isrm()) return LapNorm(c,m.transpose());
        else {
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
        }
    }
    template <> 
    double LapNorm(const char c, 
                   const GenSymBandMatrix<std::complex<double> >& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        if (m.isrm()) return LapNorm(c,m.transpose());
        else {
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
        }
    }
#endif
#ifdef INST_FLOAT
    template <>
    float LapNorm(const char c, const GenSymBandMatrix<float>& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        if (m.isrm()) return LapNorm(c,m.transpose());
        else {
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
        }
    }
    template <>
    float LapNorm(const char c, 
                  const GenSymBandMatrix<std::complex<float> >& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        if (m.isrm()) return LapNorm(c,m.transpose());
        else {
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
    SymBandMatrix<T,Upper,DiagMajor> SymTriDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2)
    {
        TMVAssert2(v2.size() == v1.size()-1);
        SymBandMatrix<T,Upper,DiagMajor> temp(v1.size(),1);
        temp.diag() = v1;
        temp.diag(1) = v2;
        return temp;
    }

    template <class T> 
    HermBandMatrix<T,Upper,DiagMajor> HermTriDiagMatrix(
        const GenVector<T>& v1, const GenVector<T>& v2, UpLoType uplo)
    {
        TMVAssert2(v2.size() == v1.size()-1);
        HermBandMatrix<T,Upper,DiagMajor> temp(v1.size(),1);
        temp.diag() = v1;
        if (uplo == Upper)
            temp.diag(1) = v2;
        else
            temp.diag(-1) = v2;
        return temp;
    }

    template <class T> 
    HermBandMatrix<std::complex<T>,Upper,DiagMajor> HermTriDiagMatrix(
        const GenVector<T>& v1, const GenVector<std::complex<T> >& v2,
        UpLoType uplo)
    {
        TMVAssert2(v2.size() == v1.size()-1);
        HermBandMatrix<std::complex<T>,Upper,DiagMajor> temp(v1.size(),1);
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

    // This bit is to workaround a bug in pgCC that was fixed in version 7.
    // I don't know if versions earlier than 6.1 had the bug, but 
    // I apply the workaround to all version before 7.
    template <class T> 
    inline T Value(const T& x) { return x; }
#ifdef PLATFORM_COMPILER_PGI
#if PLATFORM_COMPILER_VERSION < 0x070000
    inline double Value(const long double& x) { return double(x); }
    inline std::complex<double> Value(const std::complex<long double>& x) 
    { return std::complex<double>(x); }
#endif
#endif

    template <bool sym, bool conj, bool rm, bool cm, bool compact, bool th, class T>
    static void DoWrite(
        std::ostream& os, const GenSymBandMatrix<T>& m, RT thresh)
    {
        TMVAssert(m.uplo() == Lower);
        int gaplen1 = 0;
        int rowlen1 = compact ? 1 : 0;
        int rowlen2 = m.nlo()+1;
        int gaplen2 = m.size()-rowlen2;
        const int si = m.stepi();
        const int sj = rm?1:m.stepj();
        const int sd = m.diagstep();
        const T* mrowi = m.cptr();
        int k=m.nlo();

        if (compact) 
            os << (sym?"sB ":"hB ")<< m.size()<<' '<<m.nlo() << std::endl;
        else
            os << m.size() <<' '<< m.size() << std::endl;

        while (rowlen2>0) {
            os << "( ";
            if (!compact) {
                for(int j=gaplen1;j>0;--j) os << ' '<<Value(T(0))<<' ';
            }

            const T* mij = mrowi;
            for(int j=rowlen1;j>0;--j,rm?++mij:mij+=sj) {
                if (conj)
                    if (th)
                        os << ' '<<Value(TMV_ABS(*mij)<thresh ? T(0) :
                                         TMV_CONJ(*mij))<<' ';
                    else
                        os << ' '<<Value(TMV_CONJ(*mij))<<' ';
                else
                    if (th)
                        os << ' '<<Value(TMV_ABS(*mij)<thresh ? T(0) : 
                                         *mij)<<' ';
                    else
                        os << ' '<<Value(*mij)<<' ';
            }

            if (!compact) {
                for(int j=rowlen2;j>0;--j,cm?++mij:mij+=si) {
                    if (sym == conj)
                        if (th)
                            os << ' '<<Value(TMV_ABS(*mij)<thresh ? T(0) :
                                             TMV_CONJ(*mij))<<' ';
                        else
                            os << ' '<<Value(TMV_CONJ(*mij))<<' ';
                    else
                        if (th)
                            os << ' '<<Value(TMV_ABS(*mij)<thresh ? T(0) :
                                             *mij)<<' ';
                        else
                            os << ' '<<Value(*mij)<<' ';
                }

                for(int j=gaplen2;j>0;--j) os << ' '<<Value(T(0))<<' ';
            }
            os << " )\n";

            if (k>0) { --k; ++rowlen1; mrowi += si; }
            else { if(!compact) ++gaplen1; mrowi += sd; }
            if (gaplen2 > 0) --gaplen2; 
            else --rowlen2;
        }
    }

    template <bool rm, bool cm, bool compact, bool th, class T> 
    static inline void DoWrite1(
        std::ostream& os, const GenSymBandMatrix<T>& m, T thresh)
    { DoWrite<true,false,rm,cm,compact,th>(os,m,thresh); }

    template <bool rm, bool cm, bool compact, bool th, class T> 
    static void DoWrite1(
        std::ostream& os, const GenSymBandMatrix<std::complex<T> >& m, T thresh)
    {
        if (m.issym())
            if (m.isconj())
                DoWrite<true,true,rm,cm,compact,th>(os,m,thresh);
            else
                DoWrite<true,false,rm,cm,compact,th>(os,m,thresh);
        else
            if (m.isconj())
                DoWrite<false,true,rm,cm,compact,th>(os,m,thresh);
            else
                DoWrite<false,false,rm,cm,compact,th>(os,m,thresh);
    }

    template <class T> 
    void GenSymBandMatrix<T>::write(std::ostream& os) const
    {
        if (uplo() == Upper) 
            if (issym()) transpose().write(os);
            else adjoint().write(os);
        else {
            if (isrm())
                DoWrite1<true,false,false,false>(os,*this,RT(0));
            else if (iscm()) 
                DoWrite1<false,true,false,false>(os,*this,RT(0));
            else
                DoWrite1<false,false,false,false>(os,*this,RT(0));
        }
    }

    template <class T> 
    void GenSymBandMatrix<T>::write(std::ostream& os, RT thresh) const
    {
        if (uplo() == Upper) 
            if (issym()) transpose().write(os,thresh);
            else adjoint().write(os,thresh);
        else 
            if (isrm())
                DoWrite1<true,false,false,true>(os,*this,thresh);
            else if (iscm())
                DoWrite1<false,true,false,true>(os,*this,thresh);
            else
                DoWrite1<false,false,false,true>(os,*this,thresh);
    }

    template <class T> 
    void GenSymBandMatrix<T>::writeCompact(std::ostream& os) const
    {
        if (uplo() == Upper) 
            if (issym()) transpose().writeCompact(os);
            else adjoint().writeCompact(os);
        else 
            if (isrm())
                DoWrite1<true,false,true,false>(os,*this,RT(0));
            else if (iscm())
                DoWrite1<false,true,true,false>(os,*this,RT(0));
            else
                DoWrite1<false,false,true,false>(os,*this,RT(0));
    }

    template <class T> 
    void GenSymBandMatrix<T>::writeCompact(std::ostream& os, RT thresh) const
    {
        if (uplo() == Upper) 
            if (issym()) transpose().writeCompact(os);
            else adjoint().writeCompact(os);
        else 
            if (isrm())
                DoWrite1<true,false,true,true>(os,*this,thresh);
            else if (iscm())
                DoWrite1<false,true,true,true>(os,*this,thresh);
            else
                DoWrite1<false,false,true,true>(os,*this,thresh);
    }

#ifndef NOTHROW
    template <class T> 
    class SymBandMatrixReadError : public ReadError
    {
    public :
        int i,j;
        mutable auto_ptr<SymBandMatrix<T> > m;
        std::string exp,got;
        size_t s;
        int lo;
        bool is, iseof, isbad;

        SymBandMatrixReadError(
            int _i, int _j, const GenSymBandMatrix<T>& _m,
            std::istream& _is) throw() :
            ReadError("SymBandMatrix."),
            i(_i), j(_j), m(new SymBandMatrix<T>(_m)), exp(""), got(""), 
            s(_m.size()), lo(_m.nlo()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        SymBandMatrixReadError(std::istream& _is) throw() :
            ReadError("SymBandMatrix."),
            i(0), j(0), m(0), exp(""), got(""), s(0), lo(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        SymBandMatrixReadError(
            int _i, int _j, const GenSymBandMatrix<T>& _m,
            std::istream& _is, std::string _e, std::string _g) throw() :
            ReadError("SymBandMatrix."),
            i(_i), j(_j), m(new SymBandMatrix<T>(_m)), 
            exp(_e), got(_g), s(_m.size()), lo(_m.nlo()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        SymBandMatrixReadError(std::istream& _is, 
                               std::string _e, std::string _g) throw() :
            ReadError("SymBandMatrix."),
            i(0), j(0), m(0), exp(_e), got(_g), s(0), lo(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        SymBandMatrixReadError(
            const GenSymBandMatrix<T>& _m, std::istream& _is, 
            size_t _s, int _lo) throw() :
            ReadError("SymBandMatrix."),
            i(0), j(0), m(new SymBandMatrix<T>(_m)), exp(""), got(""), 
            s(_s), lo(_lo),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        SymBandMatrixReadError(const SymBandMatrixReadError<T>& rhs) :
            i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got), 
            s(rhs.s), lo(rhs.lo),
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
            if (m.get() && s != m->size()) {
                os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
            }
            if (m.get() && lo != m->nlo()) {
                os<<"Write nlo: expected "<<m->nlo()<<", got "<<lo<<".\n";
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
            if (m.get()) {
                os<<"The portion of the SymBandMatrix which was successfully read is: \n";
                ConstSymBandMatrixView<T> mm = m->view();
                int k = mm.nlo();
                int rowlen = 1;
                int rowstart = 0;
                for(int ii=0;ii<i;++ii) {
                    os<<"( ";
                    for(int jj=0;jj<rowlen;++jj)
                        os<<' '<<mm(ii,jj+rowstart)<<' ';
                    os<<" )\n";
                    if (k > 0) { ++rowlen; --k; }
                    else ++rowstart;
                }
                os<<"( ";
                for(int jj=0;jj<j;++jj)
                    os<<' '<<mm(i,jj+rowstart)<<' ';      
                os<<" )\n";
            }
        }
    };

    template <class T> 
    class HermBandMatrixReadError : public ReadError
    {
    public :
        int i,j;
        mutable auto_ptr<HermBandMatrix<T> > m;
        std::string exp,got;
        size_t s;
        int lo;
        T dv;
        bool is, iseof, isbad;

        HermBandMatrixReadError(
            int _i, int _j, const GenSymBandMatrix<T>& _m,
            std::istream& _is) throw() :
            ReadError("HermBandMatrix."),
            i(_i), j(_j), m(new HermBandMatrix<T>(_m)), exp(""), got(""), 
            s(_m.size()), lo(_m.nlo()),
            dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermBandMatrixReadError(std::istream& _is) throw() :
            ReadError("HermBandMatrix."),
            i(0), j(0), m(0), exp(""), got(""), s(0), lo(0),
            dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermBandMatrixReadError(
            int _i, int _j, const GenSymBandMatrix<T>& _m,
            std::istream& _is, std::string _e, std::string _g) throw() :
            ReadError("HermBandMatrix."),
            i(_i), j(_j), m(new HermBandMatrix<T>(_m)), exp(_e), got(_g), 
            s(_m.size()), lo(_m.nlo()),
            dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermBandMatrixReadError(
            int _i, int _j, const GenSymBandMatrix<T>& _m,
            std::istream& _is, T _dv) throw() :
            ReadError("HermBandMatrix."),
            i(_i), j(_j), m(new HermBandMatrix<T>(_m)), exp(""), got(""), 
            s(_m.size()), lo(_m.nlo()),
            dv(_dv), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermBandMatrixReadError(std::istream& _is, 
                                std::string _e, std::string _g) throw() :
            ReadError("HermBandMatrix."),
            i(0), j(0), m(0), exp(_e), got(_g), s(0), lo(0),
            dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        HermBandMatrixReadError(
            const GenSymBandMatrix<T>& _m, std::istream& _is, 
            size_t _s, int _lo) throw() :
            ReadError("HermBandMatrix."),
            i(0), j(0), m(new HermBandMatrix<T>(_m)), exp(""), got(""), 
            s(_s), lo(_lo),
            dv(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        HermBandMatrixReadError(const HermBandMatrixReadError<T>& rhs) :
            i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
            s(rhs.s), lo(rhs.lo),
            dv(rhs.dv), is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}

        virtual ~HermBandMatrixReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for HermBandMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"'";
                if (isReal(T()) && exp == "hB") os<<" (or 'sB')";
                os<<", got '"<<got<<"'.\n";
            }
            if (dv != T(0)) {
                os<<"Non-real value found on diagonal: "<<dv<<std::endl;
            }
            if (m.get() && s != m->size()) {
                os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
            }
            if (m.get() && s != m->size()) {
                os<<"Wrong nlo: expected "<<m->nlo()<<", got "<<lo<<".\n";
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
            if (m.get()) {
                os<<"The portion of the HermBandMatrix which was successfully read is: \n";
                ConstSymBandMatrixView<T> mm = m->view();
                for(int ii=0;ii<i;++ii) {
                    os<<"( ";
                    for(int jj=0;jj<(ii<j?i+1:i);++jj)
                        os<<' '<<mm(ii,jj)<<' ';
                    os<<" )\n";
                }
                os<<"( ";
                for(int jj=0;jj<j;++jj)
                    os<<' '<<mm(i,jj)<<' ';      
                os<<" )\n";
            }
        }
    };
#endif

    template <class T, IndexStyle I> 
    void SymBandMatrixView<T,I>::read(std::istream& is) const
    {
        if (uplo() == Upper) 
            if (this->issym()) transpose().read(is);
            else adjoint().read(is);
        else {
            TMVAssert(uplo() == Lower);
            std::string paren;
            T* mrowi = ptr();
            int rowlen = 1;
            int k = nlo();
            const int si = stepi();
            const int sj = stepj();
            const int sd = si+sj;
            for(int i=size();i>0;--i) {
#ifdef NOSTL
                char p;
                is >> p;
                paren = p;
#else
                is >> paren;
#endif
                if (!is || paren != "(") {
                    if (this->issym()) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix ReadError: "<<paren<<" != (\n"; 
                        exit(1); 
#else
                        throw SymBandMatrixReadError<T>(
                            size()-i,0,*this,is,"(",is?paren:"(");
#endif
                    } else {
#ifdef NOTHROW
                        std::cerr<<"HermBandMatrix ReadError: "<<paren<<" != (\n"; 
                        exit(1); 
#else
                        throw HermBandMatrixReadError<T>(
                            size()-i,0,*this,is,"(",is?paren:"(");
#endif
                    }
                }
                T* mij = mrowi;
                if (this->isrm()) {
                    for(int j=rowlen;j>0;--j,++mij) {
                        is >> *mij;
                        if (!is) {
                            if (this->issym()) {
#ifdef NOTHROW
                                std::cerr<<"SymBandMatrix ReadError: !is \n"; 
                                exit(1); 
#else
                                throw SymBandMatrixReadError<T>(
                                    size()-i,rowlen-j,*this,is);
#endif
                            } else {
#ifdef NOTHROW
                                std::cerr<<"HermBandMatrix ReadError: !is \n"; 
                                exit(1); 
#else
                                throw HermBandMatrixReadError<T>(
                                    size()-i,rowlen-j,*this,is);
#endif
                            }
                        }
                        if (this->isherm() && j==1 && TMV_IMAG(*mij) != RT(0)) {
#ifdef NOTHROW
                            std::cerr<<"HermBandMatrix ReadError: "<<*mij<" not real\n"; 
                            exit(1); 
#else
                            throw HermBandMatrixReadError<T>(
                                size()-i,rowlen-j,*this,is,*mij);
#endif
                        }
                    }
                } else {
                    for(int j=rowlen;j>0;--j,mij+=sj) {
                        is >> *mij;
                        if (!is) {
                            if (this->issym()) {
#ifdef NOTHROW
                                std::cerr<<"SymBandMatrix ReadError: !is \n"; 
                                exit(1); 
#else
                                throw SymBandMatrixReadError<T>(size()-i,rowlen-j,*this,is);
#endif
                            } else {
#ifdef NOTHROW
                                std::cerr<<"HermBandMatrix ReadError: !is \n"; 
                                exit(1); 
#else
                                throw HermBandMatrixReadError<T>(size()-i,rowlen-j,*this,is);
#endif
                            }
                        }
                        if (this->isherm() && j==1 && TMV_IMAG(*mij) != RT(0)) {
#ifdef NOTHROW
                            std::cerr<<"HermBandMatrix ReadError: "<<*mij<<" not real \n"; 
                            exit(1); 
#else
                            throw HermBandMatrixReadError<T>(size()-i,rowlen-j,*this,is,*mij);
#endif
                        }
                    }
                }
#ifdef NOSTL
                is >> p;
                paren = p;
#else
                is >> paren;
#endif
                if ((!is && i>1)  || paren != ")") {
                    if (this->issym()) {
#ifdef NOTHROW
                        std::cerr<<"SymBandMatrix ReadError: "<<paren<<" != )\n"; 
                        exit(1); 
#else
                        throw SymBandMatrixReadError<T>(size()-i,size()-i+1,*this,is,
                                                        ")",is?paren:")");
#endif
                    } else {
#ifdef NOTHROW
                        std::cerr<<"HermBandMatrix ReadError: "<<paren<<" != ) \n"; 
                        exit(1); 
#else
                        throw HermBandMatrixReadError<T>(size()-i,size()-i+1,*this,is,
                                                         ")",is?paren:")");
#endif
                    }
                }
                if (k > 0) { ++rowlen; --k; mrowi+=si; }
                else { mrowi+=sd; }
            }
        }
        if (this->isconj()) conjugateSelf();
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    std::istream& operator>>(
        std::istream& is, auto_ptr<SymBandMatrix<T,U,S,I> >& m)
    {
        std::string sh;
#ifdef NOSTL
        sh = "  ";
        is >> sh[0];
        if (is) is >> sh[1];
#else
        is >> sh;
#endif
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix ReadError: !is \n"; 
            exit(1); 
#else
            throw SymBandMatrixReadError<T>(is);
#endif
        }
        if (isReal(T())) {
            if (sh != "sB" && sh != "hB")  {
#ifdef NOTHROW
                std::cerr<<"SymBandMatrix ReadError: "<<sh<<" != sB\n"; 
                exit(1); 
#else
                throw SymBandMatrixReadError<T>(is,"sB",sh);
#endif
            }
        } else {
            if (sh != "sB")  {
#ifdef NOTHROW
                std::cerr<<"SymBandMatrix ReadError: "<<sh<<" != sB\n"; 
                exit(1); 
#else
                throw SymBandMatrixReadError<T>(is,"sB",sh);
#endif
            }
        }
        size_t size;
        int nlo;
        is >> size >> nlo;
        if (!is)  {
#ifdef NOTHROW
            std::cerr<<"SymBandMatrix ReadError: !is \n"; 
            exit(1); 
#else
            throw SymBandMatrixReadError<T>(is);
#endif
        }
        m.reset(new SymBandMatrix<T,U,S,I>(size,nlo));
        m->view().read(is); 
        return is;
    }

    template <class T, UpLoType U, StorageType S, IndexStyle I> 
    std::istream& operator>>(
        std::istream& is, auto_ptr<HermBandMatrix<T,U,S,I> >& m)
    {
        std::string sh;
#ifdef NOSTL
        sh = "  ";
        is >> sh[0];
        if (is) is >> sh[1];
#else
        is >> sh;
#endif
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"HermBandMatrix ReadError: !is \n"; 
            exit(1); 
#else
            throw HermBandMatrixReadError<T>(is);
#endif
        }
        if (isReal(T())) {
            if (sh != "sB" && sh != "hB")  {
#ifdef NOTHROW
                std::cerr<<"HermBandMatrix ReadError: "<<sh<<" != hB\n"; 
                exit(1); 
#else
                throw HermBandMatrixReadError<T>(is,"hB",sh);
#endif
            }
        } else {
            if (sh != "hB")  {
#ifdef NOTHROW
                std::cerr<<"HermBandMatrix ReadError: "<<sh<<" != hB\n"; 
                exit(1); 
#else
                throw HermBandMatrixReadError<T>(is,"hB",sh);
#endif
            }
        }
        size_t size;
        int nlo;
        is >> size >> nlo;
        if (!is)  {
#ifdef NOTHROW
            std::cerr<<"HermBandMatrix ReadError: !is \n"; 
            exit(1); 
#else
            throw HermBandMatrixReadError<T>(is);
#endif
        }
        m.reset(new HermBandMatrix<T,U,S,I>(size,nlo));
        m->view().read(is); 
        return is;
    }

    template <class T> std::istream& operator>>(
        std::istream& is, const SymBandMatrixView<T>& m)
    {
        std::string sh;
#ifdef NOSTL
        sh = "  ";
        is >> sh[0];
        if (is) is >> sh[1];
#else
        is >> sh;
#endif
        if (!is) {
            if (m.issym()) {
#ifdef NOTHROW
                std::cerr<<"SymBandMatrix ReadError: !is \n"; 
                exit(1); 
#else
                throw SymBandMatrixReadError<T>(is);
#endif
            } else {
#ifdef NOTHROW
                std::cerr<<"HermBandMatrix ReadError: !is \n"; 
                exit(1); 
#else
                throw HermBandMatrixReadError<T>(is);
#endif
            }
        }
        if (isReal(T())) {
            if (sh != "sB" && sh != "hB") {
                if (m.issym()) {
#ifdef NOTHROW
                    std::cerr<<"SymBandMatrix ReadError: "<<sh<<" != sB\n"; 
                    exit(1); 
#else
                    throw SymBandMatrixReadError<T>(is,"sB",sh);
#endif
                } else {
#ifdef NOTHROW
                    std::cerr<<"HermBandMatrix ReadError: "<<sh<<" != hB\n"; 
                    exit(1); 
#else
                    throw HermBandMatrixReadError<T>(is,"hB",sh);
#endif
                }
            }
        } else if (m.issym()) {
            if (sh != "sB")  {
#ifdef NOTHROW
                std::cerr<<"SymBandMatrix ReadError: "<<sh<<" != sB\n"; 
                exit(1); 
#else
                throw SymBandMatrixReadError<T>(is,"sB",sh);
#endif
            }
        } else {
            if (sh != "hB")  {
#ifdef NOTHROW
                std::cerr<<"HermBandMatrix ReadError: "<<sh<<" != hB\n"; 
                exit(1); 
#else
                throw HermBandMatrixReadError<T>(is,"hB",sh);
#endif
            }
        }
        size_t s;
        int nlo;
        is >> s >> nlo;
        if (!is) {
            if (m.issym()) {
#ifdef NOTHROW
                std::cerr<<"SymBandMatrix ReadError: !is \n"; 
                exit(1); 
#else
                throw SymBandMatrixReadError<T>(is);
#endif
            } else {
#ifdef NOTHROW
                std::cerr<<"HermBandMatrix ReadError: !is \n"; 
                exit(1); 
#else
                throw HermBandMatrixReadError<T>(is);
#endif
            }
        }
        if (s != m.size() || nlo != m.nlo()) {
            if (m.issym()) {
#ifdef NOTHROW
                std::cerr<<"SymBandMatrix ReadError: Wrong size\n"; 
                exit(1); 
#else
                throw SymBandMatrixReadError<T>(m,is,s,nlo);
#endif
            } else {
#ifdef NOTHROW
                std::cerr<<"HermBandMatrix ReadError: Wrong size \n"; 
                exit(1); 
#else
                throw HermBandMatrixReadError<T>(m,is,s,nlo);
#endif
            }
        }
        TMVAssert(m.size() == s);
        TMVAssert(m.nlo() == nlo);
        m.read(is);
        return is;
    }

#undef RT
#undef CT

#define InstFile "TMV_SymBandMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


