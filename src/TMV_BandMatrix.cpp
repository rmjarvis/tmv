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
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_BandLUD.h"
#include "tmv/TMV_BandQRD.h"
#include "tmv/TMV_BandSVD.h"
#include "tmv/TMV_VIt.h"
#include "tmv/TMV_BandMatrixArith.h"
#include "tmv/TMV_DiagMatrix.h"
#include <iostream>
#include "portable_platform.h"

namespace tmv {

#define RT TMV_RealType(T)
#define CT TMV_ComplexType(T)

    // 
    // Access
    //

    template <class T> 
    T GenBandMatrix<T>::cref(int i, int j) const
    {
        const T* mi = cptr()+int(i)*stepi()+int(j)*stepj();
        return isconj() ? TMV_CONJ(*mi) : *mi;
    }

    template <class T, IndexStyle I> 
    TMV_RefType(T) BandMatrixView<T,I>::ref(int i, int j) const
    {
        T* mi = ptr()+int(i)*stepi()+int(j)*stepj();
        return TMV_REF(mi,ct());
    }

    template <class T> 
    void GenBandMatrix<T>::newDivider() const
    {
        switch (this->getDivType()) {
          case LU : 
               this->setDiv(new BandLUDiv<T>(
                       *this,this->isDivInPlace())); break;
          case QR : 
               this->setDiv(new BandQRDiv<T>(
                       *this,this->isDivInPlace())); break;
          case SV : 
               this->setDiv(new BandSVDiv<T>(*this)); break;
          default : TMVAssert(TMV_FALSE);
        }
    }

    size_t BandStorageLength(
        StorageType s, size_t cs, size_t rs, int lo, int hi)
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

    //
    // OK? (SubMatrix, etc.)
    //

    template <class T> 
    bool GenBandMatrix<T>::hasSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
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
            std::cout<<"row range ("<<j2-j1<<") must be multiple of jstep (";
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
        if (!okij(i1,j2-jstep)) {
            ok = false;
            std::cout<<"Upper right corner ("<<i1<<','<<j2-jstep;
            std::cout<<") must be in band\n";
        }
        if (!okij(i2-istep,j1)) {
            ok = false;
            std::cout<<"Lower left corner ("<<i2-istep<<','<<j1;
            std::cout<<") must be in band\n";
        }
        if (!okij(i2-istep,j2-jstep)) {
            ok = false;
            std::cout<<"Lower right corner ("<<i2-istep<<','<<j2-jstep;
            std::cout<<") must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool GenBandMatrix<T>::hasSubVector(
        int i, int j, int istep, int jstep, int size) const
    {
        if (size==0) return true;
        bool ok = true;
        if (istep == 0 && jstep == 0) {
            ok = false;
            std::cout<<"istep ("<<istep<<") and jstep ("<<jstep;
            std::cout<<") can not both be 0\n";
        }
        if (i<0 || i >= int(colsize())) {
            ok = false;
            std::cout<<"i ("<<i<<") must be in 0 -- "<<colsize()-1<<std::endl;
        }
        if (j<0 || j >= int(rowsize())) {
            ok = false;
            std::cout<<"j ("<<j<<") must be in 0 -- "<<rowsize()-1<<std::endl;
        }
        int i2 = int(i)+istep*int(size-1);
        int j2 = int(j)+jstep*int(size-1);
        if (i2 < 0 || i2 >= int(colsize())) {
            ok = false;
            std::cout<<"last element's i ("<<i2<<") must be in 0 -- ";
            std::cout<<colsize()-1<<std::endl;
        }
        if (j2 < 0 || j2 >= int(rowsize())) {
            ok = false;
            std::cout<<"last element's j ("<<j2<<") must be in 0 -- ";
            std::cout<<rowsize()-1<<std::endl;
        }
        if (!okij(i,j)) {
            ok = false;
            std::cout<<"First element ("<<i<<','<<j<<") must be in band\n";
        }
        if (!okij(i2,j2)) {
            ok = false;
            std::cout<<"Last element ("<<i2<<','<<j2<<") must be in band\n";
        }
        return ok;
    } 
    template <class T> 
    bool GenBandMatrix<T>::hasSubBandMatrix(
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
        return ok;
    } 

    template <class T> 
    bool ConstBandMatrixView<T,FortranStyle>::hasSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
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
        if (j1 < 0 || j1 >= int(this->rowsize())) {
            ok = false;
            std::cout<<"first row element ("<<j1<<") must be in 1 -- ";
            std::cout<<this->rowsize()<<std::endl;
        }
        if (j2 < 0 || j2 >= int(this->rowsize())) {
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
        if (!this->okij(i1-1,j2-1)) {
            ok = false;
            std::cout<<"Upper right corner ("<<i1<<','<<j2<<") must be in band\n";
        }
        if (!this->okij(i2-1,j1-1)) {
            ok = false;
            std::cout<<"Lower left corner ("<<i2<<','<<j1<<") must be in band\n";
        }
        if (!this->okij(i2-1,j2-1)) {
            ok = false;
            std::cout<<"Lower right corner ("<<i2<<','<<j2<<") must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstBandMatrixView<T,FortranStyle>::hasSubVector(
        int i, int j, int istep, int jstep, int size) const
    {
        if (size==0) return true;
        bool ok = true;
        if (istep == 0 && jstep == 0) {
            ok = false;
            std::cout<<"istep ("<<istep<<") and jstep ("<<jstep;
            std::cout<<") can not both be 0\n";
        }
        if (i < 1 || i > int(this->colsize())) {
            ok = false;
            std::cout<<"i ("<<i<<") must be in 1 -- "<<this->colsize()<<std::endl;
        }
        if (j < 1 || j > int(this->rowsize())) {
            ok = false;
            std::cout<<"j ("<<j<<") must be in 1 -- "<<this->rowsize()<<std::endl;
        }
        int i2 = int(i)+istep*int(size-1);
        int j2 = int(j)+jstep*int(size-1);
        if (i2 < 1 || i2 > int(this->colsize())) {
            ok = false;
            std::cout<<"last element's i ("<<i2<<") must be in 1 -- ";
            std::cout<<this->colsize()<<std::endl;
        }
        if (j2 < 1 || j2 > int(this->rowsize())) {
            ok = false;
            std::cout<<"last element's j ("<<j2<<") must be in 1 -- ";
            std::cout<<this->rowsize()<<std::endl;
        }
        if (!this->okij(i-1,j-1)) {
            ok = false;
            std::cout<<"First element ("<<i<<','<<j<<") must be in band\n";
        }
        if (!this->okij(i2-1,j2-1)) {
            ok = false;
            std::cout<<"Last element ("<<i2<<','<<j2<<") must be in band\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstBandMatrixView<T,FortranStyle>::hasSubBandMatrix(
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
        return ok;
    } 

    template <class T, IndexStyle I> 
    bool ConstBandMatrixView<T,I>::canLinearize() const
    {
        if (linsize == 1) {
            size_t rs = this->rowsize();
            size_t cs = this->colsize();
            if (rs > cs+this->nhi()) rs = cs+this->nhi();
            if (cs > rs+this->nlo()) cs = rs+this->nlo();
            if (rs != 1 || cs != 1) {
                if (rs == 0 || cs == 0) linsize = 0;
                else {
                    TMVAssert(this->isrm() || this->iscm());
                    if (this->isrm()) linsize = rs + (cs-1)*this->stepi();
                    else linsize = cs + (rs-1)*this->stepj();
                }
            }
            // If linsize is only 1, then no big savings in the linear version.
            // However, doing this segment over and over would be slow.
            // So set linsize to 0 to make canLinearize() faster when linsize
            // would be legitimately equal to 1.
            // Also, the Assert in BandMatrix.h about (linsize != 1 || 
            // rowsize == 1 && colsize == 1) can fail, since really it is
            // rs == 1 && cs == 1 which lead to linsize == 1.
            // Probably a better solution would be to have the guard be 
            // linsize == -1, but I haven't done that yet. 
            if (linsize == 1) linsize = 0;
        }
        return linsize > 0;
    }

    template <class T, IndexStyle I> 
    bool BandMatrixView<T,I>::canLinearize() const
    {
        if (linsize == 1) {
            size_t rs = this->rowsize();
            size_t cs = this->colsize();
            if (rs > cs+this->nhi()) rs = cs+this->nhi();
            if (cs > rs+this->nlo()) cs = rs+this->nlo();
            if (rs != 1 || cs != 1) {
                if (rs == 0 || cs == 0) linsize = 0;
                else {
                    TMVAssert(this->isrm() || this->iscm());
                    if (this->isrm()) linsize = rs + (cs-1)*this->stepi();
                    else linsize = cs + (rs-1)*this->stepj();
                }
            }
            if (linsize == 1) linsize = 0;
        }
        return linsize > 0;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::newCopy() const
    {
        if (isrm()) return auto_ptr<BaseMatrix<T> >(
            new BandMatrix<T,RowMajor>(*this));
        else if (iscm()) return auto_ptr<BaseMatrix<T> >(
            new BandMatrix<T,ColMajor>(*this));
        else return auto_ptr<BaseMatrix<T> >(
            new BandMatrix<T,DiagMajor>(*this));
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::newView() const
    {
        return auto_ptr<BaseMatrix<T> >(new ConstBandMatrixView<T>(view()));
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::newTranspose() const
    {
        return auto_ptr<BaseMatrix<T> >(new ConstBandMatrixView<T>(transpose()));
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::newConjugate() const
    {
        return auto_ptr<BaseMatrix<T> >(new ConstBandMatrixView<T>(conjugate()));
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::newAdjoint() const
    {
        return auto_ptr<BaseMatrix<T> >(new ConstBandMatrixView<T>(adjoint()));
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenBandMatrix<T>::newInverse() const
    {
        auto_ptr<Matrix<T,ColMajor> > minv(
            new Matrix<T,ColMajor>(rowsize(),colsize()));
        makeInverse(minv->view());
        BaseMatrix<T>* ret1 = minv.release();
        auto_ptr<BaseMatrix<T> > ret(ret1);
        return ret;
    }

    template <class T> QuotXB<T,T> GenBandMatrix<T>::QInverse() const
    { return QuotXB<T,T>(T(1),*this); }

    //
    // Norms
    //
    
    template <class T>
    T GenBandMatrix<T>::sumElements() const
    {
        const int M = colsize();
        const int N = rowsize();
        T sum = 0;
        if (M > 0 && N > 0) {
            if (isrm()) {
                int j1=0;
                int j2=nhi()+1;
                int k=nlo();
                for(int i=0;i<M;++i) {
                    sum += row(i,j1,j2).sumElements();
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (iscm()) {
                int i1=0;
                int i2=nlo()+1;
                int k=nhi();
                for(int j=0;j<N;++j) {
                    sum += col(j,i1,i2).sumElements();
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(int i=-nlo();i<=nhi();++i) sum += diag(i).sumElements();
            }
        }
        return sum;
    }

    template <class T>
    RT GenBandMatrix<T>::sumAbsElements() const
    {
        const int M = colsize();
        const int N = rowsize();
        RT sum = 0;
        if (M > 0 && N > 0) {
            if (isrm()) {
                int j1=0;
                int j2=nhi()+1;
                int k=nlo();
                for(int i=0;i<M;++i) {
                    sum += row(i,j1,j2).sumAbsElements();
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (iscm()) {
                int i1=0;
                int i2=nlo()+1;
                int k=nhi();
                for(int j=0;j<N;++j) {
                    sum += col(j,i1,i2).sumAbsElements();
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(int i=-nlo();i<=nhi();++i) sum += diag(i).sumAbsElements();
            }
        }
        return sum;
    }
 
    template <class T> 
    RT GenBandMatrix<T>::normSq(RT scale) const
    {
        const int M = colsize();
        const int N = rowsize();
        RT sum = 0;
        if (M > 0 && N > 0) {
            if (isrm()) {
                int j1=0;
                int j2=nhi()+1;
                int k=nlo();
                for(int i=0;i<M;++i) {
                    sum += row(i,j1,j2).normSq(scale);
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (iscm()) {
                int i1=0;
                int i2=nlo()+1;
                int k=nhi();
                for(int j=0;j<N;++j) {
                    sum += col(j,i1,i2).normSq(scale);
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(int i=-nlo();i<=nhi();++i) sum += diag(i).normSq(scale);
            }
        }
        return sum;
    }

    template <class T> 
    static RT NonLapMaxAbsElement(const GenBandMatrix<T>& m)
    {
        RT max = 0;
        const int M = m.colsize();
        const int N = m.rowsize();
        if (M > 0 && N > 0) {
            if (m.isrm()) {
                int j1=0;
                int j2=m.nhi()+1;
                int k=m.nlo();
                for(int i=0;i<M;++i) {
                    RT temp = m.row(i,j1,j2).maxAbsElement();
                    if (temp > max) max = temp;
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (m.iscm()) {
                int i1=0;
                int i2=m.nlo()+1;
                int k=m.nhi();
                for(int j=0;j<N;++j) {
                    RT temp = m.col(j,i1,i2).maxAbsElement();
                    if (temp > max) max = temp;
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(int i=-m.nlo();i<=m.nhi();++i) {
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
        const int M = m.colsize();
        const int N = m.rowsize();
        if (M > 0 && N > 0) {
            if (m.isrm()) {
                int j1=0;
                int j2=m.nhi()+1;
                int k=m.nlo();
                for(int i=0;i<M;++i) {
                    RT temp = m.row(i,j1,j2).maxAbs2Element();
                    if (temp > max) max = temp;
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (m.iscm()) {
                int i1=0;
                int i2=m.nlo()+1;
                int k=m.nhi();
                for(int j=0;j<N;++j) {
                    RT temp = m.col(j,i1,i2).maxAbs2Element();
                    if (temp > max) max = temp;
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(int i=-m.nlo();i<=m.nhi();++i) {
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
        const int M = m.colsize();
        const int N = m.rowsize();
        RT max = 0;
        if (M > 0 && N > 0) {
            int i1=0;
            int i2=m.nlo()+1;
            int k=m.nhi();
            for(int j=0;j<N;++j) {
                RT temp = m.col(j,i1,i2).norm1();
                if (temp > max) max = temp;
                if (k>0) --k; else ++i1;
                if (i2<M) ++i2;
                else if (i1==M) break;
            }
        }
        return max;
    }

#ifdef INST_INT
    static int NonLapNormF(const GenBandMatrix<int>& m)
    { return TMV_SQRT(m.normSq()); }
    static int NonLapNormF(const GenBandMatrix<std::complex<int> >& m)
    { return TMV_SQRT(m.normSq()); }
#endif

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
        if (m.iscm()) {
            int kl = m.nlo();
            int ku = m.nhi();
            int lda = m.stepj()+1;
#ifndef LAPNOWORK
            AlignedArray<double> work(c == 'I' ? n : 0);
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
        if (m.iscm()) {
            int kl = m.nlo();
            int ku = m.nhi();
            int lda = m.stepj()+1;
#ifndef LAPNOWORK
            AlignedArray<double> work(c == 'I' ? n : 0);
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
        if (m.iscm()) {
            int kl = m.nlo();
            int ku = m.nhi();
            int lda = m.stepj()+1;
#ifndef LAPNOWORK
            AlignedArray<float> work(c == 'I' ? n : 0);
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
        if (m.iscm()) {
            int kl = m.nlo();
            int ku = m.nhi();
            int lda = m.stepj()+1;
#ifndef LAPNOWORK
            AlignedArray<float> work(c == 'I' ? n : 0);
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

    template <class T> 
    RT GenBandMatrix<T>::maxAbsElement() const
    {
#ifdef XLAP
        if (isrm() && this->isSquare()) return LapNorm('M',transpose());
        else if (iscm() && this->isSquare()) return LapNorm('M',*this);
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
        else if (isrm() && this->isSquare()) return LapNorm('M',transpose());
        else if (iscm() && this->isSquare()) return LapNorm('M',*this);
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
        if (isrm() && this->isSquare()) return LapNorm('I',transpose());
        else if (iscm() && this->isSquare()) return LapNorm('1',*this);
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
        if (isrm() && this->isSquare()) return LapNorm('F',transpose());
        else if (iscm() && this->isSquare()) return LapNorm('F',*this);
        else if (isdm() && this->isSquare() && nlo()==1 && nhi()==1)
            return LapNorm('F',*this);
        else
#endif
            return NonLapNormF(*this);
    }

    template <class T> 
    RT GenBandMatrix<T>::doNorm2() const
    {
        if (this->colsize() < this->rowsize()) return transpose().doNorm2();
        if (this->rowsize() == 0) return RT(0);
        DiagMatrix<RT> S(this->rowsize());
        SV_Decompose(*this,S.view());
        return S(0);
    }

    template <class T> 
    RT GenBandMatrix<T>::doCondition() const
    {
        if (this->colsize() < this->rowsize()) return transpose().doNorm2();
        if (this->rowsize() == 0) return RT(1);
        DiagMatrix<RT> S(this->rowsize());
        SV_Decompose(*this,S.view());
        return S(0)/S(S.size()-1);
    }

    //
    // Modifying Functions
    //

    template <class T, IndexStyle I>
    const BandMatrixView<T,I>& BandMatrixView<T,I>::clip(RT thresh) const
    {
        if (this->canLinearize()) linearView().clip(thresh);
        else {
            const int M = colsize();
            const int N = rowsize();
            if (M > 0 && N > 0) {
                if (isrm()) {
                    int j1=0;
                    int j2=nhi()+1;
                    int k=nlo();
                    for(int i=0;i<M;++i) {
                        row(i,j1,j2).clip(thresh);
                        if (k>0) --k; else ++j1;
                        if (j2<N) ++j2;
                        else if (j1==N) break;
                    }
                } else if (iscm()) {
                    int i1=0;
                    int i2=nlo()+1;
                    int k=nhi();
                    for(int j=0;j<N;++j) {
                        col(j,i1,i2).clip(thresh);
                        if (k>0) --k; else ++i1;
                        if (i2<M) ++i2;
                        else if (i1==M) break;
                    }
                } else {
                    for(int i=-nlo();i<=nhi();++i) diag(i).clip(thresh);
                }
            }
        }
        return *this;
    }

    template <class T, IndexStyle I> 
    const BandMatrixView<T,I>& BandMatrixView<T,I>::setZero() const
    {
        if (this->canLinearize()) linearView().setZero();
        else {
            const int M = colsize();
            const int N = rowsize();
            if (M > 0 && N > 0) {
                if (isrm()) {
                    int j1=0;
                    int j2=nhi()+1;
                    int k=nlo();
                    for(int i=0;i<M;++i) {
                        row(i,j1,j2).setZero();
                        if (k>0) --k; else ++j1;
                        if (j2<N) ++j2;
                        else if (j1==N) break;
                    }
                } else if (iscm()) {
                    int i1=0;
                    int i2=nlo()+1;
                    int k=nhi();
                    for(int j=0;j<N;++j) {
                        col(j,i1,i2).setZero();
                        if (k>0) --k; else ++i1;
                        if (i2<M) ++i2;
                        else if (i1==M) break;
                    }
                } else {
                    for(int i=-nlo();i<=nhi();++i) diag(i).setZero();
                }
            }
        }
        return *this;
    }

    template <class T, IndexStyle I> 
    const BandMatrixView<T,I>& BandMatrixView<T,I>::setAllTo(const T& x) const
    {
        if (this->canLinearize()) linearView().setAllTo(x);
        else {
            const int M = colsize();
            const int N = rowsize();
            if (M > 0 && N > 0) {
                if (isrm()) {
                    int j1=0;
                    int j2=nhi()+1;
                    int k=nlo();
                    for(int i=0;i<M;++i) {
                        row(i,j1,j2).setAllTo(x);
                        if (k>0) --k; else ++j1;
                        if (j2<N) ++j2;
                        else if (j1==N) break;
                    }
                } else if (iscm()) {
                    int i1=0;
                    int i2=nlo()+1;
                    int k=nhi();
                    for(int j=0;j<N;++j) {
                        col(j,i1,i2).setAllTo(x);
                        if (k>0) --k; else ++i1;
                        if (i2<M) ++i2;
                        else if (i1==M) break;
                    }
                } else {
                    for(int i=-nlo();i<=nhi();++i) diag(i).setAllTo(x);
                }
            }
        }
        return *this;
    }

    template <class T, IndexStyle I> 
    void BandMatrixView<T,I>::doTransposeSelf() const
    {
        TMVAssert(colsize() == rowsize());
        TMVAssert(nlo() == nhi());
        for(int i=1;i<=nhi();++i) Swap(diag(-i),diag(i));
    }

    template <class T, IndexStyle I> 
    const BandMatrixView<T,I>& BandMatrixView<T,I>::conjugateSelf() const
    {
        if (this->canLinearize()) linearView().conjugateSelf();
        else {
            const int M = colsize();
            const int N = rowsize();
            if (isComplex(T()) && M > 0 && N > 0) {
                if (isrm()) {
                    int j1=0;
                    int j2=nhi()+1;
                    int k=nlo();
                    for(int i=0;i<M;++i) {
                        row(i,j1,j2).conjugateSelf();
                        if (k>0) --k; else ++j1;
                        if (j2<N) ++j2;
                        else if (j1==N) break;
                    }
                } else if (iscm()) {
                    int i1=0;
                    int i2=nlo()+1;
                    int k=nhi();
                    for(int j=0;j<N;++j) {
                        col(j,i1,i2).conjugateSelf();
                        if (k>0) --k; else ++i1;
                        if (i2<M) ++i2;
                        else if (i1==M) break;
                    }
                } else {
                    for(int i=-nlo();i<=nhi();++i) diag(i).conjugateSelf();
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

    template <class T> 
    ConstBandMatrixView<T> BandMatrixViewOf(
        const T* m, size_t cs, size_t rs, int nlo, int nhi,
        StorageType stor)
    {
        TMVAssert2(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
        TMVAssert2(cs>0);
        TMVAssert2(rs>0);
        TMVAssert2(nlo<int(cs));
        TMVAssert2(nhi<int(rs));
        if (stor == DiagMajor) {
            int stepi = rs >= cs ? -int(cs)+1 : -int(rs);
            int stepj = rs >= cs ? int(cs) : int(rs)+1;
            return ConstBandMatrixView<T>(
                m-nlo*stepi,cs,rs,nlo,nhi,stepi,stepj,1,
                DiagMajor,NonConj);
        } else {
            int lohi = nlo+nhi;
            if (stor == RowMajor)
                return ConstBandMatrixView<T>(
                    m,cs,rs,nlo,nhi,lohi,1,lohi+1,RowMajor,NonConj);
            else 
                return ConstBandMatrixView<T>(
                    m,cs,rs,nlo,nhi,1,lohi,lohi+1,ColMajor,NonConj);
        }
    }

    template <class T> 
    BandMatrixView<T> BandMatrixViewOf(
        T* m, size_t cs, size_t rs, int nlo, int nhi, StorageType stor)
    {
        TMVAssert2(stor==RowMajor || stor==ColMajor || stor==DiagMajor);
        TMVAssert2(cs>0);
        TMVAssert2(rs>0);
        TMVAssert2(nlo<int(cs));
        TMVAssert2(nhi<int(rs));
        if (stor == DiagMajor) {
            int stepi = rs >= cs ? -int(cs)+1 : -int(rs);
            int stepj = rs >= cs ? int(cs) : int(rs)+1;
            return BandMatrixView<T>(
                m-nlo*stepi,cs,rs,nlo,nhi,stepi,stepj,1,
                DiagMajor,NonConj
                TMV_FIRSTLAST1(
                    m,m+BandStorageLength(DiagMajor,cs,rs,nlo,nhi)));
        } else {
            int lohi = nlo+nhi;
            if (stor == RowMajor)
                return BandMatrixView<T>(
                    m,cs,rs,nlo,nhi,lohi,1,lohi+1,RowMajor,NonConj
                    TMV_FIRSTLAST1(
                        m,m+BandStorageLength(RowMajor,cs,rs,nlo,nhi)));
            else
                return BandMatrixView<T>(
                    m,cs,rs,nlo,nhi,1,lohi,lohi+1,ColMajor,NonConj
                    TMV_FIRSTLAST1(
                        m,m+BandStorageLength(ColMajor,cs,rs,nlo,nhi)));
        }
    }



    //
    // Swap
    //

    template <class T> 
    void Swap(const BandMatrixView<T>& m1, const BandMatrixView<T>& m2)
    {
        TMVAssert2(m1.colsize() == m2.colsize());
        TMVAssert2(m1.rowsize() == m2.rowsize());
        TMVAssert2(m1.nlo() == m2.nlo());
        TMVAssert2(m1.nhi() == m2.nhi());
        const int M = m1.colsize();
        const int N = m1.rowsize();
        if (M > 0 && N > 0) {
            if (m1.isrm() && m2.isrm()) {
                int j1=0;
                int j2=m1.nhi()+1;
                int k=m1.nlo();
                for(int i=0;i<M;++i) {
                    Swap(m1.row(i,j1,j2),m2.row(i,j1,j2));
                    if (k>0) --k; else ++j1;
                    if (j2<N) ++j2;
                    else if (j1==N) break;
                }
            } else if (m1.iscm() && m2.iscm()) {
                int i1=0;
                int i2=m1.nlo()+1;
                int k=m1.nhi();
                for(int j=0;j<N;++j) {
                    Swap(m1.col(j,i1,i2),m2.col(j,i1,i2));
                    if (k>0) --k; else ++i1;
                    if (i2<M) ++i2;
                    else if (i1==M) break;
                }
            } else {
                for(int i=-m1.nlo();i<=m1.nhi();++i) 
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
        else if (m1.nlo() != m2.nlo()) return false;
        else if (m1.nhi() != m2.nhi()) return false;
        else if (m1.isSameAs(m2)) return true;
        else {
            for(int i=-m1.nlo();i<=m1.nhi();++i) 
                if (m1.diag(i) != m2.diag(i)) return false;
            return true;
        }
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

    template <bool conj, bool rm, bool compact, bool th, class T> 
    static void DoWrite(
        std::ostream& os, const GenBandMatrix<T>& m, RT thresh)
    {
        int j1=0;
        int len=m.nhi()+1;
        const T* mrowi = m.cptr();
        const int ds = m.diagstep();
        const int si = m.stepi();
        const int sj = m.stepj();
        const int M = m.colsize();
        const int N = m.rowsize();
        int k = m.nlo();

        if (compact) {
            os << "B "<<M<<' '<<N;
            os <<' '<<m.nlo()<<' '<<m.nhi()<<std::endl;
        }
        else 
            os << M<<' '<<N<<std::endl;

        for(int i=M;i>0;--i) {
            os << "( ";
            if (!compact)
                for(int j=j1;j>0;--j) os << ' '<<Value(T(0))<<' ';

            const T* mij = mrowi;
            for(int j=len;j>0;--j,rm?++mij:mij+=sj) 
                if (conj) 
                    if (th)
                        os << ' ' <<Value(TMV_ABS(*mij)<thresh ? 
                                          T(0) : TMV_CONJ(*mij)) << ' ';
                    else
                        os << ' ' <<Value(TMV_CONJ(*mij)) << ' ';
                else 
                    if (th)
                        os << ' ' <<Value(TMV_ABS(*mij)<thresh ?
                                          T(0) : *mij) << ' ';
                    else
                        os << ' ' <<Value(*mij) << ' ';

            if (!compact)
                for(int j=N-len-j1;j>0;--j) os << ' '<<Value(T(0))<<' ';
            os << " )\n";
            if (k>0) { --k; ++len; mrowi+=si; } 
            else {
                mrowi+=ds; 
                if (j1<N) ++j1; 
            }
            if (j1+len>N) --len;
        }
    }

    template <bool rm, bool compact, bool th, class T> 
    static inline void DoWrite1(
        std::ostream& os, const GenBandMatrix<T>& m, T thresh)
    { DoWrite<false,rm,compact,th>(os,m,thresh); }

    template <bool rm, bool compact, bool th, class T> 
    static inline void DoWrite1(
        std::ostream& os, const GenBandMatrix<std::complex<T> >& m, T thresh)
    {
        if (m.isconj())
            DoWrite<true,rm,compact,th>(os,m,thresh);
        else
            DoWrite<false,rm,compact,th>(os,m,thresh);
    }

    template <class T> 
    void GenBandMatrix<T>::write(std::ostream& os) const
    {
        if (isrm())
            DoWrite1<true,false,false>(os,*this,RT(0));
        else
            DoWrite1<false,false,false>(os,*this,RT(0));
    }

    template <class T> 
    void GenBandMatrix<T>::write(std::ostream& os, RT thresh) const
    {
        if (isrm())
            DoWrite1<true,false,true>(os,*this,thresh);
        else
            DoWrite1<false,false,true>(os,*this,thresh);
    }

    template <class T> 
    void GenBandMatrix<T>::writeCompact(std::ostream& os) const
    {
        if (isrm())
            DoWrite1<true,true,false>(os,*this,RT(0));
        else
            DoWrite1<false,true,false>(os,*this,RT(0));
    }

    template <class T> 
    void GenBandMatrix<T>::writeCompact(std::ostream& os, RT thresh) const
    {
        if (isrm())
            DoWrite1<true,true,true>(os,*this,thresh);
        else
            DoWrite1<false,true,true>(os,*this,thresh);
    }

#ifndef NOTHROW
    template <class T> 
    class BandMatrixReadError : public ReadError
    {
    public :
        int i,j;
        mutable auto_ptr<BandMatrix<T> > m;
        char exp,got;
        size_t cs,rs;
        int lo,hi;
        bool is, iseof, isbad;

        BandMatrixReadError(
            int _i, int _j, const GenBandMatrix<T>& _m, std::istream& _is
        ) throw() :
            ReadError("BandMatrix."),
            i(_i), j(_j), m(new BandMatrix<T>(_m)), exp(0), got(0),
            cs(_m.colsize()), rs(_m.rowsize()), lo(_m.nlo()), hi(_m.nhi()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        BandMatrixReadError(std::istream& _is) throw() :
            ReadError("BandMatrix."),
            i(0), j(0), m(0), exp(0), got(0),
            cs(0), rs(0), lo(0), hi(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        BandMatrixReadError(
            int _i, int _j, const GenBandMatrix<T>& _m,
            std::istream& _is, char _e, char _g
        ) throw() :
            ReadError("BandMatrix."),
            i(_i), j(_j), m(new BandMatrix<T>(_m)), exp(_e), got(_g),
            cs(_m.colsize()), rs(_m.rowsize()), lo(_m.nlo()), hi(_m.nhi()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        BandMatrixReadError(std::istream& _is, char _e, char _g) throw() :
            ReadError("BandMatrix."),
            i(0), j(0), m(0), exp(_e), got(_g),
            cs(0), rs(0), lo(0), hi(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        BandMatrixReadError(
            const GenBandMatrix<T>& _m,
            std::istream& _is, size_t _cs, size_t _rs, int _lo, int _hi
        ) throw() :
            ReadError("BandMatrix."),
            i(0), j(0), m(new BandMatrix<T>(_m)), exp(0), got(0),
            cs(_cs), rs(_rs), lo(_lo), hi(_hi),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        BandMatrixReadError(const BandMatrixReadError<T>& rhs) :
            i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
            cs(rhs.cs), rs(rhs.rs), lo(rhs.lo), hi(rhs.hi),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        virtual ~BandMatrixReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for BandMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (m.get() && cs != m->colsize()) {
                os<<"Wrong colsize: expected "<<m->colsize()<<
                    ", got "<<cs<<".\n";
            }
            if (m.get() && rs != m->rowsize()) {
                os<<"Wrong rowsize: expected "<<m->rowsize()<<
                    ", got "<<rs<<".\n";
            }
            if (m.get() && lo != m->nlo()) {
                os<<"Wrong nlo: expected "<<m->nlo()<<", got "<<lo<<".\n";
            }
            if (m.get() && hi != m->nhi()) {
                os<<"Wrong nhi: expected "<<m->nhi()<<", got "<<hi<<".\n";
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
                os<<"The portion of the Bandatrix which was successfully "
                    "read is: \n";
                ConstBandMatrixView<T> mm = m->view();
                for(int ii=0;ii<i;++ii) {
                    const int N = mm.rowsize();
                    os<<"( ";
                    for(int jj=0;jj<N;++jj)
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
    void BandMatrixView<T,I>::read(std::istream& is) const
    {
        char paren;
        int j1=0;
        int len=nhi()+1;
        T* mrowi = ptr();
        const int ds = diagstep();
        const int si = stepi();
        const int sj = stepj();
        const int M = colsize();
        const int N = rowsize();
        int k = nlo();
        for(int i=M;i>0;--i) {
            is >> paren;
            if (!is || paren != '(') {
#ifdef NOTHROW
                std::cerr<<"Band Matrix Read Error "<<paren<<" != (\n"; 
                exit(1); 
#else
                throw BandMatrixReadError<T>(M-i,0,*this,is,'(',is?paren:'(');
#endif
            }
            T* mij = mrowi;
            if (isrm()) 
                for(int j=len;j>0;--j,++mij) {
                    is >> *mij;
                    if (!is) {
#ifdef NOTHROW
                        std::cerr<<"Band Matrix Read Error !is \n"; 
                        exit(1); 
#else
                        throw BandMatrixReadError<T>(M-i,j1+len-j,*this,is);
#endif
                    }
                }
            else 
                for(int j=len;j>0;--j,mij+=sj) {
                    is >> *mij;
                    if (!is) {
#ifdef NOTHROW
                        std::cerr<<"Band Matrix Read Error !is \n"; 
                        exit(1); 
#else
                        throw BandMatrixReadError<T>(M-i,j1+len-j,*this,is);
#endif
                    }
                }
            is >> paren;
            if ((!is && i>1) || paren != ')') {
#ifdef NOTHROW
                std::cerr<<"Band Matrix Read Error "<<paren<<" != )\n"; 
                exit(1); 
#else
                throw BandMatrixReadError<T>(M-i,N,*this,is,')',is?paren:')');
#endif
            }
            if (k>0) { --k; ++len; mrowi+=si; } 
            else {
                mrowi+=ds; 
                if (j1<N) ++j1; 
            }
            if (j1+len>N) --len;
        }
        if (isconj()) conjugateSelf();
    }

    template <class T, StorageType S, IndexStyle I> 
    std::istream& operator>>(
        std::istream& is, auto_ptr<BandMatrix<T,S,I> >& m)
    {
        char b;
        is >> b;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"Band Matrix Read Error !is \n"; 
            exit(1); 
#else
            throw BandMatrixReadError<T>(is);
#endif
        }
        if (b != 'B') {
#ifdef NOTHROW
            std::cerr<<"Band Matrix Read Error "<<b<<" != B\n"; 
            exit(1); 
#else
            throw BandMatrixReadError<T>(is,'B',b);
#endif
        }
        size_t cs,rs;
        int lo,hi;
        is >> cs >> rs >> lo >> hi;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"Band Matrix Read Error !is \n"; 
            exit(1); 
#else
            throw BandMatrixReadError<T>(is);
#endif
        }
        m.reset(new BandMatrix<T,S,I>(cs,rs,lo,hi));
        m->view().read(is); 
        return is;
    }

    template <class T> 
    std::istream& operator>>(std::istream& is, const BandMatrixView<T>& m)
    {
        char b;
        is >> b;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"Band Matrix Read Error !is \n"; 
            exit(1); 
#else
            throw BandMatrixReadError<T>(is);
#endif
        }
        if (b != 'B') {
#ifdef NOTHROW
            std::cerr<<"Band Matrix Read Error "<<b<<" != B\n"; 
            exit(1); 
#else
            throw BandMatrixReadError<T>(is,'B',b);
#endif
        }
        size_t cs,rs;
        int lo,hi;
        is >> cs >> rs >> lo >> hi;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"Band Matrix Read Error !is \n"; 
            exit(1); 
#else
            throw BandMatrixReadError<T>(is);
#endif
        }
        if (cs != m.colsize() || rs != m.rowsize() ||
            lo != m.nlo() || hi != m.nhi()) {
#ifdef NOTHROW
            std::cerr<<"Band Matrix Read Error wrong size \n"; 
            exit(1); 
#else
            throw BandMatrixReadError<T>(m,is,cs,rs,lo,hi);
#endif
        }
        TMVAssert(m.colsize() == cs);
        TMVAssert(m.rowsize() == rs);
        TMVAssert(m.nlo() == lo);
        TMVAssert(m.nhi() == hi);
        m.read(is);
        return is;
    }

#undef RT
#undef CT

#define InstFile "TMV_BandMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


