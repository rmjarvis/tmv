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
#include "portable_platform.h"

namespace tmv {

#define RT TMV_RealType(T)

#ifdef TMV_BLOCKSIZE
#define PERM_BLOCKSIZE TMV_BLOCKSIZE/2
#else
#define PERM_BLOCKSIZE 32
#endif

    //
    // Constructors
    //

#define NEW_SIZE(cs,rs) \
    linsize((cs)*(rs)), \
    itsm(linsize), itscs(cs), itsrs(rs) \
    TMV_DEFFIRSTLAST(itsm.get(),itsm.get()+ls())

    template <class T, StorageType S, IndexStyle I>
    Matrix<T,S,I>::Matrix(const std::vector<std::vector<T> >& vv) :
        NEW_SIZE(vv.size(),(vv.size()>0?vv[0].size():0))
    {
        TMVAssert(S==RowMajor || S==ColMajor);
        T* vi=itsm.get();
        if (S == RowMajor) {
            const int M = colsize();
            const int N = rowsize();
            for(int i=0;i<M;++i) {
                TMVAssert(vv[i].size() == rowsize());
                typename std::vector<T>::const_iterator vvi = vv[i].begin();
                for(int j=0;j<N;++j,++vi,++vvi) {
#ifdef TMVFLDEBUG
                    TMVAssert(vi >= _first);
                    TMVAssert(vi < _last);
#endif
                    *vi = *vvi;
                }
            }
        } else {
            const int M = colsize();
            const int N = rowsize();
            for(int j=0;j<N;++j) {
                for(int i=0;i<M;++i,++vi) {
                    TMVAssert(vv[i].size() == rowsize());
#ifdef TMVFLDEBUG
                    TMVAssert(vi >= _first);
                    TMVAssert(vi < _last);
#endif
                    *vi = vv[i][j];
                }
            }
        }
    }

#undef NEW_SIZE


    //
    // Access
    //

    template <class T>
    T GenMatrix<T>::cref(int i, int j) const
    {
        const T* mi = cptr() + int(i)*stepi() + int(j)*stepj();
        return isconj() ? TMV_CONJ(*mi) : *mi;
    }

    template <class T, IndexStyle I>
    typename MatrixView<T,I>::reference MatrixView<T,I>::ref(int i, int j) const
    {
        T* mi = ptr() + int(i)*itssi + int(j)*stepj();
        return TMV_REF(mi,ct());
    }

    template <class T>
    void GenMatrix<T>::newDivider() const
    {
        switch (this->getDivType()) {
          case LU : this->setDiv(
                        new LUDiv<T>(*this,this->isDivInPlace())); 
                    break;
          case QR : this->setDiv(
                        new QRDiv<T>(*this,this->isDivInPlace())); 
                    break;
          case QRP : this->setDiv(
                         new QRPDiv<T>(*this,this->isDivInPlace())); 
                     break;
          case SV : this->setDiv(
                        new SVDiv<T>(*this,this->isDivInPlace())); 
                    break;
          default : TMVAssert(TMV_FALSE);
        }
    }

#ifdef INST_INT
    template <>
    void GenMatrix<int>::newDivider() const
    { TMVAssert(TMV_FALSE); }
    template <>
    void GenMatrix<std::complex<int> >::newDivider() const
    { TMVAssert(TMV_FALSE); }
#endif

    // Note: These need to be in the .cpp file, not the .h file for
    // dynamic libraries.  Apparently the typeinfo used for dynamic_cast
    // doesn't get shread correctly by different modules, so the 
    // dynamic_cast fails when called in one module for an object that 
    // was created in a different module.
    //
    // So putting these functions here puts the dynamic cast in the shared
    // library, which is also where it is created (by newDivider above).
    template <class T>
    bool GenMatrix<T>::divIsLUDiv() const
    { return static_cast<bool>(dynamic_cast<const LUDiv<T>*>(getDiv())); }

    template <class T>
    bool GenMatrix<T>::divIsQRDiv() const
    { return static_cast<bool>(dynamic_cast<const QRDiv<T>*>(getDiv())); }

    template <class T>
    bool GenMatrix<T>::divIsQRPDiv() const
    { return static_cast<bool>(dynamic_cast<const QRPDiv<T>*>(getDiv())); }

    template <class T>
    bool GenMatrix<T>::divIsSVDiv() const
    { return static_cast<bool>(dynamic_cast<const SVDiv<T>*>(getDiv())); }


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
            std::cout<<"row range ("<<j2-j1<<") must be multiple of istep (";
            std::cout<<jstep<<")\n";
        }
        if ((j2-j1)/jstep < 0) {
            ok = false;
            std::cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
        }
        return ok;
    }

    template <class T>
    bool GenMatrix<T>::hasSubVector(
        int i, int j, int istep, int jstep, int size) const 
    {
        if (size==0) return true;
        bool ok = true;
        if (istep == 0 && jstep == 0) {
            ok = false;
            std::cout<<"istep ("<<istep<<") and jstep ("<<jstep;
            std::cout<<") can not both be 0\n";
        }
        if (i < 0 || i >= int(colsize())) {
            ok = false;
            std::cout<<"i ("<<i<<") must be in 0 -- "<<colsize()-1<<std::endl;
        }
        if (j < 0 || j >= int(rowsize())) {
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
        return ok;
    }

    template <class T>
    bool ConstMatrixView<T,FortranStyle>::hasSubMatrix(
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
        return ok;
    }

    template <class T>
    bool ConstMatrixView<T,FortranStyle>::hasSubVector(
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
                const int N = rowsize();
                for(int j=0;j<N;++j) sum += col(j).sumElements();
            } else {
                const int M = colsize();
                for(int i=0;i<M;++i) sum += row(i).sumElements();
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
                const int N = m.rowsize();
                for(int j=0;j<N;++j) sum += m.col(j).sumAbsElements();
            } else {
                const int M = m.colsize();
                for(int i=0;i<M;++i) sum += m.row(i).sumAbsElements();
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
                const int N = m.rowsize();
                for(int j=0;j<N;++j) sum += m.col(j).sumAbs2Elements();
            } else {
                const int M = m.colsize();
                for(int i=0;i<M;++i) sum += m.row(i).sumAbs2Elements();
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
                const int M = colsize();
                for(int i=0;i<M;++i) sum += row(i).normSq(scale);
            } else {
                const int N = rowsize();
                for(int j=0;j<N;++j) sum += col(j).normSq(scale);
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
                const int N = m.rowsize();
                for(int j=0;j<N;++j) {
                    RT temp = m.col(j).normInf();
                    if (temp > max) max = temp;
                }
            } else {
                const int M = m.colsize();
                for(int i=0;i<M;++i) {
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
                const int N = m.rowsize();
                for(int j=0;j<N;++j) {
                    RT temp = m.col(j).maxAbs2Element();
                    if (temp > max) max = temp;
                }
            } else {
                const int M = m.colsize();
                for(int i=0;i<M;++i) {
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
        const int N = m.rowsize();
        for(int j=0;j<N;++j) {
            RT temp = m.col(j).norm1();
            if (temp > max) max = temp;
        }
        return max;
    }

    template <class T>
    static RT NonLapNormF(const GenMatrix<T>& m)
    {
        //std::cout<<"NonLapNorm "<<std::endl;
        //std::cout<<"m = "<<m<<std::endl;
        const RT eps = TMV_Epsilon<T>();

        RT mmax = m.maxAbs2Element();
        RT norm;
        if (mmax == RT(0)) {
            // Then norm is also 0
            norm = RT(0);
            //std::cout<<"case 0: norm = "<<norm<<std::endl;
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
            //std::cout<<"case 1: norm = "<<norm<<std::endl;
        } else if (RT(1) / mmax == RT(0)) {
            // Then mmax is already inf, so no hope of making it more accurate.
            norm = mmax;
            //std::cout<<"case 2: norm = "<<norm<<std::endl;
        } else if (RT(1) / (mmax*mmax) == RT(0)) {
            // Then we have overflow, so we need to rescale:
            const RT inveps = RT(1)/eps;
            RT scale = eps;
            mmax *= scale;
            while (mmax > inveps) { scale *= eps; mmax *= eps; }
            norm = TMV_SQRT(m.normSq(scale))/scale;
            //std::cout<<"case 3: norm = "<<norm<<std::endl;
        }  else {
            norm = TMV_SQRT(m.normSq());
            //std::cout<<"case 4: norm = "<<norm<<std::endl;
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

    template <class T>
    auto_ptr<BaseMatrix<T> > GenMatrix<T>::newCopy() const
    {
        auto_ptr<BaseMatrix<T> > a;
        if (isrm()) a.reset(new Matrix<T,RowMajor>(*this));
        else a.reset(new Matrix<T,ColMajor>(*this));
        return a;
    }

    template <class T>
    auto_ptr<BaseMatrix<T> > GenMatrix<T>::newView() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstMatrixView<T>(view()));
        return a;
    }

    template <class T>
    auto_ptr<BaseMatrix<T> > GenMatrix<T>::newTranspose() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstMatrixView<T>(transpose()));
        return a;
    }

    template <class T>
    auto_ptr<BaseMatrix<T> > GenMatrix<T>::newConjugate() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstMatrixView<T>(conjugate()));
        return a;
    }

    template <class T>
    auto_ptr<BaseMatrix<T> > GenMatrix<T>::newAdjoint() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstMatrixView<T>(adjoint()));
        return a;
    }

    template <class T>
    auto_ptr<BaseMatrix<T> > GenMatrix<T>::newInverse() const
    {
        auto_ptr<Matrix<T,ColMajor> > minv(
            new Matrix<T,ColMajor>(rowsize(),colsize()));
        makeInverse(minv->view());
        BaseMatrix<T>* ret1 = minv.release();
        auto_ptr<BaseMatrix<T> > ret(ret1);
        return ret;
    }
#ifdef INST_INT
    template <>
    auto_ptr<BaseMatrix<int> > GenMatrix<int>::newInverse() const
    { TMVAssert(TMV_FALSE); return auto_ptr<BaseMatrix<int> >(); }
    template <>
    auto_ptr<BaseMatrix<std::complex<int> > > 
    GenMatrix<std::complex<int> >::newInverse() const
    { 
        TMVAssert(TMV_FALSE); 
        return auto_ptr<BaseMatrix<std::complex<int> > >(); 
    }
#endif

    //
    // Modifying Functions
    //

    template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::clip(RT thresh) const
    {
        TMVAssert(I==CStyle);
        if (this->canLinearize()) linearView().clip(thresh);
        else {
            if (this->isrm()) {
                const int M = colsize();
                for(int i=0;i<M;++i) row(i).clip(thresh);
            }
            else {
                const int N = rowsize();
                for(int j=0;j<N;++j) col(j).clip(thresh);
            }
        }
        return *this; 
    }

    template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::setZero() const
    {
        TMVAssert(I==CStyle);
        if (this->canLinearize()) linearView().setZero();
        else {
            if (this->isrm()) {
                const int M = colsize();
                for(int i=0;i<M;++i) row(i).setZero();
            }
            else {
                const int N = rowsize();
                for(int j=0;j<N;++j) col(j).setZero();
            }
        }
        return *this; 
    }

    template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::setAllTo(const T& x) const
    {
        TMVAssert(I==CStyle);
        if (this->canLinearize()) linearView().setAllTo(x);
        else {
            if (this->isrm()) {
                const int M = colsize();
                for(int i=0;i<M;++i) row(i).setAllTo(x); 
            } else  {
                const int N = rowsize();
                for(int j=0;j<N;++j) col(j).setAllTo(x); 
            }
        }
        return *this; 
    }

    template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::addToAll(const T& x) const
    {
        TMVAssert(I==CStyle);
        if (this->canLinearize()) linearView().addToAll(x);
        else {
            if (this->isrm()) {
                const int M = colsize();
                for(int i=0;i<M;++i) row(i).addToAll(x); 
            } else  {
                const int N = rowsize();
                for(int j=0;j<N;++j) col(j).addToAll(x); 
            }
        }
        return *this; 
    }

    template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::transposeSelf() const
    {
        TMVAssert(I==CStyle);
        TMVAssert(colsize() == rowsize());
        const int M = colsize();
        for(int i=1;i<M;++i) Swap(row(i,0,i),col(i,0,i));
        return *this; 
    }

    template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::conjugateSelf() const
    {
        TMVAssert(I==CStyle);
        if (isComplex(T())) {
            if (this->canLinearize()) linearView().conjugateSelf();
            else {
                if (this->isrm()) {
                    const int M = colsize();
                    for(int i=0;i<M;++i) row(i).conjugateSelf();
                } else  {
                    const int N = rowsize();
                    for(int j=0;j<N;++j) col(j).conjugateSelf();
                }
            }
        }
        return *this; 
    }

    template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::setToIdentity(const T& x) const 
    {
        TMVAssert(I==CStyle);
        TMVAssert(colsize() == rowsize());
        setZero();
        diag().setAllTo(x);
        return *this;
    }

    template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::permuteRows(
        const int* p, int i1, int i2) const
    {
        TMVAssert(I==CStyle);
        TMVAssert(i2<=int(colsize()));
        TMVAssert(i1<=i2);
        // This idea of doing the permutations a block at a time
        // is cribbed from the LAPack code.  It does speed things up 
        // quite a bit for large matrices.  On my machine where BLOCKSIZE=64
        // is optimal for most routines, blocks of 32 were optimal here,
        // so I use BLOCKSIZE/2 in general.
        const int N = rowsize();
        const int Nx = N/PERM_BLOCKSIZE*PERM_BLOCKSIZE;
        if (Nx != 0) {
            for(int j=0;j<Nx;) {
                int j2 = j+PERM_BLOCKSIZE;
                const int* pi = p+i1;
                for(int i=i1;i<i2;++i,++pi) {
                    TMVAssert(*pi < int(colsize()));
                    colRange(j,j2).swapRows(i,*pi);
                }
                j = j2;
            }
        }
        if (Nx != N) {
            const int* pi = p+i1;
            for(int i=i1;i<i2;++i,++pi) {
                TMVAssert(*pi < int(colsize()));
                colRange(Nx,N).swapRows(i,*pi);
            }
        }
        return *this;
    }

    template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::reversePermuteRows(
        const int* p, int i1, int i2) const
    {
        TMVAssert(I==CStyle);
        TMVAssert(i2<=int(colsize()));
        TMVAssert(i1<=i2);
        const int N = rowsize();
        const int Nx = N/PERM_BLOCKSIZE*PERM_BLOCKSIZE;
        if (Nx != 0) {
            for(int j=0;j<Nx;) {
                int j2 = j+PERM_BLOCKSIZE;
                const int* pi = p+i2;
                for(int i=i2;i>i1;) {
                    --i; --pi;
                    TMVAssert(*pi < int(colsize()));
                    colRange(j,j2).swapRows(i,*pi);
                }
                j = j2;
            }
        }
        if (Nx != N) {
            const int* pi = p+i2;
            for(int i=i2;i>i1;) {
                --i; --pi;
                TMVAssert(*pi < int(colsize()));
                colRange(Nx,N).swapRows(i,*pi);
            }
        }
        return *this;
    }

    //
    // Copy Matrices
    //

    template <class T>
    static void NonLapCopy(const GenMatrix<T>& m1, const MatrixView<T>& m2)
    {
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m1.ct()==NonConj);
        TMVAssert(m2.ct()==NonConj);
        TMVAssert(!(m1.isrm() && m2.isrm()));
        const int M = m2.colsize();
        const int N = m2.rowsize();

        if (m1.iscm() && m2.iscm()) {
            const T* p1 = m1.cptr();
            T* p2 = m2.ptr();
            const int s1 = m1.stepj();
            const int s2 = m2.stepj();
            for(int j=0;j<N;++j,p1+=s1,p2+=s2) {
                std::copy(p1,p1+M,p2);
            }
        } else if (M > N) {
            if (shouldReverse(m1.stepi(),m2.stepi()))
                for(int j=0;j<N;++j) 
                    DoCopySameType(m1.col(j).reverse(),m2.col(j).reverse());
            else
                for(int j=0;j<N;++j) 
                    DoCopySameType(m1.col(j),m2.col(j));
        } else {
            if (shouldReverse(m1.stepj(),m2.stepj()))
                for(int i=0;i<M;++i) 
                    DoCopySameType(m1.row(i).reverse(),m2.row(i).reverse());
            else
                for(int i=0;i<M;++i) 
                    DoCopySameType(m1.row(i),m2.row(i));
        }
    }
#ifdef ELAP
    template <class T>
    static inline void LapCopy(const GenMatrix<T>& m1, const MatrixView<T>& m2)
    { NonLapCopy(m1,m2); }
#ifdef INST_DOUBLE
    template <>
    void LapCopy(const GenMatrix<double>& m1, const MatrixView<double>& m2)
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
        const MatrixView<std::complex<double> >& m2)
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
    void LapCopy(const GenMatrix<float>& m1, const MatrixView<float>& m2)
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
        const MatrixView<std::complex<float> >& m2)
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
    void DoCopySameType(const GenMatrix<T>& m1, const MatrixView<T>& m2)
    {
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m2.rowsize() > 0);
        TMVAssert(m2.colsize() > 0);
        TMVAssert(m1.ct() == NonConj);
        TMVAssert(m2.ct() == NonConj);
        TMVAssert(!m2.isrm());
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
    void Swap(const MatrixView<T>& m1, const MatrixView<T>& m2)
    {
        TMVAssert(m1.colsize() == m2.colsize());
        TMVAssert(m1.rowsize() == m2.rowsize());
        if (m1.stor() == m2.stor() && m1.canLinearize() && m2.canLinearize()) {
            TMVAssert(m1.linearView().size() == m2.linearView().size());
            TMVAssert(m1.stepi()==m2.stepi() && m1.stepj()==m2.stepj());
            Swap(m1.linearView(),m2.linearView());
        }
        else {
            if (m1.isrm() && m2.isrm()) {
                const int M = m1.colsize();
                for(int i=0;i<M;++i) Swap(m1.row(i),m2.row(i)); 
            } else {
                const int N = m1.rowsize();
                for(int j=0;j<N;++j) Swap(m1.col(j),m2.col(j)); 
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
            const int M = m1.colsize();
            for(int i=0;i<M;++i) 
                if (m1.row(i) != m2.row(i)) return false;
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

    template <bool conj, bool rm, bool th, class T> 
    static void DoWrite(std::ostream& os, const GenMatrix<T>& m, RT thresh)
    {
        const T* mrowi = m.cptr();
        const int sj = rm ? 1 : m.stepj();
        os << m.colsize() <<"  "<<m.rowsize()<<std::endl;
        for(int i=m.colsize();i>0;--i,mrowi+=m.stepi()) {
            os << "( ";
            const T* mij = mrowi;
            for(int k=m.rowsize();k>0;--k,rm?++mij:mij+=sj) {
                if (conj) {
                    if (th) {
                        os << ' '<<Value(TMV_ABS(*mij) < thresh ? T(0) :
                                         TMV_CONJ(*mij))<<' ';
                    } else {
                        os << ' '<<Value(TMV_CONJ(*mij))<<' ';
                    }
                } else {
                    if (th) {
                        os << ' '<<Value(TMV_ABS(*mij) < thresh ? T(0) :
                                         *mij)<<' ';
                    } else {
                        os << ' '<<Value(*mij)<<' ';
                    }
                }
            }
            os << " )\n";
        }
    }

    template <bool rm, bool th, class T>
    static inline void DoWrite1(
        std::ostream& os, const GenMatrix<T>& m, T thresh)
    { DoWrite<false,rm,th>(os,m,thresh); }

    template <bool rm, bool th, class T>
    static inline void DoWrite1(
        std::ostream& os, const GenMatrix<std::complex<T> >& m, T thresh)
    {
        if (m.isconj())
            DoWrite<true,rm,th>(os,m,thresh);
        else
            DoWrite<false,rm,th>(os,m,thresh);
    }

    template <class T>
    void GenMatrix<T>::write(std::ostream& os) const
    {
        if (isrm())
            DoWrite1<true,false>(os,*this,RT(0));
        else
            DoWrite1<false,false>(os,*this,RT(0));
    }

    template <class T>
    void GenMatrix<T>::write(std::ostream& os, RT thresh) const
    {
        if (isrm())
            DoWrite1<true,true>(os,*this,thresh);
        else
            DoWrite1<false,true>(os,*this,thresh);
    }

#ifndef NOTHROW
    template <class T>
    class MatrixReadError : public ReadError
    {
    public :
        int i,j;
        mutable auto_ptr<Matrix<T> > m;
        char exp,got;
        size_t cs,rs;
        bool is,iseof,isbad;

        MatrixReadError(
            int _i, int _j, const GenMatrix<T>& _m, std::istream& _is) throw() :
            ReadError("Matrix."),
            i(_i), j(_j), m(new Matrix<T>(_m)), exp(0), got(0),
            cs(_m.colsize()), rs(_m.rowsize()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        MatrixReadError(std::istream& _is) throw() :
            ReadError("Matrix."),
            i(0), j(0), m(0), exp(0), got(0), cs(0), rs(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        MatrixReadError(
            int _i, int _j, const GenMatrix<T>& _m,
            std::istream& _is, char _e, char _g) throw() :
            ReadError("Matrix."),
            i(_i), j(_j), m(new Matrix<T>(_m)), exp(_e), got(_g),
            cs(_m.colsize()), rs(_m.rowsize()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        MatrixReadError(
            const GenMatrix<T>& _m,
            std::istream& _is, size_t _cs, size_t _rs) throw() :
            ReadError("Matrix."),
            i(0), j(0), m(new Matrix<T>(_m)), exp(0), got(0), cs(_cs), rs(_rs),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        MatrixReadError(const MatrixReadError<T>& rhs) :
            i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
            cs(rhs.cs), rs(rhs.rs),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        virtual ~MatrixReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for Matrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (m.get() && cs != m->colsize()) {
                os<<"Wrong column size: expected "<<m->colsize()<<
                    ", got "<<cs<<".\n";
            }
            if (m.get() && rs != m->rowsize()) {
                os<<"Wrong row size: expected "<<m->rowsize()<<
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
            if (m.get()) {
                const int N = m->rowsize();
                os<<"The portion of the Matrix which was successfully "
                    "read is: \n";
                for(int ii=0;ii<i;++ii) {
                    os<<"( ";
                    for(int jj=0;jj<N;++jj)
                        os<<' '<<(*m)(ii,jj)<<' ';
                    os<<" )\n";
                }
                os<<"( ";
                for(int jj=0;jj<j;++jj) os<<' '<<(*m)(i,jj)<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    template <class T, IndexStyle I>
    void MatrixView<T,I>::read(std::istream& is) const
    {
        TMVAssert(I==CStyle);
        T* mrowi = ptr();
        const int sj = stepj();
        char paren;
        const int M = colsize();
        for(int i=0;i<M;++i,mrowi+=stepi()) {
            is >> paren;
            if (!is || paren != '(') {
#ifdef NOTHROW
                std::cerr<<"Matrix ReadError: "<<paren<<" != (\n";
                exit(1); 
#else
                throw MatrixReadError<T>(i,0,*this,is,'(',is?paren:'(');
#endif
            }
            T* mij = mrowi;
            if (this->isrm()) {
                for(int k=rowsize();k>0;--k,++mij) {
                    is >> *mij;
                    if (!is) {
#ifdef NOTHROW
                        std::cerr<<"Matrix ReadError: !is\n";
                        exit(1); 
#else
                        throw MatrixReadError<T>(i,rowsize()-k,*this,is);
#endif
                    }
                }
            } else {
                for(int k=rowsize();k>0;--k,mij+=sj) {
                    is >> *mij;
                    if (!is) {
#ifdef NOTHROW
                        std::cerr<<"Matrix ReadError: !is\n";
                        exit(1); 
#else
                        throw MatrixReadError<T>(i,rowsize()-k,*this,is);
#endif
                    }
                }
            }
            is >> paren;
            if ((!is && i+1<M)  || paren != ')') {
#ifdef NOTHROW
                std::cerr<<"Matrix ReadError: "<<paren<<" != )\n";
                exit(1); 
#else
                throw MatrixReadError<T>(i,rowsize(),*this,is,')',is?paren:')');
#endif
            }
        }
        if (this->isconj()) conjugateSelf();
    }

    template <class T, StorageType S, IndexStyle I>
    std::istream& operator>>(std::istream& is, auto_ptr<Matrix<T,S,I> >& m)
    {
        size_t cs,rs;
        is >> cs >> rs;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"Matrix ReadError: !is\n";
            exit(1); 
#else
            throw MatrixReadError<T>(is);
#endif
        }
        m.reset(new Matrix<T,S,I>(cs,rs));
        m->view().read(is); 
        return is;
    }

    template <class T, StorageType S, IndexStyle I>
    std::istream& operator>>(std::istream& is, Matrix<T,S,I>& m)
    {
        size_t cs,rs;
        is >> cs >> rs;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"Matrix ReadError: !is\n";
            exit(1); 
#else
            throw MatrixReadError<T>(is);
#endif
        }
        m.resize(cs,rs);
        m.view().read(is); 
        return is;
    }

    template <class T>
    std::istream& operator>>(std::istream& is, const MatrixView<T>& m)
    {
        size_t cs,rs;
        is >> cs >> rs;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"Matrix ReadError: !is\n";
            exit(1); 
#else
            throw MatrixReadError<T>(is);
#endif
        }
        if (m.colsize() != cs || m.rowsize() != rs) {
#ifdef NOTHROW
            std::cerr<<"Matrix ReadError: Wrong size\n";
            exit(1); 
#else
            throw MatrixReadError<T>(m,is,cs,rs);
#endif
        }
        TMVAssert(m.colsize() == cs);
        TMVAssert(m.rowsize() == rs);
        m.read(is);
        return is;
    }

#undef RT

#define InstFile "TMV_Matrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


