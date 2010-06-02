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
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_VIt.h"
#include <iostream>
#include "portable_platform.h"

namespace tmv {

#define RT TMV_RealType(T)

    //
    // Access
    //

    template <class T>
    T GenUpperTriMatrix<T>::cref(int i, int j) const
    {
        const T* mi = cptr() + int(i)*stepi() + int(j)*stepj();
        return (isconj() ? TMV_CONJ(*mi) : *mi);
    }

    template <class T, IndexStyle I> 
    typename UpperTriMatrixView<T,I>::reference
    UpperTriMatrixView<T,I>::ref(int i, int j) const
    {
        T* mi = ptr() + int(i)*stepi() + int(j)*stepj();
        return reference(isunit() && i==j,*mi,ct());
    }

    template <class T>
    T GenLowerTriMatrix<T>::cref(int i, int j) const
    {
        const T* mi = cptr() + int(i)*stepi() + int(j)*stepj();
        return (isconj() ? TMV_CONJ(*mi) : *mi);
    }

    template <class T, IndexStyle I> 
    typename LowerTriMatrixView<T,I>::reference
    LowerTriMatrixView<T,I>::ref(int i, int j) const
    {
        T* mi = ptr() + int(i)*stepi() + int(j)*stepj();
        return reference(isunit() && i==j,*mi,ct());
    }

    //
    // OK? (SubMatrix, etc.)
    //

    template <class T>
    bool GenUpperTriMatrix<T>::hasSubMatrix(
        int i1, int i2, int j1, int j2, int istep, int jstep) const
    {
        if (i1==i2 || j1==j2) return true; // no elements, so whatever...
        bool ok = true;
        if (istep == 0) {
            ok = false;
            std::cout<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1 < 0 || i1 >= int(size())) {
            ok = false;
            std::cout<<"first col element ("<<i1<<") must be in 0 -- ";
            std::cout<<size()-1<<std::endl;
        }
        if (i2-istep < 0 || i2-istep >= int(size())) {
            ok = false;
            std::cout<<"last col element ("<<i2-istep<<") must be in 0 -- ";
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
            std::cout<<"row range ("<<j2-j1<<") must be multiple of istep (";
            std::cout<<jstep<<")\n";
        }
        if ((j2-j1)/jstep < 0) {
            ok = false;
            std::cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
        }
        if (!this->okij(i1,j1)) {
            ok = false;
            std::cout<<"Upper left corner ("<<i1<<','<<j1;
            std::cout<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i1,j2-jstep)) {
            ok = false;
            std::cout<<"Upper right corner ("<<i1<<','<<j2-jstep;
            std::cout<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i2-istep,j1)) {
            ok = false;
            std::cout<<"Lower left corner ("<<i2-istep<<','<<j1;
            std::cout<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i2-istep,j2-jstep)) {
            ok = false;
            std::cout<<"Lower right corner ("<<i2-istep<<','<<j2-jstep;
            std::cout<<") must be in Upper Triangle\n";
        }
        return ok;
    }

    template <class T>
    bool GenUpperTriMatrix<T>::hasSubVector(
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
        if (!this->okij(i,j)) {
            ok = false;
            std::cout<<"First element ("<<i<<','<<j<<") must be in Triangle\n";
        }
        if (!this->okij(i2,j2)) {
            ok = false;
            std::cout<<"Last element ("<<i2<<','<<j2<<") must be in Triangle\n";
        }
        return ok;
    }

    template <class T>
    bool GenUpperTriMatrix<T>::hasSubTriMatrix(
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
        return ok;
    }

    template <class T>
    bool ConstUpperTriMatrixView<T,FortranStyle>::hasSubMatrix(
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
            std::cout<<"last col element ("<<i2<<") must be in 1 -- ";
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
            std::cout<<"last row element ("<<j2<<") must be in 1 -- ";
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
        if (!this->okij(i1-1,j1-1)) {
            ok = false;
            std::cout<<"Upper left corner ("<<i1<<','<<j1;
            std::cout<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i1-1,j2-1)) {
            ok = false;
            std::cout<<"Upper right corner ("<<i1<<','<<j2;
            std::cout<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i2-1,j1-1)) {
            ok = false;
            std::cout<<"Lower left corner ("<<i2<<','<<j1;
            std::cout<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i2-1,j2-1)) {
            ok = false;
            std::cout<<"Lower right corner ("<<i2<<','<<j2;
            std::cout<<") must be in Upper Triangle\n";
        }
        return ok;
    }

    template <class T>
    bool ConstUpperTriMatrixView<T,FortranStyle>::hasSubVector(
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
        if (!this->okij(i-1,j-1)) {
            ok = false;
            std::cout<<"First element ("<<i<<','<<j<<") must be in Triangle\n";
        }
        if (!this->okij(i2-1,j2-1)) {
            ok = false;
            std::cout<<"Last element ("<<i2<<','<<j2<<") must be in Triangle\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstUpperTriMatrixView<T,FortranStyle>::hasSubTriMatrix(
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
        if (i2<1 || i2 > int(this->size())) {
            ok = false;
            std::cout<<"last diag element ("<<i2<<") must be in 1 -- ";
            std::cout<<this->size()<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cout<<"range ("<<i2-i1<<") must be multiple of istep (";
            std::cout<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cout<<"n diag elements ("<<(i2-i1)/istep+1;
            std::cout<<") must be positive\n";
        }
        return ok;
    }

    //
    // norms
    //

    template <class T>
    T GenUpperTriMatrix<T>::sumElements() const
    {
        const int N = size();
        T sum(0);
        if (isrm()) 
            if (isunit())
                for(int i=0;i<N;++i) 
                    sum += row(i,i+1,N).sumElements();
            else
                for(int i=0;i<N;++i) 
                    sum += row(i,i,N).sumElements();
        else
            if (isunit())
                for(int j=0;j<N;++j) 
                    sum += col(j,0,j).sumElements();
            else
                for(int j=0;j<N;++j) 
                    sum += col(j,0,j+1).sumElements();
        if (isunit()) sum += N;
        return sum;
    }

    template <class T>
    RT GenUpperTriMatrix<T>::sumAbsElements() const
    {
        const int N = size();
        RT sum(0);
        if (isrm()) 
            if (isunit())
                for(int i=0;i<N;++i) 
                    sum += row(i,i+1,N).sumAbsElements();
            else
                for(int i=0;i<N;++i) 
                    sum += row(i,i,N).sumAbsElements();
        else
            if (isunit())
                for(int j=0;j<N;++j) 
                    sum += col(j,0,j).sumAbsElements();
            else
                for(int j=0;j<N;++j) 
                    sum += col(j,0,j+1).sumAbsElements();
        if (isunit()) sum += N;
        return sum;
    }

    template <class T>
    RT GenUpperTriMatrix<T>::normSq(const RT scale) const
    {
        const int N = size();
        RT sum(0);
        if (isrm()) 
            if (isunit())
                for(int i=0;i<N;++i) 
                    sum += row(i,i+1,N).normSq(scale);
            else
                for(int i=0;i<N;++i) 
                    sum += row(i,i,N).normSq(scale);
        else
            if (isunit())
                for(int j=0;j<N;++j) 
                    sum += col(j,0,j).normSq(scale);
            else
                for(int j=0;j<N;++j) 
                    sum += col(j,0,j+1).normSq(scale);
        if (isunit()) {
            if (scale == RT(1)) sum += N;
            else sum += N * scale * scale;
        }
        return sum;
    }

    template <class T> 
    static RT NonLapNormF(const GenUpperTriMatrix<T>& m)
    {
        const RT eps = TMV_Epsilon<T>();

        RT mmax = m.maxAbsElement();
        if (mmax == RT(0)) return RT(0);
        else if (mmax * mmax * eps == RT(0)) {
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
    static RT NonLapMaxAbsElement(const GenUpperTriMatrix<T>& m)
    {
        const int N = m.size();
        RT max(0);
        if (m.isrm()) {
            for(int i=0;i<N;++i) {
                RT temp;
                if (m.isunit())
                    temp = m.row(i,i+1,N).normInf();
                else 
                    temp = m.row(i,i,N).normInf();
                if (temp > max) max = temp;
            }
        } else {
            for(int j=0;j<N;++j) {
                RT temp;
                if (m.isunit())
                    temp = m.col(j,0,j).normInf();
                else 
                    temp = m.col(j,0,j+1).normInf();
                if (temp > max) max = temp;
            }
        }
        if (m.isunit() && max < RT(1)) max = RT(1);
        return max;
    }

    template <class T> 
    static RT NonLapNorm1(const GenUpperTriMatrix<T>& m)
    {
        RT max(0);
        const int N = m.size();
        for(int j=0;j<N;++j) {
            RT temp;
            if (m.isunit()) {
                temp = m.col(j,0,j).norm1();
                temp += RT(1);
            } else temp = m.col(j,0,j+1).norm1();
            if (temp > max) max = temp;
        }
        return max;
    } 

    template <class T> 
    static RT NonLapNormInf(const GenUpperTriMatrix<T>& m)
    {
        RT max(0);
        const int N = m.size();
        for(int j=0;j<N;++j) {
            RT temp;
            if (m.isunit()) {
                temp = m.row(j,j+1,N).norm1();
                temp += RT(1);
            } else temp = m.row(j,j,N).norm1();
            if (temp > max) max = temp;
        }
        return max;
    }
#ifdef XLAP
    template <class T> 
    static RT LapNorm(const char c, const GenUpperTriMatrix<T>& m)
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
    double LapNorm(
        const char c, const GenUpperTriMatrix<double>& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        char cc = c;
        if (m.isrm()) {
            if (c == '1') cc = 'I';
            else if (c == 'I') cc = '1';
        }
        int M = m.size();
        int N = M;
        int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
        auto_array<double> work(cc == 'I' ? new double[M] : 0);
#endif
        double norm = LAPNAME(dlantr) (LAPCM LAPV(cc),
                                       m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
                                       LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
                                       LAP1 LAP1 LAP1);
        return norm;
    }
    template <>
    double LapNorm(
        const char c, const GenUpperTriMatrix<std::complex<double> >& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        char cc = c;
        if (m.isrm()) {
            if (c == '1') cc = 'I';
            else if (c == 'I') cc = '1';
        }
        int M = m.size();
        int N = M;
        int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
        auto_array<double> work(cc == 'I' ? new double[M] : 0);
#endif
        double norm = LAPNAME(zlantr) (LAPCM LAPV(cc),
                                       m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
                                       LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
                                       LAP1 LAP1 LAP1);
        return norm;
    }
#endif
#ifdef INST_FLOAT
    template <>
    float LapNorm(
        const char c, const GenUpperTriMatrix<float>& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        char cc = c;
        if (m.isrm()) {
            if (c == '1') cc = 'I';
            else if (c == 'I') cc = '1';
        }
        int M = m.size();
        int N = M;
        int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
        auto_array<float> work(cc == 'I' ? new float[M] : 0);
#endif
        float norm = LAPNAME(slantr) (LAPCM LAPV(cc),
                                      m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
                                      LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
                                      LAP1 LAP1 LAP1);
        return norm;
    }
    template <>
    float LapNorm(
        const char c, const GenUpperTriMatrix<std::complex<float> >& m)
    {
        TMVAssert(m.iscm() || m.isrm());
        char cc = c;
        if (m.isrm()) {
            if (c == '1') cc = 'I';
            else if (c == 'I') cc = '1';
        }
        int M = m.size();
        int N = M;
        int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
        auto_array<float> work(cc == 'I' ? new float[M] : 0);
#endif
        float norm = LAPNAME(clantr) (LAPCM LAPV(cc),
                                      m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
                                      LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
                                      LAP1 LAP1 LAP1);
        return norm;
    }
#endif
#endif // XLAP

    template <class T>
    RT GenUpperTriMatrix<T>::maxAbsElement() const
    {
#ifdef XLAP
        return LapNorm('M',*this);
#else
        return NonLapMaxAbsElement(*this);
#endif
    }
    template <class T>
    RT GenUpperTriMatrix<T>::norm1() const
    {
#ifdef XLAP
        return LapNorm('1',*this);
#else
        return NonLapNorm1(*this);
#endif
    }
    template <class T>
    RT GenUpperTriMatrix<T>::normInf() const
    {
#ifdef XLAP
        return LapNorm('I',*this);
#else
        return NonLapNormInf(*this);
#endif
    }
    template <class T>
    RT GenUpperTriMatrix<T>::normF() const
    {
#ifdef XLAP
        return LapNorm('F',*this);
#else
        return NonLapNormF(*this);
#endif
    }

    template <class T>
    RT GenUpperTriMatrix<T>::doNorm2() const
    { return Matrix<T>(*this).doNorm2(); }

    template <class T>
    RT GenUpperTriMatrix<T>::doCondition() const
    { return Matrix<T>(*this).doCondition(); }

    template <class T>
    T GenUpperTriMatrix<T>::det() const
    {
        if (isunit()) return T(1);
        else return DiagMatrixViewOf(this->diag()).det(); 
    }

    template <class T>
    RT GenUpperTriMatrix<T>::logDet(T* sign) const
    {
        if (isunit()) {
            if (sign) *sign = T(1);
            return RT(0);
        } else {
            return DiagMatrixViewOf(this->diag()).logDet(sign); 
        }
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::newCopy() const
    {
        auto_ptr<BaseMatrix<T> > a;
        if (isunit()) {
            if (isrm()) a.reset(new UpperTriMatrix<T,UnitDiag,RowMajor>(*this));
            else a.reset(new UpperTriMatrix<T,UnitDiag,ColMajor>(*this));
        } else {
            if (isrm()) a.reset(
                new UpperTriMatrix<T,NonUnitDiag,RowMajor>(*this));
            else a.reset(new UpperTriMatrix<T,NonUnitDiag,ColMajor>(*this));
        }
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::newView() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstUpperTriMatrixView<T>(view()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::newTranspose() const
    {
        auto_ptr<BaseMatrix<T> > a(
            new ConstLowerTriMatrixView<T>(transpose()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::newConjugate() const
    {
        auto_ptr<BaseMatrix<T> > a(
            new ConstUpperTriMatrixView<T>(conjugate()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::newAdjoint() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstLowerTriMatrixView<T>(adjoint()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenUpperTriMatrix<T>::newInverse() const
    {
        if (isunit()) {
            if (isrm()) {
                auto_ptr<UpperTriMatrix<T,UnitDiag,ColMajor> > minv(
                    new UpperTriMatrix<T,UnitDiag,ColMajor>(*this));
                minv->invertSelf();
                BaseMatrix<T>* ret1 = minv.release();
                auto_ptr<BaseMatrix<T> > ret(ret1);
                return ret;
            } else {
                auto_ptr<UpperTriMatrix<T,UnitDiag,ColMajor> > minv(
                    new UpperTriMatrix<T,UnitDiag,ColMajor>(*this));
                minv->invertSelf();
                BaseMatrix<T>* ret1 = minv.release();
                auto_ptr<BaseMatrix<T> > ret(ret1);
                return ret;
            }
        } else {
            if (isrm()) {
                auto_ptr<UpperTriMatrix<T,NonUnitDiag,ColMajor> > minv(
                    new UpperTriMatrix<T,NonUnitDiag,ColMajor>(*this));
                minv->invertSelf();
                BaseMatrix<T>* ret1 = minv.release();
                auto_ptr<BaseMatrix<T> > ret(ret1);
                return ret;
            } else {
                auto_ptr<UpperTriMatrix<T,NonUnitDiag,ColMajor> > minv(
                    new UpperTriMatrix<T,NonUnitDiag,ColMajor>(*this));
                minv->invertSelf();
                BaseMatrix<T>* ret1 = minv.release();
                auto_ptr<BaseMatrix<T> > ret(ret1);
                return ret;
            }
        }
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::newCopy() const
    {
        auto_ptr<BaseMatrix<T> > a;
        if (isunit()) {
            if (isrm()) a.reset(new LowerTriMatrix<T,UnitDiag,RowMajor>(*this));
            else a.reset(new LowerTriMatrix<T,UnitDiag,ColMajor>(*this));
        } else {
            if (isrm()) a.reset(
                new LowerTriMatrix<T,NonUnitDiag,RowMajor>(*this));
            else a.reset(new LowerTriMatrix<T,NonUnitDiag,ColMajor>(*this));
        }
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::newView() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstLowerTriMatrixView<T>(view()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::newTranspose() const
    {
        auto_ptr<BaseMatrix<T> > a(
            new ConstUpperTriMatrixView<T>(transpose()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::newConjugate() const
    {
        auto_ptr<BaseMatrix<T> > a(
            new ConstLowerTriMatrixView<T>(conjugate()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::newAdjoint() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstUpperTriMatrixView<T>(adjoint()));
        return a;
    }

    template <class T> 
    auto_ptr<BaseMatrix<T> > GenLowerTriMatrix<T>::newInverse() const
    {
        if (isunit()) {
            if (isrm()) {
                auto_ptr<LowerTriMatrix<T,UnitDiag,ColMajor> > minv(
                    new LowerTriMatrix<T,UnitDiag,ColMajor>(*this));
                minv->invertSelf();
                BaseMatrix<T>* ret1 = minv.release();
                auto_ptr<BaseMatrix<T> > ret(ret1);
                return ret;
            } else {
                auto_ptr<LowerTriMatrix<T,UnitDiag,ColMajor> > minv(
                    new LowerTriMatrix<T,UnitDiag,ColMajor>(*this));
                minv->invertSelf();
                BaseMatrix<T>* ret1 = minv.release();
                auto_ptr<BaseMatrix<T> > ret(ret1);
                return ret;
            }
        } else {
            if (isrm()) {
                auto_ptr<LowerTriMatrix<T,NonUnitDiag,ColMajor> > minv(
                    new LowerTriMatrix<T,NonUnitDiag,ColMajor>(*this));
                minv->invertSelf();
                BaseMatrix<T>* ret1 = minv.release();
                auto_ptr<BaseMatrix<T> > ret(ret1);
                return ret;
            } else {
                auto_ptr<LowerTriMatrix<T,NonUnitDiag,ColMajor> > minv(
                    new LowerTriMatrix<T,NonUnitDiag,ColMajor>(*this));
                minv->invertSelf();
                BaseMatrix<T>* ret1 = minv.release();
                auto_ptr<BaseMatrix<T> > ret(ret1);
                return ret;
            }
        }
    }


    //
    // Modifying Functions
    //

    template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::setZero() const
    {
        const int N = size();

        if (isrm())
            if (isunit())
                for(int i=0;i<N;++i) row(i,i+1,N).setZero();
            else
                for(int i=0;i<N;++i) row(i,i,N).setZero();
        else 
            if (isunit())
                for(int j=0;j<N;++j) col(j,0,j).setZero();
            else
                for(int j=0;j<N;++j) col(j,0,j+1).setZero();
        return *this; 
    } 

    template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::setAllTo(const T& x) const
    {
        const int N = size();

        if (isrm())
            if (isunit())
                for(int i=0;i<N;++i) row(i,i+1,N).setAllTo(x); 
            else
                for(int i=0;i<N;++i) row(i,i,N).setAllTo(x); 
        else 
            if (isunit())
                for(int j=0;j<N;++j) col(j,0,j).setAllTo(x); 
            else
                for(int j=0;j<N;++j) col(j,0,j+1).setAllTo(x); 
        return *this; 
    }

    template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::clip(
        RT thresh) const
    {
        const int N = size();

        if (isrm())
            if (isunit())
                for(int i=0;i<N;++i) row(i,i+1,N).clip(thresh);
            else
                for(int i=0;i<N;++i) row(i,i,N).clip(thresh);
        else 
            if (isunit())
                for(int j=0;j<N;++j) col(j,0,j).clip(thresh);
            else
                for(int j=0;j<N;++j) col(j,0,j+1).clip(thresh);
        return *this; 
    }

    template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::conjugateSelf() const
    {
        const int N = size();

        if (isComplex(T())) {
            if (isrm())
                if (isunit())
                    for(int i=0;i<N;++i) row(i,i+1,N).conjugateSelf();
                else
                    for(int i=0;i<N;++i) row(i,i,N).conjugateSelf();
            else
                if (isunit())
                    for(int j=0;j<N;++j) col(j,0,j).conjugateSelf();
                else
                    for(int j=0;j<N;++j) col(j,0,j+1).conjugateSelf();
        }
        return *this; 
    }

    template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::setToIdentity(
        const T& x) const 
    {
        TMVAssert(!isunit() || x==T(1));
        setZero();
        if (!isunit()) diag().setAllTo(x);
        return *this;
    }

    //
    // Swap
    //

    template <class T>
    void Swap(const UpperTriMatrixView<T>& m1, const UpperTriMatrixView<T>& m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.dt() == m2.dt());
        const int N = m1.size();

        if (m1.isrm() && m2.isrm())
            if (m1.isunit()) 
                for(int i=0;i<N;++i) 
                    Swap(m1.row(i,i+1,N), m2.row(i,i+1,N));
            else
                for(int i=0;i<N;++i) 
                    Swap(m1.row(i,i,N), m2.row(i,i,N));
        else
            if (m1.isunit()) 
                for(int j=0;j<N;++j) 
                    Swap(m1.col(j,0,j), m2.col(j,0,j));
            else
                for(int j=0;j<N;++j) 
                    Swap(m1.col(j,0,j+1), m2.col(j,0,j+1));
    }

    //
    // m1 == m2
    //

    template <class T1, class T2> bool operator==(
        const GenUpperTriMatrix<T1>& m1, const GenUpperTriMatrix<T2>& m2)
    {
        if (m1.size() != m2.size()) return false;
        else if (m1.isSameAs(m2)) return true;
        else {
            const int N = m1.size();
            for(int j=0;j<N;++j) {
                if (m1.col(j,0,j) != m2.col(j,0,j)) return false;
            }

            if (m1.isunit() && !m2.isunit()) {
                for(int i=0;i<N;++i) if (m2(i,i) != T2(1)) return false;
            } else if (m2.isunit() && !m1.isunit()) {
                for(int i=0;i<N;++i) if (m1(i,i) != T1(1)) return false;
            } else if (!m1.isunit() && !m2.isunit()) {
                if (m1.diag() != m2.diag()) return false;
            }

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
        std::ostream& os, const GenUpperTriMatrix<T>& m, RT thresh)
    {
        const T* mrowi = m.cptr();
        const int sj = rm?1:m.stepj();
        const int ds = m.stepi()+sj;
        const int N = m.size();
        int len = m.size();

        if (m.isunit()) {
            mrowi += sj;
            --len;
        }

        if (compact)
            os << "U " << N << std::endl;
        else
            os << N <<' '<< N << std::endl;

        for(int i=0;i<N;++i,--len,mrowi+=ds) {
            os << "( ";
            if (!compact) {
                for(int j=0;j<i;j++) os << ' '<<Value(T(0))<<' ';
            }
            if (m.isunit()) os << ' '<<Value(T(1))<<' ';
            const T* mij = mrowi;
            for(int k=len;k>0;--k,rm?++mij:mij+=sj) 
                if (conj) 
                    if (th)
                        os << ' '<<Value(TMV_ABS(*mij) < thresh ?
                                         T(0) : TMV_CONJ(*mij))<<' ';
                    else
                        os << ' '<<Value(TMV_CONJ(*mij))<<' ';
                else 
                    if (th)
                        os << ' '<<Value(TMV_ABS(*mij) < thresh ?
                                         T(0) : *mij)<<' ';
                    else
                        os << ' '<<Value(*mij)<<' '; 
            os << " )\n";
        }
    }

    template <bool rm, bool compact, bool th, class T> 
    static inline void DoWrite1(
        std::ostream& os, const GenUpperTriMatrix<T>& m, T thresh)
    { DoWrite<false,rm,compact,th>(os,m,thresh); }

    template <bool rm, bool compact, bool th, class T> 
    static inline void DoWrite1(
        std::ostream& os, const GenUpperTriMatrix<std::complex<T> >& m,
        T thresh)
    {
        if (m.isconj())
            DoWrite<true,rm,compact,th>(os,m,thresh); 
        else
            DoWrite<false,rm,compact,th>(os,m,thresh); 
    }

    template <class T>
    void GenUpperTriMatrix<T>::write(std::ostream& os) const
    {
        if (isrm())
            DoWrite1<true,false,false>(os,*this,RT(0));
        else
            DoWrite1<false,false,false>(os,*this,RT(0));
    }

    template <class T>
    void GenUpperTriMatrix<T>::write(std::ostream& os, RT thresh) const
    {
        if (isrm())
            DoWrite1<true,false,true>(os,*this,thresh);
        else
            DoWrite1<false,false,true>(os,*this,thresh);
    }

    template <class T>
    void GenUpperTriMatrix<T>::writeCompact(std::ostream& os) const
    {
        if (isrm())
            DoWrite1<true,true,false>(os,*this,RT(0));
        else
            DoWrite1<false,true,false>(os,*this,RT(0));
    }

    template <class T>
    void GenUpperTriMatrix<T>::writeCompact(std::ostream& os, RT thresh) const
    {
        if (isrm())
            DoWrite1<true,true,true>(os,*this,thresh);
        else
            DoWrite1<false,true,true>(os,*this,thresh);
    }

#ifndef NOTHROW
    template <class T>
    class UpperTriMatrixReadError : public ReadError
    {
    public :
        int i,j;
        mutable auto_ptr<UpperTriMatrix<T> > m;
        char exp,got;
        T unitgot;
        size_t s;
        bool is, iseof, isbad;

        UpperTriMatrixReadError(
            int _i, int _j, const GenUpperTriMatrix<T>& _m, std::istream& _is
        ) throw() :
            ReadError("UpperTriMatrix."),
            i(_i), j(_j), m(new UpperTriMatrix<T>(_m)),
            exp(0), got(0), unitgot(T(1)), s(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        UpperTriMatrixReadError(std::istream& _is) throw() :
            ReadError("UpperTriMatrix."),
            i(0), j(0), m(0), exp(0), got(0), unitgot(T(1)), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        UpperTriMatrixReadError(
            int _i, int _j, const GenUpperTriMatrix<T>& _m,
            std::istream& _is, char _e, char _g
        ) throw() :
            ReadError("UpperTriMatrix."),
            i(_i), j(_j), m(new UpperTriMatrix<T>(_m)),
            exp(_e), got(_g), unitgot(T(1)), s(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        UpperTriMatrixReadError(std::istream& _is, char _e, char _g) throw() :
            ReadError("UpperTriMatrix."),
            i(0), j(0), m(0), exp(_e), got(_g), unitgot(T(1)), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        UpperTriMatrixReadError(
            int _i, int _j, const GenUpperTriMatrix<T>& _m,
            std::istream& _is, T _u
        ) throw() :
            ReadError("UpperTriMatrix."),
            i(_i), j(_j), m(new UpperTriMatrix<T>(_m)),
            exp(0), got(0), unitgot(_u), s(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        UpperTriMatrixReadError(
            const GenUpperTriMatrix<T>& _m,
            std::istream& _is, size_t _s
        ) throw() :
            ReadError("UpperTriMatrix."),
            i(0), j(0), m(new UpperTriMatrix<T>(_m)),
            exp(0), got(0), unitgot(T(1)), s(_s),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        UpperTriMatrixReadError(const UpperTriMatrixReadError<T>& rhs) :
            i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
            unitgot(rhs.unitgot), s(rhs.s),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        virtual ~UpperTriMatrixReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for UpperTriMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (unitgot != T(1)) {
                os<<"Wrong format: expected 1 on the diagonal, got "<<
                    unitgot<<".\n";
            }
            if (m.get() && s != m->size()) {
                os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
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
                os<<"The portion of the UpperTriMatrix which was successfully "
                    "read is:\n";
                ConstUpperTriMatrixView<T> mm = m->view();
                for(int ii=0;ii<i;++ii) {
                    const int N = mm.rowsize();
                    os<<"( ";
                    for(int jj=0;jj<N;++jj)
                        os<<' '<<mm(ii,jj)<<' ';
                    os<<" )\n";
                }
                os<<"( ";
                for(int jj=0;jj<j;++jj) os<<' '<<mm(i,jj)<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    template <class T, IndexStyle I> 
    void UpperTriMatrixView<T,I>::read(std::istream& is) const
    {
        T* mrowi = ptr();
        const int sj = stepj();
        const int ds = stepi()+sj;
        const int N = size();
        int len = N;
        if (isunit()) {
            mrowi += sj;
            --len;
        }
        char paren;
        for(int i=0;i<N;++i,--len,mrowi+=ds) {
            is >> paren;
            if (!is || paren != '(') {
#ifdef NOTHROW
                std::cerr<<"UpperTriMatrix ReadError: "<<paren<<" != (\n";
                exit(1); 
#else
                throw UpperTriMatrixReadError<T>(i,0,*this,is,'(',is?paren:'(');
#endif
            }
            if (isunit()) {
                T unit;
                is >> unit;
                if (!is || unit != T(1)) {
#ifdef NOTHROW
                    std::cerr<<"UpperTriMatrix ReadError: "<<unit<<" != 1\n";
                    exit(1);
#else
                    throw UpperTriMatrixReadError<T>(i,i,*this,is,is?unit:T(1));
#endif
                }
            }
            T* mij = mrowi;
            if (isrm()) {
                for(int k=len;k>0;--k,++mij) {
                    is >> *mij;
                    if (!is) {
#ifdef NOTHROW
                        std::cerr<<"UpperTriMatrix ReadError: !is\n"; 
                        exit(1); 

#else
                        throw UpperTriMatrixReadError<T>(i,N-k,*this,is);
#endif
                    }
                }
            } else  {
                for(int k=len;k>0;--k,mij+=sj) {
                    is >> *mij;
                    if (!is) {
#ifdef NOTHROW
                        std::cerr<<"UpperTriMatrix ReadError: !is\n"; 
                        exit(1); 
#else
                        throw UpperTriMatrixReadError<T>(i,N-k,*this,is);
#endif
                    }
                }
            }
            is >> paren;
            if ((!is && i+1<N)  || paren != ')') {
#ifdef NOTHROW
                std::cerr<<"UpperTriMatrix ReadError: "<<paren<<" != )\n"; 
                exit(1); 
#else
                throw UpperTriMatrixReadError<T>(i,N,*this,is,')',is?paren:')');
#endif
            }
        }
        if (isconj()) conjugateSelf();
    }

    template <class T, DiagType D, StorageType S, IndexStyle I> 
    std::istream& operator>>(
        std::istream& is, auto_ptr<UpperTriMatrix<T,D,S,I> >& m)
    {
        char ul;
        is >> ul;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw UpperTriMatrixReadError<T>(is);
#endif
        }
        if (ul != 'U') {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix ReadError: "<<ul<<" != U\n"; 
            exit(1); 
#else
            throw UpperTriMatrixReadError<T>(is,'U',ul);
#endif
        }
        size_t size;
        is >> size;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw UpperTriMatrixReadError<T>(is);
#endif
        }
        m.reset(new UpperTriMatrix<T,D,S,I>(size));
        m->view().read(is); 
        return is;
    }

    template <class T>
    std::istream& operator>>(std::istream& is, const UpperTriMatrixView<T>& m)
    {
        char ul;
        is >> ul;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw UpperTriMatrixReadError<T>(is);
#endif
        }
        if (ul != 'U') {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix ReadError: "<<ul<<" != U\n"; 
            exit(1); 
#else
            throw UpperTriMatrixReadError<T>(is,'U',ul);
#endif
        }
        size_t s;
        is >> s;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw UpperTriMatrixReadError<T>(is);
#endif
        }
        if (m.size() != s) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix ReadError: Wrong size\n"; 
            exit(1); 
#else
            throw UpperTriMatrixReadError<T>(m,is,s);
#endif
        }
        TMVAssert(m.size() == s);
        m.read(is);
        return is;
    }

    template <bool conj, bool rm, bool compact, bool th, class T> 
    static void DoWrite(
        std::ostream& os, const GenLowerTriMatrix<T>& m, RT thresh)
    {
        const T* mrowi = m.cptr();
        const int sj = rm?1:m.stepj();
        const int si = m.stepi();
        const int N = m.size();
        int len = m.isunit() ? 0 : 1;

        if (compact)
            os << "L " << N << std::endl;
        else
            os << N <<' '<< N << std::endl;

        for(int i=N;i>0;--i,++len,mrowi+=si) {
            os << "( ";
            const T* mij = mrowi;
            for(int k=len;k>0;--k,rm?++mij:mij+=sj) 
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

            if (m.isunit()) os << ' '<<Value(T(1))<<' ';
            if (!compact)
                for(int j=m.isunit()?len+1:len;j<N;j++) 
                    os << ' '<<Value(T(0))<<' ';
            os << " )\n";
        }
    }

    template <bool rm, bool compact, bool th, class T> 
    static inline void DoWrite1(
        std::ostream& os, const GenLowerTriMatrix<T>& m, T thresh)
    { DoWrite<false,rm,compact,th>(os,m,thresh); }

    template <bool rm, bool compact, bool th, class T> 
    static inline void DoWrite1(
        std::ostream& os, const GenLowerTriMatrix<std::complex<T> >& m,
        T thresh)
    {
        if (m.isconj())
            DoWrite<true,rm,compact,th>(os,m,thresh); 
        else
            DoWrite<false,rm,compact,th>(os,m,thresh); 
    }

    template <class T>
    void GenLowerTriMatrix<T>::write(std::ostream& os) const
    {
        if (isrm())
            DoWrite1<true,false,false>(os,*this,RT(0));
        else
            DoWrite1<false,false,false>(os,*this,RT(0));
    }

    template <class T>
    void GenLowerTriMatrix<T>::write(std::ostream& os, RT thresh) const
    {
        if (isrm())
            DoWrite1<true,false,true>(os,*this,thresh);
        else
            DoWrite1<false,false,true>(os,*this,thresh);
    }

    template <class T>
    void GenLowerTriMatrix<T>::writeCompact(std::ostream& os) const
    {
        if (isrm())
            DoWrite1<true,true,false>(os,*this,RT(0));
        else
            DoWrite1<false,true,false>(os,*this,RT(0));
    }

    template <class T>
    void GenLowerTriMatrix<T>::writeCompact(std::ostream& os, RT thresh) const
    {
        if (isrm())
            DoWrite1<true,true,true>(os,*this,thresh);
        else
            DoWrite1<false,true,true>(os,*this,thresh);
    }

#ifndef NOTHROW
    template <class T>
    class LowerTriMatrixReadError : public ReadError
    {
    public :
        int i,j;
        mutable auto_ptr<LowerTriMatrix<T> > m;
        char exp,got;
        T unitgot;
        size_t s;
        bool is, iseof, isbad;

        LowerTriMatrixReadError(int _i, int _j,
                                const GenLowerTriMatrix<T>& _m, std::istream& _is) throw() :
            ReadError("LowerTriMatrix."),
            i(_i), j(_j), m(new LowerTriMatrix<T>(_m)),
            exp(0), got(0), unitgot(T(1)), s(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        LowerTriMatrixReadError(std::istream& _is) throw() :
            ReadError("LowerTriMatrix."),
            i(0), j(0), m(0), exp(0), got(0), unitgot(T(1)), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        LowerTriMatrixReadError(int _i, int _j,
                                const GenLowerTriMatrix<T>& _m, std::istream& _is,
                                char _e, char _g) throw() :
            ReadError("LowerTriMatrix."),
            i(_i), j(_j), m(new LowerTriMatrix<T>(_m)),
            exp(_e), got(_g), unitgot(T(1)), s(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        LowerTriMatrixReadError(std::istream& _is, char _e, char _g) throw() :
            ReadError("LowerTriMatrix."),
            i(0), j(0), m(0), exp(_e), got(_g), unitgot(T(1)), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        LowerTriMatrixReadError(int _i, int _j,
                                const GenLowerTriMatrix<T>& _m, std::istream& _is, T _u) throw() :
            ReadError("LowerTriMatrix."),
            i(_i), j(_j), m(new LowerTriMatrix<T>(_m)),
            exp(0), got(0), unitgot(_u), s(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        LowerTriMatrixReadError(const GenLowerTriMatrix<T>& _m,
                                std::istream& _is, size_t _s) throw() :
            ReadError("LowerTriMatrix."),
            i(0), j(0), m(new LowerTriMatrix<T>(_m)),
            exp(0), got(0), unitgot(T(1)), s(_s),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        LowerTriMatrixReadError(const LowerTriMatrixReadError<T>& rhs) throw() :
            i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got),
            unitgot(rhs.unitgot), s(rhs.s),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        virtual ~LowerTriMatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for LowerTriMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (unitgot != T(1)) {
                os<<"Wrong format: expected 1 on the diagonal, got "<<unitgot<<".\n";
            }
            if (m.get() && s != m->size()) {
                os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
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
                os<<"The portion of the LowerTriMatrix which was successfully read is:\n";
                ConstLowerTriMatrixView<T> mm = m->view();
                for(int ii=0;ii<i;++ii) {
                    const int N = mm.size();
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
    void LowerTriMatrixView<T,I>::read(std::istream& is) const
    {
        T* mrowi = ptr();
        const int sj = stepj();
        int len = isunit() ? 0 : 1;
        char paren;
        const int N = size();
        for(int i=0;i<N;++i,++len,mrowi+=stepi()) {
            is >> paren;
            if (!is || paren != '(') {
#ifdef NOTHROW
                std::cerr<<"LowerTriMatrix ReadError: "<<paren<<" != (\n"; 
                exit(1); 
#else
                throw LowerTriMatrixReadError<T>(i,0,*this,is,'(',is?paren:'(');
#endif
            }
            T* mij = mrowi;
            if (isrm()) {
                for(int k=len;k>0;--k,++mij) {
                    is >> *mij;
                    if (!is) {
#ifdef NOTHROW
                        std::cerr<<"LowerTriMatrix ReadError: !is\n"; 
                        exit(1); 
#else
                        throw LowerTriMatrixReadError<T>(i,len-k,*this,is);
#endif
                    }
                }
            } else  {
                for(int k=len;k>0;--k,mij+=sj) {
                    is >> *mij;
                    if (!is) {
#ifdef NOTHROW
                        std::cerr<<"LowerTriMatrix ReadError: !is\n"; 
                        exit(1); 
#else
                        throw LowerTriMatrixReadError<T>(i,len-k,*this,is);
#endif
                    }
                }
            }
            if (isunit()) {
                T unit;
                is >> unit;
                if (!is || unit != T(1)) {
#ifdef NOTHROW
                    std::cerr<<"LowerTriMatrix ReadError: "<<unit<<" != 1\n"; 
                    exit(1); 
#else
                    throw LowerTriMatrixReadError<T>(i,i,*this,is,is?unit:T(1));
#endif
                }
            }
            is >> paren;
            if ((!is && i+1<N)  || paren != ')') {
#ifdef NOTHROW
                std::cerr<<"LowerTriMatrix ReadError: "<<paren<<" != )\n"; 
                exit(1); 
#else
                throw LowerTriMatrixReadError<T>(i,N,*this,is,')',is?paren:')');
#endif
            }
        }
        if (isconj()) conjugateSelf();
    }

    template <class T, DiagType D, StorageType S, IndexStyle I>
    std::istream& operator>>(std::istream& is, 
                             auto_ptr<LowerTriMatrix<T,D,S,I> >& m)
    {
        char ul;
        is >> ul;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw LowerTriMatrixReadError<T>(is);
#endif
        }
        if (ul != 'L') {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix ReadError: "<<ul<<" != L\n"; 
            exit(1); 
#else
            throw LowerTriMatrixReadError<T>(is,'L',ul);
#endif
        }
        size_t size;
        is >> size;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw LowerTriMatrixReadError<T>(is);
#endif
        }
        m.reset(new LowerTriMatrix<T,D,S,I>(size));
        m->view().read(is); 
        return is;
    }

    template <class T>
    std::istream& operator>>(
        std::istream& is, const LowerTriMatrixView<T>& m)
    {
        char ul;
        is >> ul;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw LowerTriMatrixReadError<T>(is);
#endif
        }
        if (ul != 'L') {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix ReadError: "<<ul<<" != L\n"; 
            exit(1); 
#else
            throw LowerTriMatrixReadError<T>(is,'L',ul);
#endif
        }
        size_t s;
        is >> s;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw LowerTriMatrixReadError<T>(is);
#endif
        }
        if (m.size() != s) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix ReadError: Wrong size\n"; 
            exit(1); 
#else
            throw LowerTriMatrixReadError<T>(m,is,s);
#endif
        }
        TMVAssert(m.size() == s);
        m.read(is);
        return is;
    }


#define InstFile "TMV_TriMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


