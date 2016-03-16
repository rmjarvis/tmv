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
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_VIt.h"
#include <iostream>

namespace tmv {

#define RT TMV_RealType(T)

    //
    // OK? (SubMatrix, etc.)
    //

    template <class T>
    bool GenUpperTriMatrix<T>::hasSubMatrix(
        ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
    {
        if (i1==i2 || j1==j2) return true; // no elements, so whatever...
        bool ok = true;
        if (istep == 0) {
            ok = false;
            std::cerr<<"istep ("<<istep<<") can not be 0\n";
        }
        if (i1 < 0 || i1 >= size()) {
            ok = false;
            std::cerr<<"first col element ("<<i1<<") must be in 0 -- ";
            std::cerr<<size()-1<<std::endl;
        }
        if (i2-istep < 0 || i2-istep >= size()) {
            ok = false;
            std::cerr<<"last col element ("<<i2-istep<<") must be in 0 -- ";
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
        if (!this->okij(i1,j1)) {
            ok = false;
            std::cerr<<"Upper left corner ("<<i1<<','<<j1;
            std::cerr<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i1,j2-jstep)) {
            ok = false;
            std::cerr<<"Upper right corner ("<<i1<<','<<j2-jstep;
            std::cerr<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i2-istep,j1)) {
            ok = false;
            std::cerr<<"Lower left corner ("<<i2-istep<<','<<j1;
            std::cerr<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i2-istep,j2-jstep)) {
            ok = false;
            std::cerr<<"Lower right corner ("<<i2-istep<<','<<j2-jstep;
            std::cerr<<") must be in Upper Triangle\n";
        }
        return ok;
    }

    template <class T>
    bool GenUpperTriMatrix<T>::hasSubVector(
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
        if (!this->okij(i,j)) {
            ok = false;
            std::cerr<<"First element ("<<i<<','<<j<<") must be in Triangle\n";
        }
        if (!this->okij(i2,j2)) {
            ok = false;
            std::cerr<<"Last element ("<<i2<<','<<j2<<") must be in Triangle\n";
        }
        return ok;
    }

    template <class T>
    bool GenUpperTriMatrix<T>::hasSubTriMatrix(
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
        return ok;
    }

    template <class T>
    bool ConstUpperTriMatrixView<T,FortranStyle>::hasSubMatrix(
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
            std::cerr<<"last col element ("<<i2<<") must be in 1 -- ";
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
            std::cerr<<"last row element ("<<j2<<") must be in 1 -- ";
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
        if (!this->okij(i1-1,j1-1)) {
            ok = false;
            std::cerr<<"Upper left corner ("<<i1<<','<<j1;
            std::cerr<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i1-1,j2-1)) {
            ok = false;
            std::cerr<<"Upper right corner ("<<i1<<','<<j2;
            std::cerr<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i2-1,j1-1)) {
            ok = false;
            std::cerr<<"Lower left corner ("<<i2<<','<<j1;
            std::cerr<<") must be in Upper Triangle\n";
        }
        if (!this->okij(i2-1,j2-1)) {
            ok = false;
            std::cerr<<"Lower right corner ("<<i2<<','<<j2;
            std::cerr<<") must be in Upper Triangle\n";
        }
        return ok;
    }

    template <class T>
    bool ConstUpperTriMatrixView<T,FortranStyle>::hasSubVector(
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
        if (!this->okij(i-1,j-1)) {
            ok = false;
            std::cerr<<"First element ("<<i<<','<<j<<") must be in Triangle\n";
        }
        if (!this->okij(i2-1,j2-1)) {
            ok = false;
            std::cerr<<"Last element ("<<i2<<','<<j2<<") must be in Triangle\n";
        }
        return ok;
    }

    template <class T> 
    bool ConstUpperTriMatrixView<T,FortranStyle>::hasSubTriMatrix(
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
        if (i2<1 || i2 > this->size()) {
            ok = false;
            std::cerr<<"last diag element ("<<i2<<") must be in 1 -- ";
            std::cerr<<this->size()<<std::endl;
        }
        if ((i2-i1)%istep != 0) {
            ok = false;
            std::cerr<<"range ("<<i2-i1<<") must be multiple of istep (";
            std::cerr<<istep<<")\n";
        }
        if ((i2-i1)/istep < 0) {
            ok = false;
            std::cerr<<"n diag elements ("<<(i2-i1)/istep+1;
            std::cerr<<") must be positive\n";
        }
        return ok;
    }

    //
    // norms
    //

    template <class T>
    T GenUpperTriMatrix<T>::sumElements() const
    {
        const ptrdiff_t N = size();
        T sum(0);
        if (isrm()) 
            if (isunit())
                for(ptrdiff_t i=0;i<N;++i) 
                    sum += row(i,i+1,N).sumElements();
            else
                for(ptrdiff_t i=0;i<N;++i) 
                    sum += row(i,i,N).sumElements();
        else
            if (isunit())
                for(ptrdiff_t j=0;j<N;++j) 
                    sum += col(j,0,j).sumElements();
            else
                for(ptrdiff_t j=0;j<N;++j) 
                    sum += col(j,0,j+1).sumElements();
        if (isunit()) sum += RT(N);
        return sum;
    }

    template <class T>
    static RT DoSumAbsElements(const GenUpperTriMatrix<T>& m)
    {
        const ptrdiff_t N = m.size();
        RT sum(0);
        if (m.isrm()) 
            if (m.isunit())
                for(ptrdiff_t i=0;i<N;++i) 
                    sum += m.row(i,i+1,N).sumAbsElements();
            else
                for(ptrdiff_t i=0;i<N;++i) 
                    sum += m.row(i,i,N).sumAbsElements();
        else
            if (m.isunit())
                for(ptrdiff_t j=0;j<N;++j) 
                    sum += m.col(j,0,j).sumAbsElements();
            else
                for(ptrdiff_t j=0;j<N;++j) 
                    sum += m.col(j,0,j+1).sumAbsElements();
        if (m.isunit()) sum += RT(N);
        return sum;
    }

    template <class T>
    static RT DoSumAbs2Elements(const GenUpperTriMatrix<T>& m)
    {
        const ptrdiff_t N = m.size();
        RT sum(0);
        if (m.isrm()) 
            if (m.isunit())
                for(ptrdiff_t i=0;i<N;++i) 
                    sum += m.row(i,i+1,N).sumAbs2Elements();
            else
                for(ptrdiff_t i=0;i<N;++i) 
                    sum += m.row(i,i,N).sumAbs2Elements();
        else
            if (m.isunit())
                for(ptrdiff_t j=0;j<N;++j) 
                    sum += m.col(j,0,j).sumAbs2Elements();
            else
                for(ptrdiff_t j=0;j<N;++j) 
                    sum += m.col(j,0,j+1).sumAbs2Elements();
        if (m.isunit()) sum += RT(N);
        return sum;
    }

#ifdef INST_INT
    template <class T>
    static RT DoSumAbsElements(const GenMatrix<std::complex<T> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T>
    RT GenUpperTriMatrix<T>::sumAbsElements() const
    { return DoSumAbsElements(*this); }

    template <class T>
    RT GenUpperTriMatrix<T>::sumAbs2Elements() const
    { return DoSumAbs2Elements(*this); }

    template <class T>
    RT GenUpperTriMatrix<T>::normSq(const RT scale) const
    {
        const ptrdiff_t N = size();
        RT sum(0);
        if (isrm()) 
            if (isunit())
                for(ptrdiff_t i=0;i<N;++i) 
                    sum += row(i,i+1,N).normSq(scale);
            else
                for(ptrdiff_t i=0;i<N;++i) 
                    sum += row(i,i,N).normSq(scale);
        else
            if (isunit())
                for(ptrdiff_t j=0;j<N;++j) 
                    sum += col(j,0,j).normSq(scale);
            else
                for(ptrdiff_t j=0;j<N;++j) 
                    sum += col(j,0,j+1).normSq(scale);
        if (isunit()) {
            if (scale == RT(1)) sum += RT(N);
            else sum += RT(N) * scale * scale;
        }
        return sum;
    }

    template <class T> 
    static RT NonLapNormF(const GenUpperTriMatrix<T>& m)
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
    static RT NonLapMaxAbsElement(const GenUpperTriMatrix<T>& m)
    {
        const ptrdiff_t N = m.size();
        RT max(0);
        if (m.isrm()) {
            for(ptrdiff_t i=0;i<N;++i) {
                RT temp;
                if (m.isunit())
                    temp = m.row(i,i+1,N).normInf();
                else 
                    temp = m.row(i,i,N).normInf();
                if (temp > max) max = temp;
            }
        } else {
            for(ptrdiff_t j=0;j<N;++j) {
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
    static RT NonLapMaxAbs2Element(const GenUpperTriMatrix<T>& m)
    {
        const ptrdiff_t N = m.size();
        RT max(0);
        if (m.isrm()) {
            for(ptrdiff_t i=0;i<N;++i) {
                RT temp;
                if (m.isunit())
                    temp = m.row(i,i+1,N).maxAbs2Element();
                else 
                    temp = m.row(i,i,N).maxAbs2Element();
                if (temp > max) max = temp;
            }
        } else {
            for(ptrdiff_t j=0;j<N;++j) {
                RT temp;
                if (m.isunit())
                    temp = m.col(j,0,j).maxAbs2Element();
                else 
                    temp = m.col(j,0,j+1).maxAbs2Element();
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
        const ptrdiff_t N = m.size();
        for(ptrdiff_t j=0;j<N;++j) {
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
        const ptrdiff_t N = m.size();
        for(ptrdiff_t j=0;j<N;++j) {
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
        int lwork = cc=='I' ? M : 0;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        double norm = LAPNAME(dlantr) (
            LAPCM LAPV(cc),
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
        int lwork = cc=='I' ? M : 0;
        AlignedArray<double> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        double norm = LAPNAME(zlantr) (
            LAPCM LAPV(cc),
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
        int lwork = cc=='I' ? M : 0;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        float norm = LAPNAME(slantr) (
            LAPCM LAPV(cc),
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
        int lwork = cc=='I' ? M : 0;
        AlignedArray<float> work(lwork);
        VectorViewOf(work.get(),lwork).setZero();
#endif
        float norm = LAPNAME(clantr) (
            LAPCM LAPV(cc),
            m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
            LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
            LAP1 LAP1 LAP1);
        return norm;
    }
#endif
#endif // XLAP

#ifdef INST_INT
    static int NonLapNormF(const GenUpperTriMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNormF(const GenUpperTriMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNorm1(const GenUpperTriMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapNormInf(const GenUpperTriMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int NonLapMaxAbsElement(const GenUpperTriMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T>
    RT GenUpperTriMatrix<T>::maxAbsElement() const
    {
#ifdef XLAP
        if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
            return LapNorm('M',*this);
        else
#endif
            return NonLapMaxAbsElement(*this);
    }
    template <class T>
    RT GenUpperTriMatrix<T>::maxAbs2Element() const
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
    RT GenUpperTriMatrix<T>::norm1() const
    {
#ifdef XLAP
        if ((isrm() && stepi()>0) || (iscm() && stepj()>0)) 
            return LapNorm('1',*this);
        else 
#endif
            return NonLapNorm1(*this);
    }
    template <class T>
    RT GenUpperTriMatrix<T>::normInf() const
    {
#ifdef XLAP
        if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
            return LapNorm('I',*this);
        else
#endif
            return NonLapNormInf(*this);
    }
    template <class T>
    RT GenUpperTriMatrix<T>::normF() const
    {
#ifdef XLAP
        if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
            return LapNorm('F',*this);
        else
#endif
            return NonLapNormF(*this);
    }

    template <class T>
    static RT DoNorm2(const GenUpperTriMatrix<T>& m)
    { return Matrix<T>(m).doNorm2(); }

    template <class T>
    static RT DoCondition(const GenUpperTriMatrix<T>& m)
    { return Matrix<T>(m).doCondition(); }

#ifdef INST_INT
    static int DoNorm2(const GenUpperTriMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoCondition(const GenUpperTriMatrix<int>& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoNorm2(const GenUpperTriMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
    static int DoCondition(const GenUpperTriMatrix<std::complex<int> >& )
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T>
    RT GenUpperTriMatrix<T>::doNorm2() const
    { return tmv::DoNorm2(*this); }
    template <class T>
    RT GenUpperTriMatrix<T>::doCondition() const
    { return tmv::DoCondition(*this); }

    //
    // Modifying Functions
    //

    template <class T, int A>
    UpperTriMatrixView<T,A>& UpperTriMatrixView<T,A>::setZero()
    {
        const ptrdiff_t N = size();

        if (isrm())
            if (isunit())
                for(ptrdiff_t i=0;i<N;++i) row(i,i+1,N).setZero();
            else
                for(ptrdiff_t i=0;i<N;++i) row(i,i,N).setZero();
        else 
            if (isunit())
                for(ptrdiff_t j=0;j<N;++j) col(j,0,j).setZero();
            else
                for(ptrdiff_t j=0;j<N;++j) col(j,0,j+1).setZero();
        return *this; 
    } 

    template <class T, int A>
    UpperTriMatrixView<T,A>& UpperTriMatrixView<T,A>::setAllTo(const T& x)
    {
        const ptrdiff_t N = size();

        if (isrm())
            if (isunit())
                for(ptrdiff_t i=0;i<N;++i) row(i,i+1,N).setAllTo(x); 
            else
                for(ptrdiff_t i=0;i<N;++i) row(i,i,N).setAllTo(x); 
        else 
            if (isunit())
                for(ptrdiff_t j=0;j<N;++j) col(j,0,j).setAllTo(x); 
            else
                for(ptrdiff_t j=0;j<N;++j) col(j,0,j+1).setAllTo(x); 
        return *this; 
    }

    template <class T, int A>
    UpperTriMatrixView<T,A>& UpperTriMatrixView<T,A>::addToAll(const T& x) 
    {
        TMVAssert(!isunit());
        const ptrdiff_t N = size();

        if (isrm())
            for(ptrdiff_t i=0;i<N;++i) row(i,i,N).addToAll(x); 
        else 
            for(ptrdiff_t j=0;j<N;++j) col(j,0,j+1).addToAll(x); 
        return *this; 
    }

    template <class T, int A>
    UpperTriMatrixView<T,A>& UpperTriMatrixView<T,A>::clip(RT thresh) 
    {
        const ptrdiff_t N = size();

        if (isrm())
            if (isunit())
                for(ptrdiff_t i=0;i<N;++i) row(i,i+1,N).clip(thresh);
            else
                for(ptrdiff_t i=0;i<N;++i) row(i,i,N).clip(thresh);
        else 
            if (isunit())
                for(ptrdiff_t j=0;j<N;++j) col(j,0,j).clip(thresh);
            else
                for(ptrdiff_t j=0;j<N;++j) col(j,0,j+1).clip(thresh);
        return *this; 
    }

    template <class T, int A>
    UpperTriMatrixView<T,A>& UpperTriMatrixView<T,A>::conjugateSelf() 
    {
        const ptrdiff_t N = size();

        if (isComplex(T())) {
            if (isrm())
                if (isunit())
                    for(ptrdiff_t i=0;i<N;++i) row(i,i+1,N).conjugateSelf();
                else
                    for(ptrdiff_t i=0;i<N;++i) row(i,i,N).conjugateSelf();
            else
                if (isunit())
                    for(ptrdiff_t j=0;j<N;++j) col(j,0,j).conjugateSelf();
                else
                    for(ptrdiff_t j=0;j<N;++j) col(j,0,j+1).conjugateSelf();
        }
        return *this; 
    }

    template <class T, int A>
    UpperTriMatrixView<T,A>& UpperTriMatrixView<T,A>::setToIdentity(const T& x) 
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
    void Swap(UpperTriMatrixView<T> m1, UpperTriMatrixView<T> m2)
    {
        TMVAssert(m1.size() == m2.size());
        TMVAssert(m1.dt() == m2.dt());
        const ptrdiff_t N = m1.size();

        if (m1.isrm() && m2.isrm())
            if (m1.isunit()) 
                for(ptrdiff_t i=0;i<N;++i) 
                    Swap(m1.row(i,i+1,N), m2.row(i,i+1,N));
            else
                for(ptrdiff_t i=0;i<N;++i) 
                    Swap(m1.row(i,i,N), m2.row(i,i,N));
        else
            if (m1.isunit()) 
                for(ptrdiff_t j=0;j<N;++j) 
                    Swap(m1.col(j,0,j), m2.col(j,0,j));
            else
                for(ptrdiff_t j=0;j<N;++j) 
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
            const ptrdiff_t N = m1.size();
            for(ptrdiff_t j=0;j<N;++j) {
                if (m1.col(j,0,j) != m2.col(j,0,j)) return false;
            }

            if (m1.isunit() && !m2.isunit()) {
                for(ptrdiff_t i=0;i<N;++i) if (m2(i,i) != T2(1)) return false;
            } else if (m2.isunit() && !m1.isunit()) {
                for(ptrdiff_t i=0;i<N;++i) if (m1(i,i) != T1(1)) return false;
            } else if (!m1.isunit() && !m2.isunit()) {
                if (m1.diag() != m2.diag()) return false;
            }

            return true;  
        }
    }

    //
    // I/O
    //

    template <class T>
    void GenUpperTriMatrix<T>::write(const TMV_Writer& writer) const
    {
        const ptrdiff_t N = size();
        writer.begin();
        writer.writeCode("U");
        writer.writeSize(N);
        writer.writeSimpleSize(N);
        writer.writeStart();
        for(ptrdiff_t i=0;i<N;++i) {
            writer.writeLParen();
            if (!writer.isCompact()) {
                for(ptrdiff_t j=0;j<i;++j) {
                    writer.writeValue(T(0));
                    writer.writeSpace();
                }
            }
            for(ptrdiff_t j=i;j<N;++j) {
                if (j > i) writer.writeSpace();
                writer.writeValue(cref(i,j));
            }
            writer.writeRParen();
            if (i < N-1) writer.writeRowEnd();
        }
        writer.writeFinal();
        writer.end();
    }

#ifndef NOTHROW
    template <class T>
    class UpperTriMatrixReadError : public ReadError
    {
    public :
        UpperTriMatrix<T> m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t s;
        T v1;
        bool is, iseof, isbad;

        UpperTriMatrixReadError(std::istream& _is) throw() :
            ReadError("UpperTriMatrix"),
            i(0), j(0), s(0), v1(1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        UpperTriMatrixReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("UpperTriMatrix"),
            i(0), j(0), exp(_e), got(_g), s(0), v1(1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        UpperTriMatrixReadError(
            const GenUpperTriMatrix<T>& _m,
            std::istream& _is, ptrdiff_t _s) throw() :
            ReadError("UpperTriMatrix"),
            m(_m), i(0), j(0), s(_s), v1(1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        UpperTriMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenUpperTriMatrix<T>& _m,
            std::istream& _is, std::string _e, std::string _g) throw() :
            ReadError("UpperTriMatrix"),
            m(_m), i(_i), j(_j), exp(_e), got(_g), 
            s(m.size()), v1(i==j?T(1):T(0)),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        UpperTriMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenUpperTriMatrix<T>& _m,
            std::istream& _is) throw() :
            ReadError("UpperTriMatrix"),
            m(_m), i(_i), j(_j), s(m.size()), v1(i==j?T(1):T(0)),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        UpperTriMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenUpperTriMatrix<T>& _m,
            std::istream& _is, T _v1) throw() :
            ReadError("UpperTriMatrix"),
            m(_m), i(_i), j(_j), s(m.size()), v1(_v1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        UpperTriMatrixReadError(const UpperTriMatrixReadError<T>& rhs) throw():
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got),
            s(rhs.s), v1(rhs.v1),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~UpperTriMatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for UpperTriMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
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
            if (i != j && v1 != T(0)) {
                os<<"Invalid input: Expected 0, got "<<v1<<".\n";
            }
            if (i == j && v1 != T(1)) {
                os<<"Invalid input: Expected 1, got "<<v1<<".\n";
            }
            if (m.size() > 0) {
                os<<"The portion of the UpperTriMatrix which was successfully "
                    "read is:\n";
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
    static void FinishRead(const TMV_Reader& reader, UpperTriMatrixView<T> m)
    {
        const ptrdiff_t N = m.size();
        std::string exp, got;
        T temp;
        if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(0,0,m,reader.getis(),exp,got);
#endif
        }
        for(ptrdiff_t i=0;i<N;++i) {
            if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw UpperTriMatrixReadError<T>(i,0,m,reader.getis(),exp,got);
#endif
            }
            if (!reader.isCompact()) {
                for(ptrdiff_t j=0;j<i;++j) {
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"UpperTriMatrix Read Error: reading value\n";
                        exit(1);
#else
                        throw UpperTriMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    if (temp != T(0)) {
#ifdef NOTHROW
                        std::cerr<<"UpperTriMatrix Read Error: "<<temp<<" != 0\n";
                        exit(1);
#else
                        throw UpperTriMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                    }
                    if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw UpperTriMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                }
            }
            for(ptrdiff_t j=i;j<N;++j) {
                if (j>i && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw UpperTriMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                }
                if (!reader.readValue(temp)) {
#ifdef NOTHROW
                    std::cerr<<"UpperTriMatrix Read Error: reading value\n";
                    exit(1);
#else
                    throw UpperTriMatrixReadError<T>(i,j,m,reader.getis());
#endif
                }
                if (j==i && m.isunit()) {
                    if (temp != T(1)) {
#ifdef NOTHROW
                        std::cerr<<"UpperTriMatrix Read Error: "<<temp<<" != 1\n"; 
                        exit(1); 
#else
                        throw UpperTriMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                    }
                } else {
                    m.ref(i,j) = temp;
                }
            }
            if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw UpperTriMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
            if (i < N-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw UpperTriMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
        }
        if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(N,0,m,reader.getis(),exp,got);
#endif
        }
    }

    template <class T, int A>
    void UpperTriMatrix<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        if (!reader.readCode("U",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) resize(s);
        s=size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        UpperTriMatrixView<T> v = view();
        FinishRead(reader,v);
    }

    template <class T, int A>
    void UpperTriMatrixView<T,A>::read(const TMV_Reader& reader) 
    {
        std::string exp,got;
        if (!reader.readCode("U",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        s=size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw UpperTriMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        FinishRead(reader,*this);
    }

 
    template <class T>
    void GenLowerTriMatrix<T>::write(const TMV_Writer& writer) const
    {
        const ptrdiff_t N = size();
        writer.begin();
        writer.writeCode("L");
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
                    writer.writeValue(T(0));
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
    class LowerTriMatrixReadError : public ReadError
    {
    public :
        LowerTriMatrix<T> m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t s;
        T v1;
        bool is, iseof, isbad;

        LowerTriMatrixReadError(std::istream& _is) throw() :
            ReadError("LowerTriMatrix"),
            i(0), j(0), s(0), v1(1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        LowerTriMatrixReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("LowerTriMatrix"),
            i(0), j(0), exp(_e), got(_g), s(0), v1(1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        LowerTriMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenLowerTriMatrix<T>& _m,
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("LowerTriMatrix"),
            m(_m), i(_i), j(_j), exp(_e), got(_g), 
            s(m.size()), v1(i==j?T(1):T(0)),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        LowerTriMatrixReadError(
            const GenLowerTriMatrix<T>& _m,
            std::istream& _is, ptrdiff_t _s) throw() :
            ReadError("LowerTriMatrix"),
            m(_m), i(0), j(0), s(_s), v1(1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        LowerTriMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenLowerTriMatrix<T>& _m,
            std::istream& _is) throw() :
            ReadError("LowerTriMatrix"),
            m(_m), i(_i), j(_j), s(m.size()), v1(i==j?T(1):T(0)),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        LowerTriMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenLowerTriMatrix<T>& _m,
            std::istream& _is, T _v1) throw() :
            ReadError("LowerTriMatrix"),
            m(_m), i(_i), j(_j), s(m.size()), v1(_v1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        LowerTriMatrixReadError(const LowerTriMatrixReadError<T>& rhs) throw() :
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got),
            s(rhs.s), v1(rhs.v1),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        virtual ~LowerTriMatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for LowerTriMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
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
            if (i != j && v1 != T(0)) {
                os<<"Invalid input: Expected 0, got "<<v1<<".\n";
            }
            if (i == j && v1 != T(1)) {
                os<<"Invalid input: Expected 1, got "<<v1<<".\n";
            }
            if (m.size() > 0) {
                os<<"The portion of the LowerTriMatrix which was successfully "
                    "read is:\n";
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
    static void FinishRead(const TMV_Reader& reader, LowerTriMatrixView<T> m)
    {
        const ptrdiff_t N = m.size();
        std::string exp, got;
        T temp;
        if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(0,0,m,reader.getis(),exp,got);
#endif
        }
        for(ptrdiff_t i=0;i<N;++i) {
            if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw LowerTriMatrixReadError<T>(i,0,m,reader.getis(),exp,got);
#endif
            }
            for(ptrdiff_t j=0;j<i+1;++j) {
                if (j>0 && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw LowerTriMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                }
                if (!reader.readValue(temp)) {
#ifdef NOTHROW
                    std::cerr<<"LowerTriMatrix Read Error: reading value\n";
                    exit(1);
#else
                    throw LowerTriMatrixReadError<T>(i,j,m,reader.getis());
#endif
                }
                if (j==i && m.isunit()) {
                    if (temp != T(1)) {
#ifdef NOTHROW
                        std::cerr<<"LowerTriMatrix Read Error: "<<temp<<" != 1\n"; 
                        exit(1); 
#else
                        throw LowerTriMatrixReadError<T>(i,j,m,reader.getis(),temp);
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
                        std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw LowerTriMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"LowerTriMatrix Read Error: reading value\n";
                        exit(1);
#else
                        throw LowerTriMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    if (temp != T(0)) {
#ifdef NOTHROW
                        std::cerr<<"LowerTriMatrix Read Error: "<<temp<<" != 0\n";
                        exit(1);
#else
                        throw LowerTriMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                    }
                }
            }
            if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw LowerTriMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
            if (i < N-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw LowerTriMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
        }
        if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(N,0,m,reader.getis(),exp,got);
#endif
        }
    }

    template <class T, int A>
    void LowerTriMatrix<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        if (!reader.readCode("L",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) resize(s);
        s=size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        LowerTriMatrixView<T> v = view();
        FinishRead(reader,v);
    }

    template <class T, int A>
    void LowerTriMatrixView<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        if (!reader.readCode("L",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        s=size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"LowerTriMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw LowerTriMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        FinishRead(reader,*this);
    }

#undef RT

#define InstFile "TMV_TriMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


