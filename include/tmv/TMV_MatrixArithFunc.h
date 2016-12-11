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


#ifndef TMV_MatrixArithFunc_H
#define TMV_MatrixArithFunc_H

#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_Array.h"

#define CT std::complex<T>

namespace tmv {

    // y (+)= alpha * A * x
    template <bool add, typename T, typename Ta, typename Tx>
    void MultMV(
        const T alpha, const GenMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<T> y);

    // A *= alpha
    template <typename T>
    void MultXM(const T alpha, MatrixView<T> A);

    // B += alpha * A
    template <typename T, typename Ta>
    void AddMM(const T alpha, const GenMatrix<Ta>& A, MatrixView<T> B);

    // C = alpha * A + beta * B
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C);

    // C (+)= alpha * A * B
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C);

    // A (+)= alpha * x * yT
    template <bool add, typename T, typename Tx, typename Ty>
    void Rank1Update(
        const T alpha, const GenVector<Tx>& x, const GenVector<Ty>& y,
        MatrixView<T> A);

    template <bool add, typename T, typename Ta, typename Tb>
    void ElemMultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<T> C);

    template <typename T>
    class MatrixComposite : public GenMatrix<T>
    {
    public:

        inline MatrixComposite() {}
        inline MatrixComposite(const MatrixComposite<T>&) {}
        virtual inline ~MatrixComposite() {}
        const T* cptr() const;
        ptrdiff_t stepi() const;
        ptrdiff_t stepj() const;
        inline ConjType ct() const { return NonConj; }
        inline bool canLinearize() const { return true; }
        ptrdiff_t ls() const;

    private:
        mutable AlignedArray<T> itsm;
    };

    // Specialize allowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tx>
    inline void MultMV(
        const T alpha, const GenMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<CT> y)
    { MultMV<add>(CT(alpha),A,x,y); }

    template <typename T>
    inline void MultXM(const T alpha, MatrixView<CT> A)
    { MultXM(CT(alpha),A); }

    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A, MatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T>
    inline void AddMM(
        const CT alpha, const GenMatrix<CT>& A,
        const CT beta, const GenMatrix<T>& B, MatrixView<CT> C)
    { AddMM(beta,B,alpha,A,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Tx, typename Ty>
    inline void Rank1Update(
        const T alpha, const GenVector<Tx>& x,
        const GenVector<Ty>& y, MatrixView<CT> A)
    { Rank1Update<add>(CT(alpha),x,y,A); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<CT> C)
    { ElemMultMM<add>(CT(alpha),A,B,C); }

    // Specialize disallowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMV(
        const CT , const GenMatrix<Ta>& ,
        const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void MultXM(const CT , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename Ta>
    inline void AddMM(const CT , const GenMatrix<Ta>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenMatrix<Ta>& ,
        const CT , const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenMatrix<Ta>& , const GenMatrix<Tb>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const CT , const GenMatrix<Ta>& , const GenMatrix<Tb>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void Rank1Update(
        const CT , const GenVector<Ta>& ,
        const GenVector<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

} // namespace tmv

#undef CT

#endif
