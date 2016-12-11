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


#ifndef TMV_BandMatrixArithFunc_H
#define TMV_BandMatrixArithFunc_H

#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_Array.h"

#define CT std::complex<T>

namespace tmv {

    // y (+)= alpha * A * x
    template <bool add, typename T, typename Ta, typename Tx>
    void MultMV(
        const T alpha, const GenBandMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y);

    // A *= alpha
    template <typename T>
    void MultXM(const T alpha, BandMatrixView<T> A);

    // B += alpha * A
    template <typename T, typename Ta>
    void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A, BandMatrixView<T> B);
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A, MatrixView<T> B)
    { AddMM(alpha,A,BandMatrixView<T>(B,A.nlo(),A.nhi())); }
    // C = alpha * A + beta * B
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const T beta, const GenBandMatrix<Tb>& B, BandMatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C);
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenBandMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B,alpha,A,C); }

    // C (+)= alpha * A * B
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const GenMatrix<Tb>& B, MatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const GenBandMatrix<Tb>& B, BandMatrixView<T> C);
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A,
        const GenBandMatrix<Tb>& B, MatrixView<T> C)
    { MultMM<add>(alpha,B.transpose(),A.transpose(),C.transpose()); }

    template <bool add, typename T, typename Ta, typename Tb>
    void ElemMultMM(
        const T alpha, const GenBandMatrix<Ta>& A,
        const GenBandMatrix<Tb>& B, BandMatrixView<T> C);

    template <typename T>
    class BandMatrixComposite : public GenBandMatrix<T>
    {
    public:

        inline BandMatrixComposite() : itsm(0) {}
        inline BandMatrixComposite(const BandMatrixComposite<T>&) : itsm(0) {}
        virtual inline ~BandMatrixComposite() {}

        // Definitions are in TMV_MultBV.cpp
        const T* cptr() const;
        ptrdiff_t stepi() const;
        ptrdiff_t stepj() const;
        ptrdiff_t diagstep() const;

        inline ConjType ct() const { return NonConj; }
        inline bool isconj() const { return false; }
        ptrdiff_t ls() const;

        ConstVectorView<T> constLinearView() const;
        inline bool canLinearize() const { return true; }

    private:

        mutable AlignedArray<T> itsm1;
        mutable T* itsm;

        BandMatrixComposite<T>& operator=(const BandMatrixComposite<T>&);
    };

    // Specialize allowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tx>
    inline void MultMV(
        const T alpha, const GenBandMatrix<Ta>& A, const GenVector<Tx>& x,
        VectorView<CT> y)
    { MultMV<add>(CT(alpha),A,x,y); }

    template <typename T>
    inline void MultXM(
        const T x, BandMatrixView<CT> m)
    { MultXM(CT(x),m); }

    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A, BandMatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A, MatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A, const T beta,
        const GenBandMatrix<Tb>& B, BandMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenBandMatrix<Ta>& A, const T beta,
        const GenMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A, const T beta,
        const GenBandMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }
    template <typename T>
    inline void AddMM(
        const CT alpha, const GenBandMatrix<CT>& A, const CT beta,
        const GenBandMatrix<T>& B, BandMatrixView<CT> C)
    { AddMM(beta,B,alpha,A,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
        BandMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
        MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const T alpha, const GenBandMatrix<Ta>& A, const GenBandMatrix<Tb>& B,
        BandMatrixView<CT> C)
    { ElemMultMM<add>(CT(alpha),A,B,C); }

    // Specialize disallowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMV(
        const CT , const GenBandMatrix<Ta>& ,
        const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void MultXM(const CT , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename Ta>
    inline void AddMM(
        const CT , const GenBandMatrix<Ta>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta>
    inline void AddMM(
        const CT , const GenBandMatrix<Ta>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenBandMatrix<Ta>& ,
        const CT , const GenBandMatrix<Tb>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenBandMatrix<Ta>& ,
        const CT , const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenMatrix<Ta>& ,
        const CT , const GenBandMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenBandMatrix<Ta>& ,
        const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenBandMatrix<Ta>& ,
        const GenBandMatrix<Tb>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenMatrix<Ta>& ,
        const GenBandMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const CT , const GenBandMatrix<Ta>& ,
        const GenBandMatrix<Tb>& , BandMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

}

#undef CT

#endif
