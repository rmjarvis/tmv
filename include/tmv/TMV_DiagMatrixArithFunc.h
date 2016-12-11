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


#ifndef TMV_DiagMatrixArithFunc_H
#define TMV_DiagMatrixArithFunc_H

#include "tmv/TMV_BaseDiagMatrix.h"
#include "tmv/TMV_VectorArithFunc.h"

#define CT std::complex<T>

namespace tmv {

    // y (+)= alpha * A * x
    template <bool add, typename T, typename Ta,  typename Tx>
    void MultMV(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<T> y);

    // B += alpha * A
    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, DiagMatrixView<T> B)
    { AddVV(alpha,A.diag(),B.diag()); }

    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, MatrixView<T> B)
    { AddVV(alpha,A.diag(),B.diag()); }

    // C = alpha * A + beta * B
    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const T beta, const GenDiagMatrix<Tb>& B, DiagMatrixView<T> C)
    { AddVV(alpha,A.diag(),beta,B.diag(),C.diag()); }

    template <typename T, typename Ta, typename Tb>
    void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<T> C);

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenDiagMatrix<Tb>& B, MatrixView<T> C)
    { AddMM(beta,B,alpha,A,C); }

    // C (+)= alpha * A * B
    template <typename T, typename Ta>
    void MultEqMM(
        const T alpha,
        const GenDiagMatrix<Ta>& A, MatrixView<T> B);

    template <bool add, typename T, typename Ta, typename Tb>
    void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<T> C);

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        MatrixView<T> C)
    { MultMM<add>(alpha,B,A.transpose(),C.transpose()); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        DiagMatrixView<T> C)
    { MultMV<add>(alpha,A,B.diag(),C.diag()); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        DiagMatrixView<T> C)
    { ElemMultVV<add>(alpha,A.diag(),B.diag(),C.diag()); }

    template <typename T>
    class DiagMatrixComposite : public GenDiagMatrix<T>
    {
    public:
        inline DiagMatrixComposite() : inst() {}
        inline DiagMatrixComposite(const DiagMatrixComposite<T>&) : inst() {}
        virtual inline ~DiagMatrixComposite() {}

    protected:
        ConstVectorView<T> cdiag() const;

    private:
        mutable auto_ptr<const DiagMatrix<T> > inst;
    };

    // Specialize allowed complex combinations:
    template <bool add, typename T, typename Ta,  typename Tx>
    inline void MultMV(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const GenVector<Tx>& x, VectorView<CT> y)
    { MultMV<add>(CT(alpha),A,x,y); }

    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, DiagMatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }

    template <typename T, typename Ta>
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A, MatrixView<CT> B)
    { AddMM(CT(alpha),A,B); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const T beta, const GenDiagMatrix<Tb>& B, DiagMatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenDiagMatrix<Ta>& A,
        const T beta, const GenMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const T alpha, const GenMatrix<Ta>& A,
        const T beta, const GenDiagMatrix<Tb>& B, MatrixView<CT> C)
    { AddMM(CT(alpha),A,CT(beta),B,C); }

    template <typename T, typename Ta>
    inline void MultEqMM(const T alpha,
                         const GenDiagMatrix<Ta>& A, MatrixView<CT> B)
    { MultEqMM(CT(alpha),A,B); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenMatrix<Tb>& B,
        MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        MatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        DiagMatrixView<CT> C)
    { MultMM<add>(CT(alpha),A,B,C); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const T alpha, const GenDiagMatrix<Ta>& A, const GenDiagMatrix<Tb>& B,
        DiagMatrixView<CT> C)
    { ElemMultMM<add>(CT(alpha),A,B,C); }

    // Specialize disallowed complex combinations:
    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMV(
        const CT , const GenDiagMatrix<Ta>& ,
        const GenVector<Tb>& , VectorView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenDiagMatrix<Ta>& ,
        const CT , const GenDiagMatrix<Ta>& , DiagMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenDiagMatrix<Ta>& ,
        const CT , const GenMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T, typename Ta, typename Tb>
    inline void AddMM(
        const CT , const GenMatrix<Ta>& ,
        const CT , const GenDiagMatrix<Tb>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenDiagMatrix<Ta>& , const GenMatrix<Tb>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenMatrix<Ta>& , const GenDiagMatrix<Tb>& ,
        MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void MultMM(
        const CT , const GenDiagMatrix<Ta>& , const GenDiagMatrix<Tb>& ,
        DiagMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <bool add, typename T, typename Ta, typename Tb>
    inline void ElemMultMM(
        const CT , const GenDiagMatrix<Ta>& , const GenDiagMatrix<Tb>& ,
        DiagMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

}

#undef CT

#endif
