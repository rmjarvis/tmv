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

#ifndef TMV_SymBandCHDiv_H
#define TMV_SymBandCHDiv_H

#include "tmv/TMV_BaseSymBandMatrix.h"

namespace tmv {

    template <typename T> 
    void LDL_Decompose(SymBandMatrixView<T> A);

    template <typename T, typename T1> 
    void CH_LDivEq(const GenSymBandMatrix<T1>& L, MatrixView<T> m);
    template <typename T, typename T1> 
    void CH_RDivEq(const GenSymBandMatrix<T1>& L, MatrixView<T> m);
    template <typename T, typename T1> 
    void CH_Inverse(
        const GenSymBandMatrix<T1>& LLx, SymMatrixView<T> sinv);

    template <typename T, typename T1> 
    void LDL_LDivEq(const GenSymBandMatrix<T1>& L, MatrixView<T> m);
    template <typename T, typename T1> 
    void LDL_RDivEq(const GenSymBandMatrix<T1>& L, MatrixView<T> m);
    template <typename T, typename T1> 
    void LDL_Inverse(
        const GenSymBandMatrix<T1>& LLx, SymMatrixView<T> sinv);

#ifndef NOTHROW
    template <typename T> 
    class NonPosDefHermBandMatrix : public NonPosDef
    {
    public:
        HermBandMatrix<T> A;

        inline NonPosDefHermBandMatrix(const GenSymBandMatrix<T>& _A) :
            NonPosDef("HermBandMatrix Cholesky decomposition"), 
            A(_A) {}
        inline NonPosDefHermBandMatrix(const NonPosDefHermBandMatrix<T>& rhs) :
            A(rhs.A) {}
        inline ~NonPosDefHermBandMatrix() throw() {}

        inline void write(std::ostream& os) const throw()
        {
            NonPosDef::write(os);
            os<<"The partially decomposed matrix is \n"<<A<<std::endl;
        }
    };

    template <typename T> 
    class NonPosDefHermBandMatrix2 : public NonPosDefHermBandMatrix<T>
    {
    public:
        HermBandMatrix<T> A0;

        inline NonPosDefHermBandMatrix2(
            const GenSymBandMatrix<T>& _A, const GenSymBandMatrix<T>& _A0) :
            NonPosDefHermBandMatrix<T>(_A), A0(_A0) {}
        inline NonPosDefHermBandMatrix2(
            const NonPosDefHermBandMatrix2<T>& rhs) :
            NonPosDefHermBandMatrix<T>(rhs), A0(rhs.A0) {}
        inline ~NonPosDefHermBandMatrix2() throw() {}

        inline void write(std::ostream& os) const throw()
        {
            NonPosDefHermBandMatrix<T>::write(os);
            os<<"The original matrix was \n"<<A0<<std::endl;
        }
    };

    template <typename T> 
    class NonPosDefSymBandLDL : public NonPosDef
    {
    public:
        SymBandMatrix<T> A;

        inline NonPosDefSymBandLDL(const GenSymBandMatrix<T>& _A) :
            NonPosDef("SymBandMatrix LDL decomposition."), A(_A) {}
        inline NonPosDefSymBandLDL(const NonPosDefSymBandLDL<T>& rhs) :
            A(rhs.A) {}
        inline ~NonPosDefSymBandLDL() throw() {}

        inline void write(std::ostream& os) const throw()
        {
            NonPosDef::write(os);
            os<<"The partially decomposed matrix is \n"<<A<<std::endl;
        }
    };

    template <typename T> 
    class NonPosDefHermBandLDL : public NonPosDef
    {
    public:
        HermBandMatrix<T> A;

        inline NonPosDefHermBandLDL(const GenSymBandMatrix<T>& _A) :
            NonPosDef("HermBandMatrix LDL decomposition."), A(_A) {}
        inline NonPosDefHermBandLDL(const NonPosDefHermBandLDL<T>& rhs) :
            A(rhs.A) {}
        inline ~NonPosDefHermBandLDL() throw() {}

        inline void write(std::ostream& os) const throw()
        {
            NonPosDef::write(os);
            os<<"The partially decomposed matrix is \n"<<A<<std::endl;
        }
    };
#endif

    // Specialize disallowed complex combinations:
#define CT std::complex<T>

    template <typename T>
    inline void CH_LDivEq(const GenSymBandMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void CH_RDivEq(const GenSymBandMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void CH_Inverse(
        const GenSymBandMatrix<CT>& , SymMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void LDL_LDivEq(
        const GenSymBandMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void LDL_RDivEq(
        const GenSymBandMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void LDL_Inverse(
        const GenSymBandMatrix<CT>& , SymMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

#undef CT

}

#endif
