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

#ifndef TMV_SymCHDiv_H
#define TMV_SymCHDiv_H

#include "tmv/TMV_BaseSymMatrix.h"

namespace tmv {

    template <typename T, typename T1> 
    void CH_LDivEq(const GenSymMatrix<T1>& L, MatrixView<T> m);

    template <typename T, typename T1> 
    void CH_RDivEq(const GenSymMatrix<T1>& L, MatrixView<T> m);

    template <typename T, typename T1> 
    void CH_Inverse(const GenSymMatrix<T1>& LLx, SymMatrixView<T> sinv);

#ifndef NOTHROW
    template <typename T> 
    class NonPosDefHermMatrix : public NonPosDef
    {
    public:

        HermMatrix<T> A;

        NonPosDefHermMatrix(const GenSymMatrix<T>& _A) :
            NonPosDef("HermMatrix Cholesky decmposition."), 
            A(_A) {}
        NonPosDefHermMatrix(const NonPosDefHermMatrix<T>& rhs) :
            A(rhs.A) {}
        ~NonPosDefHermMatrix() throw() {}

        void write(std::ostream& os) const throw()
        {
            NonPosDef::write(os);
            os<<"The partially decomposed matrix is \n"<<A<<std::endl;
        }
    };

    template <typename T> 
    class NonPosDefHermMatrix2 : public NonPosDefHermMatrix<T>
    {
    public:
        HermMatrix<T> A0;

        NonPosDefHermMatrix2(
            const GenSymMatrix<T>& _A, const GenSymMatrix<T>& _A0) :
            NonPosDefHermMatrix<T>(_A), A0(_A0) {}
        NonPosDefHermMatrix2(const NonPosDefHermMatrix2<T>& rhs) :
            NonPosDefHermMatrix<T>(rhs), A0(rhs.A0) {}
        ~NonPosDefHermMatrix2() throw() {}

        void write(std::ostream& os) const throw()
        {
            NonPosDefHermMatrix<T>::write(os);
            os<<"The original matrix was \n"<<A0<<std::endl;
        }
    };
#endif

    // Specialize disallowed complex combinations:
#define CT std::complex<T>

    template <typename T>
    inline void CH_LDivEq(const GenSymMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void CH_RDivEq(const GenSymMatrix<CT>& , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void CH_Inverse(const GenSymMatrix<CT>& , SymMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

#undef CT


}

#endif
