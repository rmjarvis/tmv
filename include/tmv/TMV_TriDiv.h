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


//---------------------------------------------------------------------------
//
// This file contains the code for doing division of Triangular matrices.
//
// This is done using back or forward substitution.


#ifndef TMV_TriDiv_H
#define TMV_TriDiv_H

#include "tmv/TMV_BaseTriMatrix.h"

namespace tmv {

#ifndef NOTHROW
    template <typename T>
    class SingularUpperTriMatrix : public Singular
    {
    public:
        UpperTriMatrix<T> A;

        SingularUpperTriMatrix(const GenUpperTriMatrix<T>& _A) :
            Singular("UpperTriMatrix."), A(_A) {}
        ~SingularUpperTriMatrix() throw() {}
        void write(std::ostream& os) const throw()
        { Singular::write(os); os<<A<<std::endl; }
    };

    template <typename T>
    class SingularLowerTriMatrix : public Singular
    {
    public:
        LowerTriMatrix<T> A;

        SingularLowerTriMatrix(const GenLowerTriMatrix<T>& _A) :
            Singular("LowerTriMatrix."), A(_A) {}
        ~SingularLowerTriMatrix() throw() {}
        void write(std::ostream& os) const throw()
        { Singular::write(os); os<<A<<std::endl; }
    };
#endif

    template <typename T, typename T1>
    void TriLDivEq(const GenUpperTriMatrix<T1>& A, VectorView<T> v);

    template <typename T, typename T1>
    void TriLDivEq(const GenLowerTriMatrix<T1>& A, VectorView<T> v);

    template <typename T, typename T1>
    void TriLDivEq(const GenUpperTriMatrix<T1>& A, MatrixView<T> m);

    template <typename T, typename T1>
    void TriLDivEq(const GenLowerTriMatrix<T1>& A, MatrixView<T> m);

    template <typename T, typename T1>
    void TriLDivEq(
        const GenUpperTriMatrix<T1>& A, UpperTriMatrixView<T> m);

    template <typename T, typename T1>
    void TriLDivEq(
        const GenLowerTriMatrix<T1>& A, LowerTriMatrixView<T> m);

    template <typename T>
    void TriInverse(UpperTriMatrixView<T> minv);

} // namespace mv

#endif
