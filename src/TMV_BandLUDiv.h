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

#ifndef TMV_BandLUDiv_H
#define TMV_BandLUDiv_H

#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"

namespace tmv {


    template <typename T, typename T1> 
    void TriLDivEq(
        const GenBandMatrix<T1>& A, MatrixView<T> B, DiagType dt);
    template <typename T, typename T1> 
    void TriLDivEq(
        const GenBandMatrix<T1>& A, VectorView<T> v, DiagType dt);
    // Solve A X = B  where A is upper or lower band triangular

    template <typename T, typename T1> 
    void LU_PackedPL_Unpack(
        const GenBandMatrix<T1>& LUx, const ptrdiff_t* p,
        LowerTriMatrixView<T> m);

    template <typename T, typename T1> 
    void LU_PackedPL_LDivEq(
        const GenBandMatrix<T1>& LUx, const ptrdiff_t* p, MatrixView<T> m);
    template <typename T, typename T1> 
    void LU_PackedPL_RDivEq(
        const GenBandMatrix<T1>& LUx, const ptrdiff_t* p, MatrixView<T> m);

    template <typename T, typename T1> 
    void LU_LDivEq(
        const GenBandMatrix<T1>& LUx, const ptrdiff_t* p, MatrixView<T> m);
    template <typename T, typename T1> 
    void LU_RDivEq(
        const GenBandMatrix<T1>& LUx, const ptrdiff_t* p, MatrixView<T> m);

    template <typename T, typename T1> 
    void LU_Inverse(
        const GenBandMatrix<T1>& LUx, const ptrdiff_t* p, MatrixView<T> m);

    template <typename T> 
    void TriInverse(UpperTriMatrixView<T> U, ptrdiff_t nhi);

#ifndef NOTHROW
    template <typename T> 
    class SingularBandLU : public Singular
    {
    public:

        BandMatrix<T> A;

        SingularBandLU(const GenBandMatrix<T>& _A) :
            Singular("BandMatrix."), A(_A) {}
        ~SingularBandLU() throw() {}

        void write(std::ostream& os) const throw()
        {
            Singular::write(os);
            os<<"In LU Decomposed form, the matrix is \n"<<A<<std::endl;
            os<<"ie. U = "<<A.upperBand()<<std::endl;
        }
    };

    template <typename T> 
    class SingularBandLU2 : public SingularBandLU<T>
    {
    public:

        BandMatrix<T> A0;

        SingularBandLU2(
            const GenBandMatrix<T>& _A, const GenBandMatrix<T>& _A0) :
            SingularBandLU<T>(_A), A0(_A0) {}
        ~SingularBandLU2() throw() {}

        void write(std::ostream& os) const throw()
        {
            SingularBandLU<T>::write(os);
            os<<"The original BandMatrix is \n"<<A0<<std::endl;
        }
    };
#endif


    // Specialize disallowed complex combinations:
#define CT std::complex<T>
    template <typename T>
    inline void TriLDivEq(
        const GenBandMatrix<CT>& , MatrixView<T> , DiagType )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void TriLDivEq(
        const GenBandMatrix<CT>& , VectorView<T> , DiagType )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void LU_PackedPL_Unpack(
        const GenBandMatrix<CT>& , const ptrdiff_t* ,
        LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void LU_PackedPL_LDivEq(
        const GenBandMatrix<CT>& , const ptrdiff_t* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void LU_PackedPL_RDivEq(
        const GenBandMatrix<CT>& , const ptrdiff_t* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void LU_LDivEq(
        const GenBandMatrix<CT>& , const ptrdiff_t* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <typename T>
    inline void LU_RDivEq(
        const GenBandMatrix<CT>& , const ptrdiff_t* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <typename T>
    inline void LU_Inverse(
        const GenBandMatrix<CT>& , const ptrdiff_t* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

#undef CT

}

#endif
