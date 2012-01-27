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
#ifndef TMV_BandLUDiv_H
#define TMV_BandLUDiv_H

#include "tmv/TMV_BaseBandMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"

namespace tmv {


    template <class T, class T1> 
    void TriLDivEq(
        const GenBandMatrix<T1>& A, MatrixView<T> B, DiagType dt);
    template <class T, class T1> 
    void TriLDivEq(
        const GenBandMatrix<T1>& A, VectorView<T> v, DiagType dt);
    // Solve A X = B  where A is upper or lower band triangular

    template <class T, class T1> 
    void LU_PackedPL_Unpack(
        const GenBandMatrix<T1>& LUx, const int* p,
        LowerTriMatrixView<T> m);

    template <class T, class T1> 
    void LU_PackedPL_LDivEq(
        const GenBandMatrix<T1>& LUx, const int* p, MatrixView<T> m);
    template <class T, class T1> 
    void LU_PackedPL_RDivEq(
        const GenBandMatrix<T1>& LUx, const int* p, MatrixView<T> m);

    template <class T, class T1> 
    void LU_LDivEq(
        const GenBandMatrix<T1>& LUx, const int* p, MatrixView<T> m);
    template <class T, class T1> 
    void LU_RDivEq(
        const GenBandMatrix<T1>& LUx, const int* p, MatrixView<T> m);

    template <class T, class T1> 
    void LU_Inverse(
        const GenBandMatrix<T1>& LUx, const int* p, MatrixView<T> m);

    template <class T> 
    void TriInverse(UpperTriMatrixView<T> U, int nhi);

#ifndef NOTHROW
    template <class T> 
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

    template <class T> 
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
    template <class T>
    inline void TriLDivEq(
        const GenBandMatrix<CT>& , MatrixView<T> , DiagType )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void TriLDivEq(
        const GenBandMatrix<CT>& , VectorView<T> , DiagType )
    { TMVAssert(TMV_FALSE); }

    template <class T>
    inline void LU_PackedPL_Unpack(
        const GenBandMatrix<CT>& , const int* ,
        LowerTriMatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T>
    inline void LU_PackedPL_LDivEq(
        const GenBandMatrix<CT>& , const int* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void LU_PackedPL_RDivEq(
        const GenBandMatrix<CT>& , const int* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T>
    inline void LU_LDivEq(
        const GenBandMatrix<CT>& , const int* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }
    template <class T>
    inline void LU_RDivEq(
        const GenBandMatrix<CT>& , const int* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

    template <class T>
    inline void LU_Inverse(
        const GenBandMatrix<CT>& , const int* , MatrixView<T> )
    { TMVAssert(TMV_FALSE); }

#undef CT

}

#endif
