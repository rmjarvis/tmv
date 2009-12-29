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



#include "tmv/TMV_TriDiv.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_TriMatrixArith.h"

namespace tmv {

    template <class T> 
    QuotXU<T,T> GenUpperTriMatrix<T>::QInverse() const
    { return QuotXU<T,T>(T(1),*this); }

    template <class T> 
    QuotXL<T,T> GenLowerTriMatrix<T>::QInverse() const
    { return QuotXL<T,T>(T(1),*this); }

    template <class T> template <class T2> 
    void GenUpperTriMatrix<T>::doLDivEq(const VectorView<T2>& v) const
    { TriLDivEq(*this,v); }

    template <class T> template <class T2> 
    void GenUpperTriMatrix<T>::doLDivEq(const MatrixView<T2>& m) const
    { TriLDivEq(*this,m); }

    template <class T> template <class T2> 
    void GenUpperTriMatrix<T>::doLDivEq(const UpperTriMatrixView<T2>& m) const
    { TriLDivEq(*this,m); }

    template <class T> template <class T1, class T2> 
    void GenUpperTriMatrix<T>::doLDiv(
        const GenVector<T1>& v1, const VectorView<T2>& v2) const
    {
        if (SameStorage(*this,v1)) {
            Vector<T2> temp = v1;
            TriLDivEq(*this,temp.view());
            v2 = temp;
        } else {
            TriLDivEq(*this,v2=v1);
        }
    }

    template <class T> template <class T1, class T2> 
    void GenUpperTriMatrix<T>::doLDiv(
        const GenMatrix<T1>& m1, const MatrixView<T2>& m2) const
    {
        if (SameStorage(*this,m2)) {
            if (m2.isrm()) {
                Matrix<T2,RowMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            } else {
                Matrix<T2,ColMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            }
        } else {
            TriLDivEq(*this,m2=m1);
        }
    }

    template <class T> template <class T1, class T2> 
    void GenUpperTriMatrix<T>::doLDiv(
        const GenUpperTriMatrix<T1>& m1, const UpperTriMatrixView<T2>& m2) const
    {
        if (SameStorage(*this,m2)) {
            if (m2.isrm()) {
                UpperTriMatrix<T2,NonUnitDiag,RowMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            } else {
                UpperTriMatrix<T2,NonUnitDiag,ColMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            }
        } else {
            TriLDivEq(*this,m2=m1);
        }
    }

    template <class T> template <class T2> 
    void GenLowerTriMatrix<T>::doLDivEq(const VectorView<T2>& v) const
    { TriLDivEq(*this,v); }

    template <class T> template <class T2> 
    void GenLowerTriMatrix<T>::doLDivEq(const MatrixView<T2>& m) const
    { TriLDivEq(*this,m); }

    template <class T> template <class T2> 
    void GenLowerTriMatrix<T>::doLDivEq(const LowerTriMatrixView<T2>& m) const
    { TriLDivEq(*this,m); }

    template <class T> template <class T1, class T2> 
    void GenLowerTriMatrix<T>::doLDiv(
        const GenVector<T1>& v1, const VectorView<T2>& v2) const
    {
        if (SameStorage(*this,v2)) {
            Vector<T2> temp = v1;
            TriLDivEq(*this,temp.view());
            v2 = temp;
        } else {
            TriLDivEq(*this,v2=v1);
        }
    }

    template <class T> template <class T1, class T2> 
    void GenLowerTriMatrix<T>::doLDiv(
        const GenMatrix<T1>& m1, const MatrixView<T2>& m2) const
    {
        if (SameStorage(*this,m2)) {
            if (m2.isrm()) {
                Matrix<T2,RowMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            } else {
                Matrix<T2,ColMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            }
        } else {
            TriLDivEq(*this,m2=m1);
        }
    }

    template <class T> template <class T1, class T2> 
    void GenLowerTriMatrix<T>::doLDiv(
        const GenLowerTriMatrix<T1>& m1, const LowerTriMatrixView<T2>& m2) const
    {
        if (SameStorage(*this,m2)) {
            if (m2.isrm()) {
                LowerTriMatrix<T2,NonUnitDiag,RowMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            } else {
                LowerTriMatrix<T2,NonUnitDiag,ColMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            }
        } else {
            TriLDivEq(*this,m2=m1);
        }
    }

    template <class T> T GenUpperTriMatrix<T>::det() const
    { return isunit() ? T(1) : DiagMatrixViewOf(diag()).det(); }                  

    template <class T> TMV_RealType(T) GenUpperTriMatrix<T>::logDet(T* s) const
    { 
        if (isunit()) {
            if (s) *s = T(1);
            return TMV_RealType(T)(0);
        } else {
            return DiagMatrixViewOf(diag()).logDet(s); 
        }                  
    }

    template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::invertSelf() const
    {
        TriInverse(*this);
        return *this;
    }

    template <class T> template <class T1> 
    void GenUpperTriMatrix<T>::doMakeInverse(const MatrixView<T1>& minv) const
    {
        bool ss = SameStorage(*this,minv);

        if (!ss) minv.zero();

        UpperTriMatrixView<T1> U = minv.upperTri();
        U = *this;
        U.invertSelf();

        if (ss && minv.rowsize() > 1) minv.lowerTri().offDiag().zero();
    }

    template <class T> template <class T1> 
    void GenUpperTriMatrix<T>::doMakeInverse(
        const UpperTriMatrixView<T1>& minv) const
    {
        minv = *this;
        minv.invertSelf();
    }

    template <class T> 
    void GenUpperTriMatrix<T>::doMakeInverseATA(const MatrixView<T>& minv) const
    {
        TMVAssert(minv.colsize() == size());
        TMVAssert(minv.rowsize() == size());

        UpperTriMatrixView<T> U = minv.upperTri();
        U = *this;
        U.invertSelf();
        minv = U * U.adjoint();
    }

    template <class T> 
    void GenLowerTriMatrix<T>::doMakeInverseATA(const MatrixView<T>& minv) const
    {
        TMVAssert(minv.colsize() == size());
        TMVAssert(minv.rowsize() == size());

        LowerTriMatrixView<T> L = minv.lowerTri();
        L = *this;
        L.invertSelf();
        minv = L * L.adjoint();
    }

#define InstFile "TMV_TriDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


