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
    void GenUpperTriMatrix<T>::doLDivEq(VectorView<T2> v) const
    { TriLDivEq(*this,v); }

    template <class T> template <class T2> 
    void GenUpperTriMatrix<T>::doLDivEq(MatrixView<T2> m) const
    { TriLDivEq(*this,m); }

    template <class T> template <class T2> 
    void GenUpperTriMatrix<T>::doLDivEq(UpperTriMatrixView<T2> m) const
    { TriLDivEq(*this,m); }

    template <class T> template <class T1, class T2> 
    void GenUpperTriMatrix<T>::doLDiv(
        const GenVector<T1>& v1, VectorView<T2> v2) const
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
        const GenMatrix<T1>& m1, MatrixView<T2> m2) const
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
        const GenUpperTriMatrix<T1>& m1, 
        UpperTriMatrixView<T2> m2) const
    {
        if (SameStorage(*this,m2)) {
            if (m2.isrm()) {
                UpperTriMatrix<T2,NonUnitDiag|RowMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            } else {
                UpperTriMatrix<T2,NonUnitDiag|ColMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            }
        } else {
            TriLDivEq(*this,m2=m1);
        }
    }

    template <class T> template <class T2> 
    void GenLowerTriMatrix<T>::doLDivEq(VectorView<T2> v) const
    { TriLDivEq(*this,v); }

    template <class T> template <class T2> 
    void GenLowerTriMatrix<T>::doLDivEq(MatrixView<T2> m) const
    { TriLDivEq(*this,m); }

    template <class T> template <class T2> 
    void GenLowerTriMatrix<T>::doLDivEq(LowerTriMatrixView<T2> m) const
    { TriLDivEq(*this,m); }

    template <class T> template <class T1, class T2> 
    void GenLowerTriMatrix<T>::doLDiv(
        const GenVector<T1>& v1, VectorView<T2> v2) const
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
        const GenMatrix<T1>& m1, MatrixView<T2> m2) const
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
        const GenLowerTriMatrix<T1>& m1,
        LowerTriMatrixView<T2> m2) const
    {
        if (SameStorage(*this,m2)) {
            if (m2.isrm()) {
                LowerTriMatrix<T2,NonUnitDiag|RowMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            } else {
                LowerTriMatrix<T2,NonUnitDiag|ColMajor> temp = m1;
                TriLDivEq(*this,temp.view());
                m2 = temp;
            }
        } else {
            TriLDivEq(*this,m2=m1);
        }
    }

    template <class T> 
    T GenUpperTriMatrix<T>::det() const
    { return isunit() ? T(1) : DiagMatrixViewOf(diag()).det(); }

    template <class T> 
    TMV_RealType(T) GenUpperTriMatrix<T>::logDet(T* s) const
    { 
        if (isunit()) {
            if (s) *s = T(1);
            return TMV_RealType(T)(0);
        } else {
            return DiagMatrixViewOf(diag()).logDet(s); 
        }                  
    }

    template <class T, int A>
    UpperTriMatrixView<T,A>& UpperTriMatrixView<T,A>::invertSelf()
    {
        TriInverse(*this);
        return *this;
    }

    template <class T> template <class T1> 
    void GenUpperTriMatrix<T>::doMakeInverse(MatrixView<T1> minv) const
    {
        bool ss = SameStorage(*this,minv);

        if (!ss) minv.setZero();

        UpperTriMatrixView<T1> U = minv.upperTri();
        U = *this;
        U.invertSelf();

        if (ss && minv.rowsize() > 1) minv.lowerTri().offDiag().setZero();
    }

    template <class T> template <class T1> 
    void GenUpperTriMatrix<T>::doMakeInverse(UpperTriMatrixView<T1> minv) const
    {
        minv = *this;
        minv.invertSelf();
    }

    template <class T> 
    void GenUpperTriMatrix<T>::doMakeInverseATA(MatrixView<T> minv) const
    {
        TMVAssert(minv.colsize() == size());
        TMVAssert(minv.rowsize() == size());

        UpperTriMatrixView<T> U = minv.upperTri();
        U = *this;
        U.invertSelf();
        minv = U * U.adjoint();
    }

    template <class T> 
    void GenLowerTriMatrix<T>::doMakeInverseATA(MatrixView<T> minv) const
    {
        TMVAssert(minv.colsize() == size());
        TMVAssert(minv.rowsize() == size());

        LowerTriMatrixView<T> L = minv.lowerTri();
        L = *this;
        L.invertSelf();
        minv = L * L.adjoint();
    }

#ifdef INST_INT
    template <> template <class T2> 
    void GenUpperTriMatrix<int>::doLDivEq(VectorView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T2> 
    void GenUpperTriMatrix<std::complex<int> >::doLDivEq(VectorView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T2> 
    void GenUpperTriMatrix<int>::doLDivEq(MatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T2> 
    void GenUpperTriMatrix<std::complex<int> >::doLDivEq(MatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T2> 
    void GenUpperTriMatrix<int>::doLDivEq(UpperTriMatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T2> 
    void GenUpperTriMatrix<std::complex<int> >::doLDivEq(
        UpperTriMatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T1, class T2> 
    void GenUpperTriMatrix<int>::doLDiv(
        const GenVector<T1>& , VectorView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1, class T2> 
    void GenUpperTriMatrix<std::complex<int> >::doLDiv(
        const GenVector<T1>& , VectorView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T1, class T2> 
    void GenUpperTriMatrix<int>::doLDiv(
        const GenMatrix<T1>& , MatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1, class T2> 
    void GenUpperTriMatrix<std::complex<int> >::doLDiv(
        const GenMatrix<T1>& , MatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T1, class T2> 
    void GenUpperTriMatrix<int>::doLDiv(
        const GenUpperTriMatrix<T1>& , 
        UpperTriMatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1, class T2> 
    void GenUpperTriMatrix<std::complex<int> >::doLDiv(
        const GenUpperTriMatrix<T1>& , 
        UpperTriMatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T2> 
    void GenLowerTriMatrix<int>::doLDivEq(VectorView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T2> 
    void GenLowerTriMatrix<std::complex<int> >::doLDivEq(
        VectorView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T2> 
    void GenLowerTriMatrix<int>::doLDivEq(MatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T2> 
    void GenLowerTriMatrix<std::complex<int> >::doLDivEq(
        MatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T2> 
    void GenLowerTriMatrix<int>::doLDivEq(LowerTriMatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T2> 
    void GenLowerTriMatrix<std::complex<int> >::doLDivEq(
        LowerTriMatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T1, class T2> 
    void GenLowerTriMatrix<int>::doLDiv(
        const GenVector<T1>& , VectorView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1, class T2> 
    void GenLowerTriMatrix<std::complex<int> >::doLDiv(
        const GenVector<T1>& , VectorView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T1, class T2> 
    void GenLowerTriMatrix<int>::doLDiv(
        const GenMatrix<T1>& , MatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1, class T2> 
    void GenLowerTriMatrix<std::complex<int> >::doLDiv(
        const GenMatrix<T1>& , MatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T1, class T2> 
    void GenLowerTriMatrix<int>::doLDiv(
        const GenLowerTriMatrix<T1>& ,
        LowerTriMatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1, class T2> 
    void GenLowerTriMatrix<std::complex<int> >::doLDiv(
        const GenLowerTriMatrix<T1>& ,
        LowerTriMatrixView<T2> ) const
    { TMVAssert(TMV_FALSE); }

    template <> 
    int GenUpperTriMatrix<int>::logDet(int* ) const
    { TMVAssert(TMV_FALSE); return 0; }
    template <> 
    int GenUpperTriMatrix<std::complex<int> >::logDet(
        std::complex<int>* ) const
    { TMVAssert(TMV_FALSE); return 0; }

    template <>
    UpperTriMatrixView<int,CStyle>&
    UpperTriMatrixView<int,CStyle>::invertSelf()
    { TMVAssert(TMV_FALSE); return *this; }
    template <>
    UpperTriMatrixView<std::complex<int>,CStyle>&
    UpperTriMatrixView<std::complex<int>,CStyle>::invertSelf()
    { TMVAssert(TMV_FALSE); return *this; }

    template <> template <class T1> 
    void GenUpperTriMatrix<int>::doMakeInverse(MatrixView<T1> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1> 
    void GenUpperTriMatrix<std::complex<int> >::doMakeInverse(
        MatrixView<T1> ) const
    { TMVAssert(TMV_FALSE); }

    template <> template <class T1> 
    void GenUpperTriMatrix<int>::doMakeInverse(UpperTriMatrixView<T1> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1> 
    void GenUpperTriMatrix<std::complex<int> >::doMakeInverse(
        UpperTriMatrixView<T1> ) const
    { TMVAssert(TMV_FALSE); }

    template <> 
    void GenUpperTriMatrix<int>::doMakeInverseATA(MatrixView<int> ) const
    { TMVAssert(TMV_FALSE); }
    template <> 
    void GenUpperTriMatrix<std::complex<int> >::doMakeInverseATA(
        MatrixView<std::complex<int> > ) const
    { TMVAssert(TMV_FALSE); }

    template <> 
    void GenLowerTriMatrix<int>::doMakeInverseATA(MatrixView<int> ) const
    { TMVAssert(TMV_FALSE); }
    template <> 
    void GenLowerTriMatrix<std::complex<int> >::doMakeInverseATA(
        MatrixView<std::complex<int> > ) const
    { TMVAssert(TMV_FALSE); }
#endif

#define InstFile "TMV_TriDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


