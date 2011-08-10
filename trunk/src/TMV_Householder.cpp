
//#define XDEBUG_HOUSE

#include "tmv/TMV_Householder.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_MinMax.h"
#include "tmv/TMV_NormV.h"
#include "tmv/TMV_MultMM.h"
#include "tmv/TMV_MultUM.h"
#include "tmv/TMV_MultUU.h"
#include "tmv/TMV_MultUL.h"
#include "tmv/TMV_MultXM.h"

namespace tmv {

    template <class T>
    void InstHouseholderReflect(
        T& x0, VectorView<T> u, typename Traits<T>::real_type& beta)
    { InlineHouseholderReflect(x0,u,beta); }

    template <class T>
    bool InstHouseholderUnReflect(
        T& y, VectorView<T> u, typename Traits<T>::real_type& beta)
    { return InlineHouseholderUnReflect(y,u,beta); }

    template <class T>
    void InstHouseholderUnpack(
        VectorView<T> u, typename Traits<T>::real_type& beta, T& u0)
    { InlineHouseholderUnpack(u,beta,u0); }

    template <class T1, int C1, class T2>
    void InstHouseholderMultEq(
        const ConstVectorView<T1,C1>& u, typename Traits<T1>::real_type beta,
        VectorView<T2> m0, MatrixView<T2> mx, VectorView<T2> temp)
    { InlineHouseholderMultEq(u,beta,m0,mx,temp); }

    template <class T>
    void InstBlockHouseholderAugment(
        const ConstMatrixView<T>& Y, UpperTriMatrixView<T> Z,
        typename Traits<T>::real_type beta)
    {
        UpperTriMatrixView<T,NonUnitDiag> Z2 = Z;
        InlineBlockHouseholderAugment(Y,Z2,beta); 
    }

    template <class T>
    void InstBlockHouseholderMakeZ(
        const ConstMatrixView<T>& Y, UpperTriMatrixView<T> Z,
        const ConstVectorView<typename Traits<T>::real_type>& beta)
    {
        UpperTriMatrixView<T,NonUnitDiag> Z2 = Z;
        InlineBlockHouseholderMakeZ(Y,Z2,beta); 
    }

    template <class T1, int C1, class T2>
    void InstBlockHouseholderLMult(
        const ConstMatrixView<T1,C1>& Y,
        const ConstUpperTriMatrixView<T1,C1>& Z,
        MatrixView<T2> m2, MatrixView<T2> temp)
    { InlineBlockHouseholderLMult(Y,Z.viewAsNonUnitDiag(),m2,temp); }

    template <class T1, int C1, class T2>
    void InstBlockHouseholderLDiv(
        const ConstMatrixView<T1,C1>& Y,
        const ConstUpperTriMatrixView<T1,C1>& Z,
        MatrixView<T2> m2, MatrixView<T2> temp)
    { InlineBlockHouseholderLDiv(Y,Z.viewAsNonUnitDiag(),m2,temp); }

    template <class T>
    void InstBlockHouseholderUnpack(
        MatrixView<T> Y, const ConstUpperTriMatrixView<T>& Z,
        UpperTriMatrixView<T> temp)
    { 
        UpperTriMatrixView<T,NonUnitDiag> t2 = temp;
        InlineBlockHouseholderUnpack(Y,Z.viewAsNonUnitDiag(),t2);
    }

    template <class T>
    void InstBlock2HouseholderAugment(
        const ConstMatrixView<T>& Y, UpperTriMatrixView<T> Z,
        typename Traits<T>::real_type beta)
    {
        UpperTriMatrixView<T,NonUnitDiag> Z2 = Z;
        InlineBlock2HouseholderAugment(Y,Z2,beta); 
    }

    template <class T>
    void InstBlock2HouseholderMakeZ(
        const ConstMatrixView<T>& Y, UpperTriMatrixView<T> Z,
        const ConstVectorView<typename Traits<T>::real_type>& beta)
    {
        UpperTriMatrixView<T,NonUnitDiag> Z2 = Z;
        InlineBlock2HouseholderMakeZ(Y,Z2,beta); 
    }

    template <class T1, int C1, class T2, int C2>
    void InstBlock2HouseholderLMult(
        const ConstMatrixView<T1,C1>& Y,
        const ConstUpperTriMatrixView<T1,C1>& Z,
        MatrixView<T2,C2> ma, MatrixView<T2> mb, MatrixView<T2> temp)
    { InlineBlock2HouseholderLMult(Y,Z.viewAsNonUnitDiag(),ma,mb,temp); }

    template <class T1, int C1, class T2, int C2>
    void InstBlock2HouseholderLDiv(
        const ConstMatrixView<T1,C1>& Y,
        const ConstUpperTriMatrixView<T1,C1>& Z,
        MatrixView<T2,C2> ma, MatrixView<T2> mb, MatrixView<T2> temp)
    { InlineBlock2HouseholderLDiv(Y,Z.viewAsNonUnitDiag(),ma,mb,temp); }

 

#define InstFile "TMV_Householder.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


