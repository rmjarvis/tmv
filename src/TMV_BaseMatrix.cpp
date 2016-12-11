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



#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Divider.h"

namespace tmv {

    template <class T>
    DivHelper<T>::DivHelper() : divider(), divtype(tmv::XX) {}

    template <class T>
    DivHelper<T>::~DivHelper() {}

    template <class T>
    T DivHelper<T>::doDet() const
    {
        TMVAssert(colsize() == rowsize());
        setDiv();
        T det = divider->det();
        doneDiv();
        return det;
    }

    template <class T>
    TMV_RealType(T) DivHelper<T>::doLogDet(T* sign) const
    {
        TMVAssert(colsize() == rowsize());
        setDiv();
        TMV_RealType(T) logdet = divider->logDet(sign);
        doneDiv();
        return logdet;
    }

    template <class T, class T1> 
    static inline void DoMakeInverse1(
        const Divider<T>& div, MatrixView<T1> minv)
    { div.makeInverse(minv); }

    template <class T> template <class T1> 
    void DivHelper<T>::doMakeInverse(MatrixView<T1> minv) const
    {
        TMVAssert(minv.colsize() == rowsize());
        TMVAssert(minv.rowsize() == colsize());
        setDiv();
        DoMakeInverse1(*divider,minv);
        doneDiv();
    }

    template <class T>
    void DivHelper<T>::doMakeInverseATA(MatrixView<T> minv) const
    {
        TMVAssert(minv.colsize() == TMV_MIN(rowsize(),colsize()));
        TMVAssert(minv.rowsize() == TMV_MIN(rowsize(),colsize()));
        setDiv();
        divider->makeInverseATA(minv);
        doneDiv();
    }

    template <class T>
    bool DivHelper<T>::doIsSingular() const
    {
        setDiv();
        bool s = divider->isSingular();
        doneDiv();
        return s;
    }

    template <class T>
    TMV_RealType(T) DivHelper<T>::doNorm2() const
    {
        TMVAssert(divIsSet());
        TMVAssert(divtype == SV);
        return divider->norm2();
    }

    template <class T>
    TMV_RealType(T) DivHelper<T>::doCondition() const
    {
        TMVAssert(divIsSet());
        TMVAssert(divtype == SV);
        return divider->condition();
    }

    template <class T, class T1> 
    static inline void DoRDivEq1(
        const Divider<T>& div, MatrixView<T1> m0)
    { div.RDivEq(m0); }

    template <class T, class T1> 
    static inline void DoLDivEq1(
        const Divider<T>& div, MatrixView<T1> m0)
    { div.LDivEq(m0); }

    template <class T, class T1, class T0> 
    static inline void DoLDiv1(
        const Divider<T>& div,
        const GenMatrix<T1>& m1, MatrixView<T0> m0)
    { div.LDiv(m1,m0); }

    template <class T, class T1, class T0> 
    static inline void DoRDiv1(
        const Divider<T>& div,
        const GenMatrix<T1>& m1, MatrixView<T0> m0)
    { div.RDiv(m1,m0); }

    template <class T> template <class T1> 
    void DivHelper<T>::doLDivEq(VectorView<T1> v) const
    {
        TMVAssert(colsize() == rowsize());
        TMVAssert(colsize() == v.size());
        setDiv();
        DoLDivEq1(*divider,ColVectorViewOf(v));
        doneDiv();
    }

    template <class T> template <class T1> 
    void DivHelper<T>::doLDivEq(MatrixView<T1> m) const
    {
        TMVAssert(colsize() == rowsize());
        TMVAssert(colsize() == m.colsize());
        setDiv();
        DoLDivEq1(*divider,m);
        doneDiv();
    }

    template <class T> template <class T1> 
    void DivHelper<T>::doRDivEq(VectorView<T1> v) const
    {
        TMVAssert(colsize() == rowsize());
        TMVAssert(colsize() == v.size());
        setDiv();
        DoRDivEq1(*divider,RowVectorViewOf(v));
        doneDiv();
    }

    template <class T> template <class T1> 
    void DivHelper<T>::doRDivEq(MatrixView<T1> m) const
    {
        TMVAssert(colsize() == rowsize());
        TMVAssert(colsize() == m.rowsize());
        setDiv();
        DoRDivEq1(*divider,m);
        doneDiv();
    }

    template <class T> template <class T1, class T0> 
    void DivHelper<T>::doLDiv(
        const GenVector<T1>& v1, VectorView<T0> v0) const
    {
        TMVAssert(rowsize() == v0.size());
        TMVAssert(colsize() == v1.size());
        setDiv();
        DoLDiv1(*divider,ColVectorViewOf(v1),ColVectorViewOf(v0));
        doneDiv();
    }

    template <class T> template <class T1, class T0> 
    void DivHelper<T>::doLDiv(
        const GenMatrix<T1>& m1, MatrixView<T0> m0) const
    {
        TMVAssert(rowsize() == m0.colsize());
        TMVAssert(colsize() == m1.colsize());
        TMVAssert(m1.rowsize() == m0.rowsize());
        setDiv();
        DoLDiv1(*divider,m1,m0);
        doneDiv();
    }

    template <class T> template <class T1, class T0> 
    void DivHelper<T>::doRDiv(
        const GenVector<T1>& v1, VectorView<T0> v0) const
    {
        TMVAssert(rowsize() == v1.size());
        TMVAssert(colsize() == v0.size());
        setDiv();
        DoRDiv1(*divider,RowVectorViewOf(v1),RowVectorViewOf(v0));
        doneDiv();
    }

    template <class T> template <class T1, class T0> 
    void DivHelper<T>::doRDiv(
        const GenMatrix<T1>& m1, MatrixView<T0> m0) const
    {
        TMVAssert(rowsize() == m1.rowsize());
        TMVAssert(colsize() == m0.rowsize());
        TMVAssert(m1.colsize() == m0.colsize());
        setDiv();
        DoRDiv1(*divider,m1,m0);
        doneDiv();
    }

    template <class T>
    void DivHelper<T>::divideUsing(DivType dt) const
    {
        if (!(divtype & dt)) {
            unsetDiv();
            divtype &= ~tmv::DivTypeFlags;
            divtype |= dt;
        }
    }

    template <class T>
    DivType DivHelper<T>::getDivType() const 
    {
        if ((divtype & tmv::DivTypeFlags) == tmv::XX) resetDivType();
        return divtype & tmv::DivTypeFlags;
    }

    template <class T>
    void DivHelper<T>::resetDivType() const
    { divideUsing(getMatrix().isSquare() ? tmv::LU : tmv::QR); }

    template <class T>
    void DivHelper<T>::divideInPlace() const
    { divtype |= tmv::DivInPlaceFlag; saveDiv(); }

    template <class T>
    void DivHelper<T>::dontDivideInPlace() const
    { divtype &= ~tmv::DivInPlaceFlag; }

    template <class T>
    bool DivHelper<T>::divIsInPlace() const
    { return divtype & tmv::DivInPlaceFlag; }

    template <class T>
    void DivHelper<T>::saveDiv() const
    { divtype |= tmv::SaveDivFlag; }

    template <class T>
    void DivHelper<T>::dontSaveDiv() const
    { divtype &= ~tmv::SaveDivFlag; }

    template <class T>
    bool DivHelper<T>::divIsSaved() const
    { return divtype & tmv::SaveDivFlag; }

    template <class T>
    void DivHelper<T>::unsetDiv() const
    { divider.reset(0); }

    template <class T>
    void DivHelper<T>::resetDiv() const
    { unsetDiv(); setDiv(); }

    template <class T>
    bool DivHelper<T>::divIsSet() const 
    { return getDiv(); }

    template <class T>
    void DivHelper<T>::doneDiv() const
    { if (!divIsSaved()) unsetDiv(); }

    template <class T>
    const Divider<T>* DivHelper<T>::getDiv() const 
    { return divider.get(); }

    template <class T>
    bool DivHelper<T>::checkDecomp(std::ostream* fout) const
    {
        TMVAssert(divider.get());
        return divider->checkDecomp(getMatrix(),fout);
    }

    template <class T>
    bool DivHelper<T>::checkDecomp(
        const BaseMatrix<T>& m2, std::ostream* fout) const
    {
        TMVAssert(divider.get());
        return divider->checkDecomp(m2,fout);
    }

#define InstFile "TMV_BaseMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


