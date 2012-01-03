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



#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Divider.h"
#include "TMV_DivImpl.h"

namespace tmv {

    template <class T>
    DivHelper<T>::DivHelper() : pdiv(0) {}

    template <class T>
    DivHelper<T>::~DivHelper() {}

    template <class T>
    void DivHelper<T>::setupDiv() const
    { if (!pdiv.get()) pdiv.reset(new DivImpl(getMatrix())); }

    template <class T>
    T DivHelper<T>::doDet() const
    {
        TMVAssert(colsize() == rowsize());
        setDiv();
        T det = pdiv->div->det();
        doneDiv();
        return det;
    }

    template <class T>
    TMV_RealType(T) DivHelper<T>::doLogDet(T* sign) const
    {
        TMVAssert(colsize() == rowsize());
        setDiv();
        TMV_RealType(T) logdet = pdiv->div->logDet(sign);
        doneDiv();
        return logdet;
    }

    template <class T, class T1> 
    static inline void DoMakeInverse1(
        const Divider<T>& div, const MatrixView<T1>& minv)
    { div.makeInverse(minv); }

    template <class T> template <class T1> 
    void DivHelper<T>::doMakeInverse(const MatrixView<T1>& minv) const
    {
        TMVAssert(minv.colsize() == rowsize());
        TMVAssert(minv.rowsize() == colsize());
        setDiv();
        DoMakeInverse1(*pdiv->div,minv);
        doneDiv();
    }

    template <class T>
    void DivHelper<T>::doMakeInverseATA(const MatrixView<T>& minv) const
    {
        TMVAssert(minv.colsize() == TMV_MIN(rowsize(),colsize()));
        TMVAssert(minv.rowsize() == TMV_MIN(rowsize(),colsize()));
        setDiv();
        pdiv->div->makeInverseATA(minv);
        doneDiv();
    }

    template <class T>
    bool DivHelper<T>::doIsSingular() const
    {
        setDiv();
        bool s = pdiv->div->isSingular();
        doneDiv();
        return s;
    }

    template <class T>
    TMV_RealType(T) DivHelper<T>::doNorm2() const
    {
        setupDiv();
        TMVAssert(divIsSet());
        TMVAssert(pdiv->dt == SV);
        return pdiv->div->norm2();
    }

    template <class T>
    TMV_RealType(T) DivHelper<T>::doCondition() const
    {
        setupDiv();
        TMVAssert(divIsSet());
        TMVAssert(pdiv->dt == SV);
        return pdiv->div->condition();
    }

    template <class T, class T1> 
    static inline void DoRDivEq1(
        const Divider<T>& div, const MatrixView<T1>& m0)
    { div.RDivEq(m0); }

    template <class T, class T1> 
    static inline void DoLDivEq1(
        const Divider<T>& div, const MatrixView<T1>& m0)
    { div.LDivEq(m0); }

    template <class T, class T1, class T0> 
    static inline void DoLDiv1(
        const Divider<T>& div,
        const GenMatrix<T1>& m1, const MatrixView<T0>& m0)
    { div.LDiv(m1,m0); }

    template <class T, class T1, class T0> 
    static inline void DoRDiv1(
        const Divider<T>& div,
        const GenMatrix<T1>& m1, const MatrixView<T0>& m0)
    { div.RDiv(m1,m0); }

    template <class T> template <class T1> 
    void DivHelper<T>::doLDivEq(const VectorView<T1>& v) const
    {
        TMVAssert(colsize() == rowsize());
        TMVAssert(colsize() == v.size());
        setDiv();
        DoLDivEq1(*pdiv->div,ColVectorViewOf(v));
        doneDiv();
    }

    template <class T> template <class T1> 
    void DivHelper<T>::doLDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(colsize() == rowsize());
        TMVAssert(colsize() == m.colsize());
        setDiv();
        DoLDivEq1(*pdiv->div,m);
        doneDiv();
    }

    template <class T> template <class T1> 
    void DivHelper<T>::doRDivEq(const VectorView<T1>& v) const
    {
        TMVAssert(colsize() == rowsize());
        TMVAssert(colsize() == v.size());
        setDiv();
        DoRDivEq1(*pdiv->div,RowVectorViewOf(v));
        doneDiv();
    }

    template <class T> template <class T1> 
    void DivHelper<T>::doRDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(colsize() == rowsize());
        TMVAssert(colsize() == m.rowsize());
        setDiv();
        DoRDivEq1(*pdiv->div,m);
        doneDiv();
    }

    template <class T> template <class T1, class T0> 
    void DivHelper<T>::doLDiv(
        const GenVector<T1>& v1, const VectorView<T0>& v0) const
    {
        TMVAssert(rowsize() == v0.size());
        TMVAssert(colsize() == v1.size());
        setDiv();
        DoLDiv1(*pdiv->div,ColVectorViewOf(v1),ColVectorViewOf(v0));
        doneDiv();
    }

    template <class T> template <class T1, class T0> 
    void DivHelper<T>::doLDiv(
        const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
    {
        TMVAssert(rowsize() == m0.colsize());
        TMVAssert(colsize() == m1.colsize());
        TMVAssert(m1.rowsize() == m0.rowsize());
        setDiv();
        DoLDiv1(*pdiv->div,m1,m0);
        doneDiv();
    }

    template <class T> template <class T1, class T0> 
    void DivHelper<T>::doRDiv(
        const GenVector<T1>& v1, const VectorView<T0>& v0) const
    {
        TMVAssert(rowsize() == v1.size());
        TMVAssert(colsize() == v0.size());
        setDiv();
        DoRDiv1(*pdiv->div,RowVectorViewOf(v1),RowVectorViewOf(v0));
        doneDiv();
    }

    template <class T> template <class T1, class T0> 
    void DivHelper<T>::doRDiv(
        const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
    {
        TMVAssert(rowsize() == m1.rowsize());
        TMVAssert(colsize() == m0.rowsize());
        TMVAssert(m1.colsize() == m0.colsize());
        setDiv();
        DoRDiv1(*pdiv->div,m1,m0);
        doneDiv();
    }

    template <class T>
    void DivHelper<T>::divideInPlace() const
    {
        setupDiv();
        saveDiv();
        pdiv->inplace = true;
    }

    template <class T>
    bool DivHelper<T>::isDivInPlace() const
    {
        setupDiv();
        return pdiv->inplace;
    }

    template <class T>
    void DivHelper<T>::saveDiv() const
    { 
        setupDiv();
        pdiv->cache = true; 
    }

    template <class T>
    void DivHelper<T>::divideUsing(DivType dt) const
    { 
        setupDiv();
        if (dt != pdiv->dt) unsetDiv();
        pdiv->dt = dt;
    }

    template <class T>
    void DivHelper<T>::setDiv() const
    { 
        setupDiv();
        if (!pdiv->div.get()) {
            if (pdiv->dt == XX) 
                pdiv->dt = (colsize() == rowsize()) ? LU : QR;
            newDivider();
        }
    }

    template <class T>
    void DivHelper<T>::unsetDiv() const
    { 
        setupDiv();
        pdiv->div.reset(0);
    }

    template <class T>
    void DivHelper<T>::resetDiv() const
    {
        unsetDiv(); 
        setDiv(); 
    }

    template <class T>
    bool DivHelper<T>::divIsSet() const 
    {
        setupDiv();
        return pdiv->div.get(); 
    }

    template <class T>
    void DivHelper<T>::doneDiv() const
    { if (!pdiv->cache) unsetDiv(); }

    template <class T>
    const Divider<T>* DivHelper<T>::getDiv() const 
    {
        setupDiv();
        return pdiv->div.get(); 
    }

    template <class T>
    void DivHelper<T>::setDiv(Divider<T>* d) const 
    {
        setupDiv();
        pdiv->div.reset(d);
    }

    template <class T>
    DivType DivHelper<T>::getDivType() const 
    {
        setupDiv();
        return pdiv->dt; 
    }

    template <class T>
    bool DivHelper<T>::checkDecomp(std::ostream* fout) const
    {
        TMVAssert(pdiv->div.get());
        return pdiv->div->checkDecomp(getMatrix(),fout);
    }

    template <class T>
    bool DivHelper<T>::checkDecomp(
        const BaseMatrix<T>& m2, std::ostream* fout) const
    {
        TMVAssert(pdiv->div.get());
        return pdiv->div->checkDecomp(m2,fout);
    }

#define InstFile "TMV_BaseMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


