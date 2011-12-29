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



#include "tmv/TMV_BandSVD.h"
#include "TMV_BandSVDiv.h"
#include "tmv/TMV_BandMatrix.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include <ostream>

namespace tmv {

#define RT TMV_RealType(T)

    template <class T> 
    struct BandSVDiv<T>::BandSVDiv_Impl
    {
    public:
        BandSVDiv_Impl(const GenBandMatrix<T>& A);

        const bool istrans;
        Matrix<T,ColMajor> U;
        DiagMatrix<RT> S;
        Matrix<T,ColMajor> V;
        RT logdet;
        T signdet;
        mutable int kmax;
    }; 

#define M TMV_MAX(A.colsize(),A.rowsize())
#define N TMV_MIN(A.colsize(),A.rowsize())

    template <class T> 
    BandSVDiv<T>::BandSVDiv_Impl::BandSVDiv_Impl(
        const GenBandMatrix<T>& A) :
        istrans(A.colsize() < A.rowsize()),
        U(M,N), S(N), V(N,N), logdet(0), signdet(1), kmax(0) {}

#undef M
#undef N

    template <class T> 
    BandSVDiv<T>::BandSVDiv(const GenBandMatrix<T>& A) :
        pimpl(new BandSVDiv_Impl(A))
    {
        SV_Decompose<T>(
            pimpl->istrans ? A.transpose() : A.view(),
            pimpl->U.view(),pimpl->S.view(),pimpl->V.view(),
            pimpl->logdet,pimpl->signdet);
        thresh(TMV_Epsilon<T>());
    }

    template <class T> 
    BandSVDiv<T>::~BandSVDiv() {}

    template <class T> template <class T1> 
    void BandSVDiv<T>::doLDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(m.colsize() == rowsize());
        TMVAssert(m.colsize() == colsize());
        if (pimpl->istrans) 
            CallSV_RDiv(
                T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
                m.transpose(),m.transpose());
        else CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
    }

    template <class T> template <class T1> 
    void BandSVDiv<T>::doRDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(m.rowsize() == rowsize());
        TMVAssert(m.rowsize() == colsize());
        if (pimpl->istrans) 
            CallSV_LDiv(
                T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
                m.transpose(),m.transpose());
        else CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,m);
    }

    template <class T> template <class T1, class T2> 
    void BandSVDiv<T>::doLDiv(
        const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        TMVAssert(m.rowsize() == x.rowsize());
        TMVAssert(m.colsize() == colsize());
        TMVAssert(x.colsize() == rowsize());
        if (pimpl->istrans) 
            CallSV_RDiv(
                T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
                m.transpose(),x.transpose());
        else CallSV_LDiv(T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
    }

    template <class T> template <class T1, class T2> 
    void BandSVDiv<T>::doRDiv(
        const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    {
        TMVAssert(m.colsize() == x.colsize());
        TMVAssert(m.rowsize() == rowsize());
        TMVAssert(x.rowsize() == colsize());
        if (pimpl->istrans) 
            CallSV_LDiv(
                T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,
                m.transpose(),x.transpose());
        else CallSV_RDiv(T(),pimpl->U,pimpl->S,pimpl->V,pimpl->kmax,m,x);
    }

    template <class T> 
    T BandSVDiv<T>::det() const 
    { 
        if (pimpl->signdet == T(0)) return T(0);
        else return pimpl->signdet * TMV_EXP(pimpl->logdet);  
    }

    template <class T> 
    RT BandSVDiv<T>::logDet(T* sign) const 
    { 
        if (sign) *sign = pimpl->signdet;
        return pimpl->logdet;
    }

    template <class T> template <class T1> 
    void BandSVDiv<T>::doMakeInverse(const MatrixView<T1>& minv) const
    { 
        if (pimpl->istrans) {
            Matrix<T,ColMajor> SinvV = 
                pimpl->V.conjugate().rowRange(0,pimpl->kmax) /
                pimpl->S.subDiagMatrix(0,pimpl->kmax);
            minv = pimpl->U.conjugate().colRange(0,pimpl->kmax) * SinvV;
        } else  {
            Matrix<T,ColMajor> SinvUt = 
                pimpl->U.adjoint().rowRange(0,pimpl->kmax) /
                pimpl->S.subDiagMatrix(0,pimpl->kmax);
            minv = pimpl->V.adjoint().colRange(0,pimpl->kmax) * SinvUt;
        }
    }

    template <class T> 
    void BandSVDiv<T>::doMakeInverseATA(const MatrixView<T>& minv) const
    {
        Matrix<T,ColMajor> SinvV = 
            pimpl->V.rowRange(0,pimpl->kmax) /
            pimpl->S.subDiagMatrix(0,pimpl->kmax);
        if (pimpl->istrans)
            minv = SinvV.transpose() * SinvV.conjugate();
        else
            minv = SinvV.adjoint() * SinvV;
    }

    template <class T> 
    bool BandSVDiv<T>::isSingular() const 
    { return pimpl->kmax < pimpl->S.size(); }

    template <class T> 
    RT BandSVDiv<T>::norm2() const 
    { return pimpl->S.size() > 0 ? pimpl->S(0) : RT(0); }

    template <class T> 
    RT BandSVDiv<T>::condition() const 
    {
        return pimpl->S.size() > 0 ? 
            pimpl->S(0)/pimpl->S(pimpl->S.size()-1) :
            RT(1); 
    }

    template <class T> 
    void BandSVDiv<T>::thresh(RT toler, std::ostream* debugout) const
    {
        if (pimpl->S.size() == 0) pimpl->kmax = 0;
        else {
            TMVAssert(toler < RT(1));
            RT thresh = pimpl->S(0)*toler;
            for(pimpl->kmax=pimpl->S.size(); 
                pimpl->kmax>0 && pimpl->S(pimpl->kmax-1)<=thresh; 
                --pimpl->kmax);
            if(debugout) {
                (*debugout)<<"S = "<<pimpl->S<<std::endl;
                (*debugout)<<"Smax = "<<pimpl->S(0)<<
                    ", thresh = "<<thresh<<std::endl;
                (*debugout)<<"kmax = "<<pimpl->kmax;
                (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
            }
        }
    }

    template <class T> 
    void BandSVDiv<T>::top(int neigen, std::ostream* debugout) const
    {
        TMVAssert(neigen > 0);
        if (neigen < pimpl->S.size()) pimpl->kmax = neigen;
        else pimpl->kmax = pimpl->S.size();
        if(debugout) {
            (*debugout)<<"S = "<<pimpl->S<<std::endl;
            (*debugout)<<"kmax = "<<pimpl->kmax;
            (*debugout)<<" (S.size = "<<pimpl->S.size()<<")"<<std::endl;
        }
    }

    template <class T> int BandSVDiv<T>::getKMax() const 
    { return pimpl->kmax; }

    template <class T> ConstMatrixView<T> BandSVDiv<T>::getU() const
    {
        if (pimpl->istrans) return pimpl->V.transpose(); 
        else return pimpl->U.view(); 
    }

    template <class T> ConstDiagMatrixView<RT> BandSVDiv<T>::getS() const
    { return pimpl->S.view(); }

    template <class T> ConstMatrixView<T> BandSVDiv<T>::getV() const
    {
        if (pimpl->istrans) return pimpl->U.transpose(); 
        else return pimpl->V.view();
    }

    template <class T> 
    bool BandSVDiv<T>::checkDecomp(
        const BaseMatrix<T>& m, std::ostream* fout) const
    {
        Matrix<T> mm = m;
        if (fout) {
            *fout << "BandSVDiv:\n";
            *fout << "M = "<<m<<std::endl;
            *fout << "U = "<<getU()<<std::endl;
            *fout << "S = "<<getS()<<std::endl;
            *fout << "V = "<<getV()<<std::endl;
        }
        Matrix<T> usv = getU()*getS()*getV();
        RT nm = Norm(usv-mm);
        nm /= Norm(getU())*Norm(getS())*Norm(getV());
        RT cond = getS()(0) / getS()(getKMax()-1);
        if (fout) {
            *fout << "USV = "<<usv<<std::endl;
            *fout << "Norm(M-USV) = "<<Norm(mm-usv)<<std::endl;
            *fout << "Norm(M-USV)/Norm(USV) = "<<nm<<std::endl;
        }
        return nm < cond*RT(mm.colsize())*TMV_Epsilon<T>();
    }

    template <class T> 
    int BandSVDiv<T>::colsize() const
    { return pimpl->istrans ? pimpl->U.rowsize() : pimpl->U.colsize(); }

    template <class T> 
    int BandSVDiv<T>::rowsize() const
    { return pimpl->istrans ? pimpl->U.colsize() : pimpl->U.rowsize(); }

#undef RT

#ifdef INST_INT
#undef INST_INT
#endif

#define InstFile "TMV_BandSVD.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


