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


//#define XDEBUG


#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VIt.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include <ostream>
#include "portable_platform.h"

#ifdef NOTHROW
#include <iostream>
#else
#include "TMV_VectorRE.h"
#endif

#ifdef XDEBUG
#include "tmv/TMV_VectorArith.h"
#include "tmv/TMV_MatrixArith.h"
#include <iostream>
using std::cerr;
using std::endl;
#endif

namespace tmv {

#define RT TMV_RealType(T)

    template <class T>
    T GenDiagMatrix<T>::det() const
    {
        T signdet(1);
        RT logdet = logDet(&signdet);
        if (signdet == T(0)) return T(0);
        else return signdet * TMV_EXP(logdet);
    }

    template <class T>
    RT GenDiagMatrix<T>::logDet(T* sign) const
    {
        const T* di = diag().cptr();
        const int ds = diag().step();
        T s(1);
        RT logdet(0);
        if (ds == 1) {
            for(int i=size();i>0;--i,++di) {
                if (*di == T(0)) { 
                    logdet = TMV_LOG(TMV_REAL(*di));
                    if (sign) s = T(0); 
                } else {
                    RT a = TMV_ABS(*di);
                    logdet += TMV_LOG(a);
                    if (sign) {
                        if (isReal(T())) {
                            if (TMV_REAL(*di) < RT(0)) s = -s;
                        } else {
                            s *= (*di/a);
                        }
                    }
                }
            }
        } else {
            for(int i=size();i>0;--i,di+=ds) {
                if (*di == T(0)) { 
                    logdet = TMV_LOG(TMV_REAL(*di));
                    if (sign) s = T(0); 
                } else {
                    RT a = TMV_ABS(*di);
                    logdet += TMV_LOG(a);
                    if (sign) {
                        if (isReal(T())) {
                            if (TMV_REAL(*di) < RT(0)) s = -s;
                        } else {
                            s *= (*di/a);
                        }
                    }
                }
            }
        }
        if (sign) {
            if (diag().isconj()) *sign = TMV_CONJ(s);
            else *sign = s;
        }
        return logdet;
    }

#ifdef INST_INT
    template <>
    int GenDiagMatrix<int>::det() const
    {
        const int* di = diag().cptr();
        const int ds = diag().step();
        int det(1);
        if (ds == 1) {
            for(int i=size();i>0;--i,++di) det *= *di;
        } else {
            for(int i=size();i>0;--i,di+=ds) det *= *di;
        }
        return det;
    }
    template <>
    std::complex<int> GenDiagMatrix<std::complex<int> >::det() const
    {
        const std::complex<int>* di = diag().cptr();
        const int ds = diag().step();
        std::complex<int> det(1);
        if (ds == 1) {
            for(int i=size();i>0;--i,++di) det *= *di;
        } else {
            for(int i=size();i>0;--i,di+=ds) det *= *di;
        }
        return diag().isconj() ? TMV_CONJ(det) : det;
    }

    template <>
    int GenDiagMatrix<int>::logDet(int* ) const
    { TMVAssert(TMV_FALSE); return 0; }
    template <>
    int GenDiagMatrix<std::complex<int> >::logDet(std::complex<int>* ) const
    { TMVAssert(TMV_FALSE); return 0; }
#endif

    template <class T>
    auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::newCopy() const
    {
        auto_ptr<BaseMatrix<T> > a(new DiagMatrix<T>(*this));
        return a;
    }

    template <class T>
    auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::newView() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(view()));
        return a;
    }

    template <class T>
    auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::newTranspose() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(transpose()));
        return a;
    }

    template <class T>
    auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::newConjugate() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(conjugate()));
        return a;
    }

    template <class T>
    auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::newAdjoint() const
    {
        auto_ptr<BaseMatrix<T> > a(new ConstDiagMatrixView<T>(adjoint()));
        return a;
    }

    template <class T>
    auto_ptr<BaseMatrix<T> > GenDiagMatrix<T>::newInverse() const
    {
        auto_ptr<DiagMatrix<T> > minv(new DiagMatrix<T>(*this));
        minv->invertSelf();
        BaseMatrix<T>* ret1 = minv.release();
        auto_ptr<BaseMatrix<T> > ret(ret1);
        return ret;
    }


#ifndef NOTHROW
    template <class T>
    class SingularDiagMatrix : public Singular
    {
    public:
        DiagMatrix<T> A;

        SingularDiagMatrix(const GenDiagMatrix<T>& _A) :
            Singular("DiagMatrix."), A(_A) {}
        ~SingularDiagMatrix() throw() {}
        void write(std::ostream& os) const throw()
        {
            Singular::write(os);
            os<<A<<std::endl;
        }
    };
#endif

    template <class T, IndexStyle I> 
    const DiagMatrixView<T,I>& DiagMatrixView<T,I>::invertSelf() const
    {
        T* di = diag().ptr();
        const int dstep = diag().step();

        if (dstep == 1) {
            for(int i=size();i>0;--i,++di) {
                if (*di == T(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular DiagMatrix found\n";
                    exit(1); 
#else
                    throw SingularDiagMatrix<T>(*this);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(di >= itsdiag._first);
                TMVAssert(di < itsdiag._last);
#endif
                if (TMV_IMAG(*di) == RT(0))
                    *di = RT(1) / TMV_REAL(*di);
                else
                    *di = RT(1) / *di;
            }
        } else {
            for(int i=size();i>0;--i,di+=dstep) {
                if (*di == T(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular DiagMatrix found\n";
                    exit(1); 
#else
                    throw SingularDiagMatrix<T>(*this);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(di >= itsdiag._first);
                TMVAssert(di < itsdiag._last);
#endif
                if (TMV_IMAG(*di) == RT(0))
                    *di = RT(1) / TMV_REAL(*di);
                else
                    *di = RT(1) / *di;
            }
        }
        return *this;
    }

#ifdef INST_INT
    template <>
    const DiagMatrixView<int,CStyle>& 
        DiagMatrixView<int,CStyle>::invertSelf() const
    { TMVAssert(TMV_FALSE); return *this; }
    template <>
    const DiagMatrixView<std::complex<int>,CStyle>& 
        DiagMatrixView<std::complex<int>,CStyle>::invertSelf() const
    { TMVAssert(TMV_FALSE); return *this; }
#endif

    template <class T> template <class T1> 
    void GenDiagMatrix<T>::doMakeInverse(const MatrixView<T1>& minv) const
    {
        bool ss = SameStorage(diag(),minv);
        if (!ss) minv.setZero();
        (DiagMatrixViewOf(minv.diag()) = *this).invertSelf();
        if (ss && size() > 1) {
            minv.upperTri().offDiag().setZero();
            minv.lowerTri().offDiag().setZero();
        }
    }

#ifdef INST_INT
    template <> template <class T1>
    void GenDiagMatrix<int>::doMakeInverse(const MatrixView<T1>& ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1>
    void GenDiagMatrix<std::complex<int> >::doMakeInverse(
        const MatrixView<T1>& ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T> template <class T1> 
    void GenDiagMatrix<T>::doMakeInverse(const DiagMatrixView<T1>& minv) const
    { (minv = *this).invertSelf(); }

#ifdef INST_INT
    template <> template <class T1>
    void GenDiagMatrix<int>::doMakeInverse(const DiagMatrixView<T1>& ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1>
    void GenDiagMatrix<std::complex<int> >::doMakeInverse(
        const DiagMatrixView<T1>& ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T>
    void GenDiagMatrix<T>::doMakeInverseATA(const DiagMatrixView<T>& ata) const
    {
        makeInverse(ata);
        T* mi = ata.diag().ptr();
        const int ds = ata.diag().step();
        if (ds==1) {
            for(int i=size();i>0;--i,++mi) {
#ifdef TMVFLDEBUG
                TMVAssert(mi >= ata.diag()._first);
                TMVAssert(mi < ata.diag()._last);
#endif
                *mi = TMV_NORM(*mi);
            }
        } else {
            for(int i=size();i>0;--i,mi+=ds) {
#ifdef TMVFLDEBUG
                TMVAssert(mi >= ata.diag()._first);
                TMVAssert(mi < ata.diag()._last);
#endif
                *mi = TMV_NORM(*mi);
            }
        }
    }

#ifdef INST_INT
    template <> 
    void GenDiagMatrix<int>::doMakeInverseATA(
        const DiagMatrixView<int>& ) const
    { TMVAssert(TMV_FALSE); }
    template <> 
    void GenDiagMatrix<std::complex<int> >::doMakeInverseATA(
        const DiagMatrixView<std::complex<int> >& ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T>
    void GenDiagMatrix<T>::doMakeInverseATA(const MatrixView<T>& ata) const
    {
        ata.setZero();
        makeInverseATA(DiagMatrixViewOf(ata.diag()));
    }

#ifdef INST_INT
    template <> 
    void GenDiagMatrix<int>::doMakeInverseATA(const MatrixView<int>& ) const
    { TMVAssert(TMV_FALSE); }
    template <> 
    void GenDiagMatrix<std::complex<int> >::doMakeInverseATA(
        const MatrixView<std::complex<int> >& ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T>
    QuotXD<T,T> GenDiagMatrix<T>::QInverse() const
    { return QuotXD<T,T>(T(1),*this); }

#define CT std::complex<T>

    template <bool cd, class T, class Td> 
    static void DoDiagLDivEq1(
        const GenDiagMatrix<Td>& d, const VectorView<T>& v)
    {
        TMVAssert(v.size() == d.size());
        TMVAssert(v.ct()==NonConj);
        TMVAssert(v.size() > 0);

        const Td* di = d.diag().cptr();
        T* vi = v.ptr();
        const int dstep = d.diag().step();
        const int vstep = v.step();

        if (dstep == 1 && vstep == 1) {
            for(int i=v.size();i>0;--i,++di,++vi) {
#ifdef TMVFLDEBUG
                TMVAssert(vi >= v._first);
                TMVAssert(vi < v._last);
#endif
                if (*di == Td(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular DiagMatrix found\n";
                    exit(1); 
#else
                    throw SingularDiagMatrix<Td>(d);
#endif
                }
                if (TMV_IMAG(*di) == RT(0)) {
                    if (TMV_REAL(*di) != RT(1))
                        *vi /= TMV_REAL(*di);
                } else {
                    *vi /= (cd?TMV_CONJ(*di):*di);
                }
            }
        } else {
            for(int i=v.size();i>0;--i,di+=dstep,vi+=vstep) {
#ifdef TMVFLDEBUG
                TMVAssert(vi >= v._first);
                TMVAssert(vi < v._last);
#endif
                if (*di == Td(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular DiagMatrix found\n";
                    exit(1); 
#else
                    throw SingularDiagMatrix<Td>(d);
#endif
                }
                if (TMV_IMAG(*di) == RT(0)) {
                    if (TMV_REAL(*di) != RT(1)) *vi /= TMV_REAL(*di);
                } else {
                    *vi /= (cd?TMV_CONJ(*di):*di);
                }
            }
        }
    }

    template <class T, class Td> 
    static inline void DoDiagLDivEq(
        const GenDiagMatrix<Td>& d, const VectorView<T>& v)
    { 
        if (d.diag().isconj()) DoDiagLDivEq1<true>(d,v);
        else DoDiagLDivEq1<false>(d,v);
    }

    template <class T> template <class T1> 
    void GenDiagMatrix<T>::doLDivEq(const VectorView<T1>& v) const
    {
#ifdef XDEBUG
        DiagMatrix<T> d0(*this);
        Vector<T1> v0(v);
#endif

        TMVAssert(v.size() == size());

        if (v.size() > 0) {
            if (v.isconj()) DoDiagLDivEq(conjugate(),v.conjugate());
            else DoDiagLDivEq(*this,v);
        }

#ifdef XDEBUG
        Vector<T1> v1 = d0*v;
        if (Norm(v1-v0) > 0.001*Norm(v0)) {
            cerr<<"DiagLDivEq v: \n";
            cerr<<"d = "<<TMV_Text(*this)<<"  "<<d0<<endl;
            cerr<<"v = "<<TMV_Text(v)<<"  "<<v0<<endl;
            cerr<<"-> v/d = "<<v<<endl;
            cerr<<"d*(v/d) = "<<v1<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
    template <> template <class T1>
    void GenDiagMatrix<int>::doLDivEq(const VectorView<T1>& ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1>
    void GenDiagMatrix<std::complex<int> >::doLDivEq(
        const VectorView<T1>& ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T> template <class T1, class T0> 
    void GenDiagMatrix<T>::doLDiv(
        const GenVector<T1>& v1, const VectorView<T0>& v0) const
    {
        TMVAssert(v1.size() == size());
        TMVAssert(v0.size() == size());
        if (SameStorage(diag(),v0)) {
            DiagMatrix<T> temp = *this;
            temp.doLDivEq(v0=v1);
        } else {
            doLDivEq(v0=v1);
        }
    }

#ifdef INST_INT
    template <> template <class T1, class T0>
    void GenDiagMatrix<int>::doLDiv(
        const GenVector<T1>& , const VectorView<T0>& ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1, class T0>
    void GenDiagMatrix<std::complex<int> >::doLDiv(
        const GenVector<T1>& , const VectorView<T0>& ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <bool rm, bool cd, class T, class Td> 
    static void RowDiagLDivEq(
        const GenDiagMatrix<Td>& d, const MatrixView<T>& m)
    {
        TMVAssert(d.size() == m.colsize());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(m.ct() == NonConj);
        TMVAssert(rm == m.isrm());
        TMVAssert(cd == d.diag().isconj());

        const Td* di = d.diag().cptr();
        T* mrowi = m.ptr();
        const int dstep = d.diag().step();
        const int stepj = m.stepj();
        const int stepi = m.stepi();
        const int M = m.colsize();
        const int N = m.rowsize();

        for(int i=M;i>0;--i,di+=dstep,mrowi+=stepi) {
            T* mij = mrowi;
            if (*di == Td(0)) {
#ifdef NOTHROW
                std::cerr<<"Singular DiagMatrix found\n";
                exit(1); 
#else
                throw SingularDiagMatrix<Td>(d);
#endif
            } else if (TMV_IMAG(*di) == RT(0)) {
                RT invdi = RT(1)/TMV_REAL(*di);
                for(int j=N;j>0;--j,(rm?++mij:mij+=stepj)) {
#ifdef TMVFLDEBUG
                    TMVAssert(mij >= m._first);
                    TMVAssert(mij < m._last);
#endif
                    *mij *= invdi;
                }
            } else {
                Td invdi = RT(1)/(cd?TMV_CONJ(*di):*di);
                for(int j=N;j>0;--j,(rm?++mij:mij+=stepj)) {
#ifdef TMVFLDEBUG
                    TMVAssert(mij >= m._first);
                    TMVAssert(mij < m._last);
#endif
                    *mij *= invdi;
                }
            }
        }
    }

    template <bool cm, bool cd, class T, class Td> 
    static void ColDiagLDivEq(
        const GenDiagMatrix<Td>& d, const MatrixView<T>& m)
    {
        TMVAssert(d.size() == m.colsize());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(m.ct()==NonConj);
        TMVAssert(cm == m.iscm());
        TMVAssert(cd == d.diag().isconj());

        DiagMatrix<Td> invd(d.size());
        const Td* di = d.diag().cptr();
        const int step = d.diag().step();
        Td* invdi = invd.diag().ptr();

        if (step == 1) {
            for(int i=d.size();i>0;--i,++di,++invdi) {
                if (*di == Td(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular DiagMatrix found\n";
                    exit(1); 
#else
                    throw SingularDiagMatrix<Td>(d);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(invdi >= invd.diag()._first);
                TMVAssert(invdi < invd.diag()._last);
#endif
                *invdi = RT(1)/(cd?TMV_CONJ(*di):*di);
            }
        } else {
            for(int i=d.size();i>0;--i,di+=step,++invdi) {
                if (*di == Td(0)) {
#ifdef NOTHROW
                    std::cerr<<"Singular DiagMatrix found\n";
                    exit(1); 
#else
                    throw SingularDiagMatrix<Td>(d);
#endif
                }
#ifdef TMVFLDEBUG
                TMVAssert(invdi >= invd.diag()._first);
                TMVAssert(invdi < invd.diag()._last);
#endif
                *invdi = RT(1)/(cd?TMV_CONJ(*di):*di);
            }
        }
        m = invd*m;
    }

    template <class T, class Td> 
    static void DoDiagLDivEq(
        const GenDiagMatrix<Td>& d, const MatrixView<T>& m)
    {
        if (d.diag().isconj())
            if (m.isrm()) RowDiagLDivEq<true,true>(d,m);
            else if (m.iscm()) ColDiagLDivEq<true,true>(d,m);
            else if (m.colsize() > m.rowsize()) ColDiagLDivEq<false,true>(d,m);
            else RowDiagLDivEq<false,true>(d,m);
        else
            if (m.isrm()) RowDiagLDivEq<true,false>(d,m);
            else if (m.iscm()) ColDiagLDivEq<true,false>(d,m);
            else if (m.colsize() > m.rowsize()) ColDiagLDivEq<false,false>(d,m);
            else RowDiagLDivEq<false,false>(d,m);
    }

    template <class T> template <class T1> 
    void GenDiagMatrix<T>::doLDivEq(const MatrixView<T1>& m) const
    {
        TMVAssert(size() == m.colsize());

#ifdef XDEBUG
        DiagMatrix<T> d0(*this);
        Matrix<T1> m0(m);
#endif

        if (m.colsize() > 0 && m.rowsize() > 0) {
            if (m.isconj()) conjugate().doLDivEq(m.conjugate());
            else if (m.rowsize() == 1) doLDivEq(m.col(0));
            else DoDiagLDivEq(*this,m);
        }

#ifdef XDEBUG
        Matrix<T1> m1 = d0*m;
        if (Norm(m1-m0) > 0.001*Norm(m0)) {
            cerr<<"DiagLDivEq m: \n";
            cerr<<"d = "<<TMV_Text(*this)<<"  "<<d0<<endl;
            cerr<<"m = "<<TMV_Text(m)<<"  "<<m0<<endl;
            cerr<<"-> m/d = "<<m<<endl;
            cerr<<"d*(m/d) = "<<m1<<endl;
            abort();
        }
#endif
    }

#ifdef INST_INT
    template <> template <class T1>
    void GenDiagMatrix<int>::doLDivEq(const MatrixView<T1>& ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1>
    void GenDiagMatrix<std::complex<int> >::doLDivEq(
        const MatrixView<T1>& ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T> template <class T1, class T0> 
    void GenDiagMatrix<T>::doLDiv(
        const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
    {
        TMVAssert(m1.rowsize() == m0.rowsize());
        TMVAssert(m1.colsize() == size());
        TMVAssert(m0.colsize() == size());
        if (SameStorage(diag(),m0)) {
            DiagMatrix<T> temp = *this;
            temp.doLDivEq(m0=m1);
        } else {
            doLDivEq(m0=m1);
        }
    }

#ifdef INST_INT
    template <> template <class T1, class T0>
    void GenDiagMatrix<int>::doLDiv(
        const GenMatrix<T1>& , const MatrixView<T0>& ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1, class T0>
    void GenDiagMatrix<std::complex<int> >::doLDiv(
        const GenMatrix<T1>& , const MatrixView<T0>& ) const
    { TMVAssert(TMV_FALSE); }
#endif

#undef CT

    // This bit is to workaround a bug in pgCC that was fixed in version 7.
    // I don't know if versions earlier than 6.1 had the bug, but 
    // I apply the workaround to all version before 7.
    template <class T>
    inline T Value(const T& x) { return x; }
#ifdef PLATFORM_COMPILER_PGI
#if PLATFORM_COMPILER_VERSION < 0x070000
    inline double Value(const long double& x) { return double(x); }
    inline std::complex<double> Value(const std::complex<long double>& x) 
    { return std::complex<double>(x); }
#endif
#endif

    template <class T>
    void GenDiagMatrix<T>::write(std::ostream& os) const
    {
        const int sd = diag().step();
        const T* di = diag().cptr();
        os << size() <<' '<< size() << std::endl;
        for (int i=0,nmi=size();nmi>0;++i,--nmi,di+=sd) {
            os << "( ";
            for(int k=i;k>0;--k) os <<' '<<Value(T(0))<<' ';
            if (diag().isconj())
                os <<' '<<Value(TMV_CONJ(*di))<<' ';
            else
                os <<' '<<Value(*di)<<' ';
            for(int k=nmi-1;k>0;--k) os <<' '<<Value(T(0))<<' ';
            os << " )\n";
        }
    }

    template <class T>
    void GenDiagMatrix<T>::write(std::ostream& os, RT thresh) const
    {
        const int sd = diag().step();
        const T* di = diag().cptr();
        os << size() <<' '<< size() << std::endl;
        for (int i=0,nmi=size();nmi>0;++i,--nmi,di+=sd) {
            os << "( ";
            for(int k=i;k>0;--k) os <<' '<<Value(T(0))<<' ';
            if (diag().isconj())
                os <<' '<<
                    Value(TMV_ABS(*di)<thresh ? T(0) : TMV_CONJ(*di))<<' ';
            else
                os <<' '<<Value(TMV_ABS(*di)<thresh ? T(0) : *di)<<' ';
            for(int k=nmi-1;k>0;--k) os <<' '<<Value(T(0))<<' ';
            os << " )\n";
        }
    }

#ifndef NOTHROW
    template <class T>
    class DiagMatrixReadError : public ReadError
    {
    public :
        int i;
        mutable auto_ptr<DiagMatrix<T> > m;
        char exp,got;
        size_t s;
        bool is, iseof, isbad;

        DiagMatrixReadError(std::istream& _is) throw() :
            ReadError("DiagMatrix."),
            i(0), m(0), exp(0), got(0), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(
            int _i, const GenDiagMatrix<T>& _m,
            char _e, char _g, size_t _s,
            bool _is, bool _iseof, bool _isbad
        ) throw() :
            ReadError("DiagMatrix."),
            i(_i), m(new DiagMatrix<T>(_m)), exp(_e), got(_g), s(_s),
            is(_is), iseof(_iseof), isbad(_isbad) {}
        DiagMatrixReadError(
            const GenDiagMatrix<T>& _m, std::istream& _is, size_t _s
        ) throw() :
            ReadError("DiagMatrix."),
            i(0), m(new DiagMatrix<T>(_m)), exp(0), got(0), s(_s),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(std::istream& _is, char _e, char _g) throw() :
            ReadError("DiagMatrix."),
            i(0), m(0), exp(_e), got(_g), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        DiagMatrixReadError(const DiagMatrixReadError<T>& rhs) :
            i(rhs.i), m(rhs.m), exp(rhs.exp), got(rhs.got), s(rhs.s),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        virtual ~DiagMatrixReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for DiagMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (m.get() && s != m->size()) {
                os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
            }
            if (!is) {
                if (iseof) {
                    os<<"Input stream reached end-of-file prematurely.\n";
                } else if (isbad) {
                    os<<"Input stream is corrupted.\n";
                } else {
                    os<<"Input stream cannot read next character.\n";
                }
            }
            if (m.get()) {
                os<<"The portion of the DiagMatrix which was successfully "
                    "read is: \n";
                ConstDiagMatrixView<T> mm = m->view();
                os<<"( ";
                for(int ii=0;ii<i;++ii)
                    os<<' '<<mm(ii,ii)<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    template <class T, IndexStyle I> 
    std::istream& operator>>(std::istream& is, auto_ptr<DiagMatrix<T,I> >& m)
    {
        char d;
        is >> d;
        if (!is || d != 'D') {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix ReadError: "<<d<<" != D\n";
            exit(1); 
#else
            throw DiagMatrixReadError<T>(is,'D',d);
#endif
        }
        size_t size;
        is >> size;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix ReadError: !is\n";
            exit(1); 
#else
            throw DiagMatrixReadError<T>(is);
#endif
        }
        m.reset(new DiagMatrix<T,I>(size));
#ifndef NOTHROW
        try {
#endif
            m->diag().read(is);
#ifndef NOTHROW
        } catch (VectorReadError<T>& ve) {
            throw DiagMatrixReadError<T>(
                ve.i,*m,ve.exp,ve.got,ve.s,ve.is,ve.iseof,ve.isbad);
        }
#endif
        return is;
    }

    template <class T>
    std::istream& operator>>(std::istream& is, const DiagMatrixView<T>& m)
    {
        char d;
        is >> d;
        if (!is || d != 'D') {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix ReadError: "<<d<<" != D\n";
            exit(1); 
#else
            throw DiagMatrixReadError<T>(is,'D',d);
#endif
        }
        size_t s;
        is >> s;
        if (!is) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix ReadError: !is\n";
            exit(1); 
#else
            throw DiagMatrixReadError<T>(is);
#endif
        }
        if (m.size() != s) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix ReadError: Wrong size\n";
            exit(1); 
#else
            throw DiagMatrixReadError<T>(m,is,s);
#endif
        }
        TMVAssert(m.size() == s);
#ifndef NOTHROW
        try {
#endif
            m.diag().read(is);
#ifndef NOTHROW
        } catch (VectorReadError<T>& ve) {
            throw DiagMatrixReadError<T>(
                ve.i,m,ve.exp,ve.got,ve.s,ve.is,ve.iseof,ve.isbad);
        }
#endif
        return is;
    }

#define InstFile "TMV_DiagMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


