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


//#define XDEBUG


#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_VIt.h"
#include "tmv/TMV_DiagMatrixArith.h"
#include <ostream>

#ifdef NOTHROW
#include <iostream>
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
        const ptrdiff_t ds = diag().step();
        T s(1);
        RT logdet(0);
        if (ds == 1) {
            for(ptrdiff_t i=size();i>0;--i,++di) {
                if (*di == T(0)) { 
                    logdet = TMV_LOG(TMV_REAL(*di));
                    s = T(0); 
                } else {
                    RT a = TMV_ABS(*di);
                    logdet += TMV_LOG(a);
                    s *= TMV_SIGN(*di,a);
                }
            }
        } else {
            for(ptrdiff_t i=size();i>0;--i,di+=ds) {
                if (*di == T(0)) { 
                    logdet = TMV_LOG(TMV_REAL(*di));
                    s = T(0); 
                } else {
                    RT a = TMV_ABS(*di);
                    logdet += TMV_LOG(a);
                    s *= TMV_SIGN(*di,a);
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
        const ptrdiff_t ds = diag().step();
        int det(1);
        if (ds == 1) {
            for(ptrdiff_t i=size();i>0;--i,++di) det *= *di;
        } else {
            for(ptrdiff_t i=size();i>0;--i,di+=ds) det *= *di;
        }
        return det;
    }
    template <>
    std::complex<int> GenDiagMatrix<std::complex<int> >::det() const
    {
        const std::complex<int>* di = diag().cptr();
        const ptrdiff_t ds = diag().step();
        std::complex<int> det(1);
        if (ds == 1) {
            for(ptrdiff_t i=size();i>0;--i,++di) det *= *di;
        } else {
            for(ptrdiff_t i=size();i>0;--i,di+=ds) det *= *di;
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

    template <class T, int A>
    DiagMatrixView<T,A>& DiagMatrixView<T,A>::invertSelf() 
    {
        T* di = diag().ptr();
        const ptrdiff_t dstep = diag().step();

        if (dstep == 1) {
            for(ptrdiff_t i=size();i>0;--i,++di) {
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
                *di = TMV_InverseOf(*di);
            }
        } else {
            for(ptrdiff_t i=size();i>0;--i,di+=dstep) {
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
                *di = TMV_InverseOf(*di);
            }
        }
        return *this;
    }

#ifdef INST_INT
    template <>
    DiagMatrixView<int,CStyle>& DiagMatrixView<int,CStyle>::invertSelf()
    { TMVAssert(TMV_FALSE); return *this; }
    template <>
    DiagMatrixView<std::complex<int>,CStyle>&
        DiagMatrixView<std::complex<int>,CStyle>::invertSelf()
    { TMVAssert(TMV_FALSE); return *this; }
#endif

    template <class T> template <class T1> 
    void GenDiagMatrix<T>::doMakeInverse(MatrixView<T1> minv) const
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
    void GenDiagMatrix<int>::doMakeInverse(MatrixView<T1> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1>
    void GenDiagMatrix<std::complex<int> >::doMakeInverse(
        MatrixView<T1> ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T> template <class T1> 
    void GenDiagMatrix<T>::doMakeInverse(DiagMatrixView<T1> minv) const
    { (minv = *this).invertSelf(); }

#ifdef INST_INT
    template <> template <class T1>
    void GenDiagMatrix<int>::doMakeInverse(DiagMatrixView<T1> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1>
    void GenDiagMatrix<std::complex<int> >::doMakeInverse(
        DiagMatrixView<T1> ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T>
    void GenDiagMatrix<T>::doMakeInverseATA(DiagMatrixView<T> ata) const
    {
        makeInverse(ata);
        T* mi = ata.diag().ptr();
        const ptrdiff_t ds = ata.diag().step();
        if (ds==1) {
            for(ptrdiff_t i=size();i>0;--i,++mi) {
#ifdef TMVFLDEBUG
                TMVAssert(mi >= ata.diag()._first);
                TMVAssert(mi < ata.diag()._last);
#endif
                *mi = TMV_NORM(*mi);
            }
        } else {
            for(ptrdiff_t i=size();i>0;--i,mi+=ds) {
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
        DiagMatrixView<int> ) const
    { TMVAssert(TMV_FALSE); }
    template <> 
    void GenDiagMatrix<std::complex<int> >::doMakeInverseATA(
        DiagMatrixView<std::complex<int> > ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T>
    void GenDiagMatrix<T>::doMakeInverseATA(MatrixView<T> ata) const
    {
        ata.setZero();
        makeInverseATA(DiagMatrixViewOf(ata.diag()));
    }

#ifdef INST_INT
    template <> 
    void GenDiagMatrix<int>::doMakeInverseATA(MatrixView<int> ) const
    { TMVAssert(TMV_FALSE); }
    template <> 
    void GenDiagMatrix<std::complex<int> >::doMakeInverseATA(
        MatrixView<std::complex<int> > ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T>
    QuotXD<T,T> GenDiagMatrix<T>::QInverse() const
    { return QuotXD<T,T>(T(1),*this); }

#define CT std::complex<T>

    template <bool cd, class T, class Td> 
    static void DoDiagLDivEq1(
        const GenDiagMatrix<Td>& d, VectorView<T> v)
    {
        TMVAssert(v.size() == d.size());
        TMVAssert(v.ct()==NonConj);
        TMVAssert(v.size() > 0);

        const Td* di = d.diag().cptr();
        T* vi = v.ptr();
        const ptrdiff_t dstep = d.diag().step();
        const ptrdiff_t vstep = v.step();

        if (dstep == 1 && vstep == 1) {
            for(ptrdiff_t i=v.size();i>0;--i,++di,++vi) {
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
                *vi = TMV_Divide(*vi,(cd?TMV_CONJ(*di):*di));
            }
        } else {
            for(ptrdiff_t i=v.size();i>0;--i,di+=dstep,vi+=vstep) {
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
                *vi = TMV_Divide(*vi,(cd?TMV_CONJ(*di):*di));
            }
        }
    }

    template <class T, class Td> 
    static inline void DoDiagLDivEq(
        const GenDiagMatrix<Td>& d, VectorView<T> v)
    { 
        if (d.diag().isconj()) DoDiagLDivEq1<true>(d,v);
        else DoDiagLDivEq1<false>(d,v);
    }

    template <class T> template <class T1> 
    void GenDiagMatrix<T>::doLDivEq(VectorView<T1> v) const
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
    void GenDiagMatrix<int>::doLDivEq(VectorView<T1> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1>
    void GenDiagMatrix<std::complex<int> >::doLDivEq(
        VectorView<T1> ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T> template <class T1, class T0> 
    void GenDiagMatrix<T>::doLDiv(
        const GenVector<T1>& v1, VectorView<T0> v0) const
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
        const GenVector<T1>& , VectorView<T0> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1, class T0>
    void GenDiagMatrix<std::complex<int> >::doLDiv(
        const GenVector<T1>& , VectorView<T0> ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <bool rm, bool cd, class T, class Td> 
    static void RowDiagLDivEq(
        const GenDiagMatrix<Td>& d, MatrixView<T> m)
    {
        TMVAssert(d.size() == m.colsize());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(m.ct() == NonConj);
        TMVAssert(rm == m.isrm());
        TMVAssert(cd == d.diag().isconj());

        const Td* di = d.diag().cptr();
        T* mrowi = m.ptr();
        const ptrdiff_t dstep = d.diag().step();
        const ptrdiff_t stepj = m.stepj();
        const ptrdiff_t stepi = m.stepi();
        const ptrdiff_t M = m.colsize();
        const ptrdiff_t N = m.rowsize();

        for(ptrdiff_t i=M;i>0;--i,di+=dstep,mrowi+=stepi) {
            T* mij = mrowi;
            if (*di == Td(0)) {
#ifdef NOTHROW
                std::cerr<<"Singular DiagMatrix found\n";
                exit(1); 
#else
                throw SingularDiagMatrix<Td>(d);
#endif
            } else {
                Td invdi = TMV_InverseOf(cd?TMV_CONJ(*di):*di);
                for(ptrdiff_t j=N;j>0;--j,(rm?++mij:mij+=stepj)) {
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
        const GenDiagMatrix<Td>& d, MatrixView<T> m)
    {
        TMVAssert(d.size() == m.colsize());
        TMVAssert(m.rowsize() > 0);
        TMVAssert(m.colsize() > 0);
        TMVAssert(m.ct()==NonConj);
        TMVAssert(cm == m.iscm());
        TMVAssert(cd == d.diag().isconj());

        DiagMatrix<Td> invd(d.size());
        const Td* di = d.diag().cptr();
        const ptrdiff_t step = d.diag().step();
        Td* invdi = invd.diag().ptr();

        if (step == 1) {
            for(ptrdiff_t i=d.size();i>0;--i,++di,++invdi) {
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
                *invdi = TMV_InverseOf(cd?TMV_CONJ(*di):*di);
            }
        } else {
            for(ptrdiff_t i=d.size();i>0;--i,di+=step,++invdi) {
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
                *invdi = TMV_InverseOf(cd?TMV_CONJ(*di):*di);
            }
        }
        m = invd*m;
    }

    template <class T, class Td> 
    static void DoDiagLDivEq(
        const GenDiagMatrix<Td>& d, MatrixView<T> m)
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
    void GenDiagMatrix<T>::doLDivEq(MatrixView<T1> m) const
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
    void GenDiagMatrix<int>::doLDivEq(MatrixView<T1> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1>
    void GenDiagMatrix<std::complex<int> >::doLDivEq(MatrixView<T1> ) const
    { TMVAssert(TMV_FALSE); }
#endif

    template <class T> template <class T1, class T0> 
    void GenDiagMatrix<T>::doLDiv(
        const GenMatrix<T1>& m1, MatrixView<T0> m0) const
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
        const GenMatrix<T1>& , MatrixView<T0> ) const
    { TMVAssert(TMV_FALSE); }
    template <> template <class T1, class T0>
    void GenDiagMatrix<std::complex<int> >::doLDiv(
        const GenMatrix<T1>& , MatrixView<T0> ) const
    { TMVAssert(TMV_FALSE); }
#endif

#undef CT

    template <class T>
    void GenDiagMatrix<T>::write(const TMV_Writer& writer) const
    {
        const ptrdiff_t N = size();
        writer.begin();
        writer.writeCode("D");
        writer.writeSize(N);
        writer.writeSimpleSize(N);
        writer.writeStart();
        for(ptrdiff_t i=0;i<N;++i) {
            writer.writeLParen();
            if (!writer.isCompact()) {
                for(ptrdiff_t j=0;j<i;++j) {
                    if (j > 0) writer.writeSpace();
                    writer.writeValue(T(0));
                }
                if (i > 0) writer.writeSpace();
            }
            writer.writeValue(diag().cref(i));
            if (!writer.isCompact()) {
                for(ptrdiff_t j=i+1;j<N;++j) {
                    writer.writeSpace();
                    writer.writeValue(T(0));
                }
            }
            writer.writeRParen();
            if (i < N-1) writer.writeRowEnd();
        }
        writer.writeFinal();
        writer.end();
    }

#ifndef NOTHROW
    template <class T>
    class DiagMatrixReadError : public ReadError
    {
    public :
        DiagMatrix<T> m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t s;
        T v1;
        bool is, iseof, isbad;

        DiagMatrixReadError(std::istream& _is) throw() :
            ReadError("DiagMatrix."),
            i(0), j(0), s(0), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("DiagMatrix."),
            i(0), j(0), exp(_e), got(_g), s(0), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        DiagMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenDiagMatrix<T>& _m, std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("DiagMatrix."),
            m(_m), i(_i), j(_j), exp(_e), got(_g), s(m.size()), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(
            const GenDiagMatrix<T>& _m, std::istream& _is, ptrdiff_t _s) throw() :
            ReadError("DiagMatrix."),
            m(_m), i(0), j(0), s(_s), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const GenDiagMatrix<T>& _m,
            std::istream& _is, T _v1=0) throw() :
            ReadError("DiagMatrix."),
            m(_m), i(_i), j(_j), s(m.size()), v1(_v1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        DiagMatrixReadError(const DiagMatrixReadError<T>& rhs) :
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got),
            s(rhs.s), v1(rhs.v1),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~DiagMatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for DiagMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (s != m.size()) {
                os<<"Wrong size: expected "<<m.size()<<", got "<<s<<".\n";
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
            if (v1 != T(0)) {
                os<<"Invalid input: Expected 0, got "<<v1<<".\n";
            }
            if (m.size() > 0) {
                os<<"The portion of the DiagMatrix which was successfully "
                    "read is: \n";
                const ptrdiff_t N = m.size();
                for(ptrdiff_t ii=0;ii<i;++ii) {
                    os<<"( ";
                    for(ptrdiff_t jj=0;jj<N;++jj) os<<' '<<m.cref(ii,jj)<<' ';
                    os<<" )\n";
                }
                os<<"( ";
                for(ptrdiff_t jj=0;jj<j;++jj) os<<' '<<m.cref(i,jj)<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    template <class T> 
    static void FinishRead(const TMV_Reader& reader, DiagMatrixView<T> m)
    {
        const ptrdiff_t N = m.size();
        std::string exp, got;
        T temp;
        if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw DiagMatrixReadError<T>(0,0,m,reader.getis(),exp,got);
#endif
        }
        for(ptrdiff_t i=0;i<N;++i) {
            if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw DiagMatrixReadError<T>(i,0,m,reader.getis(),exp,got);
#endif
            }
            if (!reader.isCompact()) {
                for(ptrdiff_t j=0;j<i;++j) {
                    if (j>0 && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw DiagMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"DiagMatrix Read Error: reading value\n";
                        exit(1);
#else
                        throw DiagMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    if (temp != T(0)) {
#ifdef NOTHROW
                        std::cerr<<"DiagMatrix Read Error: "<<temp<<" != 0\n";
                        exit(1);
#else
                        throw DiagMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                    }
                }
                if (i>0 && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw DiagMatrixReadError<T>(i,i,m,reader.getis(),exp,got);
#endif
                }
            }
            if (!reader.readValue(temp)) {
#ifdef NOTHROW
                std::cerr<<"DiagMatrix Read Error: reading value\n";
                exit(1);
#else
                throw DiagMatrixReadError<T>(i,i,m,reader.getis());
#endif
            }
            m.diag().ref(i) = temp;
            if (!reader.isCompact()) {
                for(ptrdiff_t j=i+1;j<N;++j) {
                    if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw DiagMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"DiagMatrix Read Error: reading value\n";
                        exit(1);
#else
                        throw DiagMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    if (temp != T(0)) {
#ifdef NOTHROW
                        std::cerr<<"DiagMatrix Read Error: "<<temp<<" != 0\n";
                        exit(1);
#else
                        throw DiagMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                    }
                }
            }
            if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw DiagMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
            if (i < N-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw DiagMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
            }
        }
        if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw DiagMatrixReadError<T>(N,0,m,reader.getis(),exp,got);
#endif
        }
    }

    template <class T, int A>
    void DiagMatrix<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        if (!reader.readCode("D",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: reading size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) resize(s);
        s=size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: reading size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: Wrong size\n";
            exit(1); 
#else
            throw DiagMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        DiagMatrixView<T> v = view();
        FinishRead(reader,v);
    }

    template <class T, int A>
    void DiagMatrixView<T,A>::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        if (!reader.readCode("D",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: reading size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: Wrong size\n";
            exit(1); 
#else
            throw DiagMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        s=size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: reading size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (s != size()) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: Wrong size\n";
            exit(1); 
#else
            throw DiagMatrixReadError<T>(*this,reader.getis(),s);
#endif
        }
        FinishRead(reader,*this);
    }

#undef RT
#undef CT

#define InstFile "TMV_DiagMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


