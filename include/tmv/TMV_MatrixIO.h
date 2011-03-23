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

#ifndef TMV_MatrixIO_H
#define TMV_MatrixIO_H

#include "TMV_BaseMatrix_Rec.h"

namespace tmv {

    //
    // Write Matrix
    //

    // Defined in TMV_Matrix.cpp
    template <class T, int C>
    void InstWrite(
        std::ostream& os, const ConstMatrixView<T,C>& m);
    template <class T, int C>
    void InstWrite(
        std::ostream& os, const ConstMatrixView<T,C>& m,
        typename ConstMatrixView<T>::float_type thresh);
    template <class T>
    void InstRead(std::istream& is, MatrixView<T> m);

    template <int algo, class M>
    struct WriteM_Helper;

    template <class M>
    struct WriteM_Helper<11,M>
    {
        static void call(std::ostream& os, const M& m)
        {
            const int nrows = m.nrows();
            const int ncols = m.ncols();
            os << nrows << "  " << ncols << std::endl;
            for(int i=0;i<nrows;++i) {
                os << "( ";
                for(int j=0;j<ncols;++j) {
                    os << " " << Value(m.cref(i,j)) << " ";
                }
                os << " )\n";
            }
        }
        static void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        {
            typedef typename M::value_type T;
            const int nrows = m.nrows();
            const int ncols = m.ncols();
            os << nrows << "  " << ncols << std::endl;
            for(int i=0;i<nrows;++i) {
                os << "( ";
                for(int j=0;j<ncols;++j) {
                    T temp = m.cref(i,j);
                    os << " " << Value((TMV_ABS2(temp) < thresh ? T(0) : temp)) 
                        << " ";
                }
                os << " )\n";
            }
        }
    };

    // algo 90: Call inst
    template <class M>
    struct WriteM_Helper<90,M>
    {
        static void call(std::ostream& os, const M& m)
        { InstWrite(os,m.calc().xView()); }
        static void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { InstWrite(os,m.calc().xView(),thresh); }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M>
    struct WriteM_Helper<-3,M>
    {
        static void call(std::ostream& os, const M& m)
        { WriteM_Helper<11,M>::call(os,m); }
        static void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { WriteM_Helper<11,M>::call(os,m,thresh); }
    };
             
    // algo -2: Check for inst
    template <class M>
    struct WriteM_Helper<-2,M>
    {
        typedef typename M::value_type T;
        enum { inst = (
                (M::_colsize == UNKNOWN || M::_colsize > 16) &&
                (M::_rowsize == UNKNOWN || M::_rowsize > 16) &&
                Shape(M::_shape) == Rec &&
                Traits<T>::isinst ) };
        enum { algo = (
                inst ? 90 :
                -3 ) };
        static void call(std::ostream& os, const M& m)
        { WriteM_Helper<algo,M>::call(os,m); }
        static void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { WriteM_Helper<algo,M>::call(os,m,thresh); }
    };

    template <class M>
    struct WriteM_Helper<-1,M>
    {
        static void call(std::ostream& os, const M& m)
        { WriteM_Helper<-2,M>::call(os,m); }
        static void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { WriteM_Helper<-2,M>::call(os,m,thresh); }
    };

    template <class M>
    static inline void Write(std::ostream& os, const BaseMatrix_Calc<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteM_Helper<-2,Mv>::call(os,mv);
    }

    template <class M>
    static inline void InlineWrite(
        std::ostream& os, const BaseMatrix_Calc<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteM_Helper<-3,Mv>::call(os,mv);
    }

    template <class M>
    static inline void Write(
        std::ostream& os,
        const BaseMatrix_Calc<M>& m, typename M::float_type thresh) 
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteM_Helper<-2,Mv>::call(os,mv,thresh);
    }


    template <class M>
    static inline void InlineWrite(
        std::ostream& os,
        const BaseMatrix_Calc<M>& m, typename M::float_type thresh) 
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteM_Helper<-3,Mv>::call(os,mv,thresh);
    }


    //
    // Read Matrix
    //

#ifndef TMV_NO_THROW
    template <class T>
    class MatrixReadError : 
        public ReadError
    {
    public :
        Matrix<T,NoDivider> m;
        int i,j;
        char exp,got;
        size_t cs,rs;
        bool is, iseof, isbad;

        MatrixReadError(std::istream& _is) throw() :
            ReadError("Matrix"),
            i(0), j(0), exp(0), got(0), cs(0), rs(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        MatrixReadError(
            int _i, int _j, const BaseMatrix<M>& _m, 
            std::istream& _is) throw() :
            ReadError("Matrix"),
            m(_m), i(_i), j(_j), exp(0), got(0), 
            cs(_m.colsize()), rs(_m.rowsize()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        MatrixReadError(
            int _i, int _j, const BaseMatrix<M>& _m,
            std::istream& _is, char _e, char _g) throw() :
            ReadError("Matrix"),
            m(_m), i(_i), j(_j), exp(_e), got(_g),
            cs(_m.colsize()), rs(_m.rowsize()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        MatrixReadError(
            const BaseMatrix<M>& _m,
            std::istream& _is, size_t _cs, size_t _rs) :
            ReadError("Matrix"),
            m(_m), i(0), exp(0), got(0), cs(_cs), rs(_rs),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        MatrixReadError(const MatrixReadError<T>& rhs) :
            ReadError("Matrix"),
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got), 
            cs(rhs.cs), rs(rhs.rs),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~MatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for Matrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (cs != m.colsize()) {
                os<<"Wrong column size: expected "<<m.colsize()<<
                    ", got "<<cs<<".\n";
            }
            if (rs != m.rowsize()) {
                os<<"Wrong row size: expected "<<m.rowsize()<<
                    ", got "<<rs<<".\n";
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
            if (m.colsize() > 0 || m.rowsize() > 0) {
                const int N = m.rowsize();
                os<<"The portion of the Matrix which was successfully "
                    "read is: \n";
                for(int ii=0;ii<i;++ii) {
                    os<<"( ";
                    for(int jj=0;jj<N;++jj)
                        os<<' '<<m.cref(ii,jj)<<' ';
                    os<<" )\n";
                }
                os<<"( ";
                for(int jj=0;jj<j;++jj)
                    os<<' '<<m.cref(i,jj)<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    template <int algo, class M>
    struct ReadM_Helper;

    template <class M>
    struct ReadM_Helper<11,M>
    {
        static void call(std::istream& is, M& m)
        {
            typedef typename M::value_type T;
            char paren;
            T temp;
            const int nrows = m.nrows();
            const int ncols = m.ncols();
            for(int i=0;i<nrows;++i) {
                is >> paren;
                if (!is || paren != '(') {
#ifdef TMV_NO_THROW
                    std::cerr<<"Matrix ReadError: "<<paren<<" != (\n"; 
                    exit(1); 
#else
                    throw MatrixReadError<T>(i,0,m,is,'(',is?paren:'(');
#endif
                }
                for(int j=0;j<ncols;++j) {
                    is >> temp;
                    if (!is) {
#ifdef TMV_NO_THROW
                        std::cerr<<"Matrix ReadError: !is\n"; 
                        exit(1); 
#else
                        throw MatrixReadError<T>(i,j,m,is);
#endif
                    }
                    m.ref(i,j) = temp;
                } 
                is >> paren;
                if (!is || paren != ')') {
#ifdef TMV_NO_THROW
                    std::cerr<<"Matrix ReadError: "<<paren<<" != )\n"; 
                    exit(1); 
#else
                    throw MatrixReadError<T>(i,ncols,m,is,')',is?paren:')');
#endif
                }
            }
        }
    };

    // algo 90: Call inst
    template <class M>
    struct ReadM_Helper<90,M>
    {
        static void call(std::istream& is, M& m)
        { InstRead(is,m.xView()); }
    };
             
    // algo 97: Conjugate
    template <class M>
    struct ReadM_Helper<97,M>
    {
        static void call(std::istream& is, M& m)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            ReadM_Helper<-2,Mc>::call(is,mc); 
            mc.conjugateSelf();
        }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M>
    struct ReadM_Helper<-3,M>
    {
        static void call(std::istream& is, M& m)
        { ReadM_Helper<11,M>::call(is,m); }
    };
             
    // algo -2: Check for inst
    template <class M>
    struct ReadM_Helper<-2,M>
    {
        static void call(std::istream& is, M& m)
        {
            typedef typename M::value_type T;
            const int inst = 
                (M::_colsize == UNKNOWN || M::_colsize > 16) &&
                (M::_rowsize == UNKNOWN || M::_rowsize > 16) &&
                Traits<T>::isinst;
            const int algo = 
                M::_conj ? 97 :
                inst ? 90 :
                -3;
            ReadM_Helper<algo,M>::call(is,m); 
        }
    };

    template <class M>
    struct ReadM_Helper<-1,M>
    {
        static void call(std::istream& is, M& m)
        { ReadM_Helper<-2,M>::call(is,m); }
    };

    template <class M>
    static inline void Read(std::istream& is, BaseMatrix_Rec_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadM_Helper<-2,Mv>::call(is,mv);
    }

    template <class M>
    static inline void InlineRead(
        std::istream& is, BaseMatrix_Rec_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadM_Helper<-3,Mv>::call(is,mv);
    }


    // 
    // Operator overloads for I/O
    // is >> m
    // os << m
    //

    template <class M>
    static inline std::ostream& operator<<(
        std::ostream& os, const BaseMatrix<M>& m)
    { Write(os,m.calc()); return os; }

    template <class M>
    static inline std::istream& operator>>(
        std::istream& is, BaseMatrix_Rec_Mutable<M>& m)
    {
        typedef typename M::value_type T;
        size_t cs,rs;
        is >> cs >> rs;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"Matrix ReadError: !is \n"; 
            exit(1); 
#else
            throw MatrixReadError<T>(is);
#endif
        }
        if (cs != m.colsize() || rs != m.rowsize()) {
#ifdef TMV_NO_THROW
            std::cerr<<"Matrix ReadError: Wrong size \n"; 
            exit(1); 
#else
            throw MatrixReadError<T>(m,is,cs,rs);
#endif
        }
        Read(is,m);
        return is;
    }

    template <class T, int A0, int A1>
    static inline std::istream& operator>>(
        std::istream& is, auto_ptr<Matrix<T,A0,A1> >& m)
    {
        size_t cs,rs;
        is >> cs >> rs;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"Matrix ReadError: !is \n"; 
            exit(1); 
#else
            throw MatrixReadError<T>(is);
#endif
        }
        m.reset(new Matrix<T,A0,A1>(cs,rs));
        Read(is,*m);
        return is;
    }
} // namespace mv

#endif
