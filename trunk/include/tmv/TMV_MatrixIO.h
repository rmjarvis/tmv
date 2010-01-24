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

#include "TMV_VectorIO.h"
#include "TMV_BaseMatrix_Rec.h"

namespace tmv {

    //
    // Write Matrix
    //

    template <class M>
    inline void InlineWrite(std::ostream& os, const BaseMatrix_Calc<M>& m)
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

    // Defined in TMV_Matrix.cpp
    template <class T, bool C>
    void InstWrite(std::ostream& os,
                   const ConstMatrixView<T,UNKNOWN,UNKNOWN,C>& m);

    template <bool inst, class M>
    struct CallWriteM;

    template <class M>
    struct CallWriteM<false,M> // inst = false
    {
        static inline void call(std::ostream& os, const M& m)
        { InlineWrite(os,m); }
    };
    template <class M>
    struct CallWriteM<true,M>
    {
        static inline void call(std::ostream& os, const M& m)
        { InstWrite(os,m.calc().xView()); }
    };

    template <class M>
    inline void Write(std::ostream& os, const BaseMatrix_Calc<M>& m)
    {
        typedef typename M::value_type T;
        const int inst = 
                Traits<T>::isinst &&
                Shape(M::mshape) == Rec &&
                (M::mrowmajor || M::mcolmajor) &&
                M::mcolsize == UNKNOWN &&
                M::mrowsize == UNKNOWN;
        CallWriteM<inst,M>::call(os,m.mat());
    }

    // With thresh:
    template <class M>
    inline void InlineWrite(
        std::ostream& os, const BaseMatrix_Calc<M>& m,
        typename M::real_type thresh) 
    {
        typedef typename M::value_type T;
        const int nrows = m.nrows();
        const int ncols = m.ncols();
        os << nrows << "  " << ncols << std::endl;
        for(int i=0;i<nrows;++i) {
            os << "( ";
            for(int j=0;j<ncols;++j) {
                T temp = m.cref(i,j);
                os << " " << Value((TMV_ABS(temp) < thresh ? T(0) : temp)) 
                    << " ";
            }
            os << " )\n";
        }
    }

    // Defined in TMV_Matrix.cpp
    template <class T, bool C>
    void InstWrite(
        std::ostream& os, const ConstMatrixView<T,UNKNOWN,UNKNOWN,C>& m,
        typename Traits<T>::real_type thresh);

    template <bool inst, class M>
    struct CallWriteMThresh;

    template <class M>
    struct CallWriteMThresh<false,M> // inst = false
    {
        static inline void call(
            std::ostream& os, const M& m, typename M::real_type thresh) 
        { InlineWrite(os,m,thresh); }
    };
    template <class M>
    struct CallWriteMThresh<true,M>
    {
        static inline void call(
            std::ostream& os, const M& m, typename M::real_type thresh) 
        { InstWrite(os,m.calc().xView(),thresh); }
    };

    template <class M>
    inline void Write(
        std::ostream& os,
        const BaseMatrix_Calc<M>& m, typename M::real_type thresh) 
    {
        typedef typename M::value_type T;
        const int inst = 
                Traits<T>::isinst &&
                Shape(M::mshape) == Rec &&
                (M::mrowmajor || M::mcolmajor) &&
                M::mcolsize == UNKNOWN &&
                M::mrowsize == UNKNOWN;
        CallWriteMThresh<inst,M>::call(os,m.mat(),thresh);
    }


    //
    // Read Matrix
    //

#ifndef TMV_NO_THROW
    template <class M> 
    class MatrixReadError : public ReadError
    {
    public :
        typedef typename M::copy_type copy_type;
        int i,j;
        mutable auto_ptr<copy_type> m;
        char exp,got;
        size_t cs,rs;
        bool is, iseof, isbad;

        inline MatrixReadError(std::istream& _is) throw() :
            ReadError("Matrix"),
            i(0), j(0), m(0), exp(0), got(0), cs(0), rs(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        inline MatrixReadError(int _i, int _j, const BaseMatrix<M>& _m, 
                               std::istream& _is) throw() :
            ReadError("Matrix"),
            i(_i), j(_j), m(new copy_type(_m)), exp(0), got(0), 
            cs(_m.colsize()), rs(_m.rowsize()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        inline MatrixReadError(int _i, int _j, const BaseMatrix<M>& _m,
                               std::istream& _is, char _e, char _g) throw() :
            ReadError("Matrix"),
            i(_i), j(_j), m(new copy_type(_m)), exp(_e), got(_g),
            cs(_m.colsize()), rs(_m.rowsize()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        inline MatrixReadError(const BaseMatrix<M>& _m,
                               std::istream& _is, size_t _cs, size_t _rs) :
            ReadError("Matrix"),
            i(0), m(new copy_type(_m)), exp(0), got(0), cs(_cs), rs(_rs),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        inline MatrixReadError(const MatrixReadError<M>& rhs) :
            ReadError("Matrix"),
            i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got), 
            cs(rhs.cs), rs(rhs.rs),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        inline ~MatrixReadError() throw() {}

        inline void Write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for Matrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (m.get() && cs != m->colsize()) {
                os<<"Wrong column size: expected "<<m->colsize()<<
                    ", got "<<cs<<".\n";
            }
            if (m.get() && rs != m->rowsize()) {
                os<<"Wrong row size: expected "<<m->rowsize()<<
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
            if (m.get()) {
                const int N = m->rowsize();
                os<<"The portion of the Matrix which was successfully "
                    "read is: \n";
                for(int ii=0;ii<i;++ii) {
                    os<<"( ";
                    for(int jj=0;jj<N;++jj)
                        os<<' '<<(*m).cref(ii,jj)<<' ';
                    os<<" )\n";
                }
                os<<"( ";
                for(int jj=0;jj<j;++jj)
                    os<<' '<<(*m).cref(i,jj)<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    template <class M>
    inline void InlineRead(std::istream& is, BaseMatrix_Rec_Mutable<M>& m)
    {
        char paren;
        typename M::value_type temp;
        const int nrows = m.nrows();
        const int ncols = m.ncols();
        for(int i=0;i<nrows;++i) {
            is >> paren;
            if (!is || paren != '(') {
#ifdef TMV_NO_THROW
                std::cerr<<"Matrix ReadError: "<<paren<<" != (\n"; 
                exit(1); 
#else
                throw MatrixReadError<M>(i,0,m,is,'(',is?paren:'(');
#endif
            }
            for(int j=0;j<ncols;++j) {
                is >> temp;
                if (!is) {
#ifdef TMV_NO_THROW
                    std::cerr<<"Matrix ReadError: !is\n"; 
                    exit(1); }
#else
                    throw MatrixReadError<M>(i,j,m,is);
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
            throw MatrixReadError<M>(i,ncols,m,is,')',is?paren:')');
#endif
        }
    }
}

// Defined in TMV_Matrix.cpp
template <class T, bool C>
void InstRead(std::istream& is, MatrixView<T,UNKNOWN,UNKNOWN,C> m);

template <bool inst, class M>
struct CallReadM;

template <class M>
struct CallReadM<false,M> // inst = false
{
    static inline void call(std::istream& is, M& m)
    { InlineRead(is,m); }
};
template <class M>
struct CallReadM<true,M>
{
    static inline void call(std::istream& is, M& m)
    { InstRead(is,m.xView()); }
};

template <class M>
inline void Read(std::istream& is, BaseMatrix_Rec_Mutable<M>& m)
{
    typedef typename M::value_type T;
    const int inst = 
        Traits<T>::isinst &&
        (M::mrowmajor || M::mcolmajor) &&
        M::mcolsize == UNKNOWN &&
        M::mrowsize == UNKNOWN;
    CallReadM<inst,M>::call(is,m.mat());
}


// 
// Operator overloads for I/O
// is >> m
// os << m
//

template <class M>
inline std::ostream& operator<<(
    std::ostream& os, const BaseMatrix<M>& m)
{ Write(os,m.calc()); return os; }

template <class M>
inline std::istream& operator>>(
    std::istream& is, BaseMatrix_Rec_Mutable<M>& m)
{
    size_t cs,rs;
    is >> cs >> rs;
    if (!is) {
#ifdef TMV_NO_THROW
        std::cerr<<"Matrix ReadError: !is \n"; 
        exit(1); 
#else
        throw MatrixReadError<M>(is);
#endif
    }
    if (cs != m.colsize() || rs != m.rowsize()) {
#ifdef TMV_NO_THROW
        std::cerr<<"Matrix ReadError: Wrong size \n"; 
        exit(1); 
#else
        throw MatrixReadError<M>(m,is,cs,rs);
#endif
    }
    Read(is,m);
    return is;
}

template <class T, StorageType S, IndexStyle I>
inline std::istream& operator>>(std::istream& is, 
                                auto_ptr<Matrix<T,S,I> >& m)
{
    size_t cs,rs;
    is >> cs >> rs;
    if (!is) {
#ifdef TMV_NO_THROW
        std::cerr<<"Matrix ReadError: !is \n"; 
        exit(1); 
#else
        throw MatrixReadError<Matrix<T,S,I> >(is);
#endif
    }
    m.reset(new Matrix<T,S,I>(cs,rs));
    Read(is,*m);
    return is;
}
} // namespace mv

#endif
