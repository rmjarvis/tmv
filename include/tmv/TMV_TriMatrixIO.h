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

#ifndef TMV_TriMatrixIO_H
#define TMV_TriMatrixIO_H

#include "TMV_VectorIO.h"
#include "TMV_BaseMatrix_Tri.h"

namespace tmv {

    //
    // Write Matrix
    //

    template <class M>
    static void InlineWriteCompact(
        std::ostream& os, const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::value_type T;
        const int len = m.size();
        const bool upper = m.isupper();
        const bool unit = m.isunit();
        os << (upper ? "U " : "L ");
        os << len << std::endl;
        for(int i=0;i<len;++i) {
            os << "( ";
            if (!upper) {
                for(int j=0;j<i;++j) {
                    os << " " << Value(m.cref(i,j)) << " ";
                }
            }
            if (unit)
                os << " " << Value(T(1)) << " ";
            else
                os << " " << Value(m.cref(i,i)) << " ";
            if (upper) {
                for(int j=i+1;j<len;++j) {
                    os << " " << Value(m.cref(i,j)) << " ";
                }
            }
            os << " )\n";
        }
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T, bool C>
    void InstWriteCompact(
        std::ostream& os,
        const ConstUpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C>& m);
    template <class T, bool C>
    void InstWriteCompact(
        std::ostream& os,
        const ConstLowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C>& m);

    template <bool inst, class M>
    struct CallWriteU;

    template <class M>
    struct CallWriteU<false,M> // inst = false
    {
        static void call(std::ostream& os, const M& m)
        { InlineWriteCompact(os,m); }
    };
    template <class M>
    struct CallWriteU<true,M> // inst = true
    {
        static void call(std::ostream& os, const M& m)
        { InstWriteCompact(os,m.calc().xdView()); }
    };

    template <class M>
    static void WriteCompact(std::ostream& os, const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::value_type T;
        const bool inst = 
            (M::_size == UNKNOWN || M::_size > 16) &&
            Traits<T>::isinst;
        CallWriteU<inst,M>::call(os,m.mat());
    }

    // With thresh:
    template <class M>
    static void InlineWriteCompact(
        std::ostream& os,
        const BaseMatrix_Tri<M>& m, typename M::float_type thresh) 
    {
        typedef typename M::value_type T;
        const int len = m.size();
        const bool upper = m.isupper();
        const bool unit = m.isunit();
        os << (upper ? "U " : "L ");
        os << len << std::endl;
        for(int i=0;i<len;++i) {
            os << "( ";
            if (!upper) {
                for(int j=0;j<i;++j) {
                    T temp = m.cref(i,j);
                    os << " " << 
                        Value((TMV_ABS2(temp) < thresh ? T(0) : temp)) << " ";
                }
            }
            if (unit)
                os << " " << Value(T(1)) << " ";
            else
                os << " " << Value(m.cref(i,i)) << " ";
            if (upper) {
                for(int j=i+1;j<len;++j) {
                    T temp = m.cref(i,j);
                    os << " " << 
                        Value((TMV_ABS2(temp) < thresh ? T(0) : temp)) << " ";
                }
            }
            os << " )\n";
        }
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T, bool C>
    void InstWriteCompact(
        std::ostream& os,
        const ConstUpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C>& m,
        typename ConstUpperTriMatrixView<T>::float_type thresh);
    template <class T, bool C>
    void InstWriteCompact(
        std::ostream& os,
        const ConstLowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C>& m,
        typename ConstLowerTriMatrixView<T>::float_type thresh);

    template <bool inst, class M>
    struct CallWriteUThresh;

    template <class M>
    struct CallWriteUThresh<false,M> // inst = false
    {
        static void call(
            std::ostream& os, const M& m, typename M::float_type thresh) 
        { InlineWriteCompact(os,m,thresh); }
    };
    template <class M>
    struct CallWriteUThresh<true,M> // inst = true
    {
        static void call(
            std::ostream& os, const M& m, typename M::float_type thresh) 
        { InstWriteCompact(os,m.calc().xdView(),thresh); }
    };

    template <class M>
    static void WriteCompact(
        std::ostream& os, const BaseMatrix_Tri<M>& m,
        typename M::float_type thresh) 
    {
        typedef typename M::value_type T;
        const bool inst = 
            (M::_size == UNKNOWN || M::_size > 16) &&
            Traits<T>::isinst;
        CallWriteUThresh<inst,M>::call(os,m.mat(),thresh);
    }


    //
    // Read Matrix
    //

#ifndef TMV_NO_THROW
    template <class M> 
    class TriMatrixReadError :
        public ReadError
    {
    public :
        typedef typename M::copy_type copy_type;
        int i,j;
        mutable auto_ptr<copy_type> m;
        char exp,got;
        typedef typename M::value_type T;
        T unitgot;
        size_t s;
        bool is, iseof, isbad;

        TriMatrixReadError(std::istream& _is) throw() :
            ReadError("TriMatrix"),
            i(0), j(0), m(0), exp(0), got(0), unitgot(T(1)), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        TriMatrixReadError(
            int _i, int _j, const BaseMatrix_Tri<M>& _m, 
            std::istream& _is
        ) throw() :
            ReadError("TriMatrix"),
            i(_i), j(_j), m(new copy_type(_m)), exp(0), got(0), 
            unitgot(T(1)), s(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        TriMatrixReadError(
            int _i, int _j, const BaseMatrix_Tri<M>& _m,
            std::istream& _is, char _e, char _g
        ) throw() :
            ReadError("TriMatrix"),
            i(_i), j(_j), m(new copy_type(_m)), exp(_e), got(_g),
            unitgot(T(1)), s(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        TriMatrixReadError(
            std::istream& _is, char _e, char _g
        ) throw() :
            ReadError("TriMatrix"),
            i(0), j(0), m(0), exp(_e), got(_g), unitgot(T(1)), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        TriMatrixReadError(
            int _i, int _j, const BaseMatrix_Tri<M>& _m,
            std::istream& _is, T _u
        ) throw() :
            ReadError("TriMatrix"),
            i(_i), j(_j), m(new copy_type(_m)), exp(0), got(0),
            unitgot(_u), s(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        TriMatrixReadError(
            const BaseMatrix_Tri<M>& _m, std::istream& _is, size_t _s
        ) throw():
            ReadError("TriMatrix"),
            i(0), m(new copy_type(_m)), exp(0), got(0), unitgot(T(1)), s(_s),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        TriMatrixReadError(const TriMatrixReadError<M>& rhs) throw() :
            ReadError("TriMatrix"),
            i(rhs.i), j(rhs.j), m(rhs.m), exp(rhs.exp), got(rhs.got), 
            unitgot(T(1)), s(rhs.s),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~TriMatrixReadError() throw() {}

        void Write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for ";
            if (M::_upper) os<<"UpperTriMatrix\n";
            else os<<"LowerTriMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (unitgot != T(1)) {
                os<<"Wrong format: expected 1 on the diagonal, got '"<<
                    unitgot<<"'.\n";
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
                const int N = m->size();
                os<<"The portion of the TriMatrix which was successfully "
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
    static void InlineRead(std::istream& is, BaseMatrix_Tri_Mutable<M>& m)
    {
        char paren;
        typedef typename M::value_type T;
        T temp;
        const int len = m.size();
        const bool upper = m.isupper();
        const bool unit = m.isunit();
        for(int i=0;i<len;++i) {
            is >> paren;
            if (!is || paren != '(') {
#ifdef TMV_NO_THROW
                std::cerr<<"TriMatrix ReadError: "<<paren<<" != (\n"; 
                exit(1); 
#else
                throw TriMatrixReadError<M>(i,0,m,is,'(',is?paren:'(');
#endif
            }
            if (!upper) {
                const int i1 = unit ? i : i+1;
                for(int j=0;j<i1;++j) {
                    is >> temp;
                    if (!is) {
#ifdef TMV_NO_THROW
                        std::cerr<<"TriMatrix ReadError: !is\n"; 
                        exit(1); 
#else
                        throw TriMatrixReadError<M>(i,j,m,is);
#endif
                    }
                    m.ref(i,j) = temp;
                }
            }
            if (unit) {
                is >> temp;
                if (!is || temp != T(1)) {
#ifdef TMV_NO_THROW
                    std::cerr<<"TriMatrix ReadError: "<<temp<<" != 1\n"; 
                    exit(1); 
#else
                    throw TriMatrixReadError<M>(i,i,m,is,is?temp:T(1));
#endif
                }
            }
            if (upper) {
                const int i1 = unit ? i+1 : i;
                for(int j=i1;j<len;++j) {
                    is >> temp;
                    if (!is)  {
#ifdef TMV_NO_THROW
                        std::cerr<<"TriMatrix ReadError: !is\n"; 
                        exit(1); 
#else
                        throw TriMatrixReadError<M>(i,j,m,is);
#endif
                    }
                    m.ref(i,j) = temp;
                } 
            }
            is >> paren;
            if (!is || paren != ')') {
#ifdef TMV_NO_THROW
                std::cerr<<"TriMatrix ReadError: "<<paren<<" != )\n"; 
                exit(1); 
#else
                throw TriMatrixReadError<M>(i,len,m,is,')',is?paren:')');
#endif
            }
        }
    }

    // Defined in TMV_TriMatrix.cpp
    template <class T, bool C>
    void InstRead(
        std::istream& is,
        UpperTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> m);
    template <class T, bool C>
    void InstRead(
        std::istream& is,
        LowerTriMatrixView<T,UnknownDiag,UNKNOWN,UNKNOWN,C> m);

    template <bool inst, class M>
    struct CallReadU;

    template <class M>
    struct CallReadU<false,M> // inst = false
    {
        static void call(std::istream& is, M& m)
        { InlineRead(is,m); }
    };
    template <class M>
    struct CallReadU<true,M>
    {
        static void call(std::istream& is, M& m)
        { InstRead(is,m.xdView()); }
    };

    template <class M>
    static void Read(std::istream& is, BaseMatrix_Tri_Mutable<M>& m)
    {
        typedef typename M::value_type T;
        const bool inst = 
            (M::_size == UNKNOWN || M::_size > 16) &&
            Traits<T>::isinst;
        CallReadU<inst,M>::call(is,m.mat());
    }


    // 
    // Operator overloads for I/O
    // is >> m
    //

    template <class M>
    static std::istream& operator>>(
        std::istream& is, BaseMatrix_Tri_Mutable<M>& m)
    {
        char ul;
        char ul_exp = (M::_upper ? 'U' : 'L');
        is >> ul;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"TriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw TriMatrixReadError<M>(is);
#endif
        }
        if (ul != ul_exp) {
#ifdef TMV_NO_THROW
            std::cerr<<"TriMatrix ReadError: "<<ul<<" != "<<ul_exp; 
            exit(1); 
#else
            throw TriMatrixReadError<M>(is,ul_exp,ul);
#endif
        }

        size_t s;
        is >> s;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"TriMatrix ReadError: !is \n"; 
            exit(1); 
#else
            throw TriMatrixReadError<M>(is);
#endif
        }
        if (s != m.size()) {
#ifdef TMV_NO_THROW
            std::cerr<<"TriMatrix ReadError: Wrong size \n"; 
            exit(1); 
#else
            throw TriMatrixReadError<M>(m,is,s);
#endif
        }
        Read(is,m);
        return is;
    }

    template <class T, DiagType D, StorageType S, IndexStyle I>
    static std::istream& operator>>(
        std::istream& is, auto_ptr<UpperTriMatrix<T,D,S,I> >& m)
    {
        typedef UpperTriMatrix<T,D,S,I> M;
        char ul;
        is >> ul;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"UpperTriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw TriMatrixReadError<M>(is);
#endif
        }
        if (ul != 'U') {
#ifdef TMV_NO_THROW
            std::cerr<<"UpperTriMatrix ReadError: "<<ul<<" != U\n"; 
            exit(1); 
#else
            throw TriMatrixReadError<M>(is,'U',ul);
#endif
        }
        size_t s;
        is >> s;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"UpperTriMatrix ReadError: !is \n"; 
            exit(1); 
#else
            throw TriMatrixReadError<M>(is);
#endif
        }
        m.reset(new M(s));
        Read(is,*m);
        return is;
    }

    template <class T, DiagType D, StorageType S, IndexStyle I>
    static std::istream& operator>>(
        std::istream& is, auto_ptr<LowerTriMatrix<T,D,S,I> >& m)
    {
        typedef LowerTriMatrix<T,D,S,I> M;
        char ul;
        is >> ul;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"LowerTriMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw TriMatrixReadError<M>(is);
#endif
        }
        if (ul != 'L') {
#ifdef TMV_NO_THROW
            std::cerr<<"LowerTriMatrix ReadError: "<<ul<<" != U\n"; 
            exit(1); 
#else
            throw TriMatrixReadError<M>(is,'L',ul);
#endif
        }
        size_t s;
        is >> s;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"LowerTriMatrix ReadError: !is \n"; 
            exit(1); 
#else
            throw TriMatrixReadError<M>(is);
#endif
        }
        m.reset(new M(s));
        Read(is,*m);
        return is;
    }

} // namespace mv

#endif
