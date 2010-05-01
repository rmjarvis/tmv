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

#ifndef TMV_VectorIO_H
#define TMV_VectorIO_H

#include "../util/portable_platform.h"
#include "TMV_BaseVector.h"

namespace tmv {

    //
    // Write Vector
    //

    // This bit is to workaround a bug in pgCC that was fixed in version 7.
    // I don't know if versions earlier than 6.1 had the bug, but 
    // I apply the workaround to all version before 7.
    template <class T> inline T Value(T x) { return x; }
#ifdef PLATFORM_COMPILER_PGI
#if PLATFORM_COMPILER_VERSION < 0x070000
    inline double Value(long double x) { return double(x); }
    inline std::complex<double> Value(std::complex<long double> x)
    { return std::complex<double>(x); }
#endif
#endif

    template <class V>
    inline void InlineWrite(std::ostream& os, const BaseVector<V>& v)
    {
        const int n=v.size();
        os << n << " (";
        for(int i=0;i<n;++i) os << " " << Value(v.cref(i)) << " ";
        os << ")";
    }

    // Defined in TMV_Vector.cpp
    template <class T, bool C>
    void InstWrite(std::ostream& os, const ConstVectorView<T,UNKNOWN,C>& v);

    template <bool inst, class V>
    struct CallWrite // inst = false
    {
        static inline void call(std::ostream& os, const V& v)
        { InlineWrite(os,v.calc()); }
    };
    template <class V>
    struct CallWrite<true,V>
    {
        static inline void call(std::ostream& os, const V& v)
        { InstWrite(os,v.calc().xView()); }
    };

    template <class V>
    inline void Write(std::ostream& os, const BaseVector<V>& v)
    {
        typedef typename V::value_type T;
        const bool inst = 
            V::unknownsizes &&
            Traits<T>::isinst;
        CallWrite<inst,V>::call(os,v.vec());
    }

    // With thresh:
    template <class V>
    inline void InlineWrite(
        std::ostream& os,
        const BaseVector<V>& v, typename V::real_type thresh) 
    {
        typedef typename V::value_type T;
        const int n=v.size();
        os << n << " (";
        for(int i=0;i<n;++i) {
            T temp = v.cref(i);
            os << " " << Value((TMV_ABS(temp) < thresh ? T(0) : temp)) << " ";
        }
        os << ")";
    }

    // Defined in TMV_Vector.cpp
    template <class T, bool C>
    void InstWrite(
        std::ostream& os, const ConstVectorView<T,UNKNOWN,C>& v,
        typename Traits<T>::real_type thresh);

    template <bool inst, class V>
    struct CallWriteThresh // inst = false
    {
        static inline void call(
            std::ostream& os,
            const V& v, typename V::real_type thresh) 
        { InlineWrite(os,v.calc(),thresh); }
    };
    template <class V>
    struct CallWriteThresh<true,V>
    {
        static inline void call(
            std::ostream& os,
            const V& v, typename V::real_type thresh) 
        { InstWrite(os,v.calc().xView(),thresh); }
    };

    template <class V>
    inline void Write(
        std::ostream& os,
        const BaseVector<V>& v, typename V::real_type thresh) 
    {
        typedef typename V::value_type T;
        const bool inst = 
            V::unknownsizes &&
            Traits<T>::isinst;
        CallWriteThresh<inst,V>::call(os,v.vec(),thresh);
    }


    //
    // Read Vector
    //

#ifndef TMV_NO_THROW
    template <class V> 
    class VectorReadError : public ReadError
    {
    public :
        typedef typename V::copy_type copy_type;
        int i;
        mutable auto_ptr<copy_type> v;
        char exp,got;
        size_t s;
        bool is, iseof, isbad;

        inline VectorReadError(int _i, const BaseVector<V>& _v,
                               std::istream& _is) throw() :
            ReadError("Vector"),
            i(_i), v(new copy_type(_v)), exp(0), got(0), s(_v.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        inline VectorReadError(std::istream& _is) throw() :
            ReadError("Vector"),
            i(0), v(0), exp(0), got(0), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        inline VectorReadError(
            int _i, const BaseVector<V>& _v, std::istream& _is,
            char _e, char _g
        ) throw() :
            ReadError("Vector"),
            i(_i), v(new copy_type(_v)), exp(_e), got(_g), s(_v.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        inline VectorReadError(const BaseVector<V>& _v, std::istream& _is, 
                               size_t _s) throw() :
            ReadError("Vector"),
            i(0), v(new copy_type(_v)), exp(0), got(0), s(_s),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        inline VectorReadError(const VectorReadError<V>& rhs) throw() :
            ReadError("Vector"),
            i(rhs.i), v(rhs.v), exp(rhs.exp), got(rhs.got), s(rhs.s),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        inline ~VectorReadError() throw() {}

        inline void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for Vector\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (v.get() && s != v->size()) {
                os<<"Wrong size: expected "<<v->size()<<", got "<<s<<".\n";
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
            if (v.get()) {
                os<<"The portion of the Vector which was successfully "
                    "read is: \n(";
                for(int ii=0;ii<i;++ii)
                    os<<' '<<(*v)(ii)<<' ';
                os<<")\n";
            }
        }
    };
#endif

    template <class V>
    inline void InlineRead(std::istream& is, BaseVector_Mutable<V>& v)
    {
        char paren;
        is >> paren;
        if (!is || paren != '(') {
#ifdef TMV_NO_THROW
            std::cerr<<"Vector ReadError: "<<paren<<" != (\n";
            exit(1); 
#else
            throw VectorReadError<V>(0,v,is,'(',is?paren:'(');
#endif
        }
        const int n = v.size();
        typename V::value_type temp;
        for(int i=0;i<n;++i) {
            is >> temp;
            if (!is) {
#ifdef TMV_NO_THROW
                std::cerr<<"Vector ReadError: !is \n"; 
                exit(1); 
#else
                throw VectorReadError<V>(n-i,v,is);
#endif
            }
            v.ref(i) = temp;
        }
        is >> paren;
        if (!is || paren != ')') {
#ifdef TMV_NO_THROW
            std::cerr<<"Vector ReadError: "<<paren<<" != )\n"; 
            exit(1); 
#else
            throw VectorReadError<V>(n,v,is,')',is?paren:')');
#endif
        }
    }

    // Defined in TMV_Vector.cpp
    template <class T, bool C>
    void InstRead(std::istream& is, VectorView<T,UNKNOWN,C> v);

    template <bool inst, class V>
    struct CallRead // inst = false
    {
        static inline void call(std::istream& is, V& v)
        { InlineRead(is,v); }
    };
    template <class V>
    struct CallRead<true,V>
    {
        static inline void call(std::istream& is, V& v)
        { InstRead(is,v.xView()); }
    };

    template <class V>
    inline void Read(std::istream& is, BaseVector_Mutable<V>& v)
    {
        typedef typename V::value_type T;
        const bool inst = 
            V::unknownsizes &&
            Traits<T>::isinst;
        CallRead<inst,V>::call(is,v.vec());
    }


    // 
    // Operator overloads for I/O
    // is >> v
    // os << v
    //

    template <class V>
    inline std::ostream& operator<<(
        std::ostream& os, const BaseVector<V>& v)
    { Write(os,v); return os; }

    template <class V>
    inline std::istream& operator>>(
        std::istream& is, BaseVector_Mutable<V>& v)
    {
        size_t n;
        is >> n;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"Vector ReadError: !is \n"; 
            exit(1); 
#else
            throw VectorReadError<V>(is);
#endif
        }
        if (n != v.size()) {
#ifdef TMV_NO_THROW
            std::cerr<<"Vector ReadError: Wrong size \n"; 
            exit(1); 
#else
            throw VectorReadError<V>(v,is,n);
#endif
        }
        v.read(is);
        return is;
    }

    template <class T, IndexStyle I>
    inline std::istream& operator>>(
        std::istream& is, auto_ptr<Vector<T,I> >& v)
    {
        size_t n;
        is >> n;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"Vector ReadError: !is \n"; 
            exit(1); 
#else
            throw VectorReadError<Vector<T,I> >(is);
#endif
        }
        v.reset(new Vector<T,I>(n));
        v->read(is);
        return is;
    }
} // namespace mv

#endif
