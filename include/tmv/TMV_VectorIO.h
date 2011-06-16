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

#include "TMV_BaseVector.h"

namespace tmv {

    //
    // Write Vector
    //

    // Defined in TMV_Vector.cpp
    template <class T, int C>
    void InstWrite(
        std::ostream& os, const ConstVectorView<T,C>& v);
    template <class T, int C>
    void InstWrite(
        std::ostream& os, const ConstVectorView<T,C>& v,
        typename ConstVectorView<T>::float_type thresh);
    template <class T>
    void InstRead(std::istream& is, VectorView<T> v);

    template <int algo, class V>
    struct WriteV_Helper;

    template <class V>
    struct WriteV_Helper<11,V>
    {
        static void call(std::ostream& os, const V& v)
        {
            const int n=v.size();
            os << n << " (";
            for(int i=0;i<n;++i) os << " " << Value(v.cref(i)) << " ";
            os << ")";
        }
        static void call(
            std::ostream& os, const V& v, typename V::float_type thresh)
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
    };

    // algo 90: Call inst
    template <class V>
    struct WriteV_Helper<90,V>
    {
        static TMV_INLINE void call(std::ostream& os, const V& v)
        { InstWrite(os,v.calc().xView()); }
        static TMV_INLINE void call(
            std::ostream& os, const V& v, typename V::float_type thresh)
        { InstWrite(os,v.calc().xView(),thresh); }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class V>
    struct WriteV_Helper<-3,V>
    {
        static TMV_INLINE void call(std::ostream& os, const V& v)
        { WriteV_Helper<11,V>::call(os,v); }
        static TMV_INLINE void call(
            std::ostream& os, const V& v, typename V::float_type thresh)
        { WriteV_Helper<11,V>::call(os,v,thresh); }
    };
             
    // algo -2: Check for inst
    template <class V>
    struct WriteV_Helper<-2,V>
    {
        typedef typename V::value_type T;
        enum { inst = (
                (V::_size == UNKNOWN || V::_size > 16) &&
                Traits<T>::isinst ) };
        enum { algo = (
                inst ? 90 :
                -3 ) };
        static TMV_INLINE void call(std::ostream& os, const V& v)
        { WriteV_Helper<algo,V>::call(os,v); }
        static TMV_INLINE void call(
            std::ostream& os, const V& v, typename V::float_type thresh)
        { WriteV_Helper<algo,V>::call(os,v,thresh); }
    };

    template <class V>
    struct WriteV_Helper<-1,V>
    {
        static TMV_INLINE void call(std::ostream& os, const V& v)
        { WriteV_Helper<-2,V>::call(os,v); }
        static TMV_INLINE void call(
            std::ostream& os, const V& v, typename V::float_type thresh)
        { WriteV_Helper<-2,V>::call(os,v,thresh); }
    };

    template <class V>
    static inline void Write(std::ostream& os, const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        WriteV_Helper<-2,Vv>::call(os,vv);
    }

    template <class V>
    static inline void InlineWrite(
        std::ostream& os, const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        WriteV_Helper<-3,Vv>::call(os,vv);
    }

    template <class V>
    static inline void Write(
        std::ostream& os,
        const BaseVector_Calc<V>& v, typename V::float_type thresh) 
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        WriteV_Helper<-2,Vv>::call(os,vv,thresh);
    }

    template <class V>
    static inline void InlineWrite(
        std::ostream& os,
        const BaseVector_Calc<V>& v, typename V::float_type thresh) 
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        WriteV_Helper<-3,Vv>::call(os,vv,thresh);
    }


    //
    // Read Vector
    //

#ifndef TMV_NO_THROW
    template <class T>
    class VectorReadError : 
        public ReadError
    {
    public :
        Vector<T> v;
        int i;
        char exp,got;
        size_t s;
        bool is, iseof, isbad;

        VectorReadError(std::istream& _is) throw() :
            ReadError("Vector"),
            i(0), exp(0), got(0), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class V>
        VectorReadError(
            int _i, const BaseVector<V>& _v, std::istream& _is) throw() :
            ReadError("Vector"),
            v(_v), i(_i), exp(0), got(0), s(_v.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class V>
        VectorReadError(
            int _i, const BaseVector<V>& _v, std::istream& _is,
            char _e, char _g) throw() :
            ReadError("Vector"),
            v(_v), i(_i), exp(_e), got(_g), s(_v.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class V>
        VectorReadError(
            const BaseVector<V>& _v, std::istream& _is, size_t _s) throw() :
            ReadError("Vector"),
            v(_v), i(0), exp(0), got(0), s(_s),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        VectorReadError(const VectorReadError<T>& rhs) throw() :
            ReadError("Vector"),
            v(rhs.v), i(rhs.i), exp(rhs.exp), got(rhs.got), s(rhs.s),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~VectorReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for Vector\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (s != v.size()) {
                os<<"Wrong size: expected "<<v.size()<<", got "<<s<<".\n";
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
            if (v.size() > 0) {
                os<<"The portion of the Vector which was successfully "
                    "read is: \n(";
                for(int ii=0;ii<i;++ii)
                    os<<' '<<v(ii)<<' ';
                os<<")\n";
            }
        }
    };
#endif

    template <int algo, class V>
    struct ReadV_Helper;

    template <class V>
    struct ReadV_Helper<11,V>
    {
        static void call(std::istream& is, V& v)
        {
            typedef typename V::value_type T;
            char paren;
            is >> paren;
            if (!is || paren != '(') {
#ifdef TMV_NO_THROW
                std::cerr<<"Vector ReadError: "<<paren<<" != (\n";
                exit(1);
#else
                throw VectorReadError<T>(0,v,is,'(',is?paren:'(');
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
                    throw VectorReadError<T>(n-i,v,is);
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
                throw VectorReadError<T>(n,v,is,')',is?paren:')');
#endif
            }
        }
    };

    // algo 90: Call inst
    template <class V>
    struct ReadV_Helper<90,V>
    {
        static TMV_INLINE void call(std::istream& is, V& v)
        { InstRead(is,v.xView()); }
    };
             
    // algo 97: Conjugate
    template <class V>
    struct ReadV_Helper<97,V>
    {
        static TMV_INLINE void call(std::istream& is, V& v)
        {
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            ReadV_Helper<-2,Vc>::call(is,vc); 
            vc.conjugateSelf();
        }
    };
             
    // algo -3: Only one algorithv, so call it.
    template <class V>
    struct ReadV_Helper<-3,V>
    {
        static TMV_INLINE void call(std::istream& is, V& v)
        { ReadV_Helper<11,V>::call(is,v); }
    };
             
    // algo -2: Check for inst
    template <class V>
    struct ReadV_Helper<-2,V>
    {
        static TMV_INLINE void call(std::istream& is, V& v)
        {
            typedef typename V::value_type T;
            const int inst = 
                (V::_size == UNKNOWN || V::_size > 16) &&
                Traits<T>::isinst;
            const int algo = 
                V::_conj ? 97 :
                inst ? 90 :
                -3;
            ReadV_Helper<algo,V>::call(is,v); 
        }
    };

    template <class V>
    struct ReadV_Helper<-1,V>
    {
        static TMV_INLINE void call(std::istream& is, V& v)
        { ReadV_Helper<-2,V>::call(is,v); }
    };

    template <class V>
    static inline void Read(std::istream& is, BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ReadV_Helper<-2,Vv>::call(is,vv);
    }

    template <class V>
    static inline void InlineRead(
        std::istream& is, BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ReadV_Helper<-3,Vv>::call(is,vv);
    }


    // 
    // Operator overloads for I/O
    // is >> v
    // os << v
    //

    template <class V>
    static inline std::ostream& operator<<(
        std::ostream& os, const BaseVector<V>& v)
    { Write(os,v.calc()); return os; }

    template <class V>
    static std::istream& operator>>(
        std::istream& is, BaseVector_Mutable<V>& v)
    {
        typedef typename V::value_type T;
        size_t n;
        is >> n;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"Vector ReadError: !is \n";
            exit(1);
#else
            throw VectorReadError<T>(is);
#endif
        }
        if (n != v.size()) {
#ifdef TMV_NO_THROW
            std::cerr<<"Vector ReadError: Wrong size \n";
            exit(1);
#else
            throw VectorReadError<T>(v,is,n);
#endif
        }
        v.read(is);
        return is;
    }

    template <class T, int A>
    static std::istream& operator>>(
        std::istream& is, auto_ptr<Vector<T,A> >& v)
    {
        size_t n;
        is >> n;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"Vector ReadError: !is \n";
            exit(1);
#else
            throw VectorReadError<T>(is);
#endif
        }
        v.reset(new Vector<T,A>(n));
        v->read(is);
        return is;
    }

} // namespace mv

#endif
