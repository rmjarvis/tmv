
#ifndef TMV_VectorIO_H
#define TMV_VectorIO_H

#include "TMV_BaseVector.h"
#include "TMV_IOStyle.h"

namespace tmv {

    //
    // Write Vector
    //

    // Defined in TMV_Vector.cpp
    template <class T, int C>
    void InstWrite(const TMV_Writer& writer, const ConstVectorView<T,C>& v);
    template <class T>
    void InstRead(const TMV_Reader& reader, VectorView<T> v);

    template <int algo, class V>
    struct WriteV_Helper;

    template <class V>
    struct WriteV_Helper<11,V>
    {
        static void call(const TMV_Writer& writer, const V& v)
        {
            const ptrdiff_t N = v.size();
            writer.begin();
            writer.writeCode("V");
            writer.writeSize(N);
            writer.writeLParen();
            for(ptrdiff_t i=0;i<N;++i) {
                if (i > 0) writer.writeSpace();
                writer.writeValue(v.cref(i));
            }
            writer.writeRParen();
            writer.end();
        }
    };

    // algo 90: Call inst
    template <class V>
    struct WriteV_Helper<90,V>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const V& v)
        { InstWrite(writer,v.calc().xView()); }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class V>
    struct WriteV_Helper<-3,V>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const V& v)
        { WriteV_Helper<11,V>::call(writer,v); }
    };
             
    // algo -2: Check for inst
    template <class V>
    struct WriteV_Helper<-2,V>
    {
        typedef typename V::value_type T;
        enum { inst = (
                (V::_size == Unknown || V::_size > 16) &&
                Traits<T>::isinst ) };
        enum { algo = (
                inst ? 90 :
                -3 ) };
        static TMV_INLINE void call(const TMV_Writer& writer, const V& v)
        { WriteV_Helper<algo,V>::call(writer,v); }
    };

    template <class V>
    struct WriteV_Helper<-1,V>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const V& v)
        { WriteV_Helper<-2,V>::call(writer,v); }
    };

    template <class V>
    inline void Write(const TMV_Writer& writer, const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        WriteV_Helper<-2,Vv>::call(writer,vv);
    }

    template <class V>
    inline void InlineWrite(
        const TMV_Writer& writer, const BaseVector_Calc<V>& v)
    {
        typedef typename V::const_cview_type Vv;
        TMV_MAYBE_CREF(V,Vv) vv = v.cView();
        WriteV_Helper<-3,Vv>::call(writer,vv);
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
        ptrdiff_t i;
        std::string exp,got;
        ptrdiff_t s;
        bool is, iseof, isbad;

        VectorReadError(std::istream& _is) throw() :
            ReadError("Vector"),
            i(0), s(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        VectorReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("Vector"),
            i(0), exp(_e), got(_g), s(0), 
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        template <class V>
        VectorReadError(
            ptrdiff_t _i, const BaseVector<V>& _v, std::istream& _is) throw() :
            ReadError("Vector"),
            v(_v), i(_i), s(_v.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class V>
        VectorReadError(
            ptrdiff_t _i, const BaseVector<V>& _v, std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("Vector"),
            v(_v), i(_i), exp(_e), got(_g), s(_v.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class V>
        VectorReadError(
            const BaseVector<V>& _v, std::istream& _is, ptrdiff_t _s) throw() :
            ReadError("Vector"),
            v(_v), i(0), s(_s),
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
                for(ptrdiff_t ii=0;ii<i;++ii) os<<' '<<v(ii)<<' ';
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
        static void call(const TMV_Reader& reader, V& v)
        {
            const ptrdiff_t n = v.size();
            std::string exp, got;
            typedef typename V::value_type T;
            T temp;
            if (!reader.readLParen(exp,got)) {
#ifdef TMV_NO_THROW
                std::cerr<<"Vector ReadError: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw VectorReadError<T>(0,v,reader.getis(),exp,got);
#endif
            }
            for(ptrdiff_t i=0;i<n;++i) {
                if (i>0) {
                    if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"Vector Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw VectorReadError<T>(i,v,reader.getis(),exp,got);
#endif
                    }
                }
                if (!reader.readValue(temp)) {
#ifdef NOTHROW
                    std::cerr<<"Vector Read Error: reading value\n";
                    exit(1);
#else
                    throw VectorReadError<T>(i,v,reader.getis());
#endif
                }
                v.ref(i) = temp;
            }
            if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"Vector Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw VectorReadError<T>(n,v,reader.getis(),exp,got);
#endif
            }
        }
    };

    // algo 90: Call inst
    template <class V>
    struct ReadV_Helper<90,V>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, V& v)
        { InstRead(reader,v.xView()); }
    };
             
    // algo 97: Conjugate
    template <class V>
    struct ReadV_Helper<97,V>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, V& v)
        {
            typedef typename V::conjugate_type Vc;
            Vc vc = v.conjugate();
            ReadV_Helper<-2,Vc>::call(reader,vc); 
            vc.conjugateSelf();
        }
    };
             
    // algo -3: Only one algorithv, so call it.
    template <class V>
    struct ReadV_Helper<-3,V>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, V& v)
        { ReadV_Helper<11,V>::call(reader,v); }
    };
             
    // algo -2: Check for inst
    template <class V>
    struct ReadV_Helper<-2,V>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, V& v)
        {
            typedef typename V::value_type T;
            const bool inst = 
                (V::_size == Unknown || V::_size > 16) &&
                Traits<T>::isinst;
            const int algo = 
                V::_conj ? 97 :
                inst ? 90 :
                -3;
            ReadV_Helper<algo,V>::call(reader,v); 
        }
    };

    template <class V>
    struct ReadV_Helper<-1,V>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, V& v)
        { ReadV_Helper<-2,V>::call(reader,v); }
    };

    template <class V>
    inline void Read(const TMV_Reader& reader, BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ReadV_Helper<-2,Vv>::call(reader,vv);
    }

    template <class V>
    inline void InlineRead(const TMV_Reader& reader, BaseVector_Mutable<V>& v)
    {
        typedef typename V::cview_type Vv;
        TMV_MAYBE_REF(V,Vv) vv = v.cView();
        ReadV_Helper<-3,Vv>::call(reader,vv);
    }


    // 
    // Operator overloads for I/O
    // is >> v
    // os << v
    //

    template <class V>
    inline std::ostream& operator<<(
        const TMV_Writer& writer, const BaseVector<V>& v)
    { Write(writer,v.calc()); return writer.getos(); }

    template <class V>
    static std::istream& operator>>(
        const TMV_Reader& reader, BaseVector_Mutable<V>& v)
    {
        typedef typename V::value_type T;
        std::string exp,got;
        if (!reader.readCode("V",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw VectorReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t n=v.size();
        if (!reader.readSize(n,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: reading size\n";
            exit(1);
#else
            throw VectorReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (n != v.size()) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: wrong size\n";
            exit(1);
#else
            throw VectorReadError<T>(v,reader.getis(),n);
#endif
        }
        Read(reader,v);
        return reader.getis();
    }

    template <class T, int A>
    static std::istream& operator>>(const TMV_Reader& reader, Vector<T,A>& v)
    {
        std::string exp,got;
        if (!reader.readCode("V",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw VectorReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t n=v.size();
        if (!reader.readSize(n,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Vector Read Error: reading size\n";
            exit(1);
#else
            throw VectorReadError<T>(reader.getis(),exp,got);
#endif
        }
        if (n != v.size()) v.resize(n);
        Read(reader,v);
        return reader.getis();
    }

    template <class V>
    inline std::ostream& operator<<(std::ostream& os, const BaseVector<V>& v)
    { return os << IOStyle() << v.vec(); }

    template <class V>
    inline std::istream& operator>>(std::istream& is, BaseVector_Mutable<V>& v)
    { return is >> IOStyle() >> v.vec(); }

    template <class T, int A>
    inline std::istream& operator>>(std::istream& is, Vector<T,A>& v)
    { return is >> IOStyle() >> v; }

    template <class T, int A>
    inline std::istream& operator>>(
        const TMV_Reader& reader, VectorView<T,A> v)
    { 
        return reader >> 
            static_cast<BaseVector_Mutable<VectorView<T,A> >&>(v); 
    }

    template <class T, ptrdiff_t N, ptrdiff_t S, int A>
    inline std::istream& operator>>(
        const TMV_Reader& reader, SmallVectorView<T,N,S,A> v)
    {
        return reader >> 
            static_cast<BaseVector_Mutable<SmallVectorView<T,N,S,A> >&>(v); 
    }

    template <class T, int A>
    inline std::istream& operator>>(std::istream& is, VectorView<T,A> v)
    {
        return is >> 
            static_cast<BaseVector_Mutable<VectorView<T,A> >&>(v); 
    }

    template <class T, ptrdiff_t N, ptrdiff_t S, int A>
    inline std::istream& operator>>(
        std::istream& is, SmallVectorView<T,N,S,A> v)
    { 
        return is >> 
            static_cast<BaseVector_Mutable<SmallVectorView<T,N,S,A> >&>(v); 
    }


} // namespace mv

#endif
