
#ifndef TMV_TriMatrixIO_H
#define TMV_TriMatrixIO_H

#include "TMV_BaseMatrix_Tri.h"
#include "TMV_IOStyle.h"

namespace tmv {

    //
    // Write TriMatrix
    //

    // Defined in TMV_TriMatrix.cpp
    template <class T, int C>
    void InstWrite(
        const TMV_Writer& writer, const ConstUpperTriMatrixView<T,C>& m);
    template <class T, int C>
    void InstWrite(
        const TMV_Writer& writer, const ConstLowerTriMatrixView<T,C>& m);
    template <class T>
    void InstRead(const TMV_Reader& reader, UpperTriMatrixView<T> m);
    template <class T>
    void InstRead(const TMV_Reader& reader, LowerTriMatrixView<T> m);

    template <int algo, class M>
    struct WriteU_Helper;

    // algo 11: UpperTriMatrix
    template <class M>
    struct WriteU_Helper<11,M>
    {
        static void call(const TMV_Writer& writer, const M& m)
        {
            typedef typename M::value_type T;
            const ptrdiff_t N = m.size();
            writer.begin();
            writer.writeCode("U");
            writer.writeSize(N);
            writer.writeSimpleSize(N);
            writer.writeStart();
            for(ptrdiff_t i=0;i<N;++i) {
                writer.writeLParen();
                if (!writer.isCompact()) {
                    for(ptrdiff_t j=0;j<i;++j) {
                        writer.writeValue(T(0));
                        writer.writeSpace();
                    }
                }
                for(ptrdiff_t j=i;j<N;++j) {
                    if (j > i) writer.writeSpace();
                    writer.writeValue(m.cref(i,j));
                }
                writer.writeRParen();
                if (i < N-1) writer.writeRowEnd();
            }
            writer.writeFinal();
            writer.end();
        }
    };

    // algo 12: LowerTriMatrix
    template <class M>
    struct WriteU_Helper<12,M>
    {
        static void call(const TMV_Writer& writer, const M& m)
        {
            typedef typename M::value_type T;
            const ptrdiff_t N = m.size();
            writer.begin();
            writer.writeCode("L");
            writer.writeSize(N);
            writer.writeSimpleSize(N);
            writer.writeStart();
            for(ptrdiff_t i=0;i<N;++i) {
                writer.writeLParen();
                for(ptrdiff_t j=0;j<i+1;++j) {
                    if (j > 0) writer.writeSpace();
                    writer.writeValue(m.cref(i,j));
                }
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
    };

    // algo 90: Call inst
    template <class M>
    struct WriteU_Helper<90,M>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M& m)
        { InstWrite(writer,m.calc().xView()); }
    };
             
    // algo -3: Determine which algorithm to use
    template <class M>
    struct WriteU_Helper<-3,M>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M& m)
        { 
            const int algo = M::_upper ? 11 : 12;
            WriteU_Helper<algo,M>::call(writer,m); 
        }
    };
             
    // algo -2: Check for inst
    template <class M>
    struct WriteU_Helper<-2,M>
    {
        typedef typename M::value_type T;
        enum { inst = (
                (M::_colsize == Unknown || M::_colsize > 16) &&
                (M::_rowsize == Unknown || M::_rowsize > 16) &&
                Traits<T>::isinst ) };
        enum { algo = (
                inst ? 90 :
                -3 ) };
        static TMV_INLINE void call(const TMV_Writer& writer, const M& m)
        { WriteU_Helper<algo,M>::call(writer,m); }
    };

    template <class M>
    struct WriteU_Helper<-1,M>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M& m)
        { WriteU_Helper<-2,M>::call(writer,m); }
    };

    template <class M>
    inline void Write(const TMV_Writer& writer, const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteU_Helper<-2,Mv>::call(writer,mv);
    }

    template <class M>
    inline void InlineWrite(
        const TMV_Writer& writer, const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteU_Helper<-3,Mv>::call(writer,mv);
    }


    //
    // Read TriMatrix
    //

#ifndef TMV_NO_THROW
    template <bool upper, class T>
    class TriMatrixReadError : 
        public ReadError
    {
    public :
        typedef typename 
            TypeSelect<upper,UpperTriMatrix<T>,LowerTriMatrix<T> >::type M;
#define TAG (upper ? "UpperTriMatrix" : "LowerTriMatrix")

        M m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t s;
        T v1;
        bool is, iseof, isbad;

        TriMatrixReadError(std::istream& _is) throw() :
            ReadError(TAG),
            i(0), j(0), s(0), v1(1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        TriMatrixReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError(TAG),
            i(0), j(0), exp(_e), got(_g), s(0), v1(1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        template <class M>
        TriMatrixReadError(
            const BaseMatrix_Tri<M>& _m, std::istream& _is, ptrdiff_t _s) throw() :
            ReadError(TAG),
            m(_m), i(0), j(0), s(_s), v1(1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        TriMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const BaseMatrix_Tri<M>& _m,
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError(TAG),
            m(_m), i(_i), j(_j), exp(_e), got(_g), 
            s(_m.size()), v1(i==j?T(1):T(0)),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        TriMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const BaseMatrix_Tri<M>& _m, 
            std::istream& _is) throw() :
            ReadError(TAG),
            m(_m), i(_i), j(_j), s(_m.size()), v1(i==j?T(1):T(0)),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        TriMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const BaseMatrix_Tri<M>& _m, 
            std::istream& _is, T _v1) throw() :
            ReadError(TAG),
            m(_m), i(_i), j(_j), s(_m.size()), v1(_v1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        TriMatrixReadError(const TriMatrixReadError<upper,T>& rhs) throw() :
            ReadError(TAG),
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got), 
            s(rhs.s), v1(rhs.v1),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~TriMatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for "<<TAG<<"\n";
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
            if (i != j && v1 != T(0)) {
                os<<"Invalid input: Expected 0, got "<<v1<<".\n";
            }
            if (i == j && v1 != T(1)) {
                os<<"Invalid input: Expected 1, got "<<v1<<".\n";
            }
            if (m.size() > 0) {
                os<<"The portion of the "<<TAG<<" which was successfully "
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
#undef TAG

    };
#endif

    template <int algo, class M>
    struct ReadU_Helper;

    // algo 11: UpperTriMatrix
    template <class M>
    struct ReadU_Helper<11,M>
    {
        static void call(const TMV_Reader& reader, M& m)
        {
            typedef typename M::value_type T;
            const ptrdiff_t N = m.size();
            std::string exp, got;
            T temp;
            if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw TriMatrixReadError<true,T>(0,0,m,reader.getis(),exp,got);
#endif
            }
            for(ptrdiff_t i=0;i<N;++i) {
                if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw TriMatrixReadError<true,T>(i,0,m,reader.getis(),exp,got);
#endif
                }
                if (!reader.isCompact()) {
                    for(ptrdiff_t j=0;j<i;++j) {
                        if (!reader.readValue(temp)) {
#ifdef NOTHROW
                            std::cerr<<"UpperTriMatrix Read Error: reading value\n";
                            exit(1);
#else
                            throw TriMatrixReadError<true,T>(i,j,m,reader.getis());
#endif
                        }
                        if (temp != T(0)) {
#ifdef NOTHROW
                            std::cerr<<"UpperTriMatrix Read Error: "<<temp<<" != 0\n";
                            exit(1);
#else
                            throw TriMatrixReadError<true,T>(i,j,m,reader.getis(),temp);
#endif
                        }
                        if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                            std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                            exit(1);
#else
                            throw TriMatrixReadError<true,T>(i,j,m,reader.getis(),exp,got);
#endif
                        }
                    }
                }
                for(ptrdiff_t j=i;j<N;++j) {
                    if (j>i && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw TriMatrixReadError<true,T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"UpperTriMatrix Read Error: reading value\n";
                        exit(1);
#else
                        throw TriMatrixReadError<true,T>(i,j,m,reader.getis());
#endif
                    }
                    if (j==i && m.isunit()) {
                        if (temp != T(1)) {
#ifdef NOTHROW
                            std::cerr<<"UpperTriMatrix Read Error: "<<temp<<" != 1\n";
                            exit(1);
#else
                            throw TriMatrixReadError<true,T>(i,j,m,reader.getis(),temp);
#endif
                        }
                    } else {
                        m.ref(i,j) = temp;
                    }
                }
                if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw TriMatrixReadError<true,T>(i,N,m,reader.getis(),exp,got);
#endif
                }
                if (i < N-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw TriMatrixReadError<true,T>(i,N,m,reader.getis(),exp,got);
#endif
                }
            }
            if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw TriMatrixReadError<true,T>(N,0,m,reader.getis(),exp,got);
#endif
            }
        }
    };

    // algo 12: LowerTriMatrix
    template <class M>
    struct ReadU_Helper<12,M>
    {
        static void call(const TMV_Reader& reader, M& m)
        {
            typedef typename M::value_type T;
            const ptrdiff_t N = m.size();
            std::string exp, got;
            T temp;
            if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw TriMatrixReadError<false,T>(0,0,m,reader.getis(),exp,got);
#endif
            }
            for(ptrdiff_t i=0;i<N;++i) {
                if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw TriMatrixReadError<false,T>(i,0,m,reader.getis(),exp,got);
#endif
                }
                for(ptrdiff_t j=0;j<i+1;++j) {
                    if (j>0 && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw TriMatrixReadError<false,T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"LowerTriMatrix Read Error: reading value\n";
                        exit(1);
#else
                        throw TriMatrixReadError<false,T>(i,j,m,reader.getis());
#endif
                    }
                    if (j==i && m.isunit()) {
                        if (temp != T(1)) {
#ifdef NOTHROW
                            std::cerr<<"LowerTriMatrix Read Error: "<<temp<<" != 0\n";
                            exit(1);
#else
                            throw TriMatrixReadError<false,T>(i,j,m,reader.getis(),temp);
#endif
                        }
                    } else {
                        m.ref(i,j) = temp;
                    }
                }
                if (!reader.isCompact()) {
                    for(ptrdiff_t j=i+1;j<N;++j) {
                        if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                            std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                            exit(1);
#else
                            throw TriMatrixReadError<false,T>(i,j,m,reader.getis(),exp,got);
#endif
                        }
                        if (!reader.readValue(temp)) {
#ifdef NOTHROW
                            std::cerr<<"LowerTriMatrix Read Error: reading value\n";
                            exit(1);
#else
                            throw TriMatrixReadError<false,T>(i,j,m,reader.getis());
#endif
                        }
                        if (temp != T(0)) {
#ifdef NOTHROW
                            std::cerr<<"LowerTriMatrix Read Error: "<<temp<<" != 0\n";
                            exit(1);
#else
                            throw TriMatrixReadError<false,T>(i,j,m,reader.getis(),temp);
#endif
                        }
                    }
                }
                if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw TriMatrixReadError<false,T>(i,N,m,reader.getis(),exp,got);
#endif
                }
                if (i < N-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw TriMatrixReadError<false,T>(i,N,m,reader.getis(),exp,got);
#endif
                }
            }
            if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"LowerTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw TriMatrixReadError<false,T>(N,0,m,reader.getis(),exp,got);
#endif
            }
        }
    };

    // algo 90: Call inst
    template <class M>
    struct ReadU_Helper<90,M>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M& m)
        { InstRead(reader,m.xView()); }
    };
             
    // algo 97: Conjugate
    template <class M>
    struct ReadU_Helper<97,M>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M& m)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            ReadU_Helper<-2,Mc>::call(reader,mc); 
            mc.conjugateSelf();
        }
    };
             
    // algo -3: Determine which algorithm to use
    template <class M>
    struct ReadU_Helper<-3,M>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M& m)
        {
            const int algo = M::_upper ? 11 : 12;
            ReadU_Helper<algo,M>::call(reader,m); 
        }
    };
             
    // algo -2: Check for inst
    template <class M>
    struct ReadU_Helper<-2,M>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M& m)
        {
            typedef typename M::value_type T;
            const bool inst = 
                (M::_colsize == Unknown || M::_colsize > 16) &&
                (M::_rowsize == Unknown || M::_rowsize > 16) &&
                Traits<T>::isinst;
            const int algo = 
                M::_conj ? 97 :
                inst ? 90 :
                -3;
            ReadU_Helper<algo,M>::call(reader,m); 
        }
    };

    template <class M>
    struct ReadU_Helper<-1,M>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M& m)
        { ReadU_Helper<-2,M>::call(reader,m); }
    };

    template <class M>
    inline void Read(const TMV_Reader& reader, BaseMatrix_Tri_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadU_Helper<-2,Mv>::call(reader,mv);
    }

    template <class M>
    inline void InlineRead(
        const TMV_Reader& reader, BaseMatrix_Tri_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadU_Helper<-3,Mv>::call(reader,mv);
    }


    // 
    // Operator overloads for I/O
    // is >> m
    // os << m
    //

    template <class M>
    static std::istream& operator>>(
        const TMV_Reader& reader, BaseMatrix_Tri_Mutable<M>& m)
    {
        typedef typename M::value_type T;
        std::string exp,got;
        const bool upper = M::_upper;
        const char* code = upper ? "U" : "L";
#ifdef NOTHROW
        const char* tag = upper ? "UpperTriMatrix" : "LowerTriMatrix";
#endif
        if (!reader.readCode(code,exp,got)) {
#ifdef NOTHROW
            std::cerr<<tag<<" Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw TriMatrixReadError<upper,T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=m.size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<tag<<" Read Error: reading size\n";
            exit(1);
#else
            throw TriMatrixReadError<upper,T>(reader.getis(),exp,got);
#endif
        }
        if (s != m.size()) {
#ifdef NOTHROW
            std::cerr<<tag<<" Read Error: wrong size\n";
            exit(1);
#else
            throw TriMatrixReadError<upper,T>(m,reader.getis(),s);
#endif
        }
        s=m.size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<tag<<" Read Error: reading size\n";
            exit(1);
#else
            throw TriMatrixReadError<upper,T>(reader.getis(),exp,got);
#endif
        }
        if (s != m.size()) {
#ifdef NOTHROW
            std::cerr<<tag<<" Read Error: Wrong size\n";
            exit(1);
#else
            throw TriMatrixReadError<upper,T>(m,reader.getis(),s);
#endif
        }
        Read(reader,m);
        return reader.getis();
    }

    template <class T, int A>
    static std::istream& operator>>(
        const TMV_Reader& reader, UpperTriMatrix<T,A>& m)
    {
        std::string exp,got;
        if (!reader.readCode("U",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw TriMatrixReadError<true,T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=m.size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw TriMatrixReadError<true,T>(reader.getis(),exp,got);
#endif
        }
        if (s != m.size()) m.resize(s);
        s=m.size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw TriMatrixReadError<true,T>(reader.getis(),exp,got);
#endif
        }
        if (s != m.size()) {
#ifdef NOTHROW
            std::cerr<<"UpperTriMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw TriMatrixReadError<true,T>(m,reader.getis(),s);
#endif
        }
        Read(reader,m);
        return reader.getis();
    }

    template <class T, int A>
    static std::istream& operator>>(
        const TMV_Reader& reader, LowerTriMatrix<T,A>& m)
    {
        std::string exp,got;
        if (!reader.readCode("L",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"TriMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw TriMatrixReadError<false,T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=m.size();
        if (!reader.readSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"TriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw TriMatrixReadError<false,T>(reader.getis(),exp,got);
#endif
        }
        if (s != m.size()) m.resize(s);
        s=m.size();
        if (!reader.readSimpleSize(s,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"TriMatrix Read Error: reading size\n";
            exit(1);
#else
            throw TriMatrixReadError<false,T>(reader.getis(),exp,got);
#endif
        }
        if (s != m.size()) {
#ifdef NOTHROW
            std::cerr<<"TriMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw TriMatrixReadError<false,T>(m,reader.getis(),s);
#endif
        }
        Read(reader,m);
        return reader.getis();
    }

    template <class T, int A>
    std::istream& operator>>(std::istream& is, UpperTriMatrix<T,A>& m)
    { return is >> IOStyle() >> m; }

    template <class T, int A>
    std::istream& operator>>(std::istream& is, LowerTriMatrix<T,A>& m)
    { return is >> IOStyle() >> m; }

    template <class T, int A>
    std::istream& operator>>(
        const TMV_Reader& reader, UpperTriMatrixView<T,A> m)
    {
        return reader >> 
            static_cast<BaseMatrix_Tri_Mutable<UpperTriMatrixView<T,A> >&>(m);
    }

    template <class T, int A>
    std::istream& operator>>(
        const TMV_Reader& reader, LowerTriMatrixView<T,A> m)
    {
        return reader >> 
            static_cast<BaseMatrix_Tri_Mutable<LowerTriMatrixView<T,A> >&>(m);
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    std::istream& operator>>(
        const TMV_Reader& reader, SmallUpperTriMatrixView<T,N,Si,Sj,A> m)
    {
        return reader >> 
            static_cast<BaseMatrix_Tri_Mutable<
            SmallUpperTriMatrixView<T,N,Si,Sj,A> >&>(m);
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    std::istream& operator>>(
        const TMV_Reader& reader, SmallLowerTriMatrixView<T,N,Si,Sj,A> m)
    {
        return reader >> 
            static_cast<BaseMatrix_Tri_Mutable<
            SmallLowerTriMatrixView<T,N,Si,Sj,A> >&>(m);
    }

    template <class T, int A>
    std::istream& operator>>(std::istream& is, UpperTriMatrixView<T,A> m)
    {
        return is >> 
            static_cast<BaseMatrix_Tri_Mutable<UpperTriMatrixView<T,A> >&>(m);
    }

    template <class T, int A>
    std::istream& operator>>(std::istream& is, LowerTriMatrixView<T,A> m)
    {
        return is >> 
            static_cast<BaseMatrix_Tri_Mutable<LowerTriMatrixView<T,A> >&>(m);
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    std::istream& operator>>(
        std::istream& is, SmallUpperTriMatrixView<T,N,Si,Sj,A> m)
    {
        return is >> 
            static_cast<BaseMatrix_Tri_Mutable<
            SmallUpperTriMatrixView<T,N,Si,Sj,A> >&>(m);
    }

    template <class T, ptrdiff_t N, ptrdiff_t Si, ptrdiff_t Sj, int A>
    std::istream& operator>>(
        std::istream& is, SmallLowerTriMatrixView<T,N,Si,Sj,A> m)
    {
        return is >> 
            static_cast<BaseMatrix_Tri_Mutable<
            SmallLowerTriMatrixView<T,N,Si,Sj,A> >&>(m);
    }


} // namespace mv

#endif
