
#ifndef TMV_DiagMatrixIO_H
#define TMV_DiagMatrixIO_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_IOStyle.h"

namespace tmv {

    //
    // Write DiagMatrix
    //

    // Defined in TMV_DiagMatrix.cpp
    template <class T, int C>
    void InstWrite(
        const TMV_Writer& writer, const ConstDiagMatrixView<T,C>& m);
    template <class T>
    void InstRead(const TMV_Reader& reader, DiagMatrixView<T> m);

    template <int algo, class M>
    struct WriteD_Helper;

    template <class M>
    struct WriteD_Helper<11,M>
    {
        static void call(const TMV_Writer& writer, const M& m)
        {
            typedef typename M::value_type T;
            const ptrdiff_t N = m.size();
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
                writer.writeValue(m.cref(i));
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
    struct WriteD_Helper<90,M>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M& m)
        { InstWrite(writer,m.calc().xView()); }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M>
    struct WriteD_Helper<-3,M>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M& m)
        { WriteD_Helper<11,M>::call(writer,m); }
    };
             
    // algo -2: Check for inst
    template <class M>
    struct WriteD_Helper<-2,M>
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
        { WriteD_Helper<algo,M>::call(writer,m); }
    };

    template <class M>
    struct WriteD_Helper<-1,M>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M& m)
        { WriteD_Helper<-2,M>::call(writer,m); }
    };

    template <class M>
    inline void Write(const TMV_Writer& writer, const BaseMatrix_Diag<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteD_Helper<-2,Mv>::call(writer,mv);
    }

    template <class M>
    inline void InlineWrite(
        const TMV_Writer& writer, const BaseMatrix_Diag<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteD_Helper<-3,Mv>::call(writer,mv);
    }


    //
    // Read DiagMatrix
    //

#ifndef TMV_NO_THROW
    template <class T>
    class DiagMatrixReadError : 
        public ReadError
    {
    public :
        DiagMatrix<T> m;
        ptrdiff_t i,j;
        std::string exp,got;
        ptrdiff_t s;
        T v1;
        bool is, iseof, isbad;

        DiagMatrixReadError(std::istream& _is) throw() :
            ReadError("DiagMatrix"),
            i(0), j(0), s(0), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("DiagMatrix"),
            i(0), j(0), exp(_e), got(_g), s(0), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        template <class M>
        DiagMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const BaseMatrix_Diag<M>& _m,
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("DiagMatrix"),
            m(_m), i(_i), j(_j), exp(_e), got(_g), s(_m.size()), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        DiagMatrixReadError(
            const BaseMatrix_Diag<M>& _m,
            std::istream& _is, ptrdiff_t _s) throw() :
            ReadError("DiagMatrix"),
            m(_m), i(0), j(0), s(_s), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        DiagMatrixReadError(
            ptrdiff_t _i, ptrdiff_t _j, const BaseMatrix_Diag<M>& _m, 
            std::istream& _is, T _v1=0) throw() :
            ReadError("DiagMatrix"),
            m(_m), i(_i), j(_j), s(_m.size()), v1(_v1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        DiagMatrixReadError(const DiagMatrixReadError<T>& rhs) throw() :
            ReadError("DiagMatrix"),
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

    template <int algo, class M>
    struct ReadD_Helper;

    template <class M>
    struct ReadD_Helper<11,M>
    {
        static void call(const TMV_Reader& reader, M& m)
        {
            typedef typename M::value_type T;
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
    };

    // algo 90: Call inst
    template <class M>
    struct ReadD_Helper<90,M>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M& m)
        { InstRead(reader,m.xView()); }
    };
             
    // algo 97: Conjugate
    template <class M>
    struct ReadD_Helper<97,M>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M& m)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            ReadD_Helper<-2,Mc>::call(reader,mc); 
            mc.conjugateSelf();
        }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M>
    struct ReadD_Helper<-3,M>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M& m)
        { ReadD_Helper<11,M>::call(reader,m); }
    };
             
    // algo -2: Check for inst
    template <class M>
    struct ReadD_Helper<-2,M>
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
            ReadD_Helper<algo,M>::call(reader,m); 
        }
    };

    template <class M>
    struct ReadD_Helper<-1,M>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M& m)
        { ReadD_Helper<-2,M>::call(reader,m); }
    };

    template <class M>
    inline void Read(const TMV_Reader& reader, BaseMatrix_Diag_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadD_Helper<-2,Mv>::call(reader,mv);
    }

    template <class M>
    inline void InlineRead(
        const TMV_Reader& reader, BaseMatrix_Diag_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadD_Helper<-3,Mv>::call(reader,mv);
    }


    // 
    // Operator overloads for I/O
    // is >> m
    // os << m
    //

    template <class M>
    static std::istream& operator>>(
        const TMV_Reader& reader, BaseMatrix_Diag_Mutable<M>& m)
    {
        typedef typename M::value_type T;
        std::string exp,got;
        if (!reader.readCode("D",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t s=m.size();
        if (!reader.readSize(s)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: reading size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis());
#endif
        }
        if (s != m.size()) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: wrong size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(m,reader.getis(),s);
#endif
        }
        s=m.size();
        if (!reader.readSimpleSize(s)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: reading size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis());
#endif
        }
        if (s != m.size()) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(m,reader.getis(),s);
#endif
        }
        Read(reader,m);
        return reader.getis();
    }

    template <class T, int A>
    static std::istream& operator>>(
        const TMV_Reader& reader, DiagMatrix<T,A>& m)
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
        ptrdiff_t s=m.size();
        if (!reader.readSize(s)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: reading size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis());
#endif
        }
        if (s != m.size()) m.resize(s);
        s=m.size();
        if (!reader.readSimpleSize(s)) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: reading size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(reader.getis());
#endif
        }
        if (s != m.size()) {
#ifdef NOTHROW
            std::cerr<<"DiagMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw DiagMatrixReadError<T>(m,reader.getis(),s);
#endif
        }
        Read(reader,m);
        return reader.getis();
    }

    template <class T, int A>
    std::istream& operator>>(std::istream& is, DiagMatrix<T,A>& m)
    { return is >> IOStyle() >> m; }

    template <class T, int A>
    std::istream& operator>>(
        const TMV_Reader& reader, DiagMatrixView<T,A> m)
    {
        return reader >> 
            static_cast<BaseMatrix_Diag_Mutable<DiagMatrixView<T,A> >&>(m); 
    }


    template <class T, ptrdiff_t N, ptrdiff_t S, int A>
    std::istream& operator>>(
        const TMV_Reader& reader, SmallDiagMatrixView<T,N,S,A> m)
    {
        return reader >> 
            static_cast<BaseMatrix_Diag_Mutable<
            SmallDiagMatrixView<T,N,S,A> >&>(m); 
    }


    template <class T, int A>
    std::istream& operator>>(std::istream& is, DiagMatrixView<T,A> m)
    {
        return is >> 
            static_cast<BaseMatrix_Diag_Mutable<DiagMatrixView<T,A> >&>(m); 
    }


    template <class T, ptrdiff_t N, ptrdiff_t S, int A>
    std::istream& operator>>(std::istream& is, SmallDiagMatrixView<T,N,S,A> m)
    {
        return is >> 
            static_cast<BaseMatrix_Diag_Mutable<
            SmallDiagMatrixView<T,N,S,A> >&>(m); 
    }


} // namespace mv

#endif
