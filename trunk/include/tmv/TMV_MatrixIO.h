
#ifndef TMV_MatrixIO_H
#define TMV_MatrixIO_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_IOStyle.h"

namespace tmv {

    //
    // Write Matrix
    //

    // Defined in TMV_Matrix.cpp
    template <class T, int C>
    void InstWrite(
        const TMV_Writer& writer, const ConstMatrixView<T,C>& m);
    template <class T>
    void InstRead(const TMV_Reader& reader, MatrixView<T> m);

    template <int algo, class M1>
    struct WriteM_Helper;

    template <class M1>
    struct WriteM_Helper<11,M1>
    {
        static void call(const TMV_Writer& writer, const M1& m)
        {
            const int M = m.colsize();
            const int N = m.rowsize();
            writer.begin();
            writer.writeCode("M");
            writer.writeSize(M);
            writer.writeSize(N);
            writer.writeStart();
            for(int i=0;i<M;++i) {
                writer.writeLParen();
                for(int j=0;j<N;++j) {
                    if (j > 0) writer.writeSpace();
                    writer.writeValue(m.cref(i,j));
                }
                writer.writeRParen();
                if (i < M-1) writer.writeRowEnd();
            }
            writer.writeFinal();
            writer.end();
        }
    };

    // algo 90: Call inst
    template <class M1>
    struct WriteM_Helper<90,M1>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M1& m)
        { InstWrite(writer,m.calc().xView()); }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M1>
    struct WriteM_Helper<-3,M1>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M1& m)
        { WriteM_Helper<11,M1>::call(writer,m); }
    };
             
    // algo -2: Check for inst
    template <class M1>
    struct WriteM_Helper<-2,M1>
    {
        typedef typename M1::value_type T;
        enum { inst = (
                (M1::_colsize == Unknown || M1::_colsize > 16) &&
                (M1::_rowsize == Unknown || M1::_rowsize > 16) &&
                Traits<T>::isinst ) };
        enum { algo = (
                inst ? 90 :
                -3 ) };
        static TMV_INLINE void call(const TMV_Writer& writer, const M1& m)
        { WriteM_Helper<algo,M1>::call(writer,m); }
    };

    template <class M1>
    struct WriteM_Helper<-1,M1>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M1& m)
        { WriteM_Helper<-2,M1>::call(writer,m); }
    };

    template <class M>
    inline void Write(const TMV_Writer& writer, const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteM_Helper<-2,Mv>::call(writer,mv);
    }

    template <class M>
    inline void InlineWrite(
        const TMV_Writer& writer, const BaseMatrix_Rec<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteM_Helper<-3,Mv>::call(writer,mv);
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
        std::string exp,got;
        int cs,rs;
        bool is, iseof, isbad;

        MatrixReadError(std::istream& _is) throw() :
            ReadError("Matrix"),
            i(0), j(0), cs(0), rs(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        MatrixReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("Matrix"),
            i(0), j(0), exp(_e), got(_g), cs(0), rs(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        template <class M>
        MatrixReadError(
            int _i, int _j, const BaseMatrix_Rec<M>& _m, 
            std::istream& _is) throw() :
            ReadError("Matrix"),
            m(_m), i(_i), j(_j), 
            cs(_m.colsize()), rs(_m.rowsize()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        MatrixReadError(
            int _i, int _j, const BaseMatrix_Rec<M>& _m,
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("Matrix"),
            m(_m), i(_i), j(_j), exp(_e), got(_g),
            cs(_m.colsize()), rs(_m.rowsize()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        MatrixReadError(
            const BaseMatrix_Rec<M>& _m,
            std::istream& _is, int _cs, int _rs) throw() :
            ReadError("Matrix"),
            m(_m), i(0), cs(_cs), rs(_rs),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        MatrixReadError(const MatrixReadError<T>& rhs) throw() :
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
                    for(int jj=0;jj<N;++jj) os<<' '<<m.cref(ii,jj)<<' ';
                    os<<" )\n";
                }
                os<<"( ";
                for(int jj=0;jj<j;++jj) os<<' '<<m.cref(i,jj)<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    template <int algo, class M1>
    struct ReadM_Helper;

    template <class M1>
    struct ReadM_Helper<11,M1>
    {
        static void call(const TMV_Reader& reader, M1& m)
        {
            typedef typename M1::value_type T;
            const int M = m.colsize();
            const int N = m.rowsize();
            std::string exp, got;
            T temp;
            if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw MatrixReadError<T>(0,0,m,reader.getis(),exp,got);
#endif
            }
            for(int i=0;i<M;++i) {
                if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw MatrixReadError<T>(i,0,m,reader.getis(),exp,got);
#endif
                }
                for(int j=0;j<N;++j) {
                    if (j>0) {
                        if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                            std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
                            exit(1);
#else
                            throw MatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                        }
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"Matrix Read Error: reading value\n";
                        exit(1);
#else
                        throw MatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    m.ref(i,j) = temp;
                }
                if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw MatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
                }
                if (i < M-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw MatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
                }
            }
            if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw MatrixReadError<T>(M,0,m,reader.getis(),exp,got);
#endif
            }
        }
    };

    // algo 90: Call inst
    template <class M1>
    struct ReadM_Helper<90,M1>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M1& m)
        { InstRead(reader,m.xView()); }
    };
             
    // algo 97: Conjugate
    template <class M1>
    struct ReadM_Helper<97,M1>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M1& m)
        {
            typedef typename M1::conjugate_type Mc;
            Mc mc = m.conjugate();
            ReadM_Helper<-2,Mc>::call(reader,mc); 
            mc.conjugateSelf();
        }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M1>
    struct ReadM_Helper<-3,M1>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M1& m)
        { ReadM_Helper<11,M1>::call(reader,m); }
    };
             
    // algo -2: Check for inst
    template <class M1>
    struct ReadM_Helper<-2,M1>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M1& m)
        {
            typedef typename M1::value_type T;
            const int inst = 
                (M1::_colsize == Unknown || M1::_colsize > 16) &&
                (M1::_rowsize == Unknown || M1::_rowsize > 16) &&
                Traits<T>::isinst;
            const int algo = 
                M1::_conj ? 97 :
                inst ? 90 :
                -3;
            ReadM_Helper<algo,M1>::call(reader,m); 
        }
    };

    template <class M1>
    struct ReadM_Helper<-1,M1>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M1& m)
        { ReadM_Helper<-2,M1>::call(reader,m); }
    };

    template <class M>
    inline void Read(const TMV_Reader& reader, BaseMatrix_Rec_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadM_Helper<-2,Mv>::call(reader,mv);
    }

    template <class M>
    inline void InlineRead(
        const TMV_Reader& reader, BaseMatrix_Rec_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadM_Helper<-3,Mv>::call(reader,mv);
    }


    // 
    // Operator overloads for I/O
    // is >> m
    // os << m
    //

    template <class M>
    inline std::ostream& operator<<(
        const TMV_Writer& writer, const BaseMatrix<M>& m)
    { Write(writer,m.calc()); return writer.getos(); }

    template <class M>
    static std::istream& operator>>(
        const TMV_Reader& reader, BaseMatrix_Rec_Mutable<M>& m)
    {
        typedef typename M::value_type T;
        std::string exp,got;
        if (!reader.readCode("M",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw MatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        int cs=m.colsize(), rs=m.rowsize();
        if (!reader.readSize(cs) || !reader.readSize(rs)) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: reading size\n";
            exit(1);
#else
            throw MatrixReadError<T>(reader.getis());
#endif
        }
        if (cs != m.colsize() || rs != m.rowsize()) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: wrong size\n";
            exit(1);
#else
            throw MatrixReadError<T>(m,reader.getis(),cs,rs);
#endif
        }
        Read(reader,m);
        return reader.getis();
    }

    template <class T, int A>
    static std::istream& operator>>(
        const TMV_Reader& reader, Matrix<T,A>& m)
    {
        std::string exp,got;
        if (!reader.readCode("M",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw MatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        int cs=m.colsize(), rs=m.rowsize();
        if (!reader.readSize(cs) || !reader.readSize(rs)) {
#ifdef NOTHROW
            std::cerr<<"Matrix Read Error: reading size\n";
            exit(1);
#else
            throw MatrixReadError<T>(reader.getis());
#endif
        }
        if (cs != m.colsize() || rs != m.rowsize()) m.resize(cs,rs);
        Read(reader,m);
        return reader.getis();
    }

    template <class M>
    std::ostream& operator<<(std::ostream& os, const BaseMatrix<M>& m)
    { return os << IOStyle() << m.mat(); }

    template <class M>
    std::istream& operator>>(std::istream& is, BaseMatrix_Mutable<M>& m)
    { return is >> IOStyle() >> m.mat(); }

    template <class T, int A>
    std::istream& operator>>(std::istream& is, Matrix<T,A>& m)
    { return is >> IOStyle() >> m; }

    template <class T, int A>
    inline std::istream& operator>>(
        const TMV_Reader& reader, MatrixView<T,A> m)
    { 
        return reader >> 
            static_cast<BaseMatrix_Rec_Mutable<MatrixView<T,A> >&>(m); 
    }

    template <class T, int M, int N, int Si, int Sj, int A>
    std::istream& operator>>(
        const TMV_Reader& reader, SmallMatrixView<T,M,N,Si,Sj,A> m)
    { 
        return reader >>
            static_cast<BaseMatrix_Rec_Mutable<
            SmallMatrixView<T,M,N,Si,Sj,A> >&>(m); 
    }

    template <class T, int A>
    std::istream& operator>>(std::istream& is, MatrixView<T,A> m)
    { 
        return is >> 
            static_cast<BaseMatrix_Rec_Mutable<MatrixView<T,A> >&>(m); 
    }

    template <class T, int M, int N, int Si, int Sj, int A>
    std::istream& operator>>(std::istream& is, SmallMatrixView<T,M,N,Si,Sj,A> m)
    { 
        return is >>
            static_cast<BaseMatrix_Rec_Mutable<
            SmallMatrixView<T,M,N,Si,Sj,A> >&>(m); 
    }

} // namespace mv

#endif
