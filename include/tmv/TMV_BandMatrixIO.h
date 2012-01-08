
#ifndef TMV_BandMatrixIO_H
#define TMV_BandMatrixIO_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_IOStyle.h"

namespace tmv {

    //
    // Write BandMatrix
    //

    // Defined in TMV_BandMatrix.cpp
    template <class T, int C>
    void InstWrite(
        const TMV_Writer& writer, const ConstBandMatrixView<T,C>& m);
    template <class T>
    void InstRead(const TMV_Reader& reader, BandMatrixView<T> m);

    template <int algo, class M1>
    struct WriteB_Helper;

    template <class M1>
    struct WriteB_Helper<11,M1>
    {
        static void call(const TMV_Writer& writer, const M1& m)
        {
            typedef typename M1::value_type T;
            const int M = m.colsize();
            const int N = m.rowsize();
            int j1=0;
            int j2=m.nhi()+1;

            writer.begin();
            writer.writeCode("B");
            writer.writeSize(M);
            writer.writeSize(N);
            writer.writeFullSize(m.nlo());
            writer.writeFullSize(m.nhi());
            writer.writeStart();

            for(int i=0;i<M;++i) {
                writer.writeLParen();
                if (!writer.isCompact()) {
                    for(int j=0;j<j1;++j) {
                        writer.writeValue(T(0));
                        if (j < N-1) writer.writeSpace();
                    }
                }
                for(int j=j1;j<j2;++j) {
                    if (j > j1) writer.writeSpace();
                    writer.writeValue(m.cref(i,j));
                }
                if (!writer.isCompact()) {
                    for(int j=j2;j<N;++j) {
                        writer.writeSpace();
                        writer.writeValue(T(0));
                    }
                }
                writer.writeRParen();
                if (i < M-1) writer.writeRowEnd();
                if (j2 < N) ++j2;
                if (i >= m.nlo() && j1 < N) ++j1;
            }
            writer.writeFinal();
            writer.end();
        }
    };

    // algo 90: Call inst
    template <class M1>
    struct WriteB_Helper<90,M1>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M1& m)
        { InstWrite(writer,m.calc().xView()); }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M1>
    struct WriteB_Helper<-3,M1>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M1& m)
        { WriteB_Helper<11,M1>::call(writer,m); }
    };
             
    // algo -2: Check for inst
    template <class M1>
    struct WriteB_Helper<-2,M1>
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
        { WriteB_Helper<algo,M1>::call(writer,m); }
    };

    template <class M1>
    struct WriteB_Helper<-1,M1>
    {
        static TMV_INLINE void call(const TMV_Writer& writer, const M1& m)
        { WriteB_Helper<-2,M1>::call(writer,m); }
    };

    template <class M>
    inline void Write(const TMV_Writer& writer, const BaseMatrix_Band<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteB_Helper<-2,Mv>::call(writer,mv);
    }

    template <class M>
    inline void InlineWrite(
        const TMV_Writer& writer, const BaseMatrix_Band<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteB_Helper<-3,Mv>::call(writer,mv);
    }


    //
    // Read BandMatrix
    //

#ifndef TMV_NO_THROW
    template <class T>
    class BandMatrixReadError : 
        public ReadError
    {
    public :
        BandMatrix<T,NoDivider> m;
        int i,j;
        std::string exp,got;
        int cs,rs;
        int lo,hi;
        T v1;
        bool is, iseof, isbad;

        BandMatrixReadError(std::istream& _is) throw() :
            ReadError("BandMatrix"),
            i(0), j(0), cs(0), rs(0), lo(0), hi(0), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        BandMatrixReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("BandMatrix"),
            i(0), j(0), exp(_e), got(_g), cs(0), rs(0), lo(0), hi(0), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        template <class M>
        BandMatrixReadError(
            const BaseMatrix_Band<M>& _m, std::istream& _is,
            int _cs, int _rs, int _lo, int _hi) throw() :
            ReadError("BandMatrix"),
            m(_m), i(0), j(0), cs(_cs), rs(_rs), lo(_lo), hi(_hi), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        BandMatrixReadError(
            int _i, int _j, const BaseMatrix_Band<M>& _m,
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("BandMatrix"),
            m(_m), i(_i), j(_j), exp(_e), got(_g), 
            cs(m.colsize()), rs(m.rowsize()), lo(m.nlo()), hi(m.nhi()), v1(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        BandMatrixReadError(
            int _i, int _j, const BaseMatrix_Band<M>& _m, 
            std::istream& _is, T _v1=0) throw() :
            ReadError("BandMatrix"),
            m(_m), i(_i), j(_j), 
            cs(m.colsize()), rs(m.rowsize()), lo(m.nlo()), hi(m.nhi()), v1(_v1),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        BandMatrixReadError(const BandMatrixReadError<T>& rhs) throw() :
            ReadError("BandMatrix"),
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got), 
            cs(rhs.cs), rs(rhs.rs), lo(rhs.lo), hi(rhs.hi), v1(rhs.v1),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~BandMatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for BandMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (cs != m.colsize()) {
                os<<"Wrong colsize: expected "<<m.colsize()<<
                    ", got "<<cs<<".\n";
            }
            if (rs != m.rowsize()) {
                os<<"Wrong rowsize: expected "<<m.rowsize()<<
                    ", got "<<rs<<".\n";
            }
            if (lo != m.nlo()) {
                os<<"Wrong nlo: expected "<<m.nlo()<<", got "<<lo<<".\n";
            }
            if (hi != m.nhi()) {
                os<<"Wrong nhi: expected "<<m.nhi()<<", got "<<hi<<".\n";
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
            if (m.colsize() > 0 || m.rowsize() > 0) {
                os<<"The portion of the BandMatrix which was successfully "
                    "read is: \n";
                const int N = m.rowsize();
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
    struct ReadB_Helper;

    template <class M1>
    struct ReadB_Helper<11,M1>
    {
        static void call(const TMV_Reader& reader, M1& m)
        {
            typedef typename M1::value_type T;
            const int M = m.colsize();
            const int N = m.rowsize();
            std::string exp, got;
            T temp;
            int j1=0;
            int j2=m.nhi()+1;
            if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw BandMatrixReadError<T>(0,0,m,reader.getis(),exp,got);
#endif
            }
            for(int i=0;i<M;++i) {
                if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw BandMatrixReadError<T>(i,0,m,reader.getis(),exp,got);
#endif
                }
                if (!reader.isCompact()) {
                    for(int j=0;j<j1;++j) {
                        if (!reader.readValue(temp)) {
#ifdef NOTHROW
                            std::cerr<<"BandMatrix Read Error: reading value\n";
                            exit(1);
#else
                            throw BandMatrixReadError<T>(i,j,m,reader.getis());
#endif
                        }
                        if (temp != T(0)) {
#ifdef NOTHROW
                            std::cerr<<"BandMatrix Read Error: "<<temp<<" != 0\n";
                            exit(1);
#else
                            throw BandMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                        }
                        if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                            std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                            exit(1);
#else
                            throw BandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                        }
                    }
                }
                for(int j=j1;j<j2;++j) {
                    if (j>j1 && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                        std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                        exit(1);
#else
                        throw BandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                    }
                    if (!reader.readValue(temp)) {
#ifdef NOTHROW
                        std::cerr<<"BandMatrix Read Error: reading value\n";
                        exit(1);
#else
                        throw BandMatrixReadError<T>(i,j,m,reader.getis());
#endif
                    }
                    m.ref(i,j) = temp;
                }
                if (!reader.isCompact()) {
                    for(int j=j2;j<N;++j) {
                        if (!reader.readSpace(exp,got)) {
#ifdef NOTHROW
                            std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                            exit(1);
#else
                            throw BandMatrixReadError<T>(i,j,m,reader.getis(),exp,got);
#endif
                        }
                        if (!reader.readValue(temp)) {
#ifdef NOTHROW
                            std::cerr<<"BandMatrix Read Error: reading value\n";
                            exit(1);
#else
                            throw BandMatrixReadError<T>(i,j,m,reader.getis());
#endif
                        }
                        if (temp != T(0)) {
#ifdef NOTHROW
                            std::cerr<<"BandMatrix Read Error: "<<temp<<" != 0\n";
                            exit(1);
#else
                            throw BandMatrixReadError<T>(i,j,m,reader.getis(),temp);
#endif
                        }
                    }
                }
                if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw BandMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
                }
                if (i < M-1 && !reader.readRowEnd(exp,got)) {
#ifdef NOTHROW
                    std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                    exit(1);
#else
                    throw BandMatrixReadError<T>(i,N,m,reader.getis(),exp,got);
#endif
                }
                if (j2 < N) ++j2;
                if (i >= m.nlo() && j1 < N) ++j1;
            }
            if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw BandMatrixReadError<T>(N,0,m,reader.getis(),exp,got);
#endif
            }
        }
    };

    // algo 90: Call inst
    template <class M1>
    struct ReadB_Helper<90,M1>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M1& m)
        { InstRead(reader,m.xView()); }
    };
             
    // algo 97: Conjugate
    template <class M1>
    struct ReadB_Helper<97,M1>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M1& m)
        {
            typedef typename M1::conjugate_type Mc;
            Mc mc = m.conjugate();
            ReadB_Helper<-2,Mc>::call(reader,mc); 
            mc.conjugateSelf();
        }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M1>
    struct ReadB_Helper<-3,M1>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M1& m)
        { ReadB_Helper<11,M1>::call(reader,m); }
    };
             
    // algo -2: Check for inst
    template <class M1>
    struct ReadB_Helper<-2,M1>
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
            ReadB_Helper<algo,M1>::call(reader,m); 
        }
    };

    template <class M1>
    struct ReadB_Helper<-1,M1>
    {
        static TMV_INLINE void call(const TMV_Reader& reader, M1& m)
        { ReadB_Helper<-2,M1>::call(reader,m); }
    };

    template <class M>
    inline void Read(const TMV_Reader& reader, BaseMatrix_Band_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadB_Helper<-2,Mv>::call(reader,mv);
    }

    template <class M>
    inline void InlineRead(
        const TMV_Reader& reader, BaseMatrix_Band_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadB_Helper<-3,Mv>::call(reader,mv);
    }


    // 
    // Operator overloads for I/O
    // is >> m
    // os << m
    //

    template <class M>
    static std::istream& operator>>(
        const TMV_Reader& reader, BaseMatrix_Band_Mutable<M>& m)
    {
        typedef typename M::value_type T;
        std::string exp,got;
        if (!reader.readCode("B",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw BandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        int cs=m.colsize(), rs=m.rowsize(), lo=m.nlo(), hi=m.nhi();
        if (!reader.readSize(cs) || !reader.readSize(rs) ||
            !reader.readFullSize(lo) || !reader.readFullSize(hi)) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: reading size\n";
            exit(1);
#else
            throw BandMatrixReadError<T>(reader.getis());
#endif
        }
        if (cs != m.colsize() || rs != m.rowsize() || 
            lo != m.nlo() || hi != m.nhi()) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: Wrong size\n";
            exit(1);
#else
            throw BandMatrixReadError<T>(m,reader.getis(),cs,rs,lo,hi);
#endif
        }
        Read(reader,m);
        return reader.getis();
    }

    template <class T, int A0, int A1>
    static std::istream& operator>>(
        const TMV_Reader& reader, BandMatrix<T,A0,A1>& m)
    {
        std::string exp,got;
        if (!reader.readCode("B",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw BandMatrixReadError<T>(reader.getis(),exp,got);
#endif
        }
        int cs=m.colsize(), rs=m.rowsize(), lo=m.nlo(), hi=m.nhi();
        if (!reader.readSize(cs) || !reader.readSize(rs) ||
            !reader.readFullSize(lo) || !reader.readFullSize(hi)) {
#ifdef NOTHROW
            std::cerr<<"BandMatrix Read Error: reading size\n";
            exit(1);
#else
            throw BandMatrixReadError<T>(reader.getis());
#endif
        }
        if (cs != m.colsize() || rs != m.rowsize() || 
            lo != m.nlo() || hi != m.nhi()) {
            m.resize(cs,rs,lo,hi);
        }
        Read(reader,m);
        return reader.getis();
    }

    template <class T, int A0, int A1>
    std::istream& operator>>(std::istream& is, BandMatrix<T,A0,A1>& m)
    { return is >> IOStyle() >> m; }

    template <class T, int A>
    std::istream& operator>>(const TMV_Reader& reader, BandMatrixView<T,A> m)
    {
        return reader >> 
            static_cast<BaseMatrix_Band_Mutable<BandMatrixView<T,A> >&>(m);
    }

    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A>
    std::istream& operator>>(
        const TMV_Reader& reader, SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> m)
    {
        return reader >> 
            static_cast<BaseMatrix_Band_Mutable<
            SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> >&>(m);
    }


    template <class T, int A>
    std::istream& operator>>(std::istream& is, BandMatrixView<T,A> m)
    {
        return is >> 
            static_cast<BaseMatrix_Band_Mutable<BandMatrixView<T,A> >&>(m);
    }

    template <class T, int M, int N, int LO, int HI, int Si, int Sj, int A>
    std::istream& operator>>(
        std::istream& is, SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> m)
    {
        return is >> 
            static_cast<BaseMatrix_Band_Mutable<
            SmallBandMatrixView<T,M,N,LO,HI,Si,Sj,A> >&>(m);
    }


} // namespace mv

#endif
