
#ifndef TMV_BandMatrixIO_H
#define TMV_BandMatrixIO_H

#include "TMV_BaseMatrix_Band.h"
#include "TMV_MatrixIO.h"

namespace tmv {

    //
    // Write BandMatrix
    //

    // Defined in TMV_BandMatrix.cpp
    template <class T, int C>
    void InstWriteCompact(
        std::ostream& os, const ConstBandMatrixView<T,C>& m);
    template <class T, int C>
    void InstWriteCompact(
        std::ostream& os, const ConstBandMatrixView<T,C>& m,
        typename ConstBandMatrixView<T>::float_type thresh);
    template <class T>
    void InstRead(std::istream& is, BandMatrixView<T> m);

    template <int algo, class M>
    struct WriteB_Helper;

    template <class M>
    struct WriteB_Helper<11,M>
    {
        static void call(std::ostream& os, const M& m)
        {
            const int nrows = m.nrows();
            const int ncols = m.ncols();
            os << "B " << nrows << "  " << ncols << " " <<
                m.nlo() << " " << m.nhi()<<std::endl;
            for(int i=0;i<nrows;++i) {
                os << "( ";
                for(int j=m.rowstart(i);j<m.rowend(i);++j) {
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
            os << "B " << nrows << "  " << ncols << " " <<
                m.nlo() << " " << m.nhi()<<std::endl;
            for(int i=0;i<nrows;++i) {
                os << "( ";
                for(int j=m.rowstart(i);j<m.rowend(i);++j) {
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
    struct WriteB_Helper<90,M>
    {
        static TMV_INLINE void call(std::ostream& os, const M& m)
        { InstWriteCompact(os,m.calc().xView()); }
        static TMV_INLINE void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { InstWriteCompact(os,m.calc().xView(),thresh); }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M>
    struct WriteB_Helper<-3,M>
    {
        static TMV_INLINE void call(std::ostream& os, const M& m)
        { WriteB_Helper<11,M>::call(os,m); }
        static TMV_INLINE void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { WriteB_Helper<11,M>::call(os,m,thresh); }
    };
             
    // algo -2: Check for inst
    template <class M>
    struct WriteB_Helper<-2,M>
    {
        typedef typename M::value_type T;
        enum { inst = (
                (M::_colsize == TMV_UNKNOWN || M::_colsize > 16) &&
                (M::_rowsize == TMV_UNKNOWN || M::_rowsize > 16) &&
                Shape(M::_shape) == Band &&
                Traits<T>::isinst ) };
        enum { algo = (
                inst ? 90 :
                -3 ) };
        static TMV_INLINE void call(std::ostream& os, const M& m)
        { WriteB_Helper<algo,M>::call(os,m); }
        static TMV_INLINE void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { WriteB_Helper<algo,M>::call(os,m,thresh); }
    };

    template <class M>
    struct WriteB_Helper<-1,M>
    {
        static TMV_INLINE void call(std::ostream& os, const M& m)
        { WriteB_Helper<-2,M>::call(os,m); }
        static TMV_INLINE void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { WriteB_Helper<-2,M>::call(os,m,thresh); }
    };

    template <class M>
    static inline void WriteCompact(
        std::ostream& os, const BaseMatrix_Band<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteB_Helper<-2,Mv>::call(os,mv);
    }

    template <class M>
    static inline void InlineWriteCompact(
        std::ostream& os, const BaseMatrix_Band<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteB_Helper<-3,Mv>::call(os,mv);
    }

    template <class M>
    static inline void WriteCompact(
        std::ostream& os,
        const BaseMatrix_Band<M>& m, typename M::float_type thresh) 
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteB_Helper<-2,Mv>::call(os,mv,thresh);
    }


    template <class M>
    static inline void InlineWriteCompact(
        std::ostream& os,
        const BaseMatrix_Band<M>& m, typename M::float_type thresh) 
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteB_Helper<-3,Mv>::call(os,mv,thresh);
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
        char exp,got;
        size_t cs,rs;
        int lo,hi;
        bool is, iseof, isbad;

        BandMatrixReadError(std::istream& _is) throw() :
            ReadError("BandMatrix"),
            i(0), j(0), exp(0), got(0), cs(0), rs(0), lo(0), hi(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        BandMatrixReadError(
            int _i, int _j, const BaseMatrix_Band<M>& _m, 
            std::istream& _is) throw() :
            ReadError("BandMatrix"),
            m(_m), i(_i), j(_j), exp(0), got(0), 
            cs(_m.colsize()), rs(_m.rowsize()), lo(_m.nlo()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        BandMatrixReadError(
            int _i, int _j, const BaseMatrix_Band<M>& _m,
            std::istream& _is, char _e, char _g) throw() :
            ReadError("BandMatrix"),
            m(_m), i(_i), j(_j), exp(_e), got(_g),
            cs(_m.colsize()), rs(_m.rowsize()), lo(_m.nlo()), hi(_m.nhi()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        BandMatrixReadError(std::istream& _is, char _e, char _g) throw() :
            ReadError("BandMatrix"),
            i(0), j(0), exp(_e), got(_g),
            cs(0), rs(0), lo(0), hi(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        BandMatrixReadError(
            const BaseMatrix_Band<M>& _m, std::istream& _is,
            size_t _cs, size_t _rs, int _lo, int _hi) throw() :
            ReadError("BandMatrix"),
            m(_m), i(0), exp(0), got(0), cs(_cs), rs(_rs), lo(_lo), hi(_hi),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        BandMatrixReadError(const BandMatrixReadError<T>& rhs) throw() :
            ReadError("BandMatrix"),
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got), 
            cs(rhs.cs), rs(rhs.rs), lo(rhs.lo), hi(rhs.hi),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~BandMatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for BandMatrix\n";
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
            if (m.colsize() > 0 || m.rowsize() > 0) {
                const int N = m.rowsize();
                os<<"The portion of the BandMatrix which was successfully "
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
    struct ReadB_Helper;

    template <class M>
    struct ReadB_Helper<11,M>
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
                    std::cerr<<"BandMatrix ReadError: "<<paren<<" != (\n"; 
                    exit(1); 
#else
                    throw BandMatrixReadError<T>(i,0,m,is,'(',is?paren:'(');
#endif
                }
                for(int j=m.rowstart(i);j<m.rowend(i);++j) {
                    is >> temp;
                    if (!is) {
#ifdef TMV_NO_THROW
                        std::cerr<<"BandMatrix ReadError: !is\n"; 
                        exit(1); 
#else
                        throw BandMatrixReadError<T>(i,j,m,is);
#endif
                    }
                    m.ref(i,j) = temp;
                } 
                is >> paren;
                if (!is || paren != ')') {
#ifdef TMV_NO_THROW
                    std::cerr<<"BandMatrix ReadError: "<<paren<<" != )\n"; 
                    exit(1); 
#else
                    throw BandMatrixReadError<T>(i,ncols,m,is,')',is?paren:')');
#endif
                }
            }
        }
    };

    // algo 90: Call inst
    template <class M>
    struct ReadB_Helper<90,M>
    {
        static TMV_INLINE void call(std::istream& is, M& m)
        { InstRead(is,m.xView()); }
    };
             
    // algo 97: Conjugate
    template <class M>
    struct ReadB_Helper<97,M>
    {
        static TMV_INLINE void call(std::istream& is, M& m)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            ReadB_Helper<-2,Mc>::call(is,mc); 
            mc.conjugateSelf();
        }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M>
    struct ReadB_Helper<-3,M>
    {
        static TMV_INLINE void call(std::istream& is, M& m)
        { ReadB_Helper<11,M>::call(is,m); }
    };
             
    // algo -2: Check for inst
    template <class M>
    struct ReadB_Helper<-2,M>
    {
        static TMV_INLINE void call(std::istream& is, M& m)
        {
            typedef typename M::value_type T;
            const int inst = 
                (M::_colsize == TMV_UNKNOWN || M::_colsize > 16) &&
                (M::_rowsize == TMV_UNKNOWN || M::_rowsize > 16) &&
                Traits<T>::isinst;
            const int algo = 
                M::_conj ? 97 :
                inst ? 90 :
                -3;
            ReadB_Helper<algo,M>::call(is,m); 
        }
    };

    template <class M>
    struct ReadB_Helper<-1,M>
    {
        static TMV_INLINE void call(std::istream& is, M& m)
        { ReadB_Helper<-2,M>::call(is,m); }
    };

    template <class M>
    static inline void Read(std::istream& is, BaseMatrix_Band_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadB_Helper<-2,Mv>::call(is,mv);
    }

    template <class M>
    static inline void InlineRead(
        std::istream& is, BaseMatrix_Band_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadB_Helper<-3,Mv>::call(is,mv);
    }


    // 
    // Operator overloads for I/O
    // is >> m
    // os << m
    //

    template <class M>
    static std::istream& operator>>(
        std::istream& is, BaseMatrix_Band_Mutable<M>& m)
    {
        typedef typename M::value_type T;
        char b;
        is >> b;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"BandMatrix ReadError: !is\n";
            exit(1);
#else
            throw BandMatrixReadError<T>(is);
#endif
        }
        if (b != 'B') {
#ifdef TMV_NO_THROW
            std::cerr<<"BandMatrix ReadError: "<<b<<" != B\n";
            exit(1);
#else
            throw BandMatrixReadError<T>(is,'B',b);
#endif
        }

        size_t cs,rs;
        int lo,hi;
        is >> cs >> rs >> lo >> hi;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"BandMatrix ReadError: !is \n"; 
            exit(1); 
#else
            throw BandMatrixReadError<T>(is);
#endif
        }
        if (cs != m.colsize() || rs != m.rowsize()) {
#ifdef TMV_NO_THROW
            std::cerr<<"BandMatrix ReadError: Wrong size \n"; 
            exit(1); 
#else
            throw BandMatrixReadError<T>(m,is,cs,rs,lo,hi);
#endif
        }
        Read(is,m);
        return is;
    }

    template <class T, int A0, int A1>
    static std::istream& operator>>(std::istream& is, BandMatrix<T,A0,A1>& m)
    {
        char b;
        is >> b;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"BandMatrix ReadError: !is\n";
            exit(1);
#else
            throw BandMatrixReadError<T>(is);
#endif
        }
        if (b != 'B') {
#ifdef TMV_NO_THROW
            std::cerr<<"BandMatrix ReadError: "<<b<<" != B\n";
            exit(1);
#else
            throw BandMatrixReadError<T>(is,'B',b); 
#endif
        }
        size_t cs,rs;
        int lo,hi;
        is >> cs >> rs >> lo >> hi;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"BandMatrix ReadError: !is \n"; 
            exit(1); 
#else
            throw BandMatrixReadError<T>(is);
#endif
        }
        m.resize(cs,rs,lo,hi);
        Read(is,m);
        return is;
    }

} // namespace mv

#endif
