
#ifndef TMV_TriMatrixIO_H
#define TMV_TriMatrixIO_H

#include "TMV_BaseMatrix_Tri.h"

namespace tmv {

    //
    // Write Matrix
    //
    // TODO: Have the non-compact versions instantiated as well.

    // Defined in TMV_TriMatrix.cpp
    template <class T, int C>
    void InstWriteCompact(
        std::ostream& os, const ConstUpperTriMatrixView<T,C>& m);
    template <class T, int C>
    void InstWriteCompact(
        std::ostream& os, const ConstLowerTriMatrixView<T,C>& m);
    template <class T, int C>
    void InstWriteCompact(
        std::ostream& os, const ConstUpperTriMatrixView<T,C>& m,
        typename ConstUpperTriMatrixView<T>::float_type thresh);
    template <class T, int C>
    void InstWriteCompact(
        std::ostream& os, const ConstLowerTriMatrixView<T,C>& m,
        typename ConstLowerTriMatrixView<T>::float_type thresh);
    template <class T>
    void InstRead(std::istream& is, UpperTriMatrixView<T> m);
    template <class T>
    void InstRead(std::istream& is, LowerTriMatrixView<T> m);


    template <int algo, class M>
    struct WriteU_Helper;

    template <class M>
    struct WriteU_Helper<1,M>
    {
        static void call(std::ostream& os, const M& m)
        {
            typedef typename M::value_type T;
            const int len = m.size();
            const bool upper = m.isupper();
            const bool unit = m.isunit();
            os << (upper ? "U " : "L ");
            os << len << '\n';
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
        static void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        {
            typedef typename M::value_type T;
            const int len = m.size();
            const bool upper = m.isupper();
            const bool unit = m.isunit();
            os << (upper ? "U " : "L ");
            os << len << '\n';
            for(int i=0;i<len;++i) {
                os << "( ";
                if (!upper) {
                    for(int j=0;j<i;++j) {
                        T temp = m.cref(i,j);
                        os << " " <<
                            Value((TMV_ABS2(temp)<thresh ? T(0) : temp)) << " ";
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
                            Value((TMV_ABS2(temp)<thresh ? T(0) : temp)) << " ";
                    }
                }
                os << " )\n";
            }
        }
    };

    // algo 90: Call inst
    template <class M>
    struct WriteU_Helper<90,M>
    {
        static TMV_INLINE void call(std::ostream& os, const M& m)
        { InstWriteCompact(os,m.calc().xView()); }
        static TMV_INLINE void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { InstWriteCompact(os,m.calc().xView(),thresh); }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M>
    struct WriteU_Helper<-3,M>
    {
        static TMV_INLINE void call(std::ostream& os, const M& m)
        { WriteU_Helper<1,M>::call(os,m); }
        static TMV_INLINE void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { WriteU_Helper<1,M>::call(os,m,thresh); }
    };
             
    // algo -2: Check for inst
    template <class M>
    struct WriteU_Helper<-2,M>
    {
        typedef typename M::value_type T;
        enum { inst = (
                (M::_size == TMV_UNKNOWN || M::_size > 16) &&
                Traits<T>::isinst ) };
        enum { algo = (
                inst ? 90 :
                -3 ) };
        static TMV_INLINE void call(std::ostream& os, const M& m)
        { WriteU_Helper<algo,M>::call(os,m); }
        static TMV_INLINE void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { WriteU_Helper<algo,M>::call(os,m,thresh); }
    };

    template <class M>
    struct WriteU_Helper<-1,M>
    {
        static TMV_INLINE void call(std::ostream& os, const M& m)
        { WriteU_Helper<-2,M>::call(os,m); }
        static TMV_INLINE void call(
            std::ostream& os, const M& m, typename M::float_type thresh)
        { WriteU_Helper<-2,M>::call(os,m,thresh); }
    };

    template <class M>
    inline void WriteCompact(
        std::ostream& os, const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteU_Helper<-2,Mv>::call(os,mv);
    }

    template <class M>
    inline void InlineWriteCompact(
        std::ostream& os, const BaseMatrix_Tri<M>& m)
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteU_Helper<-3,Mv>::call(os,mv);
    }

    template <class M>
    inline void WriteCompact(
        std::ostream& os,
        const BaseMatrix_Tri<M>& m, typename M::float_type thresh) 
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteU_Helper<-2,Mv>::call(os,mv,thresh);
    }


    template <class M>
    inline void InlineWriteCompact(
        std::ostream& os,
        const BaseMatrix_Tri<M>& m, typename M::float_type thresh) 
    {
        typedef typename M::const_cview_type Mv;
        TMV_MAYBE_CREF(M,Mv) mv = m.cView();
        WriteU_Helper<-3,Mv>::call(os,mv,thresh);
    }


    //
    // Read Matrix
    //

#ifndef TMV_NO_THROW
    template <class T>
    class TriMatrixReadError : 
        public ReadError
    {
    public :
        Matrix<T,NoDivider> m;
        int i,j;
        char exp,got;
        T unitgot;
        size_t s;
        bool up;
        bool is, iseof, isbad;

        TriMatrixReadError(std::istream& _is, bool _up) throw() :
            ReadError("TriMatrix"),
            i(0), j(0), exp(0), got(0), 
            unitgot(T(1)), s(0), up(_up),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        TriMatrixReadError(
            int _i, int _j, const BaseMatrix_Tri<M>& _m, 
            std::istream& _is) throw() :
            ReadError("TriMatrix"),
            m(_m), i(_i), j(_j), exp(0), got(0), 
            unitgot(T(1)), s(_m.size()), up(_m.isupper()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        TriMatrixReadError(
            int _i, int _j, const BaseMatrix_Tri<M>& _m,
            std::istream& _is, char _e, char _g) throw() :
            ReadError("TriMatrix"),
            m(_m), i(_i), j(_j), exp(_e), got(_g),
            unitgot(T(1)), s(_m.size()),  up(_m.isupper()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        TriMatrixReadError(
            std::istream& _is, char _e, char _g, bool _up) throw() :
            ReadError("TriMatrix"),
            i(0), j(0), exp(_e), got(_g),
            unitgot(T(1)), s(0),  up(_up),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        TriMatrixReadError(
            int _i, int _j, const BaseMatrix_Tri<M>& _m,
            std::istream& _is, T _u) throw() :
            ReadError("TriMatrix"),
            m(_m), i(_i), j(_j), exp(0), got(0),
            unitgot(_u), s(_m.size()),  up(_m.isupper()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        TriMatrixReadError(
            const BaseMatrix_Tri<M>& _m, std::istream& _is, size_t _s) throw():
            ReadError("TriMatrix"),
            m(_m), i(0), exp(0), got(0), 
            unitgot(T(1)), s(_s), up(_m.isupper()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        TriMatrixReadError(const TriMatrixReadError<T>& rhs) :
            ReadError("TriMatrix"),
            m(rhs.m), i(rhs.i), j(rhs.j), exp(rhs.exp), got(rhs.got), 
            unitgot(rhs.unitgot), s(rhs.s), up(rhs.up),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~TriMatrixReadError() throw() {}

        void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for ";
            if (up) os<<"UpperTriMatrix\n";
            else os<<"LowerTriMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (unitgot != T(1)) {
                os<<"Wrong format: expected 1 on the diagonal, got '"<<
                    unitgot<<"'.\n";
            }
            if (s != m.colsize()) {
                os<<"Wrong size: expected "<<m.colsize()<<", got "<<s<<".\n";
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
            if (m.colsize() > 0) {
                const int N = m.colsize();
                os<<"The portion of the TriMatrix which was successfully "
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
    struct ReadU_Helper;

    template <class M>
    struct ReadU_Helper<1,M>
    {
        static void call(std::istream& is, M& m)
        {
            typedef typename M::value_type T;
            char paren;
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
                    throw TriMatrixReadError<T>(i,0,m,is,'(',is?paren:'(');
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
                            throw TriMatrixReadError<T>(i,j,m,is);
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
                        throw TriMatrixReadError<T>(i,i,m,is,is?temp:T(1));
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
                            throw TriMatrixReadError<T>(i,j,m,is);
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
                    throw TriMatrixReadError<T>(i,len,m,is,')',is?paren:')');
#endif
                }
            }
        }
    };

    // algo 90: Call inst
    template <class M>
    struct ReadU_Helper<90,M>
    {
        static TMV_INLINE void call(std::istream& is, M& m)
        { InstRead(is,m.xView()); }
    };
             
    // algo 97: Conjugate
    template <class M>
    struct ReadU_Helper<97,M>
    {
        static TMV_INLINE void call(std::istream& is, M& m)
        {
            typedef typename M::conjugate_type Mc;
            Mc mc = m.conjugate();
            ReadU_Helper<-2,Mc>::call(is,mc); 
            mc.conjugateSelf();
        }
    };
             
    // algo -3: Only one algorithm, so call it.
    template <class M>
    struct ReadU_Helper<-3,M>
    {
        static TMV_INLINE void call(std::istream& is, M& m)
        { ReadU_Helper<1,M>::call(is,m); }
    };
             
    // algo -2: Check for inst
    template <class M>
    struct ReadU_Helper<-2,M>
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
            ReadU_Helper<algo,M>::call(is,m); 
        }
    };

    template <class M>
    struct ReadU_Helper<-1,M>
    {
        static TMV_INLINE void call(std::istream& is, M& m)
        { ReadU_Helper<-2,M>::call(is,m); }
    };

    template <class M>
    inline void Read(std::istream& is, BaseMatrix_Tri_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadU_Helper<-2,Mv>::call(is,mv);
    }

    template <class M>
    inline void InlineRead(
        std::istream& is, BaseMatrix_Tri_Mutable<M>& m)
    {
        typedef typename M::cview_type Mv;
        TMV_MAYBE_REF(M,Mv) mv = m.cView();
        ReadU_Helper<-3,Mv>::call(is,mv);
    }


    // 
    // Operator overloads for I/O
    // is >> m
    // os << m
    //

    template <class M>
    static std::istream& operator>>(
        std::istream& is, BaseMatrix_Tri_Mutable<M>& m)
    {
        typedef typename M::value_type T;
        char ul;
        char ul_exp = (M::_upper ? 'U' : 'L');
        is >> ul;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"TriMatrix ReadError: !is\n";
            exit(1);
#else
            throw TriMatrixReadError<T>(is,M::_upper);
#endif
        }
        if (ul != ul_exp) {
#ifdef TMV_NO_THROW
            std::cerr<<"TriMatrix ReadError: "<<ul<<" != "<<ul_exp<<std::endl;
            exit(1);
#else
            throw TriMatrixReadError<T>(is,ul_exp,ul,M::_upper);
#endif
        }

        size_t s;
        is >> s;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"TriMatrix ReadError: !is \n";
            exit(1);
#else
            throw TriMatrixReadError<T>(is,M::_upper);
#endif
        }
        if (s != m.size()) {
#ifdef TMV_NO_THROW
            std::cerr<<"TriMatrix ReadError: Wrong size \n";
            exit(1);
#else
            throw TriMatrixReadError<T>(m,is,s);
#endif
        }
        Read(is,m);
        return is;
    }

    template <class T, int A0, int A1, int A2>
    static std::istream& operator>>(
        std::istream& is, UpperTriMatrix<T,A0,A1,A2>& m)
    {
        char ul;
        is >> ul;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"UpperTriMatrix ReadError: !is\n";
            exit(1);
#else
            throw TriMatrixReadError<T>(is,true);
#endif
        }
        if (ul != 'U') {
#ifdef TMV_NO_THROW
            std::cerr<<"UpperTriMatrix ReadError: "<<ul<<" != U\n";
            exit(1);
#else
            throw TriMatrixReadError<T>(is,'U',ul,true);
#endif
        }
        size_t s;
        is >> s;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"UpperTriMatrix ReadError: !is \n";
            exit(1);
#else
            throw TriMatrixReadError<T>(is,true);
#endif
        }
        m.resize(s);
        Read(is,m);
        return is;
    }

    template <class T, int A0, int A1, int A2>
    static std::istream& operator>>(
        std::istream& is, LowerTriMatrix<T,A0,A1,A2>& m)
    {
        char ul;
        is >> ul;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"LowerTriMatrix ReadError: !is\n";
            exit(1);
#else
            throw TriMatrixReadError<T>(is,false);
#endif
        }
        if (ul != 'L') {
#ifdef TMV_NO_THROW
            std::cerr<<"LowerTriMatrix ReadError: "<<ul<<" != U\n";
            exit(1);
#else
            throw TriMatrixReadError<T>(is,'L',ul,false);
#endif
        }
        size_t s;
        is >> s;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"LowerTriMatrix ReadError: !is \n";
            exit(1);
#else
            throw TriMatrixReadError<T>(is,false);
#endif
        }
        m.resize(s);
        Read(is,m);
        return is;
    }

} // namespace mv

#endif
