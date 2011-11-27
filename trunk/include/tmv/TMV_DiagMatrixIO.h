

#ifndef TMV_DiagMatrixIO_H
#define TMV_DiagMatrixIO_H

#include "TMV_BaseMatrix_Diag.h"
#include "TMV_VectorIO.h"

namespace tmv {


    //
    // Write DiagMatrix
    //

    template <class M>
    static void InlineWrite(std::ostream& os, const BaseMatrix_Diag<M>& m)
    {
        typedef typename M::value_type T;
        const int n = m.size();
        os << n << "  " << n << '\n';
        for(int i=0;i<n;++i) {
            os << "( ";
            for(int j=0;j<i;++j) 
                os << " " << Value(T(0)) << " ";
            os << " " << Value(m.cref(i)) << " ";
            for(int j=i+1;j<n;++j) 
                os << " " << Value(T(0)) << " ";
            os << " )\n";
        }
    }

    // Defined in TMV_Diag.cpp
    template <class T, int C>
    void InstWrite(
        std::ostream& os, const ConstDiagMatrixView<T,C>& m);

    // With thresh:
    template <class M>
    static void InlineWrite(
        std::ostream& os, const BaseMatrix_Diag<M>& m,
        typename M::real_type thresh) 
    {
        typedef typename M::value_type T;
        const int n = m.size();
        os << n << "  " << n << '\n';
        for(int i=0;i<n;++i) {
            os << "( ";
            for(int j=0;j<i;++j) 
                os << " " << Value(T(0)) << " ";
            T temp = m.cref(i);
            os << " " << Value((TMV_ABS2(temp) < thresh ? T(0) : temp)) << " ";
            for(int j=i+1;j<n;++j) 
                os << " " << Value(T(0)) << " ";
            os << " )\n";
        }
    }

    // Defined in TMV_Diag.cpp
    template <class T, int C>
    void InstWrite(
        std::ostream& os, const ConstDiagMatrixView<T,C>& m,
        typename Traits<T>::real_type thresh);

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
        int i;
        char exp,got;
        size_t s;
        bool is, iseof, isbad;

        DiagMatrixReadError(std::istream& _is) throw() :
            ReadError("DiagMatrix"),
            i(0), exp(0), got(0), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        template <class M>
        DiagMatrixReadError(
            int _i, const BaseMatrix_Diag<M>& _m, char _e, char _g, size_t _s,
            bool _is, bool _iseof, bool _isbad) throw() :
            ReadError("DiagMatrix"),
            m(_m), i(_i), exp(_e), got(_g), s(_s),
            is(_is), iseof(_iseof), isbad(_isbad) {}
        template <class M>
        DiagMatrixReadError(
            const BaseMatrix_Diag<M>& _m,
            std::istream& _is, size_t _s) throw() :
            ReadError("DiagMatrix"),
            m(_m), i(0), exp(0), got(0), s(_s),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(std::istream& _is, char _e, char _g) throw() :
            ReadError("DiagMatrix"),
            i(0), exp(_e), got(_g), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(const DiagMatrixReadError<T>& rhs) throw() :
            ReadError("DiagMatrix"),
            m(rhs.m), i(rhs.i), exp(rhs.exp), got(rhs.got), s(rhs.s),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}

        ~DiagMatrixReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for DiagMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (m.size() > 0 && s != m.size()) {
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
            if (m.size() > 0) {
                os<<"The portion of the DiagMatrix which was successfully "
                    "read is: \n";
                os<<"( ";
                for(int ii=0;ii<i;++ii)
                    os<<' '<<m(ii,ii)<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    template <class T, int A>
    static std::istream& operator>>(std::istream& is, DiagMatrix<T,A>& m)
    {
        char d;
        is >> d;
        if (!is || d != 'D') {
#ifdef TMV_NO_THROW
            std::cerr<<"DiagMatrix ReadError: "<<d<<" != D\n"; 
            exit(1); 
#else
            throw DiagMatrixReadError<T>(is,'D',d);
#endif
        }
        size_t size;
        is >> size;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"DiagMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw DiagMatrixReadError<T>(is);
#endif
        }
        m.resize(size);
#ifndef TMV_NO_THROW
        try {
#endif
            m.diag().read(is);
#ifndef TMV_NO_THROW
        } catch (VectorReadError<T>& ve) {
            throw DiagMatrixReadError<T>(
                ve.i,m,ve.exp,ve.got,ve.s,ve.is,ve.iseof,ve.isbad);
        }
#endif
        return is;
    }

    template <class M>
    static std::istream& operator>>(
        std::istream& is, BaseMatrix_Diag_Mutable<M>& m)
    {
        typedef typename M::value_type T;
        char d;
        is >> d;
        if (!is || d != 'D') {
#ifdef TMV_NO_THROW
            std::cerr<<"DiagMatrix ReadError: "<<d<<" != D\n"; 
            exit(1); 
#else
            throw DiagMatrixReadError<T>(is,'D',d);
#endif
        }
        size_t s;
        is >> s;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"DiagMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw DiagMatrixReadError<T>(is);
#endif
        }
        if (m.size() != s) {
#ifdef TMV_NO_THROW
            std::cerr<<"DiagMatrix ReadError: Wrong size\n"; 
            exit(1); 
#else
            throw DiagMatrixReadError<T>(m,is,s);
#endif
        }
        TMVAssert(m.size() == s);
#ifndef TMV_NO_THROW
        try {
#endif
            m.diag().read(is);
#ifndef TMV_NO_THROW
        } catch (VectorReadError<T>& ve) {
            throw DiagMatrixReadError<T>(
                ve.i,m,ve.exp,ve.got,ve.s,ve.is,ve.iseof,ve.isbad);
        }
#endif
        return is;
    }

} // namespace tmv

#endif
