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


#ifndef TMV_DiagMatrixIO_H
#define TMV_DiagMatrixIO_H

#include "TMV_DiagMatrix.h"
#include "TMV_VectorIO.h"
#include <iostream>

namespace tmv {


    //
    // Write DiagMatrix
    //

    template <class M>
    static void InlineWrite(std::ostream& os, const BaseMatrix_Diag<M>& m)
    {
        typedef typename M::value_type T;
        const int n = m.size();
        os << n << "  " << n << std::endl;
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
    template <class T, bool C>
    void InstWrite(
        std::ostream& os, const ConstDiagMatrixView<T,UNKNOWN,C>& m);

    // With thresh:
    template <class M>
    static void InlineWrite(
        std::ostream& os, const BaseMatrix_Diag<M>& m,
        typename M::real_type thresh) 
    {
        typedef typename M::value_type T;
        const int n = m.size();
        os << n << "  " << n << std::endl;
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
    template <class T, bool C>
    void InstWrite(
        std::ostream& os, const ConstDiagMatrixView<T,UNKNOWN,C>& m,
        typename Traits<T>::real_type thresh);

    //
    // Read DiagMatrix
    //

#ifndef TMV_NO_THROW
    template <class M> 
    class DiagMatrixReadError : public ReadError
    {
    public :
        typedef typename M::copy_type copy_type;
        int i;
        mutable auto_ptr<copy_type> m;
        char exp,got;
        size_t s;
        bool is, iseof, isbad;

        DiagMatrixReadError(std::istream& _is) throw() :
            ReadError("DiagMatrix"),
            i(0), m(0), exp(0), got(0), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(
            int _i, const BaseMatrix_Diag<M>& _m, char _e, char _g, size_t _s,
            bool _is, bool _iseof, bool _isbad
        ) throw() :
            ReadError("DiagMatrix"),
            i(_i), m(new copy_type(_m)), exp(_e), got(_g), s(_s),
            is(_is), iseof(_iseof), isbad(_isbad) {}
        DiagMatrixReadError(
            const BaseMatrix_Diag<M>& _m, std::istream& _is, size_t _s
        ) throw() :
            ReadError("DiagMatrix"),
            i(0), m(new copy_type(_m)), exp(0), got(0), s(_s),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(
            std::istream& _is, char _e, char _g
        ) throw() :
            ReadError("DiagMatrix"),
            i(0), m(0), exp(_e), got(_g), s(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        DiagMatrixReadError(const DiagMatrixReadError<M>& rhs) :
            ReadError("DiagMatrix"),
            i(rhs.i), m(rhs.m), exp(rhs.exp), got(rhs.got), s(rhs.s),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}

        ~DiagMatrixReadError() throw() {}

        virtual void Write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for DiagMatrix\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (m.get() && s != m->size()) {
                os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
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
            if (m.get()) {
                os<<"The portion of the DiagMatrix which was successfully "
                    "read is: \n";
                os<<"( ";
                for(int ii=0;ii<i;++ii)
                    os<<' '<<(*m)(ii,ii)<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    template <class T, IndexStyle I> 
    std::istream& operator>>(std::istream& is, auto_ptr<DiagMatrix<T,I> >& m)
    {
        typedef DiagMatrix<T,I> type;
        char d;
        is >> d;
        if (!is || d != 'D') {
#ifdef TMV_NO_THROW
            std::cerr<<"DiagMatrix ReadError: "<<d<<" != D\n"; 
            exit(1); 
#else
            throw DiagMatrixReadError<type>(is,'D',d);
#endif
        }
        size_t size;
        is >> size;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"DiagMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw DiagMatrixReadError<type>(is);
#endif
        }
        m.reset(new type(size));
#ifndef TMV_NO_THROW
        try {
#endif
            m->diag().read(is);
#ifndef TMV_NO_THROW
        } catch (VectorReadError<Vector<T,I> >& ve) {
            throw DiagMatrixReadError<type>(
                ve.i,*m,ve.exp,ve.got,ve.s,ve.is,ve.iseof,ve.isbad);
        }
#endif
        return is;
    }

    template <class M> 
    std::istream& operator>>(
        std::istream& is, BaseMatrix_Diag_Mutable<M>& m)
    {
        char d;
        is >> d;
        if (!is || d != 'D') {
#ifdef TMV_NO_THROW
            std::cerr<<"DiagMatrix ReadError: "<<d<<" != D\n"; 
            exit(1); 
#else
            throw DiagMatrixReadError<M>(is,'D',d);
#endif
        }
        size_t s;
        is >> s;
        if (!is) {
#ifdef TMV_NO_THROW
            std::cerr<<"DiagMatrix ReadError: !is\n"; 
            exit(1); 
#else
            throw DiagMatrixReadError<M>(is);
#endif
        }
        if (m.size() != s) {
#ifdef TMV_NO_THROW
            std::cerr<<"DiagMatrix ReadError: Wrong size\n"; 
            exit(1); 
#else
            throw DiagMatrixReadError<M>(m,is,s);
#endif
        }
        TMVAssert(m.size() == s);
#ifndef TMV_NO_THROW
        try {
#endif
            m.diag().read(is);
#ifndef TMV_NO_THROW
        } catch (VectorReadError<typename M::diag_type>& ve) {
            throw DiagMatrixReadError<M>(
                ve.i,m,ve.exp,ve.got,ve.s,ve.is,ve.iseof,ve.isbad);
        }
#endif
        return is;
    }

} // namespace tmv

#endif
