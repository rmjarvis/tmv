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


//---------------------------------------------------------------------------
//
// This file defines the TMV SimpleMatrix class.
//
// It is used when I want some of the simpler functionality of 
// the Matrix class, but I want to make sure it works even if the
// library doesn't have the template instantiations for that type.
//
// The main place that I use it so for is in the determinant calculation
// for int (or complex<int>) matrices.  I do the calculation on a 
// long double matrix.  But I don't want to require that the 
// Matrix<long double> class is compiled into the library.
//
// So this class is completely header-only.  It doesn't have any 
// any of the arithmetic operations defined, so everything has to 
// be done on the elements directly.  Similary, it doesn't have
// all of the normal functions and methods.  Just a few that are easy
// to implement inline.
//
// Constructors:
//
//    SimpleMatrix<T,M,N,stor,I>()
//        Makes a SimpleMatrix with column size = M and row size = N
//        with _uninitialized_ values
//
//    SimpleMatrix<T,M,N,stor,I>(const GenMatrix<T2>& m)
//        Make a SimpleMatrix which copies the elements of m.
//        T2 may be different than T.
//
// SimpleMatrix doesn't have views like a regular Matrix.
// All the normal viewing kinds of routines just return a regular MatrixView.
// It is mostly useful for fast element access and simple tasks
// like multiplication and addition.  For most division routines,
// SimpleMatrix just sends the task to a regular matrix.  The exception
// is 2x2, which has specialization for det, inverse, etc.


#ifndef TMV_SimpleMatrix_H
#define TMV_SimpleMatrix_H

#include "tmv/TMV_BaseMatrix_Rec.h"
#include "tmv/TMV_Array.h"

namespace tmv {

    template <class T, StorageType S=ColMajor>
    class SimpleMatrix;

    template <class T, StorageType S>
    struct Traits<SimpleMatrix<T,S> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef SimpleMatrix<T,S> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        typedef InvalidType inverse_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _shape = Rec };
        enum { _fort = false };
        enum { _calc = true };
        enum { _rowmajor = (S == RowMajor) };
        enum { _colmajor = (S == ColMajor) };
        enum { _stor = S };
        enum { _stepi = (S==ColMajor ? 1 : UNKNOWN) };
        enum { _stepj = (S==RowMajor ? 1 : UNKNOWN) };
        enum { _diagstep = UNKNOWN };
        enum { _conj = false };
    };

    template <class T, StorageType S>
    class SimpleMatrix :
        public BaseMatrix<SimpleMatrix<T,S> >
    {
        typedef SimpleMatrix<T,S> type;
        typedef BaseMatrix<type> base;
        typedef typename Traits<T>::real_type real_type;

    public:

        //
        // Constructors
        //

        SimpleMatrix(int cs, int rs) :
            itscs(cs), itsrs(rs),
            linsize((cs)*(rs)), itsm(linsize)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMV_DEBUG
            setAllTo(T(888));
#endif
        }

        SimpleMatrix(const type& m2) :
            itscs(m2.itscs), itsrs(m2.itsrs),
            linsize(m2.linsize), itsm(linsize)
        { 
            for(int i=0;i<linsize;++i) itsm[i] = m2.itsm[i];
        }

        template <class M2>
        SimpleMatrix(const M2& m2) :
            itscs(m2.colsize()), itsrs(m2.rowsize()),
            linsize((m2.colsize())*(m2.rowsize())), itsm(linsize)
        { 
            for(int i=0;i<itscs;++i) for(int j=0;j<itsrs;++j) 
                ref(i,j) = m2.cref(i,j);
        }

        ~SimpleMatrix()
        {
#ifdef TMV_DEBUG
            setAllTo(T(999));
#endif
        }


        //
        // Op=
        //

        type& operator=(const type& m2)
        { 
            if (&m2 != this) {
                for(int i=0;i<linsize;++i) itsm[i] = m2.itsm[i];
            }
            return *this; 
        }

        template <class M2>
        type& operator=(const M2& m2)
        { 
            for(int i=0;i<itscs;++i) for(int j=0;j<itsrs;++j) 
                ref(i,j) = m2.cref(i,j);
            return *this; 
        }

        type& operator=(const T& x) 
        { return setToIdentity(x); }


        //
        // Access
        //

        T operator()(int i,int j) const
        { 
            TMVAssert(i>=0 && i<itscs);
            TMVAssert(j>=0 && j<itsrs);
            return cref(i,j); 
        }

        T& operator()(int i,int j) 
        { 
            TMVAssert(i>=0 && i<itscs);
            TMVAssert(j>=0 && j<itsrs);
            return ref(i,j);
        }

        //
        // Functions of Matrix
        //

        T trace() const
        {
            TMVAssert(itscs == itsrs);
            T sum(0);
            for(int i=0; i<itsrs; ++i) sum += cref(i,i);
            return sum;
        }

        T sumElements() const
        {
            T sum(0);
            for(int i=0;i<linsize;++i) sum += itsm[i];
            return sum;
        }

        real_type sumAbsElements() const
        {
            real_type sum(0);
            for(int i=0;i<linsize;++i) sum += TMV_ABS(itsm[i]);
            return sum;
        }

        real_type sumAbs2Elements() const
        {
            real_type sum(0);
            for(int i=0;i<linsize;++i) sum += TMV_ABS2(itsm[i]);
            return sum;
        }

        real_type norm() const 
        { return normF(); }

        real_type normF() const
        {
            TMVAssert(!std::numeric_limits<T>::is_integer);
            return TMV_SQRT(normSq()); 
        }

        real_type normSq(real_type scale=real_type(1)) const
        { 
            real_type sum(0);
            if (scale == real_type(1))
                for(int i=0;i<linsize;++i) sum += TMV_NORM(itsm[i]);
            else
                for(int i=0;i<linsize;++i) sum += TMV_NORM(itsm[i]*scale);
            return sum;
        }

        real_type norm1() const
        {
            real_type max(0);
            for(int j=0;j<itsrs;++j) {
                real_type temp(0);
                for(int i=0;i<itscs;++i) temp += TMV_ABS(cref(i,j));
                if (temp > max) max = temp;
            }
            return max;
        }

        // inf-norm = max_i (sum_j |a_ij|)
        real_type normInf() const
        {
            real_type max(0);
            for(int i=0;i<itscs;++i) {
                real_type temp(0);
                for(int j=0;j<itsrs;++j) temp += TMV_ABS(cref(i,j));
                if (temp > max) max = temp;
            }
            return max;
        }

        // = max_i,j (|a_ij|)
        real_type maxAbsElement() const
        {
            real_type max(0);
            for(int i=0;i<linsize;++i) {
                real_type temp = TMV_ABS(itsm[i]);
                if (temp > max) max = temp;
            }
            return max;
        }

        // = max_i,j (|real(a_ij)|+|imag(a_ij)|)
        real_type maxAbs2Element() const
        {
            real_type max(0);
            for(int i=0;i<linsize;++i) {
                real_type temp = TMV_ABS2(itsm[i]);
                if (temp > max) max = temp;
            }
            return max;
        }

        //
        // Modifying Functions
        //

        type& setZero() 
        { 
            for(int i=0;i<linsize;++i) itsm[i] = T(0);
            return *this;
        }

        type& clip(real_type thresh)
        { 
            for(int i=0;i<linsize;++i) 
                if (TMV_ABS(itsm[i]) < thresh) itsm[i] = T(0);
            return *this;
        }

        type& setAllTo(const T& x) 
        {
            for(int i=0;i<linsize;++i) itsm[i] = x;
            return *this;
        }

        type& addToAll(const T& x) 
        {
            for(int i=0;i<linsize;++i) itsm[i] += x;
            return *this;
        }

        type& transposeSelf() 
        {
            TMVAssert(itscs == itsrs);
            for(int i=1; i<itscs; ++i) 
                for(int j=0; j<i; ++j) TMV_SWAP(ref(i,j),ref(j,i));
            return *this;
        }

        type& conjugateSelf() 
        {
            if (isComplex(T())) {
                for(int i=0;i<linsize;++i) itsm[i] = TMV_CONJ(itsm[i]);
            }
            return *this;
        }

        type& setToIdentity(const T& x=T(1)) 
        { 
            TMVAssert(itscs == itsrs);
            setZero();
            for(int i=0; i<itsrs; ++i) ref(i,i) = T(1);
            return *this;
        }

        type& swapRows(int i1, int i2)
        {
            TMVAssert(i1 >= 0 && i1 < itscs);
            TMVAssert(i2 >= 0 && i2 < itscs);
            if (i1 != i2)
                for(int j=0; j<itsrs; ++j) TMV_SWAP(ref(i1,j),ref(i2,j));
            return *this;
        }

        type& swapCols(int j1, int j2)
        {
            TMVAssert(j1 >= 0 && j1 < itsrs);
            TMVAssert(j2 >= 0 && j2 < itsrs);
            if (j1 != j2)
                for(int i=0; i<itscs; ++i) TMV_SWAP(ref(i,j1),ref(i,j2));
            return *this;
        }

        //
        // I/O
        //

        void write(std::ostream& os) const
        { 
            os << itscs <<"  "<< itsrs <<std::endl;
            for(int i=0;i<itscs;++i) {
                os << "( ";
                for(int j=0;j<itsrs;++j) os << ' '<<cref(i,j)<<' ';
                os << " )\n";
            }
        }

        size_t colsize() const { return itscs; }
        size_t rowsize() const { return itsrs; }
        int nElements() const { return itscs*itsrs; }
        bool isSquare() const { return itscs == itsrs; }
        int stepi() const { return S == RowMajor ? itsrs : 1; }
        int stepj() const { return S == RowMajor ? 1 : itscs; }
        bool isrm() const { return S == RowMajor; }
        bool iscm() const { return S == ColMajor; }
        bool isconj() const { return false; }
        StorageType stor() const { return S; }
        size_t ls() const { return linsize; }
        const T* cptr() const { return itsm; }
        T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return S == RowMajor ? itsm[i*itsrs+j] : itsm[j*itscs+i]; }

        T& ref(int i, int j)
        { return S == RowMajor ? itsm[i*itsrs+j] : itsm[j*itscs+i]; }

        void resize(size_t cs, size_t rs)
        {
            linsize = cs*rs;
            itsm.resize(linsize);
            itscs = cs;
            itsrs = rs;
#ifdef TMV_DEBUG
            setAllTo(T(888));
#endif
        }

    protected :

        int itscs;
        int itsrs;
        int linsize;
        AlignedArray<T> itsm;

    }; // SimpleMatrix

    template <class T, StorageType S> 
    static std::ostream& operator<<(
        std::ostream& os, const SimpleMatrix<T,S>& m)
    { m.write(os); return os; }

    template <class T, int M, int N, StorageType S, IndexStyle I> 
    static std::string TMV_Text(const SimpleMatrix<T,S>& )
    { 
        std::ostringstream s;
        s << std::string("SimpleMatrix<")<<TMV_Text(T());
        s << ','<<TMV_Text(S)<<">"; 
        return s.str();
    }

} // namespace tmv

#endif
