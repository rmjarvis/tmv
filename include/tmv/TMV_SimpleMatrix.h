///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
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
//    SimpleMatrix<T,M,N,A>()
//        Makes a SimpleMatrix with column size = M and row size = N
//        with _uninitialized_ values
//
//    SimpleMatrix<T,M,N,A>(const GenMatrix<T2>& m)
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

#include "tmv/TMV_Base.h"
#include "tmv/TMV_Array.h"

namespace tmv {

    template <typename T>
    class SimpleMatrix
    {
        typedef SimpleMatrix<T> type;
        typedef T value_type;
        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        typedef real_type RT;

    public:

        //
        // Constructors
        //

        inline SimpleMatrix() : linsize(0), itsm(0), itscs(0), itsrs(0) {}

        inline SimpleMatrix(ptrdiff_t cs, ptrdiff_t rs) :
            linsize((cs)*(rs)),
            itsm(linsize), itscs(cs), itsrs(rs)
        {
            TMVAssert(cs >= 0 && rs >= 0);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline SimpleMatrix(const type& m2) :
            linsize(m2.linsize), itsm(linsize),
            itscs(m2.itscs), itsrs(m2.itsrs)
        {
            for(ptrdiff_t i=0;i<linsize;++i) itsm[i] = m2.itsm[i];
        }

        template <class M2>
        inline SimpleMatrix(const M2& m2) :
            linsize((m2.colsize())*(m2.rowsize())),
            itsm(linsize), itscs(m2.colsize()), itsrs(m2.rowsize())
        {
            for(ptrdiff_t i=0;i<itscs;++i) for(ptrdiff_t j=0;j<itsrs;++j)
                ref(i,j) = m2.cref(i,j);
        }

        inline ~SimpleMatrix()
        {
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(999));
#endif
        }


        //
        // Op=
        //

        inline type& operator=(const type& m2)
        {
            for(ptrdiff_t i=0;i<linsize;++i) itsm[i] = m2.itsm[i];
            return *this;
        }

        template <class M2>
        inline type& operator=(const M2& m2)
        {
            for(ptrdiff_t i=0;i<itscs;++i) for(ptrdiff_t j=0;j<itsrs;++j)
                ref(i,j) = m2.cref(i,j);
            return *this;
        }

        inline type& operator=(const T& x)
        { return setToIdentity(x); }


        //
        // Access
        //

        inline T operator()(ptrdiff_t i,ptrdiff_t j) const
        {
            TMVAssert(i>=0 && i<itscs);
            TMVAssert(j>=0 && j<itsrs);
            return cref(i,j);
        }

        inline T& operator()(ptrdiff_t i,ptrdiff_t j)
        {
            TMVAssert(i>=0 && i<itscs);
            TMVAssert(j>=0 && j<itsrs);
            return ref(i,j);
        }

        //
        // Functions of Matrix
        //

        inline T trace() const
        {
            TMVAssert(itscs == itsrs);
            T sum(0);
            for(ptrdiff_t i=0; i<itsrs; ++i) sum += cref(i,i);
            return sum;
        }

        inline T sumElements() const
        {
            T sum(0);
            for(ptrdiff_t i=0;i<linsize;++i) sum += itsm[i];
            return sum;
        }

        inline RT sumAbsElements() const
        {
            RT sum(0);
            for(ptrdiff_t i=0;i<linsize;++i) sum += TMV_ABS(itsm[i]);
            return sum;
        }

        inline RT sumAbs2Elements() const
        {
            RT sum(0);
            for(ptrdiff_t i=0;i<linsize;++i) sum += TMV_ABS2(itsm[i]);
            return sum;
        }

        inline RT norm() const
        { return normF(); }

        inline RT normF() const
        {
            TMVAssert(!std::numeric_limits<T>::is_integer);
            return TMV_SQRT(normSq());
        }

        inline RT normSq(RT scale=RT(1)) const
        {
            RT sum(0);
            if (scale == RT(1))
                for(ptrdiff_t i=0;i<linsize;++i) sum += TMV_NORM(itsm[i]);
            else
                for(ptrdiff_t i=0;i<linsize;++i) sum += TMV_NORM(itsm[i]*scale);
            return sum;
        }

        inline RT norm1() const
        {
            RT max(0);
            for(ptrdiff_t j=0;j<itsrs;++j) {
                RT temp(0);
                for(ptrdiff_t i=0;i<itscs;++i) temp += TMV_ABS(cref(i,j));
                if (temp > max) max = temp;
            }
            return max;
        }

        // inf-norm = max_i (sum_j |a_ij|)
        inline RT normInf() const
        {
            RT max(0);
            for(ptrdiff_t i=0;i<itscs;++i) {
                RT temp(0);
                for(ptrdiff_t j=0;j<itsrs;++j) temp += TMV_ABS(cref(i,j));
                if (temp > max) max = temp;
            }
            return max;
        }

        // = max_i,j (|a_ij|)
        inline RT maxAbsElement() const
        {
            RT max(0);
            for(ptrdiff_t i=0;i<linsize;++i) {
                RT temp = TMV_ABS(itsm[i]);
                if (temp > max) max = temp;
            }
            return max;
        }

        // = max_i,j (|real(a_ij)|+|imag(a_ij)|)
        inline RT maxAbs2Element() const
        {
            RT max(0);
            for(ptrdiff_t i=0;i<linsize;++i) {
                RT temp = TMV_ABS2(itsm[i]);
                if (temp > max) max = temp;
            }
            return max;
        }


        //
        // Modifying Functions
        //

        inline type& setZero()
        {
            for(ptrdiff_t i=0;i<linsize;++i) itsm[i] = T(0);
            return *this;
        }

        inline type& clip(RT thresh)
        {
            for(ptrdiff_t i=0;i<linsize;++i)
                if (TMV_ABS(itsm[i]) < thresh) itsm[i] = T(0);
            return *this;
        }

        inline type& setAllTo(const T& x)
        {
            for(ptrdiff_t i=0;i<linsize;++i) itsm[i] = x;
            return *this;
        }

        inline type& addToAll(const T& x)
        {
            for(ptrdiff_t i=0;i<linsize;++i) itsm[i] += x;
            return *this;
        }

        inline type& transposeSelf()
        {
            TMVAssert(itscs == itsrs);
            for(ptrdiff_t i=1; i<itscs; ++i)
                for(ptrdiff_t j=0; j<i; ++j) TMV_SWAP(ref(i,j),ref(j,i));
            return *this;
        }

        inline type& conjugateSelf()
        {
            if (isComplex(T())) {
                for(ptrdiff_t i=0;i<linsize;++i) itsm[i] = TMV_CONJ(itsm[i]);
            }
            return *this;
        }

        inline type& setToIdentity(const T& x=T(1))
        {
            TMVAssert(itscs == itsrs);
            setZero();
            for(ptrdiff_t i=0; i<itsrs; ++i) ref(i,i) = T(1);
            return *this;
        }

        inline type& swapRows(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1 >= 0 && i1 < itscs);
            TMVAssert(i2 >= 0 && i2 < itscs);
            if (i1 != i2)
                for(ptrdiff_t j=0; j<itsrs; ++j) TMV_SWAP(ref(i1,j),ref(i2,j));
            return *this;
        }

        inline type& swapCols(ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1 >= 0 && j1 < itsrs);
            TMVAssert(j2 >= 0 && j2 < itsrs);
            if (j1 != j2)
                for(ptrdiff_t i=0; i<itscs; ++i) TMV_SWAP(ref(i,j1),ref(i,j2));
            return *this;
        }

        //
        // I/O
        //

        inline void write(std::ostream& os) const
        {
            os << itscs <<"  "<< itsrs <<std::endl;
            for(ptrdiff_t i=0;i<itscs;++i) {
                os << "( ";
                for(ptrdiff_t j=0;j<itsrs;++j) os << ' '<<cref(i,j)<<' ';
                os << " )\n";
            }
        }

        inline ptrdiff_t colsize() const { return itscs; }
        inline ptrdiff_t rowsize() const { return itsrs; }
        inline bool isSquare() const { return itscs == itsrs; }
        inline ptrdiff_t stepi() const { return 1; }
        inline ptrdiff_t stepj() const { return itscs; }
        inline bool isrm() const { return false; }
        inline bool iscm() const { return true; }
        inline bool isconj() const { return false; }
        inline ConjType ct() const { return NonConj; }
        inline ptrdiff_t ls() const { return linsize; }
        inline const T* cptr() const { return itsm; }
        inline T* ptr() { return itsm; }

        inline T cref(ptrdiff_t i, ptrdiff_t j) const { return itsm[j*itscs+i]; }

        inline T& ref(ptrdiff_t i, ptrdiff_t j) { return itsm[j*itscs+i]; }

        inline void resize(ptrdiff_t cs, ptrdiff_t rs)
        {
            TMVAssert(cs >= 0 && rs >= 0);
            linsize = cs*rs;
            itsm.resize(linsize);
            itscs = cs;
            itsrs = rs;
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

    protected :

        ptrdiff_t linsize;
        AlignedArray<T> itsm;
        ptrdiff_t itscs;
        ptrdiff_t itsrs;

    }; // SimpleMatrix

    template <typename T>
    inline std::string TMV_Text(const SimpleMatrix<T>& )
    {
        return std::string("SimpleMatrix<")+TMV_Text(T())+">";
    }

} // namespace tmv

#endif
