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


//-----------------------------------------------------------------------------
//
// This file defines the TMV SmallVector class.
//
// Constructors:
//
//    explicit SmallVector<T,N,I>()  
//        Makes a Vector of size N with _uninitialized_ values
//
//    SmallVector<T,N,I>(T x)
//        Makes a Vector of size N with all values = x
//
//    SmallVector<T,N,I>(const T* vv)
//    SmallVector<T,N,I>(const vector<T>& vv)
//    SmallVector<T,N,I>(const GenVector<T>& vv)
//        Makes a SmallVector which copies the elements of vv.
//
// 
// SmallVector doesn't have views like a regular Vector.
// All the normal viewing kinds of routines just return a regular VectorView.
// It is mostly useful for fast element access and simple tasks
// like multiplication and addition.  All the calculations are done
// inline, so the compiler can optimize the calculation for the particular
// value of N.  It is ually faster, but may take longer to compile
// depending on what calculations you are doing with them.
// 

#ifndef TMV_SmallVector_H
#define TMV_SmallVector_H

#include "tmv/TMV_Vector.h"
#include "tmv/TMV_VIt.h"
#include <sstream>
#include <algorithm>

namespace tmv {

    template <class T, int N> 
    class SmallVectorComposite;

    template <class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    inline void Copy(const SmallVector<T1,N,I1>& v1, SmallVector<T2,N,I2>& v2);

#ifdef XTEST
#ifdef TMVDEBUG
#define XTEST_DEBUG
#endif
#endif

    template <class T, int N, IndexStyle I> 
    class SmallVector 
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef SmallVector<T,N,I> type;
        typedef type copy_type;
        typedef ConstVectorView<T,I> const_view_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_reverse_type;
        typedef VectorView<T,I> view_type;
        typedef view_type conjugate_type;
        typedef view_type reverse_type;
        typedef ConstVectorView<RT,I> const_real_type;
        typedef VectorView<RT,I> real_type;
        typedef T value_type;
        typedef VIt<T,Unit,NonConj> iterator;
        typedef CVIt<T,Unit,NonConj> const_iterator;
        typedef VIt<T,Step,NonConj> reverse_iterator;
        typedef CVIt<T,Step,NonConj> const_reverse_iterator;
        typedef T& reference;

        //
        // Constructors
        //

        inline SmallVector() 
        {
            TMVAssert(N > 0);
#ifdef TMVDEBUG
            setAllTo(T(888));
#endif
        }

        explicit inline SmallVector(const T& x) 
        {
            TMVAssert(N > 0);
            if (x == T(0)) setZero();
            else setAllTo(x); 
        }

        explicit inline SmallVector(const T* vv) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(N > 0);
            for(int i=0;i<N;++i) itsv[i] = vv[i];
        }

        explicit inline SmallVector(const std::vector<T>& vv) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(N > 0);
            for(int i=0;i<N;++i) itsv[i] = vv[i];
        }

        inline SmallVector(const type& v2) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(N > 0);
            for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
        }

        template <IndexStyle I2> 
        inline SmallVector(const SmallVector<T,N,I2>& v2) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(N > 0);
            TMVAssert(v2.size() == N);
            for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
        }

        template <class T2, IndexStyle I2> 
        inline SmallVector(const SmallVector<T2,N,I2>& v2) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(N > 0);
            TMVAssert(v2.size() == N);
            Copy(v2,*this); 
        }

        template <IndexStyle I2> 
        inline SmallVector(const Vector<T,I2>& v2) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(N > 0);
            TMVAssert(v2.size() == N);
            for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
        }

        template <class T2> 
        inline SmallVector(const GenVector<T2>& v2) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(N > 0);
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(v2.size() == N);
            view() = v2;
        }

        inline SmallVector(const AssignableToVector<RT>& v2) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(N > 0);
            TMVAssert(v2.size() == N);
            view() = v2;
        }

        inline SmallVector(const AssignableToVector<CT>& v2) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(N > 0);
            TMVAssert(isComplex(T()));
            TMVAssert(v2.size() == N);
            view() = v2;
        }

        inline SmallVector(const SmallVectorComposite<RT,N>& v2) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(N > 0);
            v2.assignTov(*this);
        }

        inline SmallVector(const SmallVectorComposite<CT,N>& v2) 
        {
#ifdef XTEST_DEBUG
            setAllTo(T(888));
#endif
            TMVAssert(isComplex(T()));
            TMVAssert(N > 0);
            v2.assignTov(*this);
        }

        virtual inline ~SmallVector()
        {
#ifdef TMVDEBUG
            setAllTo(T(999));
#endif
        }

        //
        // Op =
        //

        inline type& operator=(type& v2)
        { 
            if (&v2 != this) 
                for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
            return *this; 
        }

        template <IndexStyle I2> 
        inline type& operator=(SmallVector<T,N,I2>& v2)
        { 
            if (&v2 != this) 
                for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
            return *this; 
        }

        template <class T2, IndexStyle I2> 
        inline type& operator=(const SmallVector<T2,N,I2>& v2) 
        { 
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(v2.size() == N);
            Copy(v2,*this);
            return *this; 
        }

        inline type& operator=(const Vector<T>& v2) 
        {
            TMVAssert(v2.size() == N);
            for(int i=0;i<N;++i) itsv[i] = v2.cref(i);
            return *this; 
        }

        template <class T2> 
        inline type& operator=(const GenVector<T2>& v2) 
        {
            TMVAssert(isReal(T2()) || isComplex(T()));
            TMVAssert(v2.size() == N);
            view() = v2;
            return *this; 
        }

        inline type& operator=(const AssignableToVector<RT>& v2) 
        {
            TMVAssert(v2.size() == N);
            view() = v2;
            return *this; 
        }

        inline type& operator=(const AssignableToVector<CT>& v2) 
        {
            TMVAssert(isComplex(T()));
            TMVAssert(v2.size() == N);
            view() = v2;
            return *this; 
        }

        inline type& operator=(const SmallVectorComposite<RT,N>& v2) 
        {
            v2.assignTov(*this);
            return *this; 
        }

        inline type& operator=(const SmallVectorComposite<CT,N>& v2) 
        {
            TMVAssert(isComplex(T()));
            v2.assignTov(*this);
            return *this; 
        }


        //
        // Access Functions
        //

        inline const_iterator begin() const 
        { return const_iterator(itsv,1); }
        inline const_iterator end() const 
        { return begin()+N; }
        inline const_reverse_iterator rbegin() const 
        { return const_reverse_iterator(itsv+(N-1),-1); }
        inline const_reverse_iterator rend() const 
        { return rbegin()+N; }

        typedef ListAssigner<T,iterator> MyListAssigner;
        TMV_DEPRECATED(inline MyListAssigner operator=(ListInitClass))
        { return MyListAssigner(begin(),size()); }

        inline MyListAssigner operator<<(const T& x)
        { return MyListAssigner(begin(),size(),x); }

        inline T operator[](int i) const 
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<N);
                return cref(i); 
            } else {
                TMVAssert(i>=1 && i<=N);
                return cref(i-1); 
            }
        }
        inline T operator()(int i) const 
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<N);
                return cref(i); 
            } else {
                TMVAssert(i>=1 && i<=N);
                return cref(i-1); 
            }
        }

        inline iterator begin() 
        { return iterator(itsv,1 TMV_FIRSTLAST ); }
        inline iterator end() 
        { return begin() + N; }

        inline reverse_iterator rbegin() 
        { return reverse_iterator(itsv+N-1,-1 TMV_FIRSTLAST ); }
        inline reverse_iterator rend() 
        { return rbegin()+N; }

        inline T& operator[](int i) 
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<N);
                return ref(i); 
            } else {
                TMVAssert(i>=1 && i<=N);
                return ref(i-1); 
            }
        }
        inline T& operator()(int i) 
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<N);
                return ref(i); 
            } else {
                TMVAssert(i>=1 && i<=N);
                return ref(i-1); 
            }
        }

        //
        // Functions of Vector
        //

        inline RT norm() const // = norm2
        { return norm2(); }

        inline RT normSq(RT scale = RT(1)) const
        {
            RT sum(0);
            if (scale == RT(1))
                for(int i=0;i<N;++i) sum += TMV_NORM(itsv[i]);
            else
                for(int i=0;i<N;++i) sum += TMV_NORM(scale*itsv[i]);
            return sum;
        }

        inline RT norm1() const // sum_i |v_i|
        { return sumAbsElements(); }

        inline RT norm2() const // sqrt( sum_i |v_i|^2 )
        { return TMV_SQRT(normSq()); }

        inline RT normInf() const // max_i |v_i|
        { return size() > 0 ? maxAbsElement() : RT(0); }

        inline T sumElements() const
        {
            T sum(0);
            for(int i=0;i<N;++i) sum += itsv[i];
            return sum;
        }

        inline RT sumAbsElements() const
        {
            RT sum(0);
            for(int i=0;i<N;++i) sum += TMV_ABS(itsv[i]);
            return sum;
        }

        inline T minElement(int* iminout=0) const
        {
            T min = N>0 ? itsv[0] : T(0);
            if (iminout) *iminout = 0;
            for(int i=1;i<N;++i) {
                if (TMV_REAL(itsv[i]) < TMV_REAL(min)) {
                    min = itsv[i];
                    if (iminout) *iminout = i;
                }
            }
            if (I == FortranStyle && iminout) ++(*iminout);
            return min;
        }

        inline T maxElement(int* imaxout=0) const
        {
            T max = N>0 ? itsv[0] : T(0);
            if (imaxout) *imaxout = 0;
            for(int i=1;i<N;++i) {
                if (TMV_REAL(itsv[i]) > TMV_REAL(max)) {
                    max = itsv[i];
                    if (imaxout) *imaxout = i;
                }
            }
            if (I == FortranStyle && imaxout) ++(*imaxout);
            return  max;
        }

        inline RT minAbsElement(int* iminout=0) const
        {
            RT min = N>0 ? 
                (isReal(T()) ? TMV_ABS(itsv[0]) : TMV_NORM(itsv[0])) : RT(0);
            if (iminout) *iminout = 0;
            for(int i=1;i<N;++i) {
                RT absvi = isReal(T()) ? TMV_ABS(itsv[i]) : TMV_NORM(itsv[i]);
                if (absvi < min) {
                    min = absvi;
                    if (iminout) *iminout = i;
                }
            }
            if (I == FortranStyle && iminout) ++(*iminout);
            return isReal(T()) ? min : TMV_SQRT(min);
        }

        inline RT maxAbsElement(int* imaxout=0) const
        {
            RT max = N>0 ? 
                (isReal(T()) ? TMV_ABS(itsv[0]) : TMV_NORM(itsv[0])) : RT(0);
            if (imaxout) *imaxout = 0;
            for(int i=1;i<N;++i) {
                RT absvi = isReal(T()) ? TMV_ABS(itsv[i]) : TMV_NORM(itsv[i]);
                if (absvi > max) {
                    max = absvi;
                    if (imaxout) *imaxout = i;
                }
            }
            if (I == FortranStyle && imaxout) ++(*imaxout);
            return isReal(T()) ? max : TMV_SQRT(max);
        }

        TMV_DEPRECATED(RT Norm() const)
        { return norm(); }
        TMV_DEPRECATED(RT NormSq(const RT scale = RT(1)) const)
        { return normSq(scale); }
        TMV_DEPRECATED(RT Norm1() const )
        { return norm1(); }
        TMV_DEPRECATED(RT Norm2() const)
        { return norm2(); }
        TMV_DEPRECATED(RT NormInf() const )
        { return normInf(); }
        TMV_DEPRECATED(T SumElements() const)
        { return sumElements(); }
        TMV_DEPRECATED(RT SumAbsElements() const)
        { return sumAbsElements(); }
        TMV_DEPRECATED(T MinElement(int* iminout=0) const)
        { return minElement(iminout); }
        TMV_DEPRECATED(T MaxElement(int* imaxout=0) const)
        { return maxElement(imaxout); }
        TMV_DEPRECATED(RT MinAbsElement(int* iminout=0) const)
        { return minAbsElement(iminout); }
        TMV_DEPRECATED(RT MaxAbsElement(int* imaxout=0) const)
        { return maxAbsElement(imaxout); }



        //
        // Modifying Functions
        //

        inline type& setZero() 
        { 
            for(int i=0;i<N;++i) itsv[i] = T(0);
            return *this;
        }

        inline type& clip(RT thresh)
        {
            for(int i=0; i<N; ++i)
                if (std::abs(itsv[i]) < thresh) itsv[i] = T(0);
            return *this;
        }

        inline type& setAllTo(const T& x)
        { 
            for(int i=0;i<N;++i) itsv[i] = x;
            return *this;
        }

        inline type& addToAll(const T& x)
        {
            for(int i=0;i<N;++i) itsv[i] += x;
            return *this;
        }

        inline type& conjugateSelf()
        { 
            if (isComplex(T())) {
                RT* itsvi = reinterpret_cast<RT*>(itsv)+1;
                for(int i=0;i<2*N;i+=2) itsvi[i] = -itsvi[i];
            }
            return *this; 
        }

        inline type& makeBasis(int i)
        { 
            if (I == CStyle) {
                TMVAssert(i>=0 && i<N);
                setZero(); itsv[i] = T(1);
            } else {
                TMVAssert(i>=1 && i<=N);
                setZero(); itsv[i-1] = T(1);
            }
            return *this; 
        }

        inline type& swap(int i1, int i2)
        {
            if (I == CStyle) {
                TMVAssert(i1>=0 && i1<N);
                TMVAssert(i2>=0 && i2<N);
                if (i1 != i2) TMV_SWAP(itsv[i1],itsv[i2]);
            } else {
                TMVAssert(i1>=1 && i1<=N);
                TMVAssert(i2>=1 && i2<=N);
                if (i1 != i2) TMV_SWAP(itsv[i1-1],itsv[i2-1]);
            }
            return *this;
        }

        inline type& permute(const int* p, int i1, int i2) 
        {
            if (I == CStyle) {
                TMVAssert(i1>=0 && i1<=i2 && i2<=N);
                for(int i=i1;i<i2;++i) swap(i,p[i]);
            } else {
                TMVAssert(i1>=1 && i1<=i2 && i2<=N);
                for(int i=i1-1;i<i2;++i) swap(i,p[i]);
            }
            return *this;
        }

        inline type& permute(const int* p) 
        { return permute(p,I==CStyle?0:1,N); }

        inline type& reversePermute(
            const int* p, int i1, int i2)
        {
            if (I == CStyle) {
                TMVAssert(i1>=0 && i1<=i2 && i2<=N);
                for(int i=i2-1;i>=i1;--i) swap(i,p[i]);
            } else {
                TMVAssert(i1>=1 && i1<=i2 && i2<=N);
                for(int i=i2-1;i>=i1-1;--i) swap(i,p[i]);
            }
            return *this;
        }

        inline type& reversePermute(const int* p)
        { return reversePermute(p,I==CStyle?0:1,N); }

        inline type& reverseSelf()
        {
            for(int i1=0,i2=N-1;i1<i2;++i1,--i2) TMV_SWAP(itsv[i1],itsv[i2]);
            return *this;
        }

        inline type& sort(int* P, ADType ad=Ascend, CompType comp=RealComp)
        { view().sort(P,ad,comp); return *this; }

        TMV_DEPRECATED(type& Zero())
        { return setZero(); }
        TMV_DEPRECATED(type& Clip(RT thresh))
        { return clip(thresh); }
        TMV_DEPRECATED(type& SetAllTo(const T& x))
        { return setAllTo(x); }
        TMV_DEPRECATED(type& AddToAll(const T& x))
        { return addToAll(x); }
        TMV_DEPRECATED(type& ConjugateSelf())
        { return conjugateSelf(); }
        TMV_DEPRECATED(type& MakeBasis(int i))
        { return makeBasis(i); }
        TMV_DEPRECATED(type& Swap(int i1, int i2))
        { return swap(i1,i2); }
        TMV_DEPRECATED(type& Permute(const int* p, int i1, int i2))
        { return permute(p,i1,i2); }
        TMV_DEPRECATED(type& Permute(const int* p))
        { return permute(p); }
        TMV_DEPRECATED(type& ReversePermute(const int* p, int i1, int i2))
        { return reversePermute(p,i1,i2); }
        TMV_DEPRECATED(type& ReversePermute(const int* p))
        { return reversePermute(p); }
        TMV_DEPRECATED(type& ReverseSelf())
        { return reverseSelf(); }
        TMV_DEPRECATED(type& Sort(
                int* p, OldADType ad=ASCEND, OldCompType comp=REAL_COMP))
        { return sort(p,ADType(ad),CompType(comp)); }
        TMV_DEPRECATED(type& Sort(
                OldADType ad=ASCEND, OldCompType comp=REAL_COMP))
        { return sort(ADType(ad),CompType(comp)); }



        //
        // subVector
        //

        inline const_view_type cSubVector(int i1, int i2) const
        { return view().cSubVector(i1,i2); }

        inline const_view_type subVector(int i1, int i2) const
        { return view().subVector(i1,i2); }

        inline view_type cSubVector(int i1, int i2)
        { return view().cSubVector(i1,i2); }

        inline view_type subVector(int i1, int i2)
        { return view().subVector(i1,i2); }

        inline const_view_type cSubVector(int i1, int i2, int istep) const
        { return view().cSubVector(i1,i2,istep); }

        inline const_view_type subVector(int i1, int i2, int istep) const
        { return view().subVector(i1,i2,istep); }

        inline view_type cSubVector(int i1, int i2, int istep)
        { return view().cSubVector(i1,i2,istep); }

        inline view_type subVector(int i1, int i2, int istep)
        { return view().subVector(i1,i2,istep); }

        inline const_view_type reverse() const
        { return view().reverse(); }

        inline view_type reverse()
        { return view().reverse(); }

        inline const_view_type view() const
        { return const_view_type(itsv,N,1,NonConj); }

        inline view_type view()
        { return view_type(itsv,N,1,NonConj); }

        inline const_view_type conjugate() const
        { return const_view_type(itsv,N,1,isReal(T())?NonConj:Conj); }

        inline view_type conjugate()
        { return view_type(itsv,N,1,isReal(T())?NonConj:Conj); }

        inline const_real_type realPart() const
        { return view().realPart(); }

        inline const_real_type imagPart() const
        { return view().imagPart(); }

        inline const_real_type flatten() const
        { return view().flatten(); }

        inline real_type realPart()
        { return view().realPart(); }

        inline real_type imagPart()
        { return view().imagPart(); }

        inline real_type flatten()
        { return view().flatten(); }

        TMV_DEPRECATED(const_view_type SubVector(int i1, int i2) const)
        { return subVector(i1,i2); }
        TMV_DEPRECATED(const_view_type SubVector(
                int i1, int i2, int istep) const)
        { return subVector(i1,i2,istep); }
        TMV_DEPRECATED(const_view_type Reverse() const)
        { return reverse(); }
        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_real_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(const_real_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(const_real_type Flatten() const)
        { return flatten(); }

        TMV_DEPRECATED(view_type SubVector(int i1, int i2))
        { return subVector(i1,i2); }
        TMV_DEPRECATED(view_type SubVector(int i1, int i2, int istep))
        { return subVector(i1,i2,istep); }
        TMV_DEPRECATED(view_type Reverse())
        { return reverse(); }
        TMV_DEPRECATED(view_type View())
        { return view(); }
        TMV_DEPRECATED(view_type Conjugate())
        { return conjugate(); }
        TMV_DEPRECATED(real_type Real())
        { return realPart(); }
        TMV_DEPRECATED(real_type Imag())
        { return imagPart(); }
        TMV_DEPRECATED(real_type Flatten())
        { return flatten(); }


        // 
        // I/O
        //

        inline void write(std::ostream& os) const
        { view().write(os); }

        inline void write(std::ostream& os, RT thresh) const
        { view().write(os,thresh); }

        TMV_DEPRECATED(void Write(std::ostream& os) const)
        { write(os); }

        TMV_DEPRECATED(void Write(std::ostream& os, RT thresh) const)
        { write(os,thresh); }

        inline size_t size() const { return N; }
        inline const T* cptr() const { return itsv; }
        inline T* ptr() { return itsv; }
        inline int step() const { return 1; }
        inline ConjType ct() const { return NonConj; }
        inline bool isconj() const { return false; }

        inline T cref(int i) const
        { return itsv[i]; }

        inline T& ref(int i)
        { return itsv[i]; }

    protected :

        T itsv[N];

    }; // SmallVector


    //
    // Copy SmallVectors
    //

    template <int N, class T1, class T2> 
    struct DoCopyv
    {
        DoCopyv(const T1* v1, T2* v2)
        { for(int i=0;i<N;++i) v2[i] = v1[i]; }
    };

    template <int N, class T> 
    struct DoCopyv<N,std::complex<T>,T>
    {
        DoCopyv(const std::complex<T>* , T* )
        { TMVAssert(TMV_FALSE); }
    };

    template <int N, class T1, class T2> 
    inline void DoCopy(const T1* v1, T2* v2)
    { DoCopyv<N,T1,T2>(v1,v2); }

    template <class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    inline void Copy(const SmallVector<T1,N,I1>& v1, SmallVector<T2,N,I2>& v2)
    { DoCopyv<N,T1,T2>(v1.cptr(),v2.ptr()); }

    //
    // Swap SmallVectors
    //

    template <class T, int N, IndexStyle I1, IndexStyle I2> 
    inline void Swap(SmallVector<T,N,I1>& v1, SmallVector<T,N,I2>& v2)
    { for(int i=0;i<N;++i) TMV_SWAP(v1.ref(i),v2.ref(i)); }

    //
    // Functions of Vectors
    //

    template <class T, int N, IndexStyle I> 
    inline TMV_RealType(T) Norm(const SmallVector<T,N,I>& v)
    { return v.norm(); }

    template <class T, int N, IndexStyle I> 
    inline TMV_RealType(T) Norm1( const SmallVector<T,N,I>& v)
    { return v.norm1(); }

    template <class T, int N, IndexStyle I> 
    inline TMV_RealType(T) NormSq( const SmallVector<T,N,I>& v)
    { return v.normSq(); }

    template <class T, int N, IndexStyle I> 
    inline TMV_RealType(T) Norm2( const SmallVector<T,N,I>& v)
    { return v.norm2(); }

    template <class T, int N, IndexStyle I> 
    inline TMV_RealType(T) NormInf( const SmallVector<T,N,I>& v)
    { return v.normInf(); }

    template <class T, int N, IndexStyle I> 
    inline T SumElements( const SmallVector<T,N,I>& v)
    { return v.sumElements(); }

    template <class T, int N, IndexStyle I> 
    inline TMV_RealType(T) SumAbsElements(const SmallVector<T,N,I>& v)
    { return v.sumAbsElements(); }

    template <class T, int N, IndexStyle I> 
    inline T MinElement(const SmallVector<T,N,I>& v)
    { return v.minElement(); }

    template <class T, int N, IndexStyle I> 
    inline T MaxElement(const SmallVector<T,N,I>& v)
    { return v.maxElement(); }

    template <class T, int N, IndexStyle I> 
    inline TMV_RealType(T) MinAbsElement(const SmallVector<T,N,I>& v)
    { return v.minAbsElement(); }

    template <class T, int N, IndexStyle I> 
    inline TMV_RealType(T) MaxAbsElement(const SmallVector<T,N,I>& v)
    { return v.maxAbsElement(); }

    template <class T, int N, IndexStyle I> 
    inline VectorView<T,I> Conjugate(const SmallVector<T,N,I>& v)
    { return v.conjugate(); }

    template <class T, int N, IndexStyle I> 
    inline VectorView<T,I> Conjugate(SmallVector<T,N,I>& v)
    { return v.conjugate(); }

    

    //
    // Vector ==, != Vector
    //

    template <class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    inline bool operator==(
        const SmallVector<T1,N,I1>& v1, const SmallVector<T2,N,I2>& v2)
    {
        for(int i=0;i<N;++i) if (v1.cref(i) != v2.cref(i)) return false;
        return true;
    }

    template <class T1, class T2, int N, IndexStyle I1, IndexStyle I2> 
    inline bool operator!=(
        const SmallVector<T1,N,I1>& v1, const SmallVector<T2,N,I2>& v2)
    { return !(v1 == v2); }

    template <class T1, class T2, int N, IndexStyle I> 
    inline bool operator==(
        const GenVector<T1>& v1, const SmallVector<T2,N,I>& v2)
    { return v1 == v2.view(); }

    template <class T1, class T2, int N, IndexStyle I> 
    inline bool operator==(
        const SmallVector<T1,N,I>& v1, const GenVector<T2>& v2)
    { return v1.view() == v2; }

    template <class T1, class T2, int N, IndexStyle I> 
    inline bool operator!=(
        const GenVector<T1>& v1, const SmallVector<T2,N,I>& v2)
    { return v1 != v2.view(); }

    template <class T1, class T2, int N, IndexStyle I> 
    inline bool operator!=(
        const SmallVector<T1,N,I>& v1, const GenVector<T2>& v2)
    { return v1.view() != v2; }


    //
    // I/O
    //

    template <class T, int N, IndexStyle I> 
    inline std::ostream& operator<<(
        std::ostream& os, const SmallVector<T,N,I>& v)
    { v.write(os); return os; }

    template <class T, int N, IndexStyle I> 
    inline std::istream& operator>>(
        std::istream& is, SmallVector<T,N,I>& v)
    { return is >> v.view(); }

    template <class T, int N, IndexStyle I> 
    inline std::string TMV_Text(const SmallVector<T,N,I>& )
    { 
        std::ostringstream s;
        s << "SmallVector<"<<TMV_Text(T())<<","<<N<<","<<TMV_Text(I)<<">";
        return s.str();
    }

} // namespace tmv

#endif
