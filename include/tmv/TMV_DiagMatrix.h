///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
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
// This file defines the TMV DiagMatrix class.
//
// The DiagMatrix class is provided for efficient storage of a diagonal
// matrix.  You can do most of the things that you can do with a 
// regular Matrix, but it will do them more efficiently.
//
// Constructors:
//
//    DiagMatrix<T>(size_t size)
//        Makes a DiagMatrix with column size and row size = size
//        with _uninitialized_ values
//
//    DiagMatrix<T>(size_t size, T x)
//        Makes a DiagMatrix of size n with all values = x
//
//    DiagMatrix<T>(size_t size, T* vv)
//    DiagMatrix<T>(size_t size, const std::vector<T>& vv)
//        Makes a DiagMatrix of size n which copies the values is vv
//
//    DiagMatrix<T>(const Vector<T>& vv)
//        Make a DiagMatrix which copies the elements of vv.
//
//    ConstDiagMatrixView<T>(const Vector<T>& v)
//        Make a constant DiagMatrix view with v as the diagonal.
//        While this view cannon be modified, changing the original v or m
//        will cause corresponding changes in this view.
//
//    DiagMatrixView<T>(Vector<T>& v)
//        Make a mutable DiagMatrix view with v as the diagonal.
//
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//        Return the dimensions of the DiagMatrix
//
//    T& operator()(int i)
//    T operator()(int i) const
//    T& operator()(int i, int j)
//    T operator()(int i, int j) const
//        Return the (i,j) element of the DiagMatrix
//        For the single paramter version, j=i
//
//    VectorView& diag()
//    ConstVectorView& diag() const
//        Return the diagonal of the DiagMatrix as a VectorView
//
//
// Modifying Functions - The same as the regular Matrix counterparts
//
//    DiagMatrix& Zero()
//    DiagMatrix& SetAllTo(T x)
//    DiagMatrix<T>& TransposeSelf() 
//        (Does nothing.)
//    DiagMatrix& ConjugateSelf()
//    DiagMatrix& SetToIdentity(x = 1)
//    void Swap(DiagMatrix& m1, DiagMatrix& m2)
//
//
// SubDiagMatrix:
//
//    SubDiagMatrix(int i1, int i2, int istep=1)
//        Returns a Sub-DiagMatrix which extends from i1 to i2 (step istep)
//        which refers to the same physical elements as the original.
//        As usual, i2 is the "one past the end" element.
//
//
// Functions of DiagMatrices - Same as for regular Matrices:
//
//    Det(m)
//    LogDet(m)
//    Trace(m)
//    Norm(m) or NormF(m)
//    NormSq(m)
//    Norm1(m) 
//    Norm2(m) 
//    NormInf(m) 
//    MaxAbsElement(m) 
//        (Note - for diagonal matrices, 
//        Norm1 = Norm2 = NormInf = MaxAbsElement.)
//    Transpose(m)
//        (same as the original).
//    Conjugate(m)
//    Adjoint(m)
//
//    m.Inverse()
//    Inverse(m)
//    m.InvertSelf()
//    m.Inverse(minv) (takes either Matrix or DiagMatrix argument)
//    m.InverseATA(invata) (takes either Matrix or DiagMatrix argument)
//
// I/O: 
//
//    os << d 
//        Writes d to ostream os as a full matrix
//
//    d.WriteCompact(os)
//        Writes only the diagonal Vector to os
//
//    is >> d
//        Reads in d in the compact format
//
//

#ifndef TMV_DiagMatrix_H
#define TMV_DiagMatrix_H

#include "tmv/TMV_BaseMatrix_Diag.h"
#include "tmv/TMV_Shape.h"
#include "tmv/TMV_VIt.h"
#include <vector>

namespace tmv {

  template <class T, IndexStyle I>
  struct Traits<DiagMatrix<T,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef DiagMatrix<T,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef type copy_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { msize = UNKNOWN };
    enum { mfort = (I == FortranStyle) };
    enum { mshape = Diag };
    enum { mrowmajor = false };
    enum { mcolmajor = false };
    enum { mstor = ColMajor }; // arbitrary
    enum { mcalc = true };
    enum { mstep = 1 };
    enum { mconj = false };
    enum { twoS = misreal ? 1 : 2 };
    enum { notC = miscomplex };

    typedef ConstVectorView<T,mstep,false,I> const_diag_type;

    typedef ConstDiagMatrixView<T,mstep,false,I> const_subdiagmatrix_type;
    typedef ConstDiagMatrixView<T,UNKNOWN,false,I> const_subdiagmatrix_step_type;
    typedef ConstDiagMatrixView<T,mstep,false,I> const_view_type;
    typedef ConstDiagMatrixView<T,mstep,false,CStyle> const_cview_type;
    typedef ConstDiagMatrixView<T,mstep,false,FortranStyle> const_fview_type;
    typedef ConstDiagMatrixView<T> const_xview_type;
    typedef ConstDiagMatrixView<T,1,false> const_rmview_type;
    typedef ConstDiagMatrixView<T,1,false> const_cmview_type;
    typedef ConstDiagMatrixView<T,mstep,notC,I> const_conjugate_type;
    typedef ConstDiagMatrixView<T,mstep,false,I> const_transpose_type;
    typedef ConstDiagMatrixView<T,mstep,notC,I> const_adjoint_type;
    typedef ConstDiagMatrixView<real_type,twoS,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstDiagMatrixView<T,mstep,false,I> const_nonconj_type;

    typedef InvalidType inverse_type;
    typedef DiagMatrixView<T,mstep,false,I> nonconst_type;

    typedef T& reference;

    typedef VectorView<T,mstep,false,I> diag_type;

    typedef DiagMatrixView<T,mstep,false,I> subdiagmatrix_type;
    typedef DiagMatrixView<T,UNKNOWN,false,I> subdiagmatrix_step_type;
    typedef DiagMatrixView<T,mstep,false,I> view_type;
    typedef DiagMatrixView<T,mstep,false,CStyle> cview_type;
    typedef DiagMatrixView<T,mstep,false,FortranStyle> fview_type;
    typedef DiagMatrixView<T> xview_type;
    typedef DiagMatrixView<T,1,false> rmview_type;
    typedef DiagMatrixView<T,1,false> cmview_type;
    typedef DiagMatrixView<T,mstep,notC,I> conjugate_type;
    typedef DiagMatrixView<T,mstep,false,I> transpose_type;
    typedef DiagMatrixView<T,mstep,notC,I> adjoint_type;
    typedef DiagMatrixView<real_type,twoS,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef DiagMatrixView<T,mstep,false,I> nonconj_type;
  };

#ifdef XTEST
#ifdef TMV_DEBUG
#define XTEST_DEBUG
#endif
#endif

  // A Helper class to make a DiagMatrix from a BaseMatrix
  // If the BaseMatrix is assignable to Diag then we do the 
  // assign.  But if not, then we allow the construction - we
  // just make the DiagMatrix from the diagonal of the matrix.
  template <bool assignable_to_diag>
  struct DiagCopy // true
  {
    template <class M2, class T, IndexStyle I>
    static void copy(const BaseMatrix<M2>& m2, DiagMatrix<T,I>& m1)
    { m2.AssignTo(m1); }
    template <class M2, class T, int N, IndexStyle I>
    static void copy(const BaseMatrix<M2>& m2, SmallDiagMatrix<T,N,I>& m1)
    { m2.AssignTo(m1); }
  };
  template <>
  struct DiagCopy<false>
  {
    template <class M2, class T, IndexStyle I>
    static void copy(const BaseMatrix<M2>& m2, DiagMatrix<T,I>& m1)
    { m1.diag() = m2.calc().diag(); }
    template <class M2, class T, int N, IndexStyle I>
    static void copy(const BaseMatrix<M2>& m2, SmallDiagMatrix<T,N,I>& m1)
    { m1.diag() = m2.calc().diag(); }
  };

  template <class T, IndexStyle I>
  class DiagMatrix :
    public BaseMatrix_Diag_Mutable<DiagMatrix<T,I> >
  {
  public:

    typedef DiagMatrix<T,I> type;
    typedef BaseMatrix_Diag_Mutable<type> base_mut;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };
    enum { mshape = Traits<type>::mshape };

    //
    // Constructors
    //

#define NEW_SIZE(cs,rs) \

    inline DiagMatrix(size_t s) :
      itssize(s), itsm(new T[itssize])
    {
#ifdef TMV_DEBUG
      this->SetAllTo(T(888));
#endif
    }

    inline DiagMatrix(size_t s, T x) :
      itssize(s), itsm(new T[itssize])
    {
      this->SetAllTo(x);
    }

    inline DiagMatrix(size_t s, const T* vv) :
      itssize(s), itsm(new T[itssize])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      std::copy(vv,vv+itssize,itsm.get());
    }

    inline DiagMatrix(const std::vector<T>& vv) : 
      itssize(vv.size()), itsm(new T[itssize])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVAssert(vv.size() == itssize);
      std::copy(vv.begin(),vv.end(),itsm.get());
    }

    inline DiagMatrix(const type& m2) :
      itssize(m2.itssize), itsm(new T[itssize])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      std::copy(m2.cptr(),m2.cptr()+itssize,itsm.get());
    }

    template <class V2>
    inline DiagMatrix(const BaseVector<V2>& v2) : 
      itssize(v2.size()), itsm(new T[itssize])
    {
      TMVStaticAssert(V2::misreal || miscomplex);
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      v2.AssignTo(this->diag());
    }

    template <class M2>
    inline DiagMatrix(const BaseMatrix<M2>& m2) :
      itssize(m2.colsize()), itsm(new T[itssize])
    {
      TMVStaticAssert((Sizes<M2::mrowsize,M2::mcolsize>::same));
      TMVAssert(m2.colsize() == m2.rowsize());
      TMVStaticAssert(M2::misreal || miscomplex);
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      DiagCopy<ShapeTraits2<M2::mshape,Diag>::assignable>::copy(m2,*this);
    }

    inline ~DiagMatrix() 
    {
#ifdef TMV_DEBUG
      this->SetAllTo(T(999));
#endif
    }
  
    //
    // Op=
    //

    inline type& operator=(const type& m2)
    {
      if (&m2 != this)
        std::copy(m2.cptr(),m2.cptr()+itssize,itsm.get());
      return *this; 
    }

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    {
      TMVStaticAssert(M2::misreal || miscomplex);
      TMVStaticAssert((ShapeTraits2<M2::mshape,Diag>::assignable));
      base_mut::operator=(m2);
      return *this; 
    }

    inline type& operator=(T x)
    {
      base_mut::operator=(x);
      return *this; 
    }


    //
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm.get(); }
    inline T* ptr() { return itsm.get(); }

    inline T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
    inline T cref(int i) const { return itsm[i]; }
    inline T& ref(int i) { return itsm[i]; }

    inline size_t size() const { return itssize; }
    inline int step() const { return 1; }
    inline bool isconj() const { return false; }
    inline bool isrm() const { return true; }
    inline bool iscm() const { return true; }

  private:

    const size_t itssize;
    auto_array<T> itsm;

  }; // DiagMatrix

  template <class T>
  class DiagMatrixF :
    public DiagMatrix<T,FortranStyle>
  {
  public:
    typedef DiagMatrixF<T> type;
    typedef DiagMatrix<T,FortranStyle> mtype;

    inline DiagMatrixF(size_t s) : mtype(s) {}
    inline DiagMatrixF(size_t s, T x) : mtype(s,x) {}
    inline DiagMatrixF(size_t s, const T* vv) : mtype(s,vv) {}
    inline DiagMatrixF(const std::vector<T>& vv) : mtype(vv) {}
    inline DiagMatrixF(const type& m2) : mtype(m2) {}
    template <class M2>
    inline DiagMatrixF(const BaseMatrix_Diag<M2>& m2) : mtype(m2) {}
    inline ~DiagMatrixF() {}

    inline type& operator=(const type& m2)
    { mtype::operator=(m2); return *this; }
    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    { mtype::operator=(m2); return *this; }
    inline type& operator=(T x)
    { mtype::operator=(x); return *this; }

  }; // DiagMatrixF

  template <class T, int S, bool C, IndexStyle I>
  struct Traits<ConstDiagMatrixView<T,S,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef ConstDiagMatrixView<T,S,C,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef DiagMatrix<T,I> copy_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { msize = UNKNOWN };
    enum { mfort = (I == FortranStyle) };
    enum { mshape = Diag };
    enum { mrowmajor = (S==1) };
    enum { mcolmajor = (S==1) };
    enum { mstor = ColMajor }; // arbitrary
    enum { mcalc = true };
    enum { mstep = S };
    enum { mconj = C };
    enum { twoS = misreal ? S : IntTraits<S>::twoS };
    enum { notC = !C && miscomplex };

    typedef ConstVectorView<T,S,C,I> const_diag_type;

    typedef ConstDiagMatrixView<T,S,C,I> const_subdiagmatrix_type;
    typedef ConstDiagMatrixView<T,UNKNOWN,C,I> const_subdiagmatrix_step_type;
    typedef ConstDiagMatrixView<T,S,C,I> const_view_type;
    typedef ConstDiagMatrixView<T,S,C,CStyle> const_cview_type;
    typedef ConstDiagMatrixView<T,S,C,FortranStyle> const_fview_type;
    typedef ConstDiagMatrixView<T,UNKNOWN,C> const_xview_type;
    typedef ConstDiagMatrixView<T,1,C> const_rmview_type;
    typedef ConstDiagMatrixView<T,1,C> const_cmview_type;
    typedef ConstDiagMatrixView<T,S,notC,I> const_conjugate_type;
    typedef ConstDiagMatrixView<T,S,C,I> const_transpose_type;
    typedef ConstDiagMatrixView<T,S,notC,I> const_adjoint_type;
    typedef ConstDiagMatrixView<real_type,twoS,C,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstDiagMatrixView<T,S,false,I> const_nonconj_type;

    typedef InvalidType inverse_type;
    typedef DiagMatrixView<T,S,C,I> nonconst_type;
  };

  template <class T, int S, bool C, IndexStyle I>
  class ConstDiagMatrixView :
    public BaseMatrix_Diag<ConstDiagMatrixView<T,S,C,I> >
  {
  public:

    typedef ConstDiagMatrixView<T,S,C,I> type;

    //
    // Constructors
    //

    inline ConstDiagMatrixView(const T* m, size_t n, int s) :
      itsm(m), itssize(n), itsstep(s) {}

    inline ConstDiagMatrixView(const T* m, size_t n) :
      itsm(m), itssize(n), itsstep(S)
    { TMVStaticAssert(S != UNKNOWN); }

    inline ConstDiagMatrixView(const type& m2) :
      itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

    template <int S2, IndexStyle I2>
    inline ConstDiagMatrixView(const ConstDiagMatrixView<T,S2,C,I2>& m2) :
      itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

    template <int S2, IndexStyle I2>
    inline ConstDiagMatrixView(const DiagMatrixView<T,S2,C,I2>& m2) :
      itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

    template <int N2, int S2, IndexStyle I2>
    inline ConstDiagMatrixView(
        const ConstSmallDiagMatrixView<T,N2,S2,C,I2>& m2) :
      itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

    template <int N2, int S2, IndexStyle I2>
    inline ConstDiagMatrixView(
        const SmallDiagMatrixView<T,N2,S2,C,I2>& m2) :
      itsm(m2.cptr()), itssize(m2.size()), itsstep(m2.step()) {}

    inline ~ConstDiagMatrixView() { 
#ifdef TMV_DEBUG
      itsm = 0; 
#endif
    }

  private :
    inline void operator=(const type& m2);
  public :

    //
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm; }

    inline T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
    inline T cref(int i) const
    { return DoConj<C>(itsm[i*step()]); }

    inline size_t size() const { return itssize; }
    inline int step() const { return itsstep; }
    inline bool isconj() const { return C; }
    inline bool isrm() const { return step()==1; }
    inline bool iscm() const { return step()==1; }

  private :

#ifdef TMV_DEBUG
    const T* itsm;
#else
    const T*const itsm;
#endif
    const size_t itssize;
    const StepInt<S> itsstep;

  }; // ConstDiagMatrixView

  template <class T, int S=1, bool C=false>
  class ConstDiagMatrixViewF :
    public ConstDiagMatrixView<T,S,C,FortranStyle>
  {
  public:

    typedef ConstDiagMatrixViewF<T,S,C> type;
    typedef ConstDiagMatrixView<T,S,C,FortranStyle> mtype;

    inline ConstDiagMatrixViewF(const T* m, size_t n, int s) : mtype(m,n,s) {}
    inline ConstDiagMatrixViewF(const T* m, size_t n) : mtype(m,n) {}
    inline ConstDiagMatrixViewF(const type& m2) : mtype(m2) {}
    template <int S2, IndexStyle I2>
    inline ConstDiagMatrixViewF(const ConstDiagMatrixView<T,S2,C,I2>& m2) :
      mtype(m2) {}
    template <int S2, IndexStyle I2>
    inline ConstDiagMatrixViewF(const DiagMatrixView<T,S2,C,I2>& m2) :
      mtype(m2) {}
    template <int N2, int S2, IndexStyle I2>
    inline ConstDiagMatrixViewF(
        const ConstSmallDiagMatrixView<T,N2,S2,C,I2>& m2) : mtype(m2) {}
    template <int N2, int S2, IndexStyle I2>
    inline ConstDiagMatrixViewF(const SmallDiagMatrixView<T,N2,S2,C,I2>& m2) :
      mtype(m2) {}
    inline ~ConstDiagMatrixViewF() {}

  private :
    inline void operator=(const type& m2);

  }; // ConstDiagMatrixViewF


  template <class T, int S, bool C, IndexStyle I>
  struct Traits<DiagMatrixView<T,S,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef DiagMatrixView<T,S,C,I> type;
    typedef const ConstDiagMatrixView<T,S,C,I> calc_type;
    typedef calc_type eval_type;
    typedef DiagMatrix<T,I> copy_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { msize = UNKNOWN };
    enum { mfort = (I == FortranStyle) };
    enum { mshape = Diag };
    enum { mrowmajor = (S==1) }; 
    enum { mcolmajor = (S==1) }; 
    enum { mstor = ColMajor }; // arbitrary
    enum { mcalc = true };
    enum { mstep = S };
    enum { mconj = C };
    enum { twoS = misreal ? S : IntTraits<S>::twoS };
    enum { notC = !C && miscomplex };

    typedef ConstVectorView<T,S,C,I> const_diag_type;

    typedef ConstDiagMatrixView<T,S,C,I> const_subdiagmatrix_type;
    typedef ConstDiagMatrixView<T,UNKNOWN,C,I> const_subdiagmatrix_step_type;
    typedef ConstDiagMatrixView<T,S,C,I> const_view_type;
    typedef ConstDiagMatrixView<T,S,C,CStyle> const_cview_type;
    typedef ConstDiagMatrixView<T,S,C,FortranStyle> const_fview_type;
    typedef ConstDiagMatrixView<T,UNKNOWN,C> const_xview_type;
    typedef ConstDiagMatrixView<T,1,C> const_rmview_type;
    typedef ConstDiagMatrixView<T,1,C> const_cmview_type;
    typedef ConstDiagMatrixView<T,S,notC,I> const_conjugate_type;
    typedef ConstDiagMatrixView<T,S,C,I> const_transpose_type;
    typedef ConstDiagMatrixView<T,S,notC,I> const_adjoint_type;
    typedef ConstDiagMatrixView<real_type,twoS,C,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstDiagMatrixView<T,S,false,I> const_nonconj_type;

    typedef InvalidType inverse_type;
    typedef DiagMatrixView<T,S,C,I> nonconst_type;

    typedef typename AuxRef<T,C>::reference reference;

    typedef VectorView<T,S,C,I> diag_type;

    typedef DiagMatrixView<T,S,C,I> subdiagmatrix_type;
    typedef DiagMatrixView<T,UNKNOWN,C,I> subdiagmatrix_step_type;
    typedef DiagMatrixView<T,S,C,I> view_type;
    typedef DiagMatrixView<T,S,C,CStyle> cview_type;
    typedef DiagMatrixView<T,S,C,FortranStyle> fview_type;
    typedef DiagMatrixView<T,UNKNOWN,C> xview_type;
    typedef DiagMatrixView<T,1,C> rmview_type;
    typedef DiagMatrixView<T,1,C> cmview_type;
    typedef DiagMatrixView<T,S,notC,I> conjugate_type;
    typedef DiagMatrixView<T,S,C,I> transpose_type;
    typedef DiagMatrixView<T,S,notC,I> adjoint_type;
    typedef DiagMatrixView<real_type,twoS,C,I> realview_type;
    typedef realview_type imagview_type;
    typedef DiagMatrixView<T,S,false,I> nonconj_type;
  };

  template <class T, int S, bool C, IndexStyle I>
  class DiagMatrixView :
    public BaseMatrix_Diag_Mutable<DiagMatrixView<T,S,C,I> >
  {
  public:

    typedef DiagMatrixView<T,S,C,I> type;
    typedef BaseMatrix_Diag_Mutable<type> base_mut;
    typedef typename base_mut::reference reference;

    //
    // Constructors
    //

    inline DiagMatrixView(T* m, size_t n, int s) :
      itsm(m), itssize(n), itsstep(s) {}

    inline DiagMatrixView(T* m, size_t n) :
      itsm(m), itssize(n), itsstep(S)
    { TMVStaticAssert(S != UNKNOWN); }

    inline DiagMatrixView(const type& m2) :
      itsm(m2.itsm), itssize(m2.size()), itsstep(m2.step()) {}

    template <int S2, IndexStyle I2>
    inline DiagMatrixView(DiagMatrixView<T,S2,C,I2> m2) :
      itsm(m2.ptr()), itssize(m2.size()), itsstep(m2.step()) {}

    template <int N2, int S2, IndexStyle I2>
    inline DiagMatrixView(SmallDiagMatrixView<T,N2,S2,C,I2> m2) :
      itsm(m2.ptr()), itssize(m2.size()), itsstep(m2.step()) {}

    inline ~DiagMatrixView() { 
#ifdef TMV_DEBUG
      itsm = 0; 
#endif
    }


    //
    // Op = 
    //

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    {
      base_mut::operator=(m2); 
      return *this;
    }

    inline type& operator=(const type& m2)
    {
      base_mut::operator=(m2); 
      return *this;
    }

    inline type& operator=(const T x)
    {
      base_mut::operator=(x); 
      return *this;
    }


    //
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm; }
    inline T* ptr() { return itsm; }

    inline T cref(int i, int j) const { return (i!=j ? T(0) : cref(i)); }
    inline T cref(int i) const
    { return DoConj<C>(itsm[i*step()]); }

    inline reference ref(int i) 
    { return reference(itsm[i*step()]); }

    inline size_t size() const { return itssize; }
    inline int step() const { return itsstep; }
    inline bool isconj() const { return C; }
    inline bool isrm() const { return step()==1; }
    inline bool iscm() const { return step()==1; }

  private :

#ifdef TMV_DEBUG
    T* itsm;
#else
    T*const itsm;
#endif
    const size_t itssize;
    const StepInt<S> itsstep;

  }; // DiagMatrixView

  template <class T, int S=1, bool C=false>
  class DiagMatrixViewF :
    public DiagMatrixView<T,S,C,FortranStyle>
  {
  public:

    typedef DiagMatrixViewF<T,S,C> type;
    typedef DiagMatrixView<T,S,C,FortranStyle> mtype;

    inline DiagMatrixViewF(T* m, size_t n, int s) : mtype(m,n,s) {}
    inline DiagMatrixViewF(T* m, size_t n) : mtype(m,n) {}
    inline DiagMatrixViewF(const type& m2) : mtype(m2) {}
    template <int S2, IndexStyle I2>
    inline DiagMatrixViewF(DiagMatrixView<T,S2,C,I2> m2) : mtype(m2) {}
    template <int N2, int S2, IndexStyle I2>
    inline DiagMatrixViewF(SmallMatrixView<T,N2,S2,C,I2> m2) : mtype(m2) {}
    inline ~DiagMatrixViewF() {}

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    { mtype::operator=(m2); return *this; }
    inline type& operator=(const type& m2)
    { mtype::operator=(m2); return *this; }
    inline type& operator=(const T x)
    { mtype::operator=(x); return *this; }

  }; // DiagMatrixViewF

  //
  // Special Constructors:
  //   DiagMatrixViewOf(v)
  //   DiagMatrixViewOf(T* v, n)
  //

#define T typename V::value_type
#define S V::vstep
#define C V::vconj
#define I (V::vfort ? FortranStyle : CStyle)
  template <class V> 
  inline ConstDiagMatrixView<T,S,C,I> DiagMatrixViewOf(const BaseVector_Calc<V>& v)
  { return ConstDiagMatrixView<T,S,C,I>(v.cptr(),v.size(),v.step()); }

  template <class V> 
  inline DiagMatrixView<T,S,C,I> DiagMatrixViewOf(BaseVector_Mutable<V>& v)
  { return DiagMatrixView<T,S,C,I>(v.ptr(),v.size(),v.step()); }
#undef T
#undef S
#undef C
#undef I

  template <class T> 
  inline ConstDiagMatrixView<T> DiagMatrixViewOf(const T* v, size_t size)
  { return ConstDiagMatrixView<T,1>(v,size); }

  template <class T> 
  inline DiagMatrixView<T> DiagMatrixViewOf(T* v, size_t size)
  { return DiagMatrixView<T,1>(v,size); }


  //
  // TypeText 
  //

  template <class T, IndexStyle I>
  inline std::string TypeText(const DiagMatrix<T,I>& )
  {
    std::ostringstream s;
    s << "DiagMatrix<"<<TypeText(T())<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, int S, bool C, IndexStyle I>
  inline std::string TypeText(const ConstDiagMatrixView<T,S,C,I>& m)
  {
    std::ostringstream s;
    s << "ConstDiagMatrixView<"<<TypeText(T());
    s << ","<<IntTraits<S>::text();
    if (S == UNKNOWN) s << "("<<m.step()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, int S, bool C, IndexStyle I>
  inline std::string TypeText(const DiagMatrixView<T,S,C,I>& m)
  {
    std::ostringstream s;
    s << "DiagMatrixView<"<<TypeText(T());
    s << ","<<IntTraits<S>::text();
    if (S == UNKNOWN) s << "("<<m.step()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }

} // namespace tmv

#endif
