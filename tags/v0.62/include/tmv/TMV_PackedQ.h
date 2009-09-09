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


#ifndef TMV_PackedQ_H
#define TMV_PackedQ_H

#include "tmv/TMV_MatrixArith.h"
#include "tmv/TMV_VectorArith.h"

namespace tmv {

  template <class T> 
  class PackedQ :
    public MatrixComposite<T>
  {
  public :

    inline PackedQ(const tmv::GenMatrix<T>& _Q,
        const tmv::GenVector<T>& _beta) : Q(_Q), beta(_beta) 
    { TMVAssert(beta.size() == Q.rowsize()); }
    inline size_t colsize() const { return Q.colsize(); }
    inline size_t rowsize() const { return Q.rowsize(); }
    inline StorageType stor() const { return ColMajor; }
    inline const GenMatrix<T>& GetQ() const { return Q; }
    inline const GenVector<T>& GetBeta() const { return beta; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize());
      TMVAssert(m0.rowsize() == rowsize());
      DoAssignToM(m0);
    }

    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize());
      TMVAssert(m0.rowsize() == rowsize());
      DoAssignToM(m0);
    }

    // v = Qt v
    template <class T1> 
    void LDivEq(const VectorView<T1>& v) const;

    // v = v Qt
    template <class T1> 
    void RDivEq(const VectorView<T1>& v) const;

    // x = Qt v
    template <class T1, class T2> 
    void LDiv(const GenVector<T1>& v, const VectorView<T2>& x) const;
    // x = v Qt
    template <class T1, class T2> 
    void RDiv(const GenVector<T1>& v, const VectorView<T2>& x) const;

    // m = Qt m
    template <class T1> 
    void LDivEq(const MatrixView<T1>& m) const;

    // m = m Qt
    template <class T1> 
    void RDivEq(const MatrixView<T1>& m) const;

    // x = Qt m
    template <class T1, class T2> 
    void LDiv(const GenMatrix<T1>& m, const MatrixView<T2>& x) const;

    // x = m Qt
    template <class T1, class T2> 
    void RDiv(const GenMatrix<T1>& m, const MatrixView<T2>& x) const;

    // v = Q v
    template <class T1> 
    inline void LMultEq(const VectorView<T1>& v) const
    { RDivEq(v.Conjugate()); }

    // v = v Q
    template <class T1> 
    inline void RMultEq(const VectorView<T1>& v) const
    { LDivEq(v.Conjugate()); }

    // x = Q v
    template <class T1, class T2> 
    inline void LMult(const GenVector<T1>& v, const VectorView<T2>& x) const
    { RDiv(v.Conjugate(),x.Conjugate()); }

    // x = v Q
    template <class T1, class T2> 
    inline void RMult(const GenVector<T1>& v, const VectorView<T2>& x) const
    { LDiv(v.Conjugate(),x.Conjugate()); }

    // m = Q m
    template <class T1> 
    inline void LMultEq(const MatrixView<T1>& m) const
    { RDivEq(m.Adjoint()); }

    // m = m Q
    template <class T1> 
    inline void RMultEq(const MatrixView<T1>& m) const
    { LDivEq(m.Adjoint()); }

    // x = Q m
    template <class T1, class T2> 
    inline void LMult(const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    { RDiv(m.Adjoint(),x.Adjoint()); }

    // x = m Q
    template <class T1, class T2> 
    inline void RMult(const GenMatrix<T1>& m, const MatrixView<T2>& x) const
    { LDiv(m.Adjoint(),x.Adjoint()); }

  private : 

    const tmv::GenMatrix<T>& Q;
    const tmv::GenVector<T>& beta;

    template <class T1> 
    void DoAssignToM(const MatrixView<T1>& m0) const;

  };

#define CT std::complex<T>

  template <class T, class Tm> 
  class ProdXpQ :
    public MatrixComposite<T>
  {
  public:
    inline ProdXpQ(const T _x, const PackedQ<Tm>& _q) : x(_x), q(_q) {}
    inline size_t colsize() const { return q.colsize(); }
    inline size_t rowsize() const { return q.rowsize(); }
    inline StorageType stor() const { return ColMajor; }
    inline T GetX() const { return x; }
    inline const GenMatrix<Tm>& GetQ() const { return q; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      TMVAssert(IsReal(T()));
      MultXM(x,m0=q);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      MultXM(x,m0=q);
    }
  private:
    const T x;
    const PackedQ<Tm>& q;
  };

#define GENMATRIX PackedQ
#define PRODXM ProdXpQ
#include "tmv/TMV_AuxProdXM.h"
#undef GENMATRIX
#undef PRODXM

  template <class T, class T1, class T2> 
  class ProdpQV :
    public VectorComposite<T>
  {   
  public:
    inline ProdpQV(const T _x, const PackedQ<T1>& _q,
        const GenVector<T2>& _v) : x(_x), q(_q), v(_v)
    { TMVAssert(v.size()==q.rowsize()); }
    inline size_t size() const { return q.colsize(); }
    inline T GetX() const { return x; }
    inline const PackedQ<T1>& GetQ() const { return q; }
    inline const GenVector<T2>& GetV() const { return v; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      TMVAssert(IsReal(T()));
      q.LMult(v,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      TMVAssert(v0.size() == size());
      q.LMult(v,v0);
    }
  private:
    const T x;
    const PackedQ<T1>& q;
    const GenVector<T2>& v;
  };

  template <class T, class T1, class T2> 
  class ProdVpQ :
    public VectorComposite<T>
  {   
  public:
    inline ProdVpQ(const T _x, const GenVector<T1>& _v,
        const PackedQ<T2>& _q) : x(_x), v(_v), q(_q)
    { TMVAssert(v.size()==q.colsize()); }
    inline size_t size() const { return q.rowsize(); }
    inline T GetX() const { return x; }
    inline const GenVector<T1>& GetV() const { return v; }
    inline const PackedQ<T2>& GetQ() const { return q; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    { 
      TMVAssert(v0.size() == size());
      TMVAssert(IsReal(T()));
      q.RMult(v,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      TMVAssert(v0.size() == size());
      q.RMult(v,v0);
    }
  private:
    const T x;
    const GenVector<T1>& v;
    const PackedQ<T2>& q;
  };

  template <class T> 
  inline const VectorView<T>& operator*=(
      const VectorView<T>& v, const PackedQ<T>& m)
  {
    TMVAssert(v.size() == m.colsize());
    MultMV<false>(T(1),m.Transpose(),v,v);
    return v;
  }

  template <class T> 
  inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v, const PackedQ<T>& m)
  {
    TMVAssert(v.size() == m.colsize());
    MultMV<false>(T(1),m.Transpose(),v,v);
    return v;
  }

  template <class T, class T2> 
  inline const VectorView<T>& operator*=(
      const VectorView<T>& v, const ProdXpQ<T,T2>& pxq)
  {
    TMVAssert(pxq.colsize()==pxq.rowsize());
    TMVAssert(v.size()==pxq.rowsize());
    pxq.GetQ().RMult(v,v);
    MultXV(pxq.GetX(),v);
    return v;
  }

  template <class T> 
  inline const VectorView<CT>& operator*=(
      const VectorView<CT>& v, const ProdXpQ<T,T>& pxq)
  {
    TMVAssert(pxq.colsize()==pxq.rowsize());
    TMVAssert(v.size()==pxq.rowsize());
    pxq.GetQ().RMult(v,v);
    MultXV(pxq.GetX(),v);
    return v;
  }

#define GENMATRIX1 PackedQ
#define GENMATRIX2 GenVector
#define PRODXM1 ProdXpQ
#define PRODXM2 ProdXV
#define PRODMM ProdpQV
#define GETM1 .GetQ()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .GetQ()
#define GETM2 .GetV()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM
#define GENMATRIX1 GenVector
#define GENMATRIX2 PackedQ
#define PRODXM1 ProdXV
#define PRODXM2 ProdXpQ
#define PRODMM ProdVpQ
#define GETM1 .GetV()
#define GETM2 .GetQ()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .GetV()
#define GETM2 .GetQ()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM


  template <class T, class T1, class T2> 
  class ProdpQM :
    public MatrixComposite<T>
  {   
  public:
    inline ProdpQM(const T _x, const PackedQ<T1>& _q,
        const GenMatrix<T2>& _m) : x(_x), q(_q), m(_m)
    { TMVAssert(m.colsize()==q.rowsize()); }
    inline size_t colsize() const { return q.colsize(); }
    inline size_t rowsize() const { return m.rowsize(); }
    inline StorageType stor() const { return ColMajor; }
    inline T GetX() const { return x; }
    inline const PackedQ<T1>& GetQ() const { return q; }
    inline const GenMatrix<T2>& GetM() const { return m; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == colsize());
      TMVAssert(m0.rowsize() == rowsize());
      TMVAssert(IsReal(T()));
      q.LMult(m,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize());
      TMVAssert(m0.rowsize() == rowsize());
      q.LMult(m,m0);
    }
  private:
    const T x;
    const PackedQ<T1>& q;
    const GenMatrix<T2>& m;
  };

  template <class T, class T1, class T2> 
  class ProdMpQ :
    public MatrixComposite<T>
  {   
  public:
    inline ProdMpQ(const T _x, const GenMatrix<T1>& _m,
        const PackedQ<T2>& _q) : x(_x), m(_m), q(_q)
    { TMVAssert(m.rowsize()==q.colsize()); }
    inline size_t colsize() const { return m.colsize(); }
    inline size_t rowsize() const { return q.rowsize(); }
    inline StorageType stor() const { return ColMajor; }
    inline T GetX() const { return x; }
    inline const GenMatrix<T1>& GetM() const { return m; }
    inline const PackedQ<T2>& GetQ() const { return q; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    { 
      TMVAssert(m0.colsize() == colsize());
      TMVAssert(m0.rowsize() == rowsize());
      TMVAssert(IsReal(T()));
      q.RMult(m,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize());
      TMVAssert(m0.rowsize() == rowsize());
      q.RMult(m,m0);
    }
  private:
    const T x;
    const GenMatrix<T1>& m;
    const PackedQ<T2>& q;
  };

  template <class T> 
  inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m, const PackedQ<T>& q)
  {
    TMVAssert(q.colsize()==q.rowsize());
    TMVAssert(m.rowsize()==q.rowsize());
    q.RMultEq(m);
    return m;
  }

  template <class T> 
  inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m, const PackedQ<T>& q)
  {
    TMVAssert(q.colsize()==q.rowsize());
    TMVAssert(m.rowsize()==q.rowsize());
    q.RMultEq(m);
    return m;
  }

  template <class T, class T2> 
  inline const MatrixView<T>& operator*=(
      const MatrixView<T>& m, const ProdXpQ<T,T2>& pxq)
  {
    TMVAssert(pxq.colsize()==pxq.rowsize());
    TMVAssert(m.rowsize()==pxq.rowsize());
    pxq.GetQ().RMult(m,m);
    MultXM(pxq.GetX(),m);
    return m;
  }

  template <class T> 
  inline const MatrixView<CT>& operator*=(
      const MatrixView<CT>& m, const ProdXpQ<T,T>& pxq)
  {
    TMVAssert(pxq.colsize()==pxq.rowsize());
    TMVAssert(m.rowsize()==pxq.rowsize());
    pxq.GetQ().RMult(m,m);
    MultXM(pxq.GetX(),m);
    return m;
  }

#define GENMATRIX1 PackedQ
#define GENMATRIX2 GenMatrix
#define PRODXM1 ProdXpQ
#define PRODXM2 ProdXM
#define PRODMM ProdpQM
#define GETM1 .GetQ()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .GetQ()
#define GETM2 .GetM()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM
#define GENMATRIX1 GenMatrix
#define GENMATRIX2 PackedQ
#define PRODXM1 ProdXM
#define PRODXM2 ProdXpQ
#define PRODMM ProdMpQ
#define GETM1 .GetM()
#define GETM2 .GetQ()
#include "tmv/TMV_AuxProdMM.h"
#define GETM1 .GetM()
#define GETM2 .GetQ()
#include "tmv/TMV_AuxProdMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef PRODMM

  template <class T, class T1, class T2> 
  class QuotVpQ :
    public VectorComposite<T>
  {
  public:
    inline QuotVpQ(const T _x, const GenVector<T1>& _v,
        const PackedQ<T2>& _q) : x(_x), v(_v), q(_q)
    { TMVAssert(v.size()==q.colsize()); }
    inline size_t size() const { return q.rowsize(); }
    inline T GetX() const { return x; }
    inline const GenVector<T1>& GetV() const { return v; }
    inline const PackedQ<T2>& GetQ() const { return q; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    {
      TMVAssert(v0.size() == size());
      TMVAssert(IsReal(T()));
      q.LDiv(v,v0);
      MultXV(x,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      TMVAssert(v0.size() == size());
      q.LDiv(v,v0);
      MultXV(x,v0);
    }
  private:
    const T x;
    const GenVector<T1>& v;
    const PackedQ<T2>& q;
  };

  template <class T, class T1, class T2> 
  class RQuotVpQ :
    public VectorComposite<T>
  {
  public:
    inline RQuotVpQ(const T _x, const GenVector<T1>& _v,
        const PackedQ<T2>& _q) : x(_x), v(_v), q(_q)
    { TMVAssert(v.size()==q.colsize()); }
    inline size_t size() const { return q.rowsize(); }
    inline T GetX() const { return x; }
    inline const GenVector<T1>& GetV() const { return v; }
    inline const PackedQ<T2>& GetQ() const { return q; }
    inline void AssignToV(const VectorView<RealType(T)>& v0) const
    {
      TMVAssert(v0.size() == size());
      TMVAssert(IsReal(T()));
      q.RDiv(v,v0);
      MultXV(x,v0);
    }
    inline void AssignToV(const VectorView<ComplexType(T)>& v0) const
    {
      TMVAssert(v0.size() == size());
      q.RDiv(v,v0);
      MultXV(x,v0);
    }
  private:
    const T x;
    const GenVector<T1>& v;
    const PackedQ<T2>& q;
  };

  template <class T> 
  inline const VectorView<T>& operator/=(
      const VectorView<T>& v, const PackedQ<T>& q)
  {
    TMVAssert(q.colsize() == q.rowsize());
    TMVAssert(q.rowsize() == v.size());
    q.LDivEq(v);
    return v;
  }

  template <class T> 
  inline const VectorView<CT>& operator/=(
      const VectorView<CT>& v, const PackedQ<T>& q)
  {
    TMVAssert(q.colsize() == q.rowsize());
    TMVAssert(q.rowsize() == v.size());
    q.LDivEq(v);
    return v;
  }

  template <class T> 
  inline const VectorView<T>& operator%=(
      const VectorView<T>& v, const PackedQ<T>& q)
  {
    TMVAssert(q.colsize() == q.rowsize());
    TMVAssert(q.rowsize() == v.size());
    q.RDivEq(v);
    return v;
  }

  template <class T> 
  inline const VectorView<CT>& operator%=(
      const VectorView<CT>& v, const PackedQ<T>& q)
  {
    TMVAssert(q.colsize() == q.rowsize());
    TMVAssert(q.rowsize() == v.size());
    q.RDivEq(v);
    return v;
  }

#define GENMATRIX1 GenVector
#define GENMATRIX2 PackedQ
#define PRODXM1 ProdXV
#define PRODXM2 ProdXpQ
#define QUOTMM QuotVpQ
#define RQUOTMM RQuotVpQ
#include "tmv/TMV_AuxQuotMM.h"
#include "tmv/TMV_AuxQuotMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTMM
#undef RQUOTMM

  template <class T, class T1, class T2> 
  class QuotMpQ :
    public MatrixComposite<T>
  {
  public:
    inline QuotMpQ(const T _x, const GenMatrix<T1>& _m,
        const PackedQ<T2>& _q) : x(_x), m(_m), q(_q)
    { TMVAssert( m.colsize() == q.colsize() ); }
    inline size_t colsize() const { return q.rowsize(); }
    inline size_t rowsize() const { return m.rowsize(); }
    inline StorageType stor() const { return ColMajor; }
    inline T GetX() const { return x; }
    inline const GenMatrix<T1>& GetM() const { return m; }
    inline const GenMatrix<T2>& GetQ() const { return q; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      TMVAssert(IsReal(T()));
      q.LDiv(m,m0);
      MultXM(x,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      q.LDiv(m,m0);
      MultXM(x,m0);
    }
  protected:
    const T x;
    const GenMatrix<T1>& m;
    const PackedQ<T2>& q;
  };


  template <class T, class T1, class T2> 
  class RQuotMpQ :
    public MatrixComposite<T>
  {
  public:
    inline RQuotMpQ(const T _x, const GenMatrix<T1>& _m,
        const PackedQ<T2>& _q) : x(_x), m(_m), q(_q)
    { TMVAssert( m.rowsize() == q.rowsize() ); }
    inline size_t colsize() const { return m.colsize(); }
    inline size_t rowsize() const { return q.colsize(); }
    inline StorageType stor() const { return ColMajor; }
    inline T GetX() const { return x; }
    inline const GenMatrix<T1>& GetM() const { return m; }
    inline const GenMatrix<T2>& GetQ() const { return q; }
    inline void AssignToM(const MatrixView<RealType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      TMVAssert(IsReal(T()));
      q.RDiv(m,m0);
      MultXM(x,m0);
    }
    inline void AssignToM(const MatrixView<ComplexType(T)>& m0) const
    {
      TMVAssert(m0.colsize() == colsize() && m0.rowsize() == rowsize());
      q.RDiv(m,m0);
      MultXM(x,m0);
    }
  protected:
    const T x;
    const GenMatrix<T1>& m;
    const PackedQ<T2>& q;
  };

  template <class T> 
  inline const MatrixView<T>& operator/=(
      const MatrixView<T>& m, const PackedQ<T>& q)
  {
    TMVAssert(q.colsize() == q.rowsize());
    TMVAssert(m.colsize() == q.rowsize());
    q.LDivEq(m);
    return m;
  }

  template <class T> 
  inline const MatrixView<CT>& operator/=(
      const MatrixView<CT>& m, const PackedQ<T>& q)
  {
    TMVAssert(q.colsize() == q.rowsize());
    TMVAssert(m.colsize() == q.rowsize());
    q.LDivEq(m);
    return m;
  }

  template <class T> 
  inline const MatrixView<T>& operator%=(
      const MatrixView<T>& m, const PackedQ<T>& q)
  {
    TMVAssert(q.colsize() == q.rowsize());
    TMVAssert(m.rowsize() == q.rowsize());
    q.RDivEq(m);
    return m;
  }

  template <class T> 
  inline const MatrixView<CT>& operator%=(
      const MatrixView<CT>& m, const PackedQ<T>& q)
  {
    TMVAssert(q.colsize() == q.rowsize());
    TMVAssert(m.rowsize() == q.rowsize());
    q.RDivEq(m);
    return m;
  }


#define GENMATRIX1 GenMatrix
#define GENMATRIX2 PackedQ
#define PRODXM1 ProdXM
#define PRODXM2 ProdXpQ
#define QUOTMM QuotMpQ
#define RQUOTMM RQuotMpQ
#include "tmv/TMV_AuxQuotMM.h"
#include "tmv/TMV_AuxQuotMMa.h"
#undef GENMATRIX1
#undef GENMATRIX2
#undef PRODXM1
#undef PRODXM2
#undef QUOTXM
#undef QUOTMM
#undef RQUOTMM

#undef CT

} // namespace mv


#endif
