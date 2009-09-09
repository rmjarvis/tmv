//---------------------------------------------------------------------------
//
// This file defines the TMV BlockDiagMatrix class.
//
// The BlockDiagMatrix class is provided for efficient storage of a block
// diagonal matrix.  You can do most of the things that you can do with a 
// regular Matrix, but it will do them more efficiently.
//
//
// Constructors:
//
//    BlockDiagMatrix<T>(vector<size_t> sizes)
//        Make an uninitialized BlockDiagMatrix with each block
//        a square matrix with size sizes[i]
//
//    BlockDiagMatrix<T>(vector<Matrix<T>*> blocks)
//        Make a BlockDiagMatrix with blocks copied from *blocks[i]
//
//    BlockDiagSparseMatrix<T>(vector<const BaseMatrix<T>*> blocks)
//        Make a ConstSubBlockDiagMatrix with each block being the 
//        (potentially, but not necessarily sparse) matrix blocks[i]
//        These matrices are assumed to be created with new, and they 
//        will be deleted when the BlockDiagSparseMatrix is deleted.
//        Also, they are all const, so functions which modify matrices
//        are invalid for BlockDiagSparseMatrices.  Modifying the 
//        constituent blocks after creation in here may cause 
//        unexpected behavior.
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//        Return the dimensions of the BlockDiagMatrix
//
//    T& operator()(size_t i, size_t j)
//    T operator()(size_t i, size_t j) const
//        Return the (i,j) element of the BlockDiagMatrix
//
//    bool IsSquare()
//        Always true.
//
//    Matrix<T>* GetBlock(size_t i)
//        Returns the ith block of the BlockDiagMatrix.
//
//
// Modifying Functions
//
//    BlockDiagMatrix& Zero()
//        Sets all elements to 0
//
//    BlockDiagMatrix<T>& TransposeSelf() 
//        Transposes the matrix in place.
//
//    BlockDiagMatrix& ConjugateSelf()
//        Sets all elements to its conjugate
//
//    BlockDiagMatrix& SetToIdentity(x = 1)
//        Set to Identity BlockDiagMatrix, or 
//        with a parameter, set to x times Identity BlockDiagMatrix
//
//    BlockDiagMatrix& SwapRows(size_t i1, size_t i2)
//    BlockDiagMatrix& SwapCols(size_t j1, size_t j2)
//    BlockDiagMatrix& SwapRowsCols(size_t i1, size_t i2)
//        Swaps rows/cols as for a normal matrix.  
//        However, the elements must be in the same block.
//
//    BlockDiagMatrix& SwapWith(BlockDiagMatrix& m2)
//        Swap values with those in matrix m2
//        The Matrices must be the same size and have the same block structure.
//
//    void Swap(BlockDiagMatrix& m1, BlockDiagMatrix& m2)
//        Swap the values of two Matrices
//        The Matrices must be the same size and have the same block structure.
//
//
// Functions of BlockDiagMatrixs:
//        (These are all both member functions and functions of a BlockDiagMatrix,
//         so Norm(m) and m.Norm() for example are equivalent.)
//
//    Det(m)
//        Returns the determinant of a BlockDiagMatrix.
//
//    Trace(m)
//        Returns the trace of a BlockDiagMatrix.
//        = sum_i ( a_ii )
//
//    Norm(m) or NormF(m)
//        Return the Frobenius norm of a BlockDiagMatrix.
//        = sqrt( sum_ij |a_ij|^2 )
//
//    NormSq(m)
//        Returns the square of Norm().
//
//    Norm1(m) 
//        Returns the 1-norm of a BlockDiagMatrix.
//        = max_j (sum_i |a_ij|)
//
//    Norm2(m) 
//        Returns the 2-norm of a BlockDiagMatrix.
//        = sqrt( Max Eigenvalue of (A.Dagger * A) )
//        = Maximum singular value
//        (Note - for diagonal matrices, Norm1 = Norm2 = NormInf.)
//
//    NormInf(m) 
//        Returns the infinity-norm of a BlockDiagMatrix.
//        = max_i (sum_j |a_ij|)
//
//    Transpose(m)
//        Returns the transpose of a BlockDiagMatrix
//
//    Conjugate(m)
//        Returns the conjugate of a BlockDiagMatrix
//
//    Adjoint(m)
//        Returns the conjugate of the transpose of a BlockDiagMatrix
//        Note: Some people define the adjoint as the cofactor matrix.
//              This is not the same as our definition of the Adjoint.
//
//    Dagger(m)
//        An alias for Adjoint(m)
//
//    Inverse(m)
//        Returns the inverse of m if it exists.
//        If m is singular, then a run-time error will occur.
//
//
// I/O: 
//
//    os << k 
//        Writes k to ostream os as a full matrix
//
//

#ifndef TMV_BlockDiagMatrix_H
#define TMV_BlockDiagMatrix_H

#include "TMV_BaseSparseMatrix.h"

namespace jmv {

  template <class T> class GenBlockDiagMatrix;
  template <class T> class ConstSubBlockDiagMatrix;
  template <class T> class ModSubBlockDiagMatrix;
  template <class T> class BlockDiagMatrix;
  template <class T> class BlockDiagMatrixComposite;

  template <class T> class BlockDiagLUDiv;
  template <class T> class BlockDiagQRDiv;
  template <class T> class BlockDiagQRPDiv;
  template <class T> class BlockDiagSVDiv;

  template <class T> class GenBlockDiagMatrix : public BaseSparseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      GenBlockDiagMatrix() : BaseSparseMatrix<T>(jmv::LU) {}
      GenBlockDiagMatrix(const GenBlockDiagMatrix<T>& rhs) : BaseSparseMatrix<T>(rhs.itsdt) {}
      ~GenBlockDiagMatrix() {}

      //
      // Access Functions
      //

      inline size_t size() const { return starts().back(); }
      inline size_t colsize() const { return size(); }
      inline size_t rowsize() const { return size(); }

      inline bool IsSquare() const { return true; }

      //
      // Functions of BlockDiagMatrix
      //

      T Det() const
      {
	T det = T(1);
	for(size_t i=0;i<cblocks().size();++i) {
	  det *= cblocks()[i]->Det();
	}
	return det;
      }

      T Trace() const
      { 
	T trace = T(0);
	for(size_t i=0;i<cblocks().size();++i) {
	  trace += cblocks()[i]->Trace();
	}
	return trace;
      }

      // NormF()^2
      RealType(T) NormSq() const 
      {
	RealType(T) normsq = T(0);
	for(size_t i=0;i<cblocks().size();++i) {
	  normsq += cblocks()[i]->NormSq();
	}
	return normsq;
      }

      // 1-Norm = max_j (sum_i |a_ij|)
      RealType(T) Norm1() const 
      {
	RealType(T) norm = T(0);
	for(size_t i=0;i<cblocks().size();++i) {
	  RealType(T) subnorm = cblocks()[i]->Norm1();
	  if (subnorm > norm) norm = subnorm;
	}
	return norm;
      }

      // 2-Norm = sqrt(max_eigenvalue(A At))
      RealType(T) Norm2() const 
      {
	RealType(T) norm = T(0);
	for(size_t i=0;i<cblocks().size();++i) {
	  RealType(T) subnorm = cblocks()[i]->Norm2();
	  if (subnorm > norm) norm = subnorm;
	}
	return norm;
      }

      // inf-Norm = max_i (sum_j |a_ij|)
      RealType(T) NormInf() const
      {
	RealType(T) norm = T(0);
	for(size_t i=0;i<cblocks()->size();++i) {
	  RealType(T) subnorm = cblocks()[i]->NormInf();
	  if (subnorm > norm) norm = subnorm;
	}
	return norm;
      }

      ConstSubBlockDiagMatrix<T> Transpose() const
      {
	vector<BaseMatrix<T>*> newblocks;
	for(size_t i=0;i<cblocks().size();++i) {
	  newblocks.push_back(cblocks()[i]->NewConstTranspose());
	}
	return ConstSubBlockDiagMatrix<T>(newblocks);
      }

      BaseMatrix<T>* NewConstTranspose() const
      { return new ConstSubBlockDiagMatrix<T>(Transpose()); }

      ConstSubBlockDiagMatrix<T> Conjugate() const
      {
	vector<BaseMatrix<T>*> newblocks;
	for(size_t i=0;i<cblocks().size();++i) {
	  newblocks.push_back(cblocks()[i]->NewConjugate());
	}
	return ConstSubBlockDiagMatrix<T>(newblocks);
      }

      BaseMatrix<T>* NewConjugate() const
      { return new ConstSubBlockDiagMatrix<T>(Conjugate()); }

      ConstSubBlockDiagMatrix<T> Adjoint() const
      {
	vector<BaseMatrix<T>*> newblocks;
	for(size_t i=0;i<cblocks().size();++i) {
	  newblocks.push_back(cblocks()[i]->NewAdjoint());
	}
	return ConstSubBlockDiagMatrix<T>(newblocks);
      }

      BaseMatrix<T>* NewAdjoint() const
      { return new ConstSubBlockDiagMatrix<T>(Adjoint()); }

      ConstSubBlockDiagMatrix<T> Dagger() const { return Adjoint(); } 

      BlockDiagMatrix<T> Inverse() const 
      {
	vector<Matrix<T>*> newblocks;
	for(size_t i=0;i<blocks().size();++i) {
	  newblocks.push_back(blocks()[i].NewInverse());
	}
	return BlockDiagMatrix<T>(newblocks);
      }

      BaseMatrix<T>* NewInverse() const
      { return new BlockDiagMatrix<T>(Inverse()); }

      BaseMatrix<T>* NewView() const
      { return new ConstSubBlockDiagMatrix<T>(blocks()); }

      //
      // Arithmetic Helpers
      //

      // Note that K in these functions is the letter for "Block".
      // B was already taken for Band Matrices.
      // In this case, KD means BlockDiag, since one can also 
      // have a BlockSparse matrix which would be KS.
      
      template <class T1, class T2> friend inline void AddEqMKD(
	  const ModSubMatrix<T1>& m, const GenBlockDiagMatrix<T2>& d);
      void AddEq(const ModSubMatrix<RealType(T)>& m) const
      { AddEqMKD(m,*this); }
      void AddEq(const ModSubMatrix<ComplexType(T)>& m) const 
      { AddEqMKD(m,*this); }
      void Add(const GenMatrix<RealType(T)>& m1, 
	  const ModSubMatrix<T>& m0) const 
      { AddEqMKD(m0=m1,*this); }
      void Add(const GenMatrix<ComplexType(T)>& m1, 
	  const ModSubMatrix<ComplexType(T)>& m0) const 
      { AddEqMKD(m0=m1,*this); }

      template <class T1, class T2> friend inline void SubtractEqMKD(
	  const ModSubMatrix<T1>& m, const GenBlockDiagMatrix<T2>& d);
      void SubtractEq(const ModSubMatrix<RealType(T)>& m) const 
      { SubtractEqMKD(m,*this); }
      void SubtractEq(const ModSubMatrix<ComplexType(T)>& m) const 
      { SubtractEqMKD(m,*this); }
      void Subtract(const GenMatrix<RealType(T)>& m1, 
	  const ModSubMatrix<T>& m0) const 
      { SubtractEqMKD(m0=m1,*this); }
      void Subtract(const GenMatrix<ComplexType(T)>& m1, 
	  const ModSubMatrix<ComplexType(T)>& m0) const 
      { SubtractEqMKD(m0=m1,*this); }

      template <class T1, class T2> friend inline void MultEqKDV(
	  const GenBlockDiagMatrix<T2>& d, const ModSubVector<T1>& v);
      void MultEq(bool trans, const ModSubVector<RealType(T)>& v) const 
      { MultEqKDV(*this,v); }
      void MultEq(bool trans, const ModSubVector<ComplexType(T)>& v) const 
      { MultEqKDV(*this,v); }
      void Mult(bool trans, const GenVector<RealType(T)>& v1, 
	  const ModSubVector<T>& v0) const 
      { MultKDV(*this,v0=v1); }
      void Mult(bool trans, const GenVector<ComplexType(T)>& v1, 
	  const ModSubVector<ComplexType(T)>& v0) const 
      { MultKDV(*this,v0=v1); }

      template <class T1, class T2> friend inline void DivEqKDV(
	  const GenBlockDiagMatrix<T2>& d, const ModSubVector<T1>& v);
      void LDivEq(const ModSubVector<RealType(T)>& v) const 
      { DivEqKDV(*this,v); }
      void LDivEq(const ModSubVector<ComplexType(T)>& v) const 
      { DivEqKDV(*this,v); }
      void LDiv(const GenVector<RealType(T)>& v1, 
	  const ModSubVector<T>& v0) const 
      { DivEqKDV(*this,v0=v1); }
      void LDiv(const GenVector<ComplexType(T)>& v1, 
	  const ModSubVector<ComplexType(T)>& v0) const 
      { DivEqKDV(*this,v0=v1); }

      template <class T1, class T2> friend inline void DivEqVKD(
	  const ModSubVector<T1>& v, const GenBlockDiagMatrix<T2>& d);
      void RDivEq(const ModSubVector<RealType(T)>& v) const 
      { DivEqVKD(v,*this); }
      void RDivEq(const ModSubVector<ComplexType(T)>& v) const 
      { DivEqVKD(v,*this); }
      void RDiv(const GenVector<RealType(T)>& v1, 
	  const ModSubVector<T>& v0) const 
      { DivEqVKD(v0=v1,*this); }
      void RDiv(const GenVector<ComplexType(T)>& v1, 
	  const ModSubVector<ComplexType(T)>& v0) const 
      { DivEqVKD(v0=v1,*this); }
      

      virtual const vector<BaseMatrix<T>*>& cblocks() const = 0;

    protected :

      const vector<size_t>& starts() const
      {
	static vector<size_t> itsstarts;
	if (itsstarts.size() == 0) {
	  const vector<BaseMa5trix<T>*>& cb = cblocks();
	  size_t k=0;
	  itsstarts.push_back(k);
	  for(size_t i=0;i<cb.size();++i) {
	    TMVAssert(cb[i]);
	    TMVAssert(cb[i]->colsize() > 0);
	    TMVAssert(cb[i]->colsize() == cb[i]->rowsize());
	    k += cb[i]->colsize();
	    itsstarts.push_back(k);
	  }
	}
	return itsstarts;
      }

      void GetBlockNum(size_t i, size_t& n, size_t& blocki) const
      {
	TMVAssert(cblocks().size() > 0);
	TMVAssert(i < size());
	vector<size_t>::const_iterator mit = 
	  upper_bound(starts().begin(),starts().end(),i);
	TMVAssert(mit == starts().end() || *mit > i);
	TMVAssert(mit != starts().begin());
	mit--;
	TMVAssert(*mit <= i);
	n = mit - starts.begin();
	blocki = i-*mit;
      }

      T cref(size_t i, size_t j) const
      {
	size_t bi, bj;
	size_t n, n2;
	GetBlockNum(i,n,bi);
	GetBlockNum(j,n2,bj);
	if (n == n2) return (*(cblocks()[n]))(bi,bj);
	else return T(0);
      }

    private :

      void operator=(const GenBlockDiagMatrix<T>&) { TMVAssert(false); }

  }; // GenBlockDiagMatrix

  template <class T1, class T2> inline void AddEqMKD(
      const ModSubMatrix<T1>& m, const BlockDiagMatrix<T2>& d)
  { 
    TMVAssert(m.colsize() == d.colsize());
    TMVAssert(m.rowsize() == d.rowsize());
    size_t k0 = 0;
    for(size_t k=0;k<d.cblocks().size();++k) {
      size_t subn = d.cblocks()[k]->colsize();
      size_t k1 = k0 + subn;
      m.SubMatrix(k0,k1,k0,k1) += *d.cblocks()[k];
      k0 = k1;
    }
  }
  template <class T> inline void AddEqMKD(
      const ModSubMatrix<T>& m, const BlockDiagMatrix<complex<T> >& d)
  { TMVAssert(false); }

  template <class T1, class T2> inline void SubtractEqMKD(
      const ModSubMatrix<T1>& m, const BlockDiagMatrix<T2>& d)
  { 
    TMVAssert(m.colsize() == d.colsize());
    TMVAssert(m.rowsize() == d.rowsize());
    size_t k0 = 0;
    for(size_t k=0;k<d.cblocks().size();++k) {
      size_t subn = d.cblocks()[k]->colsize();
      size_t k1 = k0 + subn;
      m.SubMatrix(k0,k1,k0,k1) -= *d.cblocks()[k];
      k0 = k1;
    }
  }
  template <class T> inline void SubtractEqMKD(
      const ModSubMatrix<T>& m, const BlockDiagMatrix<complex<T> >& d)
  { TMVAssert(false); }

  template <class T1, class T2> inline void MultEqVKD(
      const ModSubVector<T1>& v, const BlockDiagMatrix<T2>& d)
  { 
    TMVAssert(v.size() == d.colsize());
    size_t k0 = 0;
    for(size_t k=0;k<d.cblocks().size();++k) {
      size_t subn = d.cblocks()[k]->colsize();
      size_t k1 = k0 + subn;
      v.SubVector(k0,k1) *= *d.cblocks()[k];
      k0 = k1;
    }
  }
  template <class T> inline void MultEqVKD(
      const ModSubVector<T>& v, const BlockDiagMatrix<complex<T> >& d)
  { TMVAssert(false); }

  template <class T1, class T2> inline void MultEqKDV(
      const BlockDiagMatrix<T2>& d, const ModSubVector<T3>& v) 
  { 
    TMVAssert(v.size() == d.rowsize());
    size_t k0 = 0;
    for(size_t k=0;k<d.cblocks().size();++k) {
      size_t subn = d.cblocks()[k]->colsize();
      size_t k1 = k0 + subn;
      ModSubVector subv = v.SubVector(k0,k1);
      subv = *d.cblocks()[k] * subv;
      k0 = k1;
    }
  }
  template <class T> inline void MultEqKDV(
      const BlockDiagMatrix<complex<T> >& d, const ModSubVector<T>& v)
  { TMVAssert(false); }

  template <class T1, class T2> inline void DivEqVKD(
      const ModSubVector<T1>& v, const BlockDiagMatrix<T2>& d)
  { 
    TMVAssert(v.size() == d.colsize());
    size_t k0 = 0;
    for(size_t k=0;k<d.cblocks().size();++k) {
      size_t subn = d.cblocks()[k]->colsize();
      size_t k1 = k0 + subn;
      v.SubVector(k0,k1) /= *d.cblocks()[k];
      k0 = k1;
    }
  }
  template <class T> inline void DivEqVKD(
      const ModSubVector<T>& v, const BlockDiagMatrix<complex<T> >& d)
  { TMVAssert(false); }

  template <class T1, class T2, class T3> inline void DivEqKDV(
      const GenVector<T1>& v1, const BlockDiagMatrix<T2>& d,
      const ModSubVector<T3>& v0)
  { 
    TMVAssert(v.size() == d.rowsize());
    size_t k0 = 0;
    for(size_t k=0;k<d.cblocks().size();++k) {
      size_t subn = d.cblocks()[k]->colsize();
      size_t k1 = k0 + subn;
      v.SubVector(k0,k1) %= *d.cblocks()[k];
      k0 = k1;
    }
  }
  template <class T> inline void DivEqKDV(
      const BlockDiagMatrix<complex<T> >& d, const ModSubVector<T>& v)
  { TMVAssert(false); }

  template <class T> class ConstSubBlockDiagMatrix : public GenBlockDiagMatrix<T>
  {
    public:

      //
      // Constructors
      //
      
      ConstSubBlockDiagMatrix(const ConstSubBlockDiagMatrix<T>& rhs)
	: itsblocks(rhs.itsblocks.size())
      {
	for(size_t k=0;k<itsblocks.size();++k)
	  itsblocks[k] = rhs.itsblocks[k]->NewView();
      }
      ConstSubBlockDiagMatrix(const GenBlockDiagMatrix<T>& rhs)
	: itsblocks(rhs.cblocks().size())
      {
	for(size_t k=0;k<itsblocks.size();++k)
	  itsblocks[k] = rhs.cblocks()[k]->NewView();
      }
      ConstSubBlockDiagMatrix(const vector<BaseMatrix<T>*> _blocks)
	: itsblocks(_blocks()) {}

      ~ConstSubBlockDiagMatrix()
      { for(size_t k=0;k<itsblocks.size();++k) delete itsblocks[k]; }

    protected :

      vector<BaseMatrix*> itsblocks;

      virtual const vector<BaseMatrix<T>*>& cblocks() const
      { return itsblocks; }

    private :

      void operator=(const ConstSubBlockDiagMatrix<T>&) { TMVAssert(false); }

  }; // ConstSubBlockDiagMatrix

  // This name is a big more intuitive for making a block-diag sparse matrix.
  template <class T> inline ConstSubBlockDiagMatrix<T> BlockDiagSparseMatrix(
      const vector<BaseMatrix<T>*> blocks)
  { return ConstSubBlockDiagMatrix<T>(blocks); }

  template <class T, class IT> class MosSubBlockDiagMatrix : public GenBlockDiagMatrix<T>
  {
    public:

      //
      // Constructors
      //

      ModSubBlockDiagMatrix(const ModSubBlockDiagMatrix<T,IT>& rhs)
	: itsblocks(rhs.itsblocks.size()) 
	{ 
	  for(size_t k=0;k<itsblocks.size();++k) 
	    itsblocks[k] = new ModSubMatrix<T,IT>(rhs.itsblocks[k]);
	}

      ModSubBlockDiagMatrix(const vector<ModSubMatrix<T,IT>*> _blocks)
	: itsblocks(_blocks()) {}

      ~ModSubBlockDiagMatrix()
      { for(size_t i=0;i<itsblocks.size();++i) delete itsblocks[i]; }

    private:
      template <class T2> inline void Copy(const GenBlockDiagMatrix<T2>& rhs) const
      { 
	TMVAssert(itsblocks.size() == rhs.cblocks().size());
	for(size_t k=0;k<itsblocks.size();++k) {
	  TMVAssert(itsblocks[k]->colsize() == rhs.cblocks()[k]->colsize());
	  TMVAssert(itsblocks[k]->IsSquare());
	  TMVAssert(rhs.cblocks()[k]->IsSquare());
	  *itsblocks[k] = *rhs.cblocks()[k];
      }
    protected :

      vector<const ModSubMatrix*> itsblocks;

      virtual const vector<BaseMatrix<T>*>& cblocks() const
      { return itsblocks; }

  }; // BlockDiagModSubMatrix

  template <class T> class BlockDiagMatrix : public GenBlockDiagMatrix<T>
  {

    public:

      //
      // Constructors
      //

      BlockDiagMatrix(const vector<Matrix<T>*> _blocks) : blocks(_blocks)
      { 
	for(size_t  k=0;k<blocks.size();++k) {
	  whichblock[totsize] = k;
	  blocks[k] = new Matrix<T>(_blocks[k]);
	  totsize += blocks[k]->colsize();
	}
      }
      explicit BlockDiagMatrix(size_t _size) : itsdiag(_size) {}

      BlockDiagMatrix(size_t _size, T x) : itsdiag(_size,x) {}

      explicit BlockDiagMatrix(const Vector<T>& rhs) : itsdiag(rhs) {}

      BlockDiagMatrix(const BlockDiagMatrix<T>& rhs) : itsdiag(rhs.diag()) {}

      template <class T2> BlockDiagMatrix(const BlockDiagMatrix<T2>& rhs) :
	itsdiag(rhs) {}

      BlockDiagMatrix(const BlockDiagMatrixComposite<T>& dcomp) : itsdiag(dcomp.size())
      { dcomp.AssignTo(View()); }

      ~BlockDiagMatrix() {}


      //
      // Op=
      //

      BlockDiagMatrix<T>& operator=(const BlockDiagMatrix<T>& m2)
      {
	if (&m2 != this) View() = m2.View(); 
	return *this; 
      }

      template <class T2> BlockDiagMatrix<T>& operator=(const BlockDiagMatrix<T2>& m2)
      { View() = m2; return *this; }

      BlockDiagMatrix<T>& operator=(T x) { View() = x; return *this; }

      BlockDiagMatrix<T>& operator=(const BlockDiagMatrixComposite<T>& dcomp)
      { View() = dcomp; return *this; }


      //
      // Access
      //

      ModSubVector<T> diag() { return View().diag(); }

      T& operator()(size_t i) { return diag()(i); }
      T& operator()(size_t i, size_t j) { TMVAssert(i==j); return diag()(i); }

      ConstSubVector<T> diag() const { return cdiag(); }

      T operator()(size_t i) const { return cdiag()(i); }
      T operator()(size_t i,size_t j) const 
      { if (i==j) return cdiag()(i); else return T(0); }


      //
      // Modifying Functions
      //

      BlockDiagMatrix<T>& Zero() { return SetAllTo(0); }

      BlockDiagMatrix<T>& SetAllTo(T x) 
      { View().SetAllTo(x); return *this; }

      BlockDiagMatrix<T>& TransposeSelf() 
      { return *this; }

      BlockDiagMatrix<T>& ConjugateSelf() 
      { View().ConjugateSelf(); return *this; }

      BlockDiagMatrix<T>& SetToIdentity(T x=T(1)) 
      { View().SetToIdentity(x); return *this; }

      BlockDiagMatrix<T>& SwapRowsCols(size_t i1, size_t i2) 
      { View().SwapRowsCols(i1,i2); return *this; }

      BlockDiagMatrix<T>& SwapWith(const ModSubBlockDiagMatrix<T>& m2) 
      { diag().SwapWith(m2.diag()); return *this; }

      BlockDiagMatrix<T>& SwapWith(BlockDiagMatrix<T>& m2) 
      { diag().SwapWith(m2.diag()); return *this; }


      //
      // SubBlockDiagMatrix
      //

      ModSubBlockDiagMatrix<T> SubBlockDiagMatrix(int i1, int i2, int istep=1) 
      { return View().SubBlockDiagMatrix(i1,i2,istep); }

      ConstSubBlockDiagMatrix<T> SubBlockDiagMatrix(int i1, int i2, int istep=1) const
      { return BlockDiagMatrix<T>::SubBlockDiagMatrix(i1,i2,istep); }

      ConstSubBlockDiagMatrix<T> Transpose() const
      { return BlockDiagMatrix<T>::Transpose(); }

      ModSubBlockDiagMatrix<T> Transpose() 
      { return View().Transpose(); }

      ConstSubBlockDiagMatrix<T> View() const
      { return ConstSubBlockDiagMatrix<T>(itsdiag.View()); }

      ModSubBlockDiagMatrix<T> View() 
      { return ModSubBlockDiagMatrix<T>(itsdiag.View()); }


    protected :

      Vector<T> itsdiag;
      ConstSubVector<T> cdiag() const { return itsdiag.View(); }

      friend BlockDiagMatrixComposite<T>;

  }; // BlockDiagMatrix

//---------------------------------------------------------------------------


  //
  // Swap Matrices
  //

  template <class T> inline void Swap(const ModSubBlockDiagMatrix<T>& m1, 
      const ModSubBlockDiagMatrix<T>& m2)
  { m1.SwapWith(m2); }

  template <class T> inline void Swap(const BlockDiagMatrix<T>& m1, 
      const ModSubBlockDiagMatrix<T>& m2)
  { m1.View().SwapWith(m2); }

  template <class T> inline void Swap(const ModSubBlockDiagMatrix<T>& m1, 
      const BlockDiagMatrix<T>& m2)
  { m1.SwapWith(m2.View()); }

  template <class T> inline void Swap(const BlockDiagMatrix<T>& m1, 
      const BlockDiagMatrix<T>& m2)
  { m1.View().SwapWith(m2.View()); }

  //
  // Functions of Matrices:
  //

  template <class T> inline T Det(const BlockDiagMatrix<T>& m)
  { return m.Det(); }

  template <class T> inline T Trace(const BlockDiagMatrix<T>& m)
  { return m.Trace(); }

  template <class T> inline RealType(T) Norm(const BlockDiagMatrix<T>& m)
  { return m.Norm(); }

  template <class T> inline RealType(T) NormSq(const BlockDiagMatrix<T>& m)
  { return m.NormSq(); }
  
  template <class T> inline RealType(T) NormF(const BlockDiagMatrix<T>& m)
  { return m.NormF(); }

  template <class T> inline RealType(T) Norm1(const BlockDiagMatrix<T>& m)
  { return m.Norm1(); }

  template <class T> inline RealType(T) Norm2(const BlockDiagMatrix<T>& m)
  { return m.Norm2(); }

  template <class T> inline RealType(T) NormInf(const BlockDiagMatrix<T>& m)
  { return m.NormInf(); }

  template <class T> inline ConstSubBlockDiagMatrix<T> Transpose(const BlockDiagMatrix<T>& m)
  { return m.Transpose(); }

  template <class T> inline const BlockDiagMatrix<T>& Conjugate(const BlockDiagMatrix<T>& m)
  { return m; }

  template <class T> inline BlockDiagMatrix<complex<T> >& Conjugate(const BlockDiagMatrix<complex<T> >& m)
  { return m.Conjugate(); }

  template <class T> inline ConstSubBlockDiagMatrix<T> Adjoint(const BlockDiagMatrix<T>& m)
  { return m.Transpose(); }

  template <class T> inline BlockDiagMatrix<complex<T> > Adjoint(const BlockDiagMatrix<complex<T> >& m)
  { return m.Adjoint(); }

  template <class T> inline ConstSubBlockDiagMatrix<T> Dagger(const BlockDiagMatrix<T>& m)
  { return m.Transpose(); }

  template <class T> inline BlockDiagMatrix<complex<T> > Dagger(const BlockDiagMatrix<complex<T> >& m)
  { return m.Adjoint(); }

  template <class T> inline BlockDiagMatrix<T> Inverse(const BlockDiagMatrix<T>& m)
  { return m.Inverse(); }

  template <class T> inline BlockDiagMatrix<T> InverseATA(const BlockDiagMatrix<T>& m)
  { return m.InverseATA(); }

  template <class T> inline T SumElements(const BlockDiagMatrix<T>& m)
  { return m.SumElements(); }

  template <class T> inline T MaxElement(const BlockDiagMatrix<T>& m,
      size_t* imax=0, size_t* jmax=0, const vector<bool>* rowmask=0,
      const vector<bool>* colmask=0)
  { return m.MinMaxElement(imax,jmax,rowmask,colmask); }

  template <class T> inline T MinElement(const BlockDiagMatrix<T>& m,
      size_t* imin=0, size_t* jmin=0, const vector<bool>* rowmask=0,
      const vector<bool>* colmask=0)
  { return m.MinMaxElement(imin,jmin,rowmask,colmask); }

  template <class T> inline RealType(T) MaxAbsElement(const BlockDiagMatrix<T>& m,
      size_t* imax=0, size_t* jmax=0, 
      const vector<bool>* rowmask=0, const vector<bool>* colmask=0)
  { return m.MinMaxAbsElement(imax,jmax,rowmask,colmask); }

  template <class T> inline RealType(T) MinAbsElement(const BlockDiagMatrix<T>& m,
      size_t* imin=0, size_t* jmin=0,
      const vector<bool>* rowmask=0, const vector<bool>* colmask=0)
  { return m.MinMaxAbsElement(imin,jmin,rowmask,colmask); }


  //
  // BlockDiagMatrix ==, != BlockDiagMatrix
  //

  template <class T> inline bool operator==(
      const BlockDiagMatrix<T>& m1, const BlockDiagMatrix<T>& m2)
  { return m1.diag() == m2.diag(); }

  template <class T> inline bool operator!=(
      const BlockDiagMatrix<T>& m1, const BlockDiagMatrix<T>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T> inline ostream& operator<<(
      ostream& fout, const BlockDiagMatrix<T>& m)
  { m.Write(fout); return fout;}

}; // namespace jmv

#include "TMV_BlockDiagMatrixArith.h"

#endif
