//---------------------------------------------------------------------------
//
// This file contains the code for doing division of BandMatrices using 
// LU Decomposition.
//
// The basics of LU decomposition for band matrices are the same as 
// for regular matrices.  However, there are a few wrinkles about doing
// it efficiently.  
//
// We leave the details to the comments in TMV_BandLUDiv.cpp, but 
// the main difference for the routines in this file is that L can
// be stored in a lower band matrix with m.nlo() subdiagonals.  
// However, U needs m.nlo() + m.nhi() superdiagonals for its storage.
//
//


#ifndef TMV_BandLUDiv_H
#define TMV_BandLUDiv_H

#include "TMV_BandMatrix.h"
#include "TMV_Divider.h"
#include "TMV_LUDiv.h"

namespace tmv {

  template <class T> void BandLU_Decompose(
      const BandMatrixView<T>& LUx, vector<size_t>& p, T& det);
  // Decompose A (input as LUx) into L * U
  
  template <class T1, class T2> void BandTriLDivEq(
      const GenBandMatrix<T1>& A, const MatrixView<T2>& B, DiagType dt);
  // Solve A X = B  where A is upper or lower band triangular

#define CT complex<T>
  template <class T> inline void BandTriLDivEq(
      const GenBandMatrix<CT>& A, const MatrixView<T>& B, DiagType dt)
  { TMVAssert(false); }
#undef CT

  template <class T1, class T2> void BandLU_LDivEq(
      const GenBandMatrix<T1>& LUx,
      const vector<size_t>& p, const MatrixView<T2>& m);
  template <class T1, class T2> void BandLU_RDivEq(
      const GenBandMatrix<T1>& LUx,
      const vector<size_t>& p, const MatrixView<T2>& m);

#define CT complex<T>
  template <class T> inline void BandLU_LDivEq(const GenBandMatrix<CT>& LUx,
      const vector<size_t>& p, const MatrixView<T>& m) 
  { TMVAssert(false); }
  template <class T> inline void BandLU_RDivEq(const GenBandMatrix<CT>& LUx,
      const vector<size_t>& p, const MatrixView<T>& m) 
  { TMVAssert(false); }
#undef CT

  template <class T> class BandLUDiv : 
    virtual public BaseLUDiv<T> 
  {

    public :

      //
      // Constructors
      //

      BandLUDiv(const GenBandMatrix<T>& A) :
	istrans(A.nhi()<A.nlo()),
	LUx(A.colsize(),A.colsize(), min(A.nhi(),A.nlo()), 
	    min(A.nlo()+A.nhi(),int(A.colsize())-1),
	    istrans?TransOf(BaseStorOf(A)):BaseStorOf(A)),
	p(A.colsize()), det(T(1)), donedet(false)
      {
	TMVAssert(A.IsSquare());
	if (istrans) {
	  BandMatrixViewOf(LUx.QuickTranspose(),A.nlo(),A.nhi()) = A;
	  for(int i=A.nlo()+1;i<=LUx.nhi();++i) LUx.diag(i).Zero();
	} else {
	  BandMatrixViewOf(LUx,A.nlo(),A.nhi()) = A;
	  for(int i=A.nhi()+1;i<=LUx.nhi();++i) LUx.diag(i).Zero();
	}
	BandLU_Decompose(LUx.QuickView(),p,det);
      }

      ~BandLUDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(m.colsize() == LUx.colsize());
	if (istrans) BandLU_RDivEq(LUx,p,m.QuickTranspose());
	else BandLU_LDivEq(LUx,p,m);
      }

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(m.rowsize() == LUx.colsize());
	if (istrans) BandLU_LDivEq(LUx,p,m.QuickTranspose());
	else BandLU_RDivEq(LUx,p,m);
      }

      template <class T1, class T2> inline void DoLDiv(
	    const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
	{ 
	  TMVAssert(m.rowsize() == x.rowsize());
	  TMVAssert(m.colsize() == LUx.colsize());
	  TMVAssert(x.colsize() == LUx.colsize());
	  if (istrans) BandLU_RDivEq(LUx,p,(x=m).QuickTranspose()); 
	  else BandLU_LDivEq(LUx,p,x=m); 
	}

      template <class T1, class T2> inline void DoRDiv(
	    const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
	{ 
	  TMVAssert(m.colsize() == x.colsize());
	  TMVAssert(m.rowsize() == LUx.colsize());
	  TMVAssert(x.rowsize() == LUx.colsize());
	  if (istrans) BandLU_LDivEq(LUx,p,(x=m).QuickTranspose()); 
	  else BandLU_RDivEq(LUx,p,x=m); 
	}

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      Matrix<T,ColMajor> Inverse() const;
      Matrix<T,ColMajor> InverseATA() const;
      inline bool Singular() const { return Det() == T(0); }

      //
      // Access Decomposition
      //

      Matrix<T,ColMajor> LU_GetL() const;
      Matrix<T,ColMajor> LU_GetU() const;
      Permutation LU_GetP() const;
      Permutation LU_GetQ() const;

      inline ConstBandMatrixView<T> GetBandU() const 
      { return BandMatrixViewOf(LUx,0,LUx.nhi()); }
      const BandMatrix<T,ColMajor>& GetBandLU() const { return LUx; }
      const vector<size_t>& Getp() const { return p; }

      inline std::string Type() const
      { return std::string("BandLUDiv<") + tmv::Type(T()) + ">"; }

    private :
      const bool istrans;
      BandMatrix<T,ColMajor> LUx;
      vector<size_t> p;
      mutable T det;
      mutable bool donedet;
  };

} // namespace mv

#endif
