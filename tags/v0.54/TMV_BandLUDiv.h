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
      const BandMatrixView<T>& LUx, size_t* p, T& det, int Anhi);
  // Decompose A (input as LUx) into L * U
  
  template <class T, class Ta> void BandTriLDivEq(
      const GenBandMatrix<Ta>& A, const MatrixView<T>& B, DiagType dt);
  template <class T, class Ta> void BandTriLDivEq(
      const GenBandMatrix<Ta>& A, const VectorView<T>& v, DiagType dt);
  // Solve A X = B  where A is upper or lower band triangular

#define CT complex<T>
  template <class T> inline void BandTriLDivEq(
      const GenBandMatrix<CT>& , const VectorView<T>& , DiagType )
  { TMVAssert(FALSE); }
  template <class T> inline void BandTriLDivEq(
      const GenBandMatrix<CT>& , const MatrixView<T>& , DiagType )
  { TMVAssert(FALSE); }
#undef CT

  template <class T, class T1> void BandLU_LDivEq(
      const GenBandMatrix<T1>& LUx, const size_t* p, const MatrixView<T>& m);
  template <class T, class T1> void BandLU_RDivEq(
      const GenBandMatrix<T1>& LUx, const size_t* p, const MatrixView<T>& m);

#define CT complex<T>
  template <class T> inline void BandLU_LDivEq(const GenBandMatrix<CT>& ,
      const size_t* , const MatrixView<T>& ) 
  { TMVAssert(FALSE); }
  template <class T> inline void BandLU_RDivEq(const GenBandMatrix<CT>& ,
      const size_t* , const MatrixView<T>& ) 
  { TMVAssert(FALSE); }
#undef CT

  template <class T> class SingularBandLU :
    public Singular
    {
      public:

	const GenBandMatrix<T>& A;

	SingularBandLU(const GenBandMatrix<T>& _A) : 
	  Singular("BandMatrix"), A(_A) {}
	~SingularBandLU() throw() {}

	void Write(ostream& os) const throw()
	{
	  Singular::Write(os);
	  os<<"In LU Decomposed form, the matrix is \n"<<A<<endl;
	  os<<"ie. U = "<<A.UpperBand()<<endl;
	}
    };

  template <class T> class SingularBandLU2 :
    public SingularBandLU<T>
    {
      public:

	const GenBandMatrix<T>& A0;

	SingularBandLU2(const GenBandMatrix<T>& _A,
	    const GenBandMatrix<T>& _A0) : 
	  SingularBandLU<T>(_A), A0(_A0) {}
	~SingularBandLU2() throw() {}

	void Write(ostream& os) const throw() 
	{
	  SingularBandLU<T>::Write(os);
	  os<<"The original BandMatrix is \n"<<A0<<endl;
	}
    };

  template <class T> class BandLUDiv : 
    public Divider<T> 
  {

    public :

      //
      // Constructors
      //

      BandLUDiv(const GenBandMatrix<T>& A, bool _inplace);
      ~BandLUDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> inline void DoLDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(m.colsize() == LUx.colsize());
	if (istrans) BandLU_RDivEq(LUx,P.get(),m.Transpose());
	else BandLU_LDivEq(LUx,P.get(),m);
      }

      template <class T1> inline void DoRDivEq(const MatrixView<T1>& m) const 
      {
	TMVAssert(m.rowsize() == LUx.colsize());
	if (istrans) BandLU_LDivEq(LUx,P.get(),m.Transpose());
	else BandLU_RDivEq(LUx,P.get(),m);
      }

      template <class T1, class T2> inline void DoLDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
      { 
	TMVAssert(m.rowsize() == x.rowsize());
	TMVAssert(m.colsize() == LUx.colsize());
	TMVAssert(x.colsize() == LUx.colsize());
	if (istrans) BandLU_RDivEq(LUx,P.get(),(x=m).Transpose()); 
	else BandLU_LDivEq(LUx,P.get(),x=m); 
      }

      template <class T1, class T2> inline void DoRDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
      { 
	TMVAssert(m.colsize() == x.colsize());
	TMVAssert(m.rowsize() == LUx.colsize());
	TMVAssert(x.rowsize() == LUx.colsize());
	if (istrans) BandLU_LDivEq(LUx,P.get(),(x=m).Transpose()); 
	else BandLU_RDivEq(LUx,P.get(),x=m); 
      }

#include "TMV_AuxAllDiv.h"

      //
      // Determinant, Inverse
      //

      T Det() const;
      void Inverse(const MatrixView<T>& minv) const;
      void InverseATA(const MatrixView<T>& minv) const;
      inline bool Singular() const { return Det() == T(0); }

      //
      // Access Decomposition
      //

      inline bool IsTrans() const { return istrans; }
      inline ConstBandMatrixView<T> GetU() const 
      { return BandMatrixViewOf(LUx,0,LUx.nhi()); }
      inline LowerTriMatrix<T,UnitDiag> GetL() const;
      inline const GenBandMatrix<T>& GetLU() const { return LUx; }
      const size_t* GetP() const;

      inline std::string Type() const
      { return std::string("BandLUDiv<") + tmv::Type(T()) + ">"; }
      inline DivType GetDivType() const { return LU; }

      bool CheckDecomp(const BaseMatrix<T>& m, ostream* fout) const;

    private :
      const bool istrans;
      const bool inplace;
      auto_array<T> Aptr1;
      T* Aptr;
      BandMatrixView<T> LUx;
      auto_array<size_t> P;
      mutable T det;
      mutable bool donedet;
  };

} // namespace mv

#endif
