//---------------------------------------------------------------------------
#ifndef TMV_TriSVDiv_H
#define TMV_TriSVDiv_H

#ifdef TMVDEBUG
#define CHECKTRI
#endif

#include "TMV_TriDivider.h"
#include "TMV_SVDiv.h"

namespace tmv {

  // Everything here is identical to the regular SVDiv
  template <class T> class UpperTriSVDiv : 
    virtual public BaseSVDiv<T>, 
    virtual public UpperTriDivider<T>
  {

    public :

      //
      // Constructors
      //

      UpperTriSVDiv(const GenUpperTriMatrix<T>& A) :
	SS(A), Adt(A.dt()) 
#ifdef CHECKTRI
	  , norma(Norm(A))
#endif
      {}

      ~UpperTriSVDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const 
      { SS.DoLDivEq(m); } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      { SS.DoRDivEq(m); }

      template <class T1, class T2> void DoLDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
      { SS.DoLDiv(m,x); } 

      template <class T1, class T2> void DoRDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
      { SS.DoRDiv(m,x); }

      template <class T1> void DoLDivEq(const UpperTriMatrixView<T1>& m) const
      { 
	TMVAssert(!m.isunit() || Adt == UnitDiag);
	Matrix<T1,RowMajor> mm(m);
	DoLDivEq(mm.QuickView());
	m = UpperTriMatrixViewOf(mm.QuickView(),m.dt());
#ifdef CHECKTRI
	cerr<<"SVUTriLDivEq: "<<Norm(mm-m)<<"  "<<Epsilon<T>()*Norm(m)*norma<<endl;
	TMVAssert(Norm(mm-m)<10*Epsilon<T>()*Norm(m)*norma);
#endif
      }

      template <class T1> void DoRDivEq(const UpperTriMatrixView<T1>& m) const
      { 
	TMVAssert(!m.isunit() || Adt == UnitDiag);
	Matrix<T1,RowMajor> mm(m);
	DoRDivEq(mm.QuickView());
	m = UpperTriMatrixViewOf(mm.QuickView(),m.dt());
#ifdef CHECKTRI
	cerr<<"SVUTriRDivEq: "<<Norm(mm-m)<<"  "<<Epsilon<T>()*Norm(m)*norma<<endl;
	TMVAssert(Norm(mm-m)<10*Epsilon<T>()*Norm(m)*norma);
#endif
      }

      template <class T1, class T2> void DoLDiv(
	  const GenUpperTriMatrix<T1>& m, const UpperTriMatrixView<T2>& x) const
      { 
	TMVAssert(!x.isunit() || (Adt == UnitDiag && m.isunit()));
	Matrix<T1,RowMajor> mm(m);
	Matrix<T2,RowMajor> xx(x);
	DoLDiv(mm,xx.QuickView());
	x = UpperTriMatrixViewOf(xx.QuickView(),x.dt());
#ifdef CHECKTRI
	cerr<<"SVUTriLDiv: "<<Norm(xx-x)<<"  "<<Epsilon<T>()*Norm(x)*norma<<endl;
	TMVAssert(Norm(xx-x)<10*Epsilon<T>()*Norm(x)*norma);
#endif
      }

      template <class T1, class T2> void DoRDiv(
	  const GenUpperTriMatrix<T1>& m, const UpperTriMatrixView<T2>& x) const
      { 
	TMVAssert(!x.isunit() || (Adt == UnitDiag && m.isunit()));
	Matrix<T1,RowMajor> mm(m);
	Matrix<T2,RowMajor> xx(x);
	DoRDiv(mm,xx.QuickView());
	x = UpperTriMatrixViewOf(xx.QuickView(),x.dt());
#ifdef CHECKTRI
	cerr<<"SVUTriRDiv: "<<Norm(xx-x)<<"  "<<Epsilon<T>()*Norm(x)*norma<<endl;
	TMVAssert(Norm(xx-x)<10*Epsilon<T>()*Norm(x)*norma);
#endif
      }

#include "TMV_AuxAllDiv.h"
#include "TMV_AuxUpperTriAllDiv.h"

      inline T Det() const { return SS.Det(); }
      inline Matrix<T,ColMajor> Inverse() const { return SS.Inverse(); }
      inline bool Singular() const { return SS.Singular(); }
      inline RealType(T) Norm2() const { return SS.Norm2(); }
      inline Matrix<T,ColMajor> InverseATA() const { return SS.InverseATA(); }
      inline void Thresh(RealType(T) toler,ostream* debugout=0) const
      { SS.Thresh(toler,debugout); }
      inline void Top(size_t neigen,ostream* debugout=0) const 
      { SS.Top(neigen,debugout); }
      inline size_t GetKMax() const { return SS.GetKMax(); }

      //
      // Access Decomposition
      //

      inline Matrix<T,ColMajor> SV_GetU() const { return SS.SV_GetU(); }
      inline Vector<RealType(T)> SV_GetS() const 
      { return SS.SV_GetS(); }
      inline Matrix<T,ColMajor> SV_GetV() const { return SS.SV_GetV(); }
      inline const Matrix<T,ColMajor>& GetU() const { return SS.GetU(); }
      inline const Vector<RealType(T)>& GetS() const 
      { return SS.GetS(); }
      inline const Matrix<T,ColMajor>& GetV() const { return SS.GetV(); }

      UpperTriMatrix<T,NonUnitDiag,ColMajor> TInverse() const
      {
	Matrix<T,ColMajor> inv = Inverse();
#ifdef CHECKTRI
	cerr<<"SVUTriInv: "<<Norm(inv-UpperTriMatrixViewOf(inv,Adt))<<"  "<<Epsilon<T>()*Norm(inv)*norma<<endl;
	TMVAssert(Norm(inv-UpperTriMatrixViewOf(inv,Adt))<
	    10*Epsilon<T>()*Norm(inv)*norma);
#endif
	return UpperTriMatrix<T,NonUnitDiag,ColMajor>(inv);
      }

      inline std::string Type() const
      { return std::string("UpperTriSVDiv<")+tmv::Type(T())+">"; }

    private :

      SVDiv<T> SS;
      DiagType Adt;
#ifdef CHECKTRI
      RealType(T) norma;
#endif

  }; // UpperTriSVDiv

  template <class T> class LowerTriSVDiv : 
    virtual public BaseSVDiv<T>, 
    virtual public LowerTriDivider<T>
  {

    public :

      //
      // Constructors
      //

      LowerTriSVDiv(const GenLowerTriMatrix<T>& A) :
	SS(A), Adt(A.dt()) 
#ifdef CHECKTRI
	  , norma(Norm(A))
#endif
      {}

      ~LowerTriSVDiv() {}

      //
      // Div, DivEq
      //

      template <class T1> void DoLDivEq(const MatrixView<T1>& m) const 
      { SS.DoLDivEq(m); } 

      template <class T1> void DoRDivEq(const MatrixView<T1>& m) const 
      { SS.DoRDivEq(m); }

      template <class T1, class T2> void DoLDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
      { SS.DoLDiv(m,x); } 

      template <class T1, class T2> void DoRDiv(
	  const GenMatrix<T1>& m, const MatrixView<T2>& x) const 
      { SS.DoRDiv(m,x); }

      template <class T1> void DoLDivEq(const LowerTriMatrixView<T1>& m) const
      { 
	TMVAssert(!m.isunit() || Adt == UnitDiag);
	Matrix<T1,RowMajor> mm(m);
	DoLDivEq(mm.QuickView());
	m = LowerTriMatrixViewOf(mm.QuickView(),m.dt());
#ifdef CHECKTRI
	cerr<<"SVLTriLDivEq: "<<Norm(mm-m)<<"  "<<Epsilon<T>()*Norm(m)*norma<<endl;
	TMVAssert(Norm(mm-m)<10*Epsilon<T>()*Norm(m)*norma);
#endif
      }

      template <class T1> void DoRDivEq(const LowerTriMatrixView<T1>& m) const
      { 
	TMVAssert(!m.isunit() || Adt == UnitDiag);
	Matrix<T1,RowMajor> mm(m);
	DoRDivEq(mm.QuickView());
	m = LowerTriMatrixViewOf(mm.QuickView(),m.dt());
#ifdef CHECKTRI
	cerr<<"SVLTriRDivEq: "<<Norm(mm-m)<<"  "<<Epsilon<T>()*Norm(m)*norma<<endl;
	TMVAssert(Norm(mm-m)<10*Epsilon<T>()*Norm(m)*norma);
#endif
      }

      template <class T1, class T2> void DoLDiv(
	  const GenLowerTriMatrix<T1>& m, const LowerTriMatrixView<T2>& x) const
      { 
	TMVAssert(!x.isunit() || (Adt == UnitDiag && m.isunit()));
	Matrix<T1,RowMajor> mm(m);
	Matrix<T2,RowMajor> xx(x);
	DoLDiv(mm,xx.QuickView());
	x = LowerTriMatrixViewOf(xx.QuickView(),x.dt());
#ifdef CHECKTRI
	cerr<<"SVLTriLDiv: "<<Norm(xx-x)<<"  "<<Epsilon<T>()*Norm(x)*norma<<endl;
	TMVAssert(Norm(xx-x)<10*Epsilon<T>()*Norm(x)*norma);
#endif
      }

      template <class T1, class T2> void DoRDiv(
	  const GenLowerTriMatrix<T1>& m, const LowerTriMatrixView<T2>& x) const
      { 
	TMVAssert(!x.isunit() || (Adt == UnitDiag && m.isunit()));
	Matrix<T1,RowMajor> mm(m);
	Matrix<T2,RowMajor> xx(x);
	DoRDiv(mm,xx.QuickView());
	x = LowerTriMatrixViewOf(xx.QuickView(),x.dt());
#ifdef CHECKTRI
	cerr<<"SVLTriRDiv: "<<Norm(xx-x)<<"  "<<Epsilon<T>()*Norm(x)*norma<<endl;
	TMVAssert(Norm(xx-x)<10*Epsilon<T>()*Norm(x)*norma);
#endif
      }

#include "TMV_AuxAllDiv.h"
#include "TMV_AuxLowerTriAllDiv.h"

      inline T Det() const { return SS.Det(); }
      inline Matrix<T,ColMajor> Inverse() const { return SS.Inverse(); }
      inline bool Singular() const { return SS.Singular(); }
      inline RealType(T) Norm2() const { return SS.Norm2(); }
      inline Matrix<T,ColMajor> InverseATA() const { return SS.InverseATA(); }
      inline void Thresh(RealType(T) toler,ostream* debugout=0) const
      { SS.Thresh(toler,debugout); }
      inline void Top(size_t neigen,ostream* debugout=0) const 
      { SS.Top(neigen,debugout); }
      inline size_t GetKMax() const { return SS.GetKMax(); }

      //
      // Access Decomposition
      //

      inline Matrix<T,ColMajor> SV_GetU() const { return SS.SV_GetU(); }
      inline Vector<RealType(T)> SV_GetS() const 
      { return SS.SV_GetS(); }
      inline Matrix<T,ColMajor> SV_GetV() const { return SS.SV_GetV(); }
      inline const Matrix<T,ColMajor>& GetU() const { return SS.GetU(); }
      inline const Vector<RealType(T)>& GetS() const 
      { return SS.GetS(); }
      inline const Matrix<T,ColMajor>& GetV() const { return SS.GetV(); }

      LowerTriMatrix<T,NonUnitDiag,ColMajor> TInverse() const
      {
	Matrix<T,ColMajor> inv = Inverse();
#ifdef CHECKTRI
	cerr<<"SVLInv: "<<Norm(inv-LowerTriMatrixViewOf(inv,Adt))<<"  "<<Epsilon<T>()*Norm(inv)*norma<<endl;
	TMVAssert(Norm(inv-LowerTriMatrixViewOf(inv,Adt))<
	    10*Epsilon<T>()*Norm(inv)*norma);
#endif
	return LowerTriMatrix<T,NonUnitDiag,ColMajor>(inv);
      }

      inline std::string Type() const
      { return std::string("LowerTriSVDiv<")+tmv::Type(T())+">"; }

    private :

      SVDiv<T> SS;
      DiagType Adt;
#ifdef CHECKTRI
      RealType(T) norma;
#endif

  };

}; // namespace tmv;

#endif
