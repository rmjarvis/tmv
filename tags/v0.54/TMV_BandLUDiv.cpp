
#include "TMV_Band.h"
#include "TMV_Diag.h"
#include "TMV_Tri.h"

namespace tmv {

#define NEWLO min(A.nlo(),A.nhi())
#define NEWHI min(A.nlo()+A.nhi(),int(A.colsize())-1)
#define APTR1 inplace ? 0 : \
  new T[BandStorageLength(ColMajor,A.colsize(),A.colsize(),NEWLO,NEWHI)]
#define APTR inplace ? A.NonConst().ptr() : Aptr1.get()

#define LUX istrans ? \
  inplace ? A.NonConst().Transpose() : \
  BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nhi(), NEWHI,  \
      (A.nlo() == 1 && A.nhi() == 1) ? DiagMajor : ColMajor) : \
  inplace ? A.NonConst().View() : \
  BandMatrixViewOf(Aptr,A.colsize(),A.colsize(),A.nlo(), NEWHI, \
      (A.nlo() == 1 && A.nhi() == 1) ? DiagMajor : ColMajor)

  template <class T> BandLUDiv<T>::BandLUDiv(const GenBandMatrix<T>& A,
      bool _inplace) :
    istrans(A.nhi()<A.nlo() || (A.nhi()==A.nlo() && A.isrm())),
    inplace(NEWLO == 0 || (_inplace && 
	  ((A.isrm() && istrans) || (A.iscm() && !istrans) || 
	   (A.isdm() && A.nlo()==1 && A.nhi()==1)))),
    Aptr1(APTR1), Aptr(APTR), LUx(LUX),
    P(new size_t[A.colsize()]), det(T(1)), donedet(false)
  {
    TMVAssert(A.IsSquare());
    if (inplace) {
      // For inplace decomposition, make sure the original band matrix
      // has room for the extra upper diagonals...
      // if iscm stepj >= (2*A.nlo()+A.nhi())
      // if isdm extra diags appear at end, so can't really check
      // if isrm stepi >= (2*A.nlo()+A.nhi())
      // -- although rm is currently not allowed.
      //TMVAssert(!LUx.isrm() || LUx.stepi()>=2*LUx.nlo()+LUx.nhi());
      TMVAssert(!LUx.iscm() || LUx.stepj()>=2*LUx.nlo()+LUx.nhi());
      TMVAssert(LUx.nlo() == 0 || LUx.iscm() || 
	  (LUx.isdm() && LUx.nlo() == 1 && LUx.nhi()==2));
      TMVAssert(LUx == (istrans ? A.Transpose() : A.View()));
    } else {
      if (istrans) BandMatrixViewOf(LUx,A.nhi(),A.nlo()) = A.Transpose();
      else BandMatrixViewOf(LUx,A.nlo(),A.nhi()) = A;
    }

    if (LUx.nlo() > 0) {
      int Anhi = istrans ? A.nlo() : A.nhi();
      if (Anhi < int(A.colsize())-1)
	LUx.Diags(Anhi+1,LUx.nhi()+1).Zero();
      BandLU_Decompose(LUx,P.get(),det,Anhi);
    } else {
      P.get()[0] = 1; // A tag to indicate that P is not set yet.
    }
  }

#undef LUX
#undef APTR
#undef NEWLO
#undef NEWHI

  template <class T> bool BandLUDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenBandMatrix<T>* bm = dynamic_cast<const GenBandMatrix<T>*>(&m);
    TMVAssert(bm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*bm)<<"  ";
      *fout <<(istrans ? bm->Transpose() : bm->View()) <<endl;
      *fout << "L = "<<GetL()<<endl;
      *fout << "U = "<<GetU()<<endl;
    }
    Matrix<T> lu = GetL() * Matrix<T>(GetU());
    lu.ReversePermuteRows(GetP());
    RealType(T) nm = Norm(lu- (istrans ? bm->Transpose() : bm->View()) );
    nm /= Norm(GetL())*Norm(GetU());
    if (fout) {
      *fout << "PLU = "<<lu<<endl;
      *fout << "Norm(M-PLU)/Norm(PLU) = "<<nm<<endl;
    }
    BandMatrix<T> bm2 = *bm;
    bm2.DivideUsing(SVS);
    bm2.SetDiv();
    return nm < bm2.Condition()*bm2.colsize()*Epsilon<T>();
  }

  template <class T> T BandLUDiv<T>::Det() const
  {
    if (!donedet) {
      det *= DiagMatrixViewOf(LUx.diag()).Det();
      donedet = true;
    }         
    return det;  
  }                  

  template <class T> LowerTriMatrix<T,UnitDiag> BandLUDiv<T>::GetL() const
  {
    size_t N = LUx.colsize();
    int nlo = LUx.nlo();
    LowerTriMatrix<T,UnitDiag> L(N,T(0));
    if (nlo == 0) 
      L.SetToIdentity();
    else {
      for(size_t i=0;i<N;++i) {
	Swap(L.row(i,0,i),L.row(P.get()[i],0,i));
	size_t end = min(i+nlo+1,N);
	L.col(i,i+1,end) = LUx.col(i,i+1,end);
      }
    }
    return L;
  }

  template <class T> const size_t* BandLUDiv<T>::GetP() const
  {
    if (LUx.nlo() == 0 && LUx.colsize() > 0 && P.get()[0] == 1) {
      for(size_t i=0;i<LUx.colsize();++i) P.get()[i] = i;
    }
    return P.get();
  }

#define InstFile "TMV_BandLUDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


