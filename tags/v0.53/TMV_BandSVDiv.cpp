
#include "TMV_Band.h"
#include "TMV_SVDiv.h"
#include "TMV_Diag.h"

namespace tmv {

  template <class T> BandSVDiv<T>::BandSVDiv(const GenBandMatrix<T>& A,
      bool StoreU, bool StoreV) :
    istrans(A.colsize() < A.rowsize()),
    U(0), S(min(A.rowsize(),A.colsize())), V(0), det(T(1))
  {
    size_t M = max(A.colsize(),A.rowsize());
    size_t N = min(A.colsize(),A.rowsize());

    if (istrans) swap(StoreU,StoreV);

    if (StoreU) U = new Matrix<T,ColMajor>(M,N);
    if (StoreV) V = new Matrix<T,ColMajor>(N,N);
    MatrixView<T>* Uv = U ? new MatrixView<T>(U->View()) : 0;
    MatrixView<T>* Vv = V ? new MatrixView<T>(V->View()) : 0;

    if (istrans)
      BandSV_Decompose(A.Transpose(),Uv,S.View(),Vv,det);
    else
      BandSV_Decompose(A,Uv,S.View(),Vv,det);

    delete Uv;
    delete Vv;

    Thresh(Epsilon<T>());
  }

  template <class T> BandSVDiv<T>::~BandSVDiv() 
  { if (U) delete U; if (V) delete V; }

  template <class T> bool BandSVDiv<T>::CheckDecomp(
      const BaseMatrix<T>& m, ostream* fout) const
  {
    const GenBandMatrix<T>* bm = dynamic_cast<const GenBandMatrix<T>*>(&m);
    TMVAssert(bm);
    if (fout) {
      *fout << "M = "<<tmv::Type(*bm)<<"  "<<*bm<<endl;
      *fout << "U = "<<GetU()<<endl;
      *fout << "S = "<<GetS()<<endl;
      *fout << "V = "<<GetV()<<endl;
    }
    Matrix<T> usv = GetU()*DiagMatrixViewOf(GetS())*GetV();
    RealType(T) nm = Norm(usv-(*bm));
    nm /= Norm(GetU())*Norm(GetS())*Norm(GetV());
    if (fout) {
      *fout << "USV = "<<usv<<endl;
      *fout << "Norm(M-USV) = "<<Norm(*bm-usv)<<endl;
      *fout << "Norm(M-USV)/Norm(USV) = "<<nm<<endl;
    }
    return nm < Condition()*bm->colsize()*Epsilon<T>();
  }

  template <class T> void BandSVDiv<T>::Thresh(RealType(T) toler,
      ostream* debugout) const
  {
    TMVAssert(toler < RealType(T)(1));
    RealType(T) thresh = S(0)*toler;
    for(kmax=S.size(); kmax>0 && S(kmax-1)<=thresh; --kmax);
    if(debugout) {
      (*debugout)<<"S = "<<S<<endl;
      (*debugout)<<"Smax = "<<S(0)<<", thresh = "<<thresh<<std::endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

  template <class T> void BandSVDiv<T>::Top(size_t neigen,
      ostream* debugout) const
  {
    TMVAssert(neigen > 0);
    TMVAssert(neigen <= S.size());
    kmax = neigen;
    if(debugout) {
      (*debugout)<<"S = "<<S<<endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

#define InstFile "TMV_BandSVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


