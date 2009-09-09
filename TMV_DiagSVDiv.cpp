
#include "TMV.h"
#include "TMV_Diag.h"

namespace tmv {

#define RT RealType(T)

  template <class T> void MakeRealS(
      const GenDiagMatrix<T>& M, const DiagMatrixView<T>& S,
      const DiagMatrixView<T>& V1, T& detuv)
  {
    for(size_t i=0;i<S.size();++i) {
      if (M(i) < 0.) {
	S(i) = -M(i);
	V1(i) = T(-1);
	detuv = -detuv;
      } else {
	S(i) = M(i);
	V1(i) = T(1);
      }
    }
  }

  template <class T> void MakeRealS(
      const GenDiagMatrix<complex<T> >& M, const DiagMatrixView<T>& S, 
      const DiagMatrixView<complex<T> >& V1, complex<T>& detuv)
  {
    for(size_t i=0;i<S.size();++i) {
      S(i) = abs(M(i));
      V1(i) = SIGN(M(i),S(i)); // M(i) = V1(i) * S(i)
      detuv *= V1(i);
    }
  }

  template <class T> void DiagSV_Decompose(
      const GenDiagMatrix<T>& m, Permutation& P,
      const DiagMatrixView<RealType(T)>& S, 
      const DiagMatrixView<T>& V1, T& det)
  {
    MakeRealS(m,S,V1,det);
    P = SortPermutation(S.diag(),DESCEND);
    S.diag() = P * S.diag();
    V1.diag() = P * V1.diag();
  }

  template <class T1, class T2> void SV_LDivEq(
      const Permutation& P, const GenDiagMatrix<RealType(T1)>& S,
      const GenDiagMatrix<T1>& V1, const size_t kmax,
      const MatrixView<T2>& m)
  {
    TMVAssert(P.size() == S.size());
    TMVAssert(V1.size() == S.size());
    TMVAssert(kmax <= S.size());
    TMVAssert(m.colsize() == S.size());

    // Solve A x = m  (returning x in place of m)
    // where A = Pt S V1 P

    // Pt y = m
    m = P * m;
    // S V1 z = y
    for(size_t k=0;k<kmax;++k) {
      m.row(k) *= CONJ(V1(k))/S(k);
    }
    m.SubMatrix(kmax,m.colsize(),0,m.rowsize()).Zero();
    // P x = z;
    m /= P;
  }

  template <class T> T DiagSVDiv<T>::Det() const 
  {
    if (!calcdet) {
      const CVIt<RT,Unit,NonConj> _end = S.diag().end();
      for(CVIt<RT,Unit,NonConj> it=S.diag().begin(); it!=_end; ++it)
	det *= *it;
      calcdet = true;
    }
    return det;
  }

  template <class T> DiagMatrix<T> DiagSVDiv<T>::DInverse() const 
  {
    // A^-1 = Vt S^-1 Ut
    //      = Pt V1t S^-1 P
    DiagMatrix<T> temp(S.size());
    VIt<T,Unit,NonConj> mit = temp.diag().begin();
    CVIt<T,Unit,NonConj> vit = V1.diag().begin();
    CVIt<RealType(T),Unit,NonConj> sit = S.diag().begin();
    for (size_t k=0;k<kmax;++k,++mit,++vit,++sit) 
      *mit = CONJ(*vit)/(*sit);
    temp.diag() *= P;
    return temp;
  }

  template <class T> DiagMatrix<T> DiagSVDiv<T>::DInverseATA() const 
  { 
    // (AtA)^-1 = Vt S^-2 V
    //          = Pt V1t S^-2 V1 P
    //          = Pt S^-2 P
    DiagMatrix<T> temp(S.size());
    VIt<T,Unit,NonConj> mit = temp.diag().begin();
    CVIt<RealType(T),Unit,NonConj> sit = S.diag().begin();
    for(size_t k=0;k<kmax;++k,++sit)
      *mit = RealType(T)(1) / ((*sit)*(*sit));
    temp.diag() *= P;
    return temp;
  }
		
  template <class T> void DiagSVDiv<T>::Thresh(
      RT toler, ostream* debugout) const
  {
    TMVAssert(toler < RT(1));
    RT thresh = S(0)*toler;
    for(kmax=S.size(); kmax>0 && S(kmax-1)<=thresh; kmax--);
    if (debugout) {
      (*debugout)<<"S = "<<S.diag()<<endl;
      (*debugout)<<"Smax = "<<S(0)<<", thresh = "<<thresh<<std::endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

  template <class T> void DiagSVDiv<T>::Top(
      size_t neigen, ostream* debugout) const
  {
    TMVAssert(neigen > 0);
    TMVAssert(neigen <= S.size());
    kmax = neigen;
    if(debugout) {
      (*debugout)<<"S = "<<S.diag()<<endl;
      (*debugout)<<"kmax = "<<kmax<<" (S.size = "<<S.size()<<")"<<endl;
    }
  }

#undef RT

#define InstFile "TMV_DiagSVDiv.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


