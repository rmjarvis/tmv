
#include "TMV.h"
#include "TMV_Tri.h"
#include "TMV_Sym.h"

//#define XDEBUG

namespace tmv {

  //
  // Access
  //

  template <class T> T GenSymMatrix<T>::cref(size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    const T* mi;
    if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
      if (isrm()) mi = cptr() + int(i)*stepi() + j;
      else if (iscm()) mi = cptr() + i + int(j)*stepj();
      else mi = cptr() + int(i)*stepi() + int(j)*stepj();
      return isconj() ? CONJ(*mi) : *mi;
    } else {
      if (isrm()) mi = cptr() + int(j)*stepi() + i;
      else if (iscm()) mi = cptr() + j + int(i)*stepj();
      else mi = cptr() + int(j)*stepi() + int(i)*stepj();
      return IsReal(T()) ? *mi : 
	(ct()==Conj) != (sym()==Herm) ? CONJ(*mi) : *mi;
    }
  }

  template <class T, IndexStyle I> RefType(T) SymMatrixView<T,I>::ref(
      size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    T* mi;
    if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) {
      if (isrm()) mi = ptr() + int(i)*stepi() + j;
      else if (iscm()) mi = ptr() + i + int(j)*stepj();
      else mi = ptr() + int(i)*stepi() + int(j)*stepj();
    } else {
      if (isrm()) mi = ptr() + int(j)*stepi() + i;
      else if (iscm()) mi = ptr() + j + int(i)*stepj();
      else mi = ptr() + int(j)*stepi() + int(i)*stepj();
    }
#ifdef TMVFLDEBUG
    TMVAssert(mi>=first);
    TMVAssert(mi<last);
#endif
    if ((uplo() == Upper && i<=j) || (uplo() == Lower && i>=j)) 
      return REF(mi,ct());
    else
      return IsReal(T()) ? REF(mi,NonConj) : 
	(ct()==Conj) != (sym()==Herm) ? REF(mi,Conj) : REF(mi,NonConj);
  }

  template <class T> void GenSymMatrix<T>::NewDivider() const
  { 
    switch(this->itsdt) {
      case LU : {
		  this->itsdiv.reset(new SymLUDiv<T>(*this,this->divinplace)); 
		  break; 
		}
      case SV : {
		  if (IsReal(T()) || sym()==Herm) {
		    this->itsdiv.reset(new HermSVDiv<T>(*this,true)); 
		  } else {
		    this->itsdiv.reset(new SymSVDiv<T>(*this,true,true)); 
		  }
		  break; 
		}
      case SVU : {
		  if (IsReal(T()) || sym()==Herm) {
		    this->itsdiv.reset(new HermSVDiv<T>(*this,true)); 
		  } else {
		    this->itsdiv.reset(new SymSVDiv<T>(*this,true,false)); 
		  }
		  break; 
		}
      case SVV : {
		  if (IsReal(T()) || sym()==Herm) {
		    this->itsdiv.reset(new HermSVDiv<T>(*this,true)); 
		  } else {
		    this->itsdiv.reset(new SymSVDiv<T>(*this,false,true)); 
		  }
		  break; 
		}
      case SVS : {
		  if (IsReal(T()) || sym()==Herm) {
		    this->itsdiv.reset(new HermSVDiv<T>(*this,false)); 
		  } else {
		    this->itsdiv.reset(new SymSVDiv<T>(*this,false,false)); 
		  }
		  break; 
		 }
      case CH : {
		  this->itsdiv.reset(new HermCHDiv<T>(*this,this->divinplace));
		  break;
		}
      default : TMVAssert(FALSE); 
    }
  }

  template <class T> void GenSymMatrix<T>::Inverse(
      const SymMatrixView<T>& sinv) const
  {
    TMVAssert(issym() == sinv.issym());
    TMVAssert(isherm() == sinv.isherm());

    this->SetDiv();
    const SymDivider<T>* sdiv = dynamic_cast<const SymDivider<T>*>(
	this->GetDiv());
    TMVAssert(sdiv);
    sdiv->Inverse(sinv);
  }

  //
  // OK? (SubMatrix, etc.)
  //

  template <class T> bool GenSymMatrix<T>::OKSubMatrix(
      int i1, int i2, int j1, int j2, int istep, int jstep) const
  { 
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    int i2x = i2-istep;
    int j2x = j2-jstep;
    if (istep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 0 || i1 >= int(size())) {
      ok = false;
      cout<<"first col element ("<<i1<<") must be in 0 -- "<<size()-1<<endl;
    }
    if (i2x < 0 || i2x >= int(size())) {
      ok = false;
      cout<<"last col element ("<<i2x<<") must be in 0 -- "<<size()-1<<endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      cout<<"col range ("<<i2-i1<<") must be multiple of istep ("<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      cout<<"n col elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
    }
    if (jstep == 0) {
      ok = false;
      cout<<"jstep ("<<jstep<<") can not be 0\n";
    }
    if (j1 < 0 || j1 >= int(size())) {
      ok = false;
      cout<<"first row element ("<<j1<<") must be in 0 -- "<<size()-1<<endl;
    }
    if (j2-jstep < 0 || j2-jstep >= int(size())) {
      ok = false;
      cout<<"last row element ("<<j2-jstep<<") must be in 0 -- "<<size()-1<<endl;
    }
    if ((j2-j1)%jstep != 0) {
      ok = false;
      cout<<"row range ("<<j2-j1<<") must be multiple of istep ("<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
    }
    if ((i1<j1 && i2x>j2x) || (i1>j1 && i2x<j2x)) {
      ok = false;
      cout<<"Upper left ("<<i1<<','<<j1<<") and lower right ("<<i2x<<','<<j2x<<") corners must be in same triangle\n";
    }
    if ((i2x<j1 && i1>j2x) || (i2x>j1 && i1<j2x)) {
      ok = false;
      cout<<"Upper right ("<<i1<<','<<j2x<<") and lower left ("<<i2x<<','<<j1<<") corners must be in same triangle\n";
    }
    return ok;
  }

  template <class T> bool GenSymMatrix<T>::OKSubVector(
      size_t i, size_t j, int istep, int jstep, size_t n) const 
  {
    if (n==0) return true;
    bool ok = true;
    if (istep == 0 && jstep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") and jstep ("<<jstep<<") can not both be 0\n";
    }
    if (i >= size()) {
      ok = false;
      cout<<"i ("<<i<<") must be in 0 -- "<<size()-1<<endl;
    }
    if (j >= size()) {
      ok = false;
      cout<<"j ("<<j<<") must be in 0 -- "<<size()-1<<endl;
    }
    int i2 = int(i)+istep*int(n-1);
    int j2 = int(j)+jstep*int(n-1);
    if (i2 < 0 || i2 >= int(size())) {
      ok = false;
      cout<<"last element's i ("<<i2<<") must be in 0 -- "<<size()-1<<endl;
    }
    if (j2 < 0 || j2 >= int(size())) {
      ok = false;
      cout<<"last element's j ("<<j2<<") must be in 0 -- "<<size()-1<<endl;
    }
    if ((i<j && i2>j2) || (i>j && i2<j2)) {
      ok = false;
      cout<<"First ("<<i<<','<<j<<") and last ("<<i2<<','<<j2<<") elements must be in same triangle\n";
    }
    return ok;
  }

  template <class T> bool GenSymMatrix<T>::OKSubSymMatrix(
      int i1, int i2, int istep) const 
  {
    if (i1==i2) return true;
    bool ok=true;
    if (istep == 0) {
      ok = false; 
      cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1<0 || i1 >= int(size())) {
      ok = false;
      cout<<"first diag element ("<<i1<<") must be in 0 -- "<<size()-1<<endl;
    }
    if (i2-istep<0 || i2-istep >= int(size())) {
      ok = false;
      cout<<"last diag element ("<<i2-istep<<") must be in 0 -- "<<size()-1<<endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      cout<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      cout<<"n diag elements ("<<(i2-i1)/istep<<") must be nonnegative\n";
    }
    return ok;
  }

  template <class T> bool ConstSymMatrixView<T,FortranStyle>::OKSubMatrix(
      int i1, int i2, int j1, int j2, int istep, int jstep) const
  { 
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    if (istep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 1 || i1 > int(this->size())) {
      ok = false;
      cout<<"first col element ("<<i1<<") must be in 1 -- "<<this->size()<<endl;
    }
    if (i2 < 1 || i2 > int(this->size())) {
      ok = false;
      cout<<"last col element ("<<i2-istep<<") must be in 1 -- "<<this->size()<<endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      cout<<"col range ("<<i2-i1<<") must be multiple of istep ("<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      cout<<"n col elements ("<<(i2-i1)/istep+1<<") must be positive\n";
    }
    if (jstep == 0) {
      ok = false;
      cout<<"jstep ("<<jstep<<") can not be 0\n";
    }
    if (j1 < 1 || j1 > int(this->size())) {
      ok = false;
      cout<<"first row element ("<<j1<<") must be in 1 -- "<<this->size()<<endl;
    }
    if (j2 < 1 || j2 > int(this->size())) {
      ok = false;
      cout<<"last row element ("<<j2-jstep<<") must be in 1 -- "<<this->size()<<endl;
    }
    if ((j2-j1)%jstep != 0) {
      ok = false;
      cout<<"row range ("<<j2-j1<<") must be multiple of istep ("<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      cout<<"n row elements ("<<(j2-j1)/jstep+1<<") must be positive\n";
    }
    if ((i1<j1 && i2>j2) || (i1>j1 && i2<j2)) {
      ok = false;
      cout<<"Upper left ("<<i1<<','<<j1<<") and lower right ("<<i2<<','<<j2<<") corners must be in same triangle\n";
    }
    if ((i2<j1 && i1>j2) || (i2>j1 && i1<j2)) {
      ok = false;
      cout<<"Upper right ("<<i1<<','<<j2<<") and lower left ("<<i2<<','<<j1<<") corners must be in same triangle\n";
    }
    return ok;
  }

  template <class T> bool ConstSymMatrixView<T,FortranStyle>::OKSubVector(
      size_t i, size_t j, int istep, int jstep, size_t n) const 
  {
    if (n==0) return true;
    bool ok = true;
    if (istep == 0 && jstep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") and jstep ("<<jstep<<") can not both be 0\n";
    }
    if (i < 1 || i > this->size()) {
      ok = false;
      cout<<"i ("<<i<<") must be in 1 -- "<<this->size()<<endl;
    }
    if (i < 1 || j > this->size()) {
      ok = false;
      cout<<"j ("<<j<<") must be in 1 -- "<<this->size()<<endl;
    }
    int i2 = int(i)+istep*int(n-1);
    int j2 = int(j)+jstep*int(n-1);
    if (i2 < 1 || i2 > int(this->size())) {
      ok = false;
      cout<<"last element's i ("<<i2<<") must be in 1 -- "<<this->size()<<endl;
    }
    if (j2 < 1 || j2 > int(this->size())) {
      ok = false;
      cout<<"last element's j ("<<j2<<") must be in 1 -- "<<this->size()<<endl;
    }
    if ((i<j && i2>j2) || (i>j && i2<j2)) {
      ok = false;
      cout<<"First ("<<i<<','<<j<<") and last ("<<i2<<','<<j2<<") elements must be in same triangle\n";
    }
    return ok;
  }

  template <class T> bool ConstSymMatrixView<T,FortranStyle>::OKSubSymMatrix(
      int i1, int i2, int istep) const 
  {
    if (i1==i2) return true;
    bool ok=true;
    if (istep == 0) {
      ok = false; 
      cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1<1 || i1 > int(this->size())) {
      ok = false;
      cout<<"first diag element ("<<i1<<") must be in 1 -- "<<this->size()<<endl;
    }
    if (i2-istep<1 || i2-istep > int(this->size())) {
      ok = false;
      cout<<"last diag element ("<<i2-istep<<") must be in 1 -- "<<this->size()<<endl;
    }
    if ((i2-i1)%istep != 0) {
      ok = false;
      cout<<"range ("<<i2-i1<<") must be multiple of istep ("<<istep<<")\n";
    }
    if ((i2-i1)/istep < 0) {
      ok = false;
      cout<<"n diag elements ("<<(i2-i1)/istep+1<<") must be positive\n";
    }
    return ok;
  }

  // 
  // SwapRowsCols
  //

  template <class T, IndexStyle I> 
    const SymMatrixView<T,I>& SymMatrixView<T,I>::SwapRowsCols(
	size_t i1, size_t i2) const
    {
      TMVAssert(i1<size());
      TMVAssert(i2<size());
      if (i1 == i2) return *this;
      else {
#ifdef XDEBUG
	Matrix<T> M = *this;
	Matrix<T> M0 = *this;
	M.SwapRows(i1,i2);
	M.SwapCols(i1,i2);
#endif
	if (i1 > i2) swap(i1,i2);
	// Now i1 < i2
	if (uplo() == Upper) Transpose().SwapRowsCols(i1,i2);
	else {
	  Swap(row(i1,0,i1),row(i2,0,i1));
	  Swap(row(i2,i1+1,i2),col(i1,i1+1,i2));
	  if (!issym()) {
	    row(i2,i1,i2).ConjugateSelf(); // Conj m(i2,i1) as well
	    col(i1,i1+1,i2).ConjugateSelf();
	  }
	  Swap(col(i1,i2+1,size()),col(i2,i2+1,size()));
	  diag().Swap(i1,i2);
	}
#ifdef XDEBUG
	if (tmv::Norm(M-*this) > 1.e-5*max(RealType(T)(1),tmv::Norm(M))) {
	  cerr<<"SwapRowsCols: i1,i2 = "<<i1<<','<<i2<<endl;
	  cerr<<"M0 = "<<M0<<endl;
	  cerr<<"M = "<<M<<endl;
	  cerr<<"S = "<<*this<<endl;
	  abort();
	}
#endif
	return *this;
      }
    }

  template <class T, IndexStyle I> 
    const SymMatrixView<T,I>& SymMatrixView<T,I>::PermuteRowsCols(
	const size_t* p, size_t i1, size_t i2) const
    {
      const size_t* pi = p+i1;
      for(size_t i=i1;i<i2;++i,++pi) {
	TMVAssert(*pi < size());
	SwapRowsCols(i,*pi);
      }
      return *this;
    }
 
  template <class T, IndexStyle I> 
    const SymMatrixView<T,I>& SymMatrixView<T,I>::ReversePermuteRowsCols(
	const size_t* p, size_t i1, size_t i2) const
    {
      const size_t* pi = p+i2;
      for(size_t i=i2;i>i1;) {
	--i; --pi;
	TMVAssert(*pi < size());
	SwapRowsCols(i,*pi);
      }
      return *this;
    }
 

  //
  // Norms
  //

  template <class T> RealType(T) GenSymMatrix<T>::NormSq() const
  { 
    RealType(T) ans = diag().NormSq();
    if (size() > 0) ans += RealType(T)(2) * UpperTri().OffDiag().NormSq();
    return ans;
  }

  template <class T> inline RealType(T) NonLapNorm1(const GenSymMatrix<T>& m) 
  { 
    RealType(T) max(0);
    for(size_t j=0;j<m.size();++j) {
      RealType(T) temp = m.col(j,0,j).Norm1();
      temp += m.col(j,j,m.size()).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  } 

  template <class T> inline RealType(T) NonLapNormF(const GenSymMatrix<T>& m)
  { return SQRT(m.NormSq()); }

  template <class T> inline RealType(T) NonLapMaxAbsElement(const GenSymMatrix<T>& m)
  { return m.UpperTri().MaxAbsElement(); }

#ifdef XLAP
  template <class T> inline RealType(T) LapNorm(const char c, 
      const GenSymMatrix<T>& m)
  {
    switch(c) {
      case 'M' : return NonLapMaxAbsElement(m);
      case '1' : return NonLapNorm1(m);
      case 'F' : return NonLapNormF(m);
      default : TMVAssert(FALSE);
    }
    return RealType(T)(0);
  }
  template <> double LapNorm(const char c, const GenSymMatrix<double>& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    int N = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<double> work(cc == '1' ? new double[N] : 0);
#endif
    double norm = LAPNAMEX(dlansy) (LAPCM LAPV(cc),
	(m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
	LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1);
    return norm;
  }
  template <> double LapNorm(const char c, 
      const GenSymMatrix<complex<double> >& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    int N = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<double> work(cc == '1' ? new double[N] : 0);
#endif
    double norm = m.isherm() ?
      LAPNAMEX(zlanhe) (LAPCM LAPV(cc),
	  (m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
	  LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1) :
      LAPNAMEX(zlansy) (LAPCM LAPV(cc),
	  (m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
	  LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1);
    return norm;
  }
#ifndef NOFLOAT
  template <> float LapNorm(const char c, const GenSymMatrix<float>& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    int N = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<float> work(cc == '1' ? new float[N] : 0);
#endif
    float norm = LAPNAMEX(slansy) (LAPCM LAPV(cc),
	(m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
	LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1);
    return norm;
  }
  template <> float LapNorm(const char c, 
      const GenSymMatrix<complex<float> >& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    int N = m.size();
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<float> work(cc == '1' ? new float[N] : 0);
#endif
    float norm = m.isherm() ?
      LAPNAMEX(clanhe) (LAPCM LAPV(cc),
	  (m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
	  LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1) :
      LAPNAMEX(clansy) (LAPCM LAPV(cc),
	  (m.iscm() == (m.uplo()==Upper) ? LAPCH_UP : LAPCH_LO),
	  LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()) LAP1 LAP1);
    return norm;
  }
#endif
#endif // XLAP

  template <class T> RealType(T) GenSymMatrix<T>::MaxAbsElement() const
  {
#ifdef XLAP
    if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
      return LapNorm('M',*this);
    else
#endif
      return NonLapMaxAbsElement(*this);
  }
  template <class T> RealType(T) GenSymMatrix<T>::Norm1() const
  {
#ifdef XLAP
    if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
      return LapNorm('1',*this);
    else
#endif
      return NonLapNorm1(*this);
  }
  template <class T> RealType(T) GenSymMatrix<T>::NormF() const
  {
#ifdef XLAP
    if ((isrm() && stepi()>0) || (iscm() && stepj()>0))
      return LapNorm('F',*this);
    else
#endif
      return NonLapNormF(*this);
  }

  //
  // I/O
  //

  template <class T> void GenSymMatrix<T>::Write(ostream& os) const
  {
    if (uplo() == Upper) 
      if (issym()) Transpose().Write(os);
      else Adjoint().Write(os);
    else {
      TMVAssert(uplo() == Lower);
      size_t rowlen = 0;
      size_t collen = size();
      const int si = stepi();
      const int sj = stepj();
      const T* mrowi = cptr();
      os << size() <<' '<<size() << endl;
      for(;collen>0;++rowlen,--collen,mrowi+=si) {
	os << "( ";
	const T* mij = mrowi;
	if (issym()) {
	  if (isconj()) {
	    for(size_t k=rowlen;k>0;--k,mij+=sj) os << ' '<<PRINTCONJ(*mij)<<' ';
	    for(size_t k=collen;k>0;--k,mij+=si) os << ' '<<PRINTCONJ(*mij)<<' ';
	  } else {
	    for(size_t k=rowlen;k>0;--k,mij+=sj) os << ' '<<*mij<<' ';
	    for(size_t k=collen;k>0;--k,mij+=si) os << ' '<<*mij<<' ';
	  }
	} else {
	  if (isconj()) {
	    for(size_t k=rowlen;k>0;--k,mij+=sj) os << ' '<<PRINTCONJ(*mij)<<' ';
	    for(size_t k=collen;k>0;--k,mij+=si) os << ' '<<*mij<<' ';
	  } else {
	    for(size_t k=rowlen;k>0;--k,mij+=sj) os << ' '<<*mij<<' ';
	    for(size_t k=collen;k>0;--k,mij+=si) os << ' '<<PRINTCONJ(*mij)<<' ';
	  }
	}
	os << " )\n";
      }
    }
  }

  template <class T> void GenSymMatrix<T>::Write(ostream& os,
      RealType(T) thresh) const
  {
    if (uplo() == Upper) 
      if (issym()) Transpose().Write(os);
      else Adjoint().Write(os);
    else {
      TMVAssert(uplo() == Lower);
      size_t rowlen = 0;
      size_t collen = size();
      const int si = stepi();
      const int sj = stepj();
      const T* mrowi = cptr();
      os << size() <<' '<<size() << endl;
      for(;collen>0;++rowlen,--collen,mrowi+=si) {
	os << "( ";
	const T* mij = mrowi;
	if (issym()) {
	  if (isconj()) {
	    for(size_t k=rowlen;k>0;--k,mij+=sj) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : PRINTCONJ(*mij))<<' ';
	    for(size_t k=collen;k>0;--k,mij+=si) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : PRINTCONJ(*mij))<<' ';
	  } else {
	    for(size_t k=rowlen;k>0;--k,mij+=sj) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : *mij)<<' ';
	    for(size_t k=collen;k>0;--k,mij+=si) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : *mij)<<' ';
	  }
	} else {
	  if (isconj()) {
	    for(size_t k=rowlen;k>0;--k,mij+=sj) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : PRINTCONJ(*mij))<<' ';
	    for(size_t k=collen;k>0;--k,mij+=si) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : *mij)<<' ';
	  } else {
	    for(size_t k=rowlen;k>0;--k,mij+=sj) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : *mij)<<' ';
	    for(size_t k=collen;k>0;--k,mij+=si) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : PRINTCONJ(*mij))<<' ';
	  }
	}
	os << " )\n";
      }
    }
  }

  template <class T> void GenSymMatrix<T>::WriteCompact(ostream& os) const
  {
    if (uplo() == Upper) 
      if (issym()) Transpose().WriteCompact(os);
      else Adjoint().WriteCompact(os);
    else {
      TMVAssert(uplo() == Lower);
      const T* mrowi = cptr();
      size_t rowlen = 1;
      const int si = stepi();
      const int sj = stepj();
      os << (issym() ? "S " : "H ") << size() << endl;
      for(size_t i=size();i>0;--i,++rowlen,mrowi+=si) {
	os << "( ";
	const T* mij = mrowi;
	if (isrm())
	  if (isconj()) 
	    for(size_t k=rowlen;k>0;--k,++mij) os << ' '<<PRINTCONJ(*mij)<<' ';
	  else 
	    for(size_t k=rowlen;k>0;--k,++mij) os << ' '<<*mij<<' ';
	else
	  if (isconj()) 
	    for(size_t k=rowlen;k>0;--k,mij+=sj) os << ' '<<PRINTCONJ(*mij)<<' ';
	  else 
	    for(size_t k=rowlen;k>0;--k,mij+=sj) os << ' '<<*mij<<' ';
	os << " )\n";
      }
    }
  }

  template <class T> void GenSymMatrix<T>::WriteCompact(ostream& os,
      RealType(T) thresh) const
  {
    if (uplo() == Upper) 
      if (issym()) Transpose().WriteCompact(os);
      else Adjoint().WriteCompact(os);
    else {
      TMVAssert(uplo() == Lower);
      const T* mrowi = cptr();
      size_t rowlen = 1;
      const int si = stepi();
      const int sj = stepj();
      os << (issym() ? "S " : "H ") << size() << endl;
      for(size_t i=size();i>0;--i,++rowlen,mrowi+=si) {
	os << "( ";
	const T* mij = mrowi;
	if (isrm())
	  if (isconj()) 
	    for(size_t k=rowlen;k>0;--k,++mij) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : PRINTCONJ(*mij))<<' ';
	  else 
	    for(size_t k=rowlen;k>0;--k,++mij) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : *mij)<<' ';
	else
	  if (isconj()) 
	    for(size_t k=rowlen;k>0;--k,mij+=sj) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : PRINTCONJ(*mij))<<' ';
	  else 
	    for(size_t k=rowlen;k>0;--k,mij+=sj) 
	      os << ' '<<(abs(*mij)<thresh ? T(0) : *mij)<<' ';
	os << " )\n";
      }
    }
  }

  template <class T> void SymMatrixReadError<T>::Write(
      ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for SymMatrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"'";
      if (IsReal(T()) && exp == 'S') os<<" (or 'H')";
      os<<", got '"<<got<<"'.\n";
    }
    if (m.get() && s != m->size()) {
      os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
    }
    if (!is) {
      if (iseof) {
	os<<"Input stream reached end-of-file prematurely.\n";
      } else if (isbad) {
	os<<"Input stream is corrupted.\n";
      } else {
	os<<"Input stream cannot read next character.\n";
      }
    }
    if (m.get()) {
      os<<"The portion of the SymMatrix which was successfully read is: \n";
      ConstSymMatrixView<T> mm = m->View();
      for(size_t ii=0;ii<i;++ii) {
	os<<"( ";
	for(size_t jj=0;jj<(ii<j?i+1:i);++jj)
	  os<<' '<<mm(ii,jj)<<' ';
	os<<" )\n";
      }
      os<<"( ";
      for(size_t jj=0;jj<j;++jj)
	os<<' '<<mm(i,jj)<<' ';      
      os<<" )\n";
    }
  }

  template <class T> void HermMatrixReadError<T>::Write(
      ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for HermMatrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"'";
      if (IsReal(T()) && exp == 'H') os<<" (or 'S')";
      os<<", got '"<<got<<"'.\n";
    }
    if (dv != T(0)) {
      os<<"Non-real value found on diagonal: "<<dv<<endl;
    }
    if (m.get() && s != m->size()) {
      os<<"Wrong size: expected "<<m->size()<<", got "<<s<<".\n";
    }
    if (!is) {
      if (iseof) {
	os<<"Input stream reached end-of-file prematurely.\n";
      } else if (isbad) {
	os<<"Input stream is corrupted.\n";
      } else {
	os<<"Input stream cannot read next character.\n";
      }
    }
    if (m.get()) {
      os<<"The portion of the HermMatrix which was successfully read is: \n";
      ConstSymMatrixView<T> mm = m->View();
      for(size_t ii=0;ii<i;++ii) {
	os<<"( ";
	for(size_t jj=0;jj<(ii<j?i+1:i);++jj)
	  os<<' '<<mm(ii,jj)<<' ';
	os<<" )\n";
      }
      os<<"( ";
      for(size_t jj=0;jj<j;++jj)
	os<<' '<<mm(i,jj)<<' ';      
      os<<" )\n";
    }
  }

  template <class T, IndexStyle I> void SymMatrixView<T,I>::Read(
      istream& is) const
  {
    if (uplo() == Upper) 
      if (issym()) Transpose().Read(is);
      else Adjoint().Read(is);
    else {
      TMVAssert(uplo() == Lower);
      char paren;
      T* mrowi = ptr();
      size_t rowlen = 1;
      const int si = stepi();
      const int sj = stepj();
      for(size_t i=size();i>0;--i,++rowlen,mrowi+=si) {
	is >> paren;
	if (!is || paren != '(') {
	  if (issym())
	    throw SymMatrixReadError<T>(size()-i,0,*this,is,'(',is?paren:'(');
	  else
	    throw HermMatrixReadError<T>(size()-i,0,*this,is,'(',is?paren:'(');
	}
	T* mij = mrowi;
	if (isrm()) {
	  for(size_t k=rowlen;k>0;--k,++mij) {
	    is >> *mij;
	    if (!is)
	      if (issym())
		throw SymMatrixReadError<T>(size()-i,rowlen-k,*this,is);
	      else
		throw HermMatrixReadError<T>(size()-i,rowlen-k,*this,is);
	    if (isherm() && k==1 && IMAG(*mij) != RealType(T)(0))
	      throw HermMatrixReadError<T>(size()-i,rowlen-k,*this,is,*mij);
	  }
	}
	else {
	  for(size_t k=rowlen;k>0;--k,mij+=sj) {
	    is >> *mij;
	    if (!is)
	      if (issym())
		throw SymMatrixReadError<T>(size()-i,rowlen-k,*this,is);
	      else
		throw HermMatrixReadError<T>(size()-i,rowlen-k,*this,is);
	    if (isherm() && k==1 && IMAG(*mij) != RealType(T)(0))
	      throw HermMatrixReadError<T>(size()-i,rowlen-k,*this,is,*mij);
	  }
	}
	is >> paren;
	if ((!is && i>1)  || paren != ')') {
	  if (issym())
	    throw SymMatrixReadError<T>(size()-i,size()-i+1,*this,is,
		')',is?paren:'(');
	  else
	    throw HermMatrixReadError<T>(size()-i,size()-i+1,*this,is,
		')',is?paren:'(');
	}
      }
    }
    if (isconj()) ConjugateSelf();
  }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    istream& operator>>(istream& is, auto_ptr<SymMatrix<T,U,S,I> >& m)
    { 
      char sh;
      is >> sh;
      if (!is)
	throw SymMatrixReadError<T>(is);
      if (IsReal(T())) {
	if (sh != 'S' && sh != 'H') 
	  throw SymMatrixReadError<T>(is,'S',sh);
      } else {
	if (sh != 'S') 
	  throw SymMatrixReadError<T>(is,'S',sh);
      }
      size_t size;
      is >> size;
      if (!is) 
	throw SymMatrixReadError<T>(is);
      m.reset(new SymMatrix<T,U,S,I>(size));
      m->View().Read(is); 
      return is;
    }

  template <class T, UpLoType U, StorageType S, IndexStyle I> 
    istream& operator>>(istream& is, auto_ptr<HermMatrix<T,U,S,I> >& m)
    { 
      char sh;
      is >> sh;
      if (!is)
	throw HermMatrixReadError<T>(is);
      if (IsReal(T())) {
	if (sh != 'S' && sh != 'H') 
	  throw HermMatrixReadError<T>(is,'H',sh);
      } else {
	if (sh != 'H') 
	  throw HermMatrixReadError<T>(is,'H',sh);
      }
      size_t size;
      is >> size;
      if (!is) 
	throw HermMatrixReadError<T>(is);
      m.reset(new HermMatrix<T,U,S,I>(size));
      m->View().Read(is); 
      return is;
    }

  template <class T> istream& operator>>(
      istream& is, const SymMatrixView<T>& m)
  { 
    char sh;
    is >> sh;
    if (!is)
      if (m.issym())
	throw SymMatrixReadError<T>(is);
      else
	throw HermMatrixReadError<T>(is);
    if (IsReal(T())) {
      if (sh != 'S' && sh != 'H') 
	if (m.issym())
	  throw SymMatrixReadError<T>(is,'S',sh);
	else
	  throw HermMatrixReadError<T>(is,'H',sh);
    } else if (m.issym()) {
      if (sh != 'S') 
	throw SymMatrixReadError<T>(is,'S',sh);
    } else {
      if (sh != 'H') 
	throw HermMatrixReadError<T>(is,'H',sh);
    }
    size_t s;
    is >> s;
    if (!is) 
      if (m.issym())
	throw SymMatrixReadError<T>(is);
      else
	throw HermMatrixReadError<T>(is);
    if (s != m.size())
      if (m.issym())
	throw SymMatrixReadError<T>(m,is,s);
      else
	throw HermMatrixReadError<T>(m,is,s);
    TMVAssert(m.size() == s);
    m.Read(is);
    return is;
  }


#define InstFile "TMV_SymMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


