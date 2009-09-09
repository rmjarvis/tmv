
#include "TMV_Band.h"

namespace tmv {

  // 
  // Access
  //

  template <class T> T GenBandMatrix<T>::cref(size_t i, size_t j) const
  {
    TMVAssert(i<colsize() && j<rowsize());
    if (!okij(i,j)) return T(0);
    else {
      const T* mi;
      if (isrm()) mi = cptr() + int(i)*stepi() + j;
      else if (iscm()) mi = cptr() + i + int(j)*stepj();
      else mi = cptr()+int(i)*stepi()+int(j)*stepj();
      return isconj() ? CONJ(*mi) : *mi;
    }
  }

  template <class T, IndexStyle I> RefType(T) BandMatrixView<T,I>::ref(
      size_t i, size_t j) const
  {
    TMVAssert(this->okij(i,j));
    T* mi;
    if (isrm()) mi = ptr() + int(i)*stepi() + j;
    else if (iscm()) mi = ptr() + i + int(j)*stepj();
    else mi = ptr()+int(i)*stepi()+int(j)*stepj();
#ifdef TMVFLDEBUG
    TMVAssert(mi >= first);
    TMVAssert(mi < last);
#endif
    return REF(mi,ct());
  }

  template <class T> void GenBandMatrix<T>::NewDivider() const
  {
    switch (this->itsdt) {
      case LU : 
	this->itsdiv.reset(new BandLUDiv<T>(*this,this->divinplace)); break;
      case QR : 
	this->itsdiv.reset(new BandQRDiv<T>(*this,this->divinplace)); break;
      case SV : 
	this->itsdiv.reset(new BandSVDiv<T>(*this,true,true)); break;
      case SVS :
	this->itsdiv.reset(new BandSVDiv<T>(*this,false,false)); break;
      case SVU :
	this->itsdiv.reset(new BandSVDiv<T>(*this,true,false)); break;
      case SVV :
	this->itsdiv.reset(new BandSVDiv<T>(*this,false,true)); break;
      default : TMVAssert(FALSE);
    }
  }

  size_t BandStorageLength(StorageType s, size_t cs, size_t rs, int lo, int hi)
  {
    TMVAssert(s == RowMajor || s == ColMajor || s == DiagMajor);
    if (s == RowMajor) {
      if (cs == rs) return (cs-1)*(lo+hi+1)+1;
      else if (cs > rs) return min(rs+lo,cs)*(lo+hi)+rs;
      else return (cs-1)*(lo+hi)+min(cs+hi,rs);
    } else if (s == ColMajor) {
      if (cs == rs) return (cs-1)*(lo+hi+1)+1;
      else if (cs > rs) return (rs-1)*(lo+hi)+min(rs+lo,cs);
      else return min(cs+hi,rs)*(lo+hi)+cs;
    } else {
      if (cs == rs) return (cs-1)*(lo+hi+1)+1;
      else if (cs > rs) return rs*(lo+hi+1);
      else return min(cs+hi,rs)*(lo+hi)+cs;
    }
  }

  //
  // OK? (SubMatrix, etc.)
  //

  template <class T> bool GenBandMatrix<T>::OKSubMatrix(
      int i1, int i2, int j1, int j2, int istep, int jstep) const
  {
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    if (istep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 0 || i1 >= int(colsize())) {
      ok = false;
      cout<<"first col element ("<<i1<<") must be in 0 -- "<<colsize()-1<<endl;
    }
    if (i2-istep < 0 || i2-istep >= int(colsize())) {
      ok = false;
      cout<<"last col element ("<<i2-istep<<") must be in 0 -- "<<colsize()-1<<endl;
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
    if (j1 < 0 || j1 >= int(rowsize())) {
      ok = false;
      cout<<"first row element ("<<j1<<") must be in 0 -- "<<rowsize()-1<<endl;
    }
    if (j2-jstep < 0 || j2-jstep >= int(rowsize())) {
      ok = false;
      cout<<"last row element ("<<j2-jstep<<") must be in 0 -- "<<rowsize()-1<<endl;
    }
    if ((j2-j1)%jstep != 0) {
      ok = false;
      cout<<"row range ("<<j2-j1<<") must be multiple of istep ("<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
    }
    if (!okij(i1,j1)) {
      ok = false;
      cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in band\n";
    }
    if (!okij(i1,j2-jstep)) {
      ok = false;
      cout<<"Upper right corner ("<<i1<<','<<j2-jstep<<") must be in band\n";
    }
    if (!okij(i2-istep,j1)) {
      ok = false;
      cout<<"Lower left corner ("<<i2-istep<<','<<j1<<") must be in band\n";
    }
    if (!okij(i2-istep,j2-jstep)) {
      ok = false;
      cout<<"Lower right corner ("<<i2-istep<<','<<j2-jstep<<") must be in band\n";
    }
    return ok;
  }

  template <class T> bool GenBandMatrix<T>::OKSubVector(
      size_t i, size_t j, int istep, int jstep, size_t size) const
  {
    if (size==0) return true;
    bool ok = true;
    if (istep == 0 && jstep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") and jstep ("<<jstep<<") can not both be 0\n";
    }
    if (i >= colsize()) {
      ok = false;
      cout<<"i ("<<i<<") must be in 0 -- "<<colsize()-1<<endl;
    }
    if (j >= rowsize()) {
      ok = false;
      cout<<"j ("<<j<<") must be in 0 -- "<<rowsize()-1<<endl;
    }
    int i2 = int(i)+istep*int(size-1);
    int j2 = int(j)+jstep*int(size-1);
    if (i2 < 0 || i2 >= int(colsize())) {
      ok = false;
      cout<<"last element's i ("<<i2<<") must be in 0 -- "<<colsize()-1<<endl;
    }
    if (j2 < 0 || j2 >= int(rowsize())) {
      ok = false;
      cout<<"last element's j ("<<j2<<") must be in 0 -- "<<rowsize()-1<<endl;
    }
    if (!okij(i,j)) {
      ok = false;
      cout<<"First element ("<<i<<','<<j<<") must be in band\n";
    }
    if (!okij(i2,j2)) {
      ok = false;
      cout<<"Last element ("<<i2<<','<<j2<<") must be in band\n";
    }
    return ok;
  }
  
  template <class T> bool GenBandMatrix<T>::OKSubBandMatrix(
      int i1, int i2, int j1, int j2, int newnlo, int newnhi,
      int istep, int jstep) const
  {
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    if (istep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 0 || i1 >= int(colsize())) {
      ok = false;
      cout<<"first col element ("<<i1<<") must be in 0 -- "<<colsize()-1<<endl;
    }
    if (i2-istep < 0 || i2-istep >= int(colsize())) {
      ok = false;
      cout<<"last col element ("<<i2-istep<<") must be in 0 -- "<<colsize()-1<<endl;
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
    if (j1 < 0 || j1 >= int(rowsize())) {
      ok = false;
      cout<<"first row element ("<<j1<<") must be in 0 -- "<<rowsize()-1<<endl;
    }
    if (j2-jstep < 0 || j2-jstep >= int(rowsize())) {
      ok = false;
      cout<<"last row element ("<<j2-jstep<<") must be in 0 -- "<<rowsize()-1<<endl;
    }
    if ((j2-j1)%jstep != 0) {
      ok = false;
      cout<<"row range ("<<j2-j1<<") must be multiple of istep ("<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      cout<<"n row elements ("<<(j2-j1)/jstep<<") must be nonnegative\n";
    }
    if (!okij(i1,j1)) {
      ok = false;
      cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in band\n";
    }
    if (!okij(i1,j1+newnhi)) {
      ok = false;
      cout<<"Start of top diagonal ("<<i1<<','<<j1+newnhi<<") must be in band\n";
    }
    if (!okij(i1+newnlo,j1)) {
      ok = false;
      cout<<"Start of bottom diagonal ("<<i1+newnlo<<','<<j1<<") must be in band\n";
    }
    if (newnhi >= j2-j1) {
      ok = false;
      cout<<"new nhi ("<<newnhi<<") must be less than the new rowsize ("<<j2-j1<<")\n";
    }
    if (newnlo >= i2-i1) {
      ok = false;
      cout<<"new nlo ("<<newnlo<<") must be less than the new colsize ("<<i2-i1<<")\n";
    }
    return ok;
  } 

  template <class T> bool ConstBandMatrixView<T,FortranStyle>::OKSubMatrix(
      int i1, int i2, int j1, int j2, int istep, int jstep) const
  {
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    if (istep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 1 || i1 > int(this->colsize())) {
      ok = false;
      cout<<"first col element ("<<i1<<") must be in 1 -- "<<this->colsize()<<endl;
    }
    if (i2 < 1 || i2 > int(this->colsize())) {
      ok = false;
      cout<<"last col element ("<<i2<<") must be in 1 -- "<<this->colsize()<<endl;
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
    if (j1 < 0 || j1 >= int(this->rowsize())) {
      ok = false;
      cout<<"first row element ("<<j1<<") must be in 1 -- "<<this->rowsize()<<endl;
    }
    if (j2 < 0 || j2 >= int(this->rowsize())) {
      ok = false;
      cout<<"last row element ("<<j2<<") must be in 1 -- "<<this->rowsize()<<endl;
    }
    if ((j2-j1)%jstep != 0) {
      ok = false;
      cout<<"row range ("<<j2-j1<<") must be multiple of istep ("<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      cout<<"n row elements ("<<(j2-j1)/jstep+1<<") must be positive\n";
    }
    if (!this->okij(i1-1,j1-1)) {
      ok = false;
      cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in band\n";
    }
    if (!this->okij(i1-1,j2-1)) {
      ok = false;
      cout<<"Upper right corner ("<<i1<<','<<j2<<") must be in band\n";
    }
    if (!this->okij(i2-1,j1-1)) {
      ok = false;
      cout<<"Lower left corner ("<<i2<<','<<j1<<") must be in band\n";
    }
    if (!this->okij(i2-1,j2-1)) {
      ok = false;
      cout<<"Lower right corner ("<<i2<<','<<j2<<") must be in band\n";
    }
    return ok;
  }

  template <class T> bool ConstBandMatrixView<T,FortranStyle>::OKSubVector(
      size_t i, size_t j, int istep, int jstep, size_t size) const
  {
    if (size==0) return true;
    bool ok = true;
    if (istep == 0 && jstep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") and jstep ("<<jstep<<") can not both be 0\n";
    }
    if (i < 1 || i > this->colsize()) {
      ok = false;
      cout<<"i ("<<i<<") must be in 1 -- "<<this->colsize()<<endl;
    }
    if (j < 1 || j > this->rowsize()) {
      ok = false;
      cout<<"j ("<<j<<") must be in 1 -- "<<this->rowsize()<<endl;
    }
    int i2 = int(i)+istep*int(size-1);
    int j2 = int(j)+jstep*int(size-1);
    if (i2 < 1 || i2 > int(this->colsize())) {
      ok = false;
      cout<<"last element's i ("<<i2<<") must be in 1 -- "<<this->colsize()<<endl;
    }
    if (j2 < 1 || j2 > int(this->rowsize())) {
      ok = false;
      cout<<"last element's j ("<<j2<<") must be in 1 -- "<<this->rowsize()<<endl;
    }
    if (!this->okij(i-1,j-1)) {
      ok = false;
      cout<<"First element ("<<i<<','<<j<<") must be in band\n";
    }
    if (!this->okij(i2-1,j2-1)) {
      ok = false;
      cout<<"Last element ("<<i2<<','<<j2<<") must be in band\n";
    }
    return ok;
  }
  
  template <class T> bool ConstBandMatrixView<T,FortranStyle>::OKSubBandMatrix(
      int i1, int i2, int j1, int j2, int newnlo, int newnhi,
      int istep, int jstep) const
  {
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    if (istep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 1 || i1 > int(this->colsize())) {
      ok = false;
      cout<<"first col element ("<<i1<<") must be in 1 -- "<<this->colsize()<<endl;
    }
    if (i2 < 1 || i2 > int(this->colsize())) {
      ok = false;
      cout<<"last col element ("<<i2<<") must be in 1 -- "<<this->colsize()<<endl;
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
    if (j1 < 1 || j1 > int(this->rowsize())) {
      ok = false;
      cout<<"first row element ("<<j1<<") must be in 1 -- "<<this->rowsize()<<endl;
    }
    if (j2 < 1 || j2 > int(this->rowsize())) {
      ok = false;
      cout<<"last row element ("<<j2<<") must be in 1 -- "<<this->rowsize()<<endl;
    }
    if ((j2-j1)%jstep != 0) {
      ok = false;
      cout<<"row range ("<<j2-j1<<") must be multiple of istep ("<<jstep<<")\n";
    }
    if ((j2-j1)/jstep < 0) {
      ok = false;
      cout<<"n row elements ("<<(j2-j1)/jstep+1<<") must be positive\n";
    }
    if (!this->okij(i1-1,j1-1)) {
      ok = false;
      cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in band\n";
    }
    if (!this->okij(i1-1,j1-1+newnhi)) {
      ok = false;
      cout<<"Start of top diagonal ("<<i1<<','<<j1+newnhi<<") must be in band\n";
    }
    if (!this->okij(i1-1+newnlo,j1-1)) {
      ok = false;
      cout<<"Start of bottom diagonal ("<<i1+newnlo<<','<<j1<<") must be in band\n";
    }
    if (newnhi >= j2-j1+1) {
      ok = false;
      cout<<"new nhi ("<<newnhi<<") must be less than the new rowsize ("<<j2-j1+1<<")\n";
    }
    if (newnlo >= i2-i1+1) {
      ok = false;
      cout<<"new nlo ("<<newnlo<<") must be less than the new colsize ("<<i2-i1+1<<")\n";
    }
    return ok;
  } 

  //
  // Norms
  //

  template <class T> RealType(T) GenBandMatrix<T>::NormSq() const
  { 
    RealType(T) sum = 0;
    if (isrm()) {
      size_t j1=0;
      size_t j2=nhi()+1;
      size_t k=nlo();
      for(size_t i=0;i<colsize();++i) {
	sum += row(i,j1,j2).NormSq();
	if (k>0) --k; else ++j1;
	if (j2<rowsize()) ++j2;
	else if (j1==rowsize()) break;
      }
    } else if (iscm()) {
      size_t i1=0;
      size_t i2=nlo()+1;
      size_t k=nhi();
      for(size_t j=0;j<rowsize();++j) {
	sum += col(j,i1,i2).NormSq();
	if (k>0) --k; else ++i1;
	if (i2<colsize()) ++i2;
	else if (i1==colsize()) break;
      }
    } else {
      for(int i=-nlo();i<=nhi();++i) sum += diag(i).NormSq();
    }
    return sum;
  }

  template <class T> inline RealType(T) NonLapMaxAbsElement(
      const GenBandMatrix<T>& m)
  { 
    RealType(T) max = 0;
    if (m.isrm()) {
      size_t j1=0;
      size_t j2=m.nhi()+1;
      size_t k=m.nlo();
      for(size_t i=0;i<m.colsize();++i) {
	RealType(T) temp = m.row(i,j1,j2).MaxAbsElement();
	if (temp > max) max = temp;
	if (k>0) --k; else ++j1;
	if (j2<m.rowsize()) ++j2;
	else if (j1==m.rowsize()) break;
      }
    } else if (m.iscm()) {
      size_t i1=0;
      size_t i2=m.nlo()+1;
      size_t k=m.nhi();
      for(size_t j=0;j<m.rowsize();++j) {
	RealType(T) temp = m.col(j,i1,i2).MaxAbsElement();
	if (temp > max) max = temp;
	if (k>0) --k; else ++i1;
	if (i2<m.colsize()) ++i2;
	else if (i1==m.colsize()) break;
      }
    } else {
      for(int i=-m.nlo();i<=m.nhi();++i) {
	RealType(T) temp = m.diag(i).MaxAbsElement();
	if (temp > max) max = temp;
      }
    }
    return max;
  }

  // 1-Norm = max_j (sum_i |a_ij|)
  template <class T> inline RealType(T) NonLapNorm1(
      const GenBandMatrix<T>& m)
  { 
    RealType(T) max = 0;
    size_t i1=0;
    size_t i2=m.nlo()+1;
    size_t k=m.nhi();
    for(size_t j=0;j<m.rowsize();++j) {
      RealType(T) temp = m.col(j,i1,i2).Norm1();
      if (temp > max) max = temp;
      if (k>0) --k; else ++i1;
      if (i2<m.colsize()) ++i2;
      else if (i1==m.colsize()) break;
    }
    return max;
  }

  template <class T> inline RealType(T) NonLapNormF(
      const GenBandMatrix<T>& m)
  { return tmv::SQRT(m.NormSq()); }

  template <class T> inline RealType(T) NonLapNormInf(
      const GenBandMatrix<T>& m)
  { return NonLapNorm1(m.Transpose()); }

#ifdef XLAP
  template <class T> inline RealType(T) LapNorm(
      const char c, const GenBandMatrix<T>& m)
  { 
    switch(c) {
      case 'M' : return NonLapMaxAbsElement(m);
      case '1' : return NonLapNorm1(m);
      case 'F' : return NonLapNormF(m);
      case 'I' : return NonLapNormInf(m);
      default : TMVAssert(FALSE); 
    }
    return RealType(T)(0);
  }
  template <> inline double LapNorm(
      const char c, const GenBandMatrix<double>& m)
  { 
    TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
    TMVAssert(m.IsSquare());
    char cc = c;
    int n = m.colsize();
    double norm;
    if (m.iscm()) {
      int kl = m.nlo();
      int ku = m.nhi();
      int lda = m.stepj()+1;
#ifndef LAPNOWORK
      auto_array<double> work(c == 'I' ? new double[n] : 0);
#endif
      norm = LAPNAMEX(dlangb) (LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
	  LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
    } else {
      norm = LAPNAMEX(dlangt) (LAPCM LAPV(cc),LAPV(n),
	  LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),LAPP(m.cptr()+m.stepj()));
    }
    return norm;
  }
  template <> inline double LapNorm(
      const char c, const GenBandMatrix<complex<double> >& m)
  { 
    TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
    TMVAssert(m.IsSquare());
    char cc = c;
    int n = m.colsize();
    double norm;
    if (m.iscm()) {
      int kl = m.nlo();
      int ku = m.nhi();
      int lda = m.stepj()+1;
#ifndef LAPNOWORK
      auto_array<double> work(c == 'I' ? new double[n] : 0);
#endif
      norm = LAPNAMEX(zlangb) (LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
	  LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
    } else {
      norm = LAPNAMEX(zlangt) (LAPCM LAPV(cc),LAPV(n),
	  LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),LAPP(m.cptr()+m.stepj()));
    }
    return norm;
  }
#ifndef NOFLOAT
  template <> inline float LapNorm(
      const char c, const GenBandMatrix<float>& m)
  { 
    TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
    TMVAssert(m.IsSquare());
    char cc = c;
    int n = m.colsize();
    float norm;
    if (m.iscm()) {
      int kl = m.nlo();
      int ku = m.nhi();
      int lda = m.stepj()+1;
#ifndef LAPNOWORK
      auto_array<float> work(c == 'I' ? new float[n] : 0);
#endif
      norm = LAPNAMEX(slangb) (LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
	  LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
    } else {
      norm = LAPNAMEX(slangt) (LAPCM LAPV(cc),LAPV(n),
	  LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),LAPP(m.cptr()+m.stepj()));
    }
    return norm;
  }
  template <> inline float LapNorm(
      const char c, const GenBandMatrix<complex<float> >& m)
  { 
    TMVAssert(m.iscm() || (m.isdm() && m.nlo()==1 && m.nhi()==1));
    TMVAssert(m.IsSquare());
    char cc = c;
    int n = m.colsize();
    float norm;
    if (m.iscm()) {
      int kl = m.nlo();
      int ku = m.nhi();
      int lda = m.stepj()+1;
#ifndef LAPNOWORK
      auto_array<float> work(c == 'I' ? new float[n] : 0);
#endif
      norm = LAPNAMEX(clangb) (LAPCM LAPV(cc),LAPV(n),LAPV(kl),LAPV(ku),
	  LAPP(m.cptr()-m.nhi()),LAPV(lda) LAPWK(work.get()));
    } else {
      norm = LAPNAMEX(clangt) (LAPCM LAPV(cc),LAPV(n),
	  LAPP(m.cptr()+m.stepi()),LAPP(m.cptr()),LAPP(m.cptr()+m.stepj()));
    }
    return norm;
  }
#endif
#endif // XLAP

  template <class T> RealType(T) GenBandMatrix<T>::MaxAbsElement() const
  {
#ifdef XLAP
    if (isrm() && IsSquare()) return LapNorm('M',Transpose());
    else if (iscm() && IsSquare()) return LapNorm('M',*this);
    else if (isdm() && IsSquare() && nlo()==1 && nhi()==1)
      return LapNorm('M',*this);
    else
#endif
      return NonLapMaxAbsElement(*this);
  }
  template <class T> RealType(T) GenBandMatrix<T>::Norm1() const
  {
#ifdef XLAP
    if (isrm() && IsSquare()) return LapNorm('I',Transpose());
    else if (iscm() && IsSquare()) return LapNorm('1',*this);
    else if (isdm() && IsSquare() && nlo()==1 && nhi()==1)
      return LapNorm('1',*this);
    else
#endif
      return NonLapNorm1(*this);
  }
  template <class T> RealType(T) GenBandMatrix<T>::NormF() const
  {
#ifdef XLAP
    if (isrm() && IsSquare()) return LapNorm('F',Transpose());
    else if (iscm() && IsSquare()) return LapNorm('F',*this);
    else if (isdm() && IsSquare() && nlo()==1 && nhi()==1)
      return LapNorm('F',*this);
    else
#endif
      return NonLapNormF(*this);
  }

  //
  // Modifying Functions
  //

  template <class T, IndexStyle I>
    const BandMatrixView<T,I>& BandMatrixView<T,I>::Clip(RealType(T) thresh) const
    {
      if (isrm()) {
	size_t j1=0;
	size_t j2=nhi()+1;
	size_t k=nlo();
	for(size_t i=0;i<colsize();++i) {
	  row(i,j1,j2).Clip(thresh);
	  if (k>0) --k; else ++j1;
	  if (j2<rowsize()) ++j2;
	  else if (j1==rowsize()) break;
	}
      } else if (iscm()) {
	size_t i1=0;
	size_t i2=nlo()+1;
	size_t k=nhi();
	for(size_t j=0;j<rowsize();++j) {
	  col(j,i1,i2).Clip(thresh);
	  if (k>0) --k; else ++i1;
	  if (i2<colsize()) ++i2;
	  else if (i1==colsize()) break;
	}
      } else {
	for(int i=-nlo();i<=nhi();++i) diag(i).Clip(thresh);
      }
      return *this;
    }

  template <class T, IndexStyle I> 
    const BandMatrixView<T,I>& BandMatrixView<T,I>::SetAllTo(T x) const
    {
      if (isrm()) {
	size_t j1=0;
	size_t j2=nhi()+1;
	size_t k=nlo();
	for(size_t i=0;i<colsize();++i) {
	  row(i,j1,j2).SetAllTo(x);
	  if (k>0) --k; else ++j1;
	  if (j2<rowsize()) ++j2;
	  else if (j1==rowsize()) break;
	}
      } else if (iscm()) {
	size_t i1=0;
	size_t i2=nlo()+1;
	size_t k=nhi();
	for(size_t j=0;j<rowsize();++j) {
	  col(j,i1,i2).SetAllTo(x);
	  if (k>0) --k; else ++i1;
	  if (i2<colsize()) ++i2;
	  else if (i1==colsize()) break;
	}
      } else {
	for(int i=-nlo();i<=nhi();++i) diag(i).SetAllTo(x);
      }
      return *this;
    }

  template <class T, IndexStyle I> 
    const BandMatrixView<T,I>& BandMatrixView<T,I>::ConjugateSelf() const
    {
      if (IsComplex(T())) {
	if (isrm()) {
	  size_t j1=0;
	  size_t j2=nhi()+1;
	  size_t k=nlo();
	  for(size_t i=0;i<colsize();++i) {
	    row(i,j1,j2).ConjugateSelf();
	    if (k>0) --k; else ++j1;
	    if (j2<rowsize()) ++j2;
	    else if (j1==rowsize()) break;
	  }
	} else if (iscm()) {
	  size_t i1=0;
	  size_t i2=nlo()+1;
	  size_t k=nhi();
	  for(size_t j=0;j<rowsize();++j) {
	    col(j,i1,i2).ConjugateSelf();
	    if (k>0) --k; else ++i1;
	    if (i2<colsize()) ++i2;
	    else if (i1==colsize()) break;
	  }
	} else {
	  for(int i=-nlo();i<=nhi();++i) diag(i).ConjugateSelf();
	}
      }
      return *this;
    }

  //
  // Special Constructors
  //

  template <class T> BandMatrix<T,DiagMajor> UpperBiDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2)
  {
    if (v1.size() == v2.size()) {
      BandMatrix<T,DiagMajor> temp(v1.size(),v1.size()+1,0,1);
      temp.diag(0) = v1;
      temp.diag(1) = v2;
      return temp;
    } else {
      TMVAssert(v2.size() == v1.size()-1);
      BandMatrix<T,DiagMajor> temp(v1.size(),v1.size(),0,1);
      temp.diag(0) = v1;
      temp.diag(1) = v2;
      return temp;
    }
  }

  template <class T> BandMatrix<T,DiagMajor> LowerBiDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2)
  {
    if (v1.size() == v2.size()) {
      BandMatrix<T,DiagMajor> temp(v2.size()+1,v2.size(),1,0);
      temp.diag(-1) = v1;
      temp.diag(0) = v2;
      return temp;
    } else {
      TMVAssert(v1.size() == v2.size()-1);
      BandMatrix<T,DiagMajor> temp(v2.size(),v2.size(),1,0);
      temp.diag(-1) = v1;
      temp.diag(0) = v2;
      return temp;
    }
  }

  template <class T> BandMatrix<T,DiagMajor> TriDiagMatrix(
      const GenVector<T>& v1, const GenVector<T>& v2,
      const GenVector<T>& v3)
  {
    if (v1.size() == v2.size()) {
      TMVAssert(v3.size() == v2.size()-1);
      BandMatrix<T,DiagMajor> temp(v2.size()+1,v2.size(),1,1);
      temp.diag(-1) = v1;
      temp.diag(0) = v2;
      temp.diag(1) = v3;
      return temp;
    } else if (v2.size() == v3.size()) {
      TMVAssert(v1.size() == v2.size()-1);
      BandMatrix<T,DiagMajor> temp(v2.size(),v2.size()+1,1,1);
      temp.diag(-1) = v1;
      temp.diag(0) = v2;
      temp.diag(1) = v3;
      return temp;
    } else {
      TMVAssert(v1.size() == v2.size()-1);
      TMVAssert(v3.size() == v2.size()-1);
      BandMatrix<T,DiagMajor> temp(v2.size(),v2.size(),1,1);
      temp.diag(-1) = v1;
      temp.diag(0) = v2;
      temp.diag(1) = v3;
      return temp;
    }
  }

  //
  // Swap
  //

  template <class T> void Swap(const BandMatrixView<T>& m1,
      const BandMatrixView<T>& m2)
  {
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    TMVAssert(m1.nlo() == m2.nlo());
    TMVAssert(m1.nhi() == m2.nhi());
    if (m1.isrm() && m2.isrm()) {
      size_t j1=0;
      size_t j2=m1.nhi()+1;
      size_t k=m1.nlo();
      for(size_t i=0;i<m1.colsize();++i) {
	Swap(m1.row(i,j1,j2),m2.row(i,j1,j2));
	if (k>0) --k; else ++j1;
	if (j2<m1.rowsize()) ++j2;
	else if (j1==m1.rowsize()) break;
      }
    } else if (m1.iscm() && m2.iscm()) {
      size_t i1=0;
      size_t i2=m1.nlo()+1;
      size_t k=m1.nhi();
      for(size_t j=0;j<m1.rowsize();++j) {
	Swap(m1.col(j,i1,i2),m2.col(j,i1,i2));
	if (k>0) --k; else ++i1;
	if (i2<m1.colsize()) ++i2;
	else if (i1==m1.colsize()) break;
      }
    } else {
      for(int i=-m1.nlo();i<=m1.nhi();++i) Swap(m1.diag(i),m2.diag(i));
    }
  }

  //
  // m1 == m2
  //

  template <class T1, class T2> bool operator==(const GenBandMatrix<T1>& m1,
      const GenBandMatrix<T2>& m2)
  {
    if (m1.colsize() != m2.colsize()) return false;
    if (m1.rowsize() != m2.rowsize()) return false;
    if (m1.nlo() != m2.nlo()) return false;
    if (m1.nhi() != m2.nhi()) return false;
    if (m1.SameAs(m2)) return true;
    for(int i=-m1.nlo();i<=m1.nhi();++i) 
      if (m1.diag(i) != m2.diag(i)) return false;
    return true;
  }

  //
  // I/O
  //
 
  template <class T> void GenBandMatrix<T>::Write(ostream& os) const
  {
    size_t j1=0;
    size_t len=nhi()+1;
    const T* mrowi = cptr();
    const int ds = diagstep();
    const int si = stepi();
    const int sj = stepj();
    os << colsize()<<' '<<rowsize()<<endl;
    size_t k = nlo();
    for(size_t i=colsize();i>0;--i) {
      os << "( ";
      for(size_t j=j1;j>0;--j) os << ' '<<T(0)<<' ';
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' ' << CONJ(*mij) << ' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' ' << *mij << ' ';
      else 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj)
	    os << ' ' << CONJ(*mij) << ' ';
	else 
	  for(size_t k=len;k>0;--k,mij+=sj)
	    os << ' ' << *mij << ' ';
      for(size_t j=rowsize()-len-j1;j>0;--j) os << ' '<<T(0)<<' ';
      os << " )\n";
      if (k>0) { --k; ++len; mrowi+=si; } 
      else { 
	mrowi+=ds; 
	if (j1<rowsize()) ++j1; 
      }
      if (j1+len>rowsize()) --len;
    }
  }

  template <class T> void GenBandMatrix<T>::Write(ostream& os, 
      RealType(T) thresh) const
  {
    size_t j1=0;
    size_t len=nhi()+1;
    const T* mrowi = cptr();
    const int ds = diagstep();
    const int si = stepi();
    const int sj = stepj();
    os << colsize()<<' '<<rowsize()<<endl;
    size_t k = nlo();
    for(size_t i=colsize();i>0;--i) {
      os << "( ";
      for(size_t j=j1;j>0;--j) os << ' '<<T(0)<<' ';
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' ' <<(abs(*mij)<thresh ? T(0) : CONJ(*mij)) << ' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' ' <<(abs(*mij)<thresh ? T(0) : *mij) << ' ';
      else 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj)
	    os << ' ' <<(abs(*mij)<thresh ? T(0) : CONJ(*mij)) << ' ';
	else 
	  for(size_t k=len;k>0;--k,mij+=sj)
	    os << ' ' <<(abs(*mij)<thresh ? T(0) : *mij) << ' ';
      for(size_t j=rowsize()-len-j1;j>0;--j) os << ' '<<T(0)<<' ';
      os << " )\n";
      if (k>0) { --k; ++len; mrowi+=si; } 
      else { 
	mrowi+=ds; 
	if (j1<rowsize()) ++j1; 
      }
      if (j1+len>rowsize()) --len;
    }
  }

  template <class T> void GenBandMatrix<T>::WriteCompact(ostream& os) const
  {
    os << "B "<<colsize()<<' '<<rowsize()<<' '<<nlo()<<' '<<nhi()<<endl;
    size_t j1=0;
    size_t len=nhi()+1;
    const T* mrowi = cptr();
    const int ds = diagstep();
    const int si = stepi();
    const int sj = stepj();
    size_t k = nlo();
    for(size_t i=colsize();i>0;--i) {
      os << "( ";
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' ' << CONJ(*mij) << ' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' ' << *mij << ' ';
      else 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj)
	    os << ' ' << CONJ(*mij) << ' ';
	else 
	  for(size_t k=len;k>0;--k,mij+=sj)
	    os << ' ' << *mij << ' ';
      os << " )\n";
      if (k>0) { --k; ++len; mrowi+=si; } 
      else { 
	mrowi+=ds; 
	if (j1<rowsize()) ++j1; 
      }
      if (j1+len>rowsize()) --len;
    }
  }

  template <class T> void GenBandMatrix<T>::WriteCompact(ostream& os,
      RealType(T) thresh) const
  {
    os << "B "<<colsize()<<' '<<rowsize()<<' '<<nlo()<<' '<<nhi()<<endl;
    size_t j1=0;
    size_t len=nhi()+1;
    const T* mrowi = cptr();
    const int ds = diagstep();
    const int si = stepi();
    const int sj = stepj();
    size_t k = nlo();
    for(size_t i=colsize();i>0;--i) {
      os << "( ";
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' ' <<(abs(*mij)<thresh ? T(0) : CONJ(*mij)) << ' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' ' <<(abs(*mij)<thresh ? T(0) : *mij) << ' ';
      else 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj)
	    os << ' ' <<(abs(*mij)<thresh ? T(0) : CONJ(*mij)) << ' ';
	else 
	  for(size_t k=len;k>0;--k,mij+=sj)
	    os << ' ' <<(abs(*mij)<thresh ? T(0) : *mij) << ' ';
      os << " )\n";
      if (k>0) { --k; ++len; mrowi+=si; } 
      else { 
	mrowi+=ds; 
	if (j1<rowsize()) ++j1; 
      }
      if (j1+len>rowsize()) --len;
    }
  }

  template <class T> void BandMatrixReadError<T>::Write(
      ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for BandMatrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (m.get() && cs != m->colsize()) {
      os<<"Wrong colsize: expected "<<m->colsize()<<", got "<<cs<<".\n";
    }
    if (m.get() && rs != m->rowsize()) {
      os<<"Wrong rowsize: expected "<<m->rowsize()<<", got "<<rs<<".\n";
    }
    if (m.get() && lo != m->nlo()) {
      os<<"Wrong nlo: expected "<<m->nlo()<<", got "<<lo<<".\n";
    }
    if (m.get() && hi != m->nhi()) {
      os<<"Wrong nhi: expected "<<m->nhi()<<", got "<<hi<<".\n";
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
      os<<"The portion of the Bandatrix which was successfully read is: \n";
      ConstBandMatrixView<T> mm = m->View();
      for(size_t ii=0;ii<i;++ii) {
	os<<"( ";
	for(size_t jj=0;jj<mm.rowsize();++jj)
	  os<<' '<<mm(ii,jj)<<' ';
	os<<" )\n";
      }
      os<<"( ";
      for(size_t jj=0;jj<j;++jj)
	os<<' '<<mm(i,jj)<<' ';      
      os<<" )\n";
    }
  }


  template <class T, IndexStyle I> void BandMatrixView<T,I>::Read(
      istream& is) const
  {
    char paren;
    size_t j1=0;
    size_t len=nhi()+1;
    T* mrowi = ptr();
    const int ds = diagstep();
    const int si = stepi();
    const int sj = stepj();
    size_t k = nlo();
    for(size_t i=colsize();i>0;--i) {
      is >> paren;
      if (!is || paren != '(') 
	throw BandMatrixReadError<T>(colsize()-i,0,*this,is,
	    '(',is?paren:'(');
      T* mij = mrowi;
      if (isrm()) 
	for(size_t k=len;k>0;--k,++mij) {
	  is >> *mij;
	  if (!is)
	    throw BandMatrixReadError<T>(colsize()-i,j1+len-k,*this,is);
	}
      else 
	for(size_t k=len;k>0;--k,mij+=sj) {
	  is >> *mij;
	  if (!is)
	    throw BandMatrixReadError<T>(colsize()-i,j1+len-k,*this,is);
	}
      is >> paren;
      if ((!is && i>1) || paren != ')') 
	throw BandMatrixReadError<T>(colsize()-i,rowsize(),*this,is,
	    ')',is?paren:')');
      if (k>0) { --k; ++len; mrowi+=si; } 
      else { 
	mrowi+=ds; 
	if (j1<rowsize()) ++j1; 
      }
      if (j1+len>rowsize()) --len;
    }
    if (isconj()) ConjugateSelf();
  }

  template <class T, StorageType S, IndexStyle I> istream& operator>>(
      istream& is, auto_ptr<BandMatrix<T,S,I> >& m)
  { 
    char b;
    is >> b;
    if (!is)
      throw BandMatrixReadError<T>(is);
    if (b != 'B') 
      throw BandMatrixReadError<T>(is,'B',b);
    size_t cs,rs;
    int lo,hi;
    is >> cs >> rs >> lo >> hi;
    if (!is) 
      throw BandMatrixReadError<T>(is);
    m.reset(new BandMatrix<T,S,I>(cs,rs,lo,hi));
    m->View().Read(is); 
    return is;
  }

  template <class T> istream& operator>>(istream& is,
      const BandMatrixView<T>& m)
  { 
    char b;
    is >> b;
    if (!is)
      throw BandMatrixReadError<T>(is);
    if (b != 'B') 
      throw BandMatrixReadError<T>(is,'B',b);
    size_t cs,rs;
    int lo,hi;
    is >> cs >> rs >> lo >> hi;
    if (!is)
      throw BandMatrixReadError<T>(is);
    if (cs != m.colsize() || rs != m.rowsize() ||
	lo != m.nlo() || hi != m.nhi())
      throw BandMatrixReadError<T>(m,is,cs,rs,lo,hi);
    TMVAssert(m.colsize() == cs);
    TMVAssert(m.rowsize() == rs);
    TMVAssert(m.nlo() == lo);
    TMVAssert(m.nhi() == hi);
    m.Read(is);
    return is;
  }

#define InstFile "TMV_BandMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


