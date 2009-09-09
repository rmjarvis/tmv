
#include "TMV.h"

namespace tmv {

#ifdef TMV_BLOCKSIZE
  const size_t PERM_BLOCKSIZE = TMV_BLOCKSIZE/2;
#else
  const size_t PERM_BLOCKSIZE = 32;
#endif

  //
  // Access
  //

  template <class T> T GenMatrix<T>::cref(size_t i, size_t j) const
  {
    TMVAssert(i<colsize() && j<rowsize());
    const T* mi;
    if (isrm()) mi = cptr() + int(i)*stepi() + j;
    else if (iscm()) mi = cptr() + i + int(j)*stepj();
    else mi = cptr() + int(i)*stepi() + int(j)*stepj();
    return itsct==Conj ? CONJ(*mi) : *mi;
  }

  template <class T, IndexStyle I> RefType(T) MatrixView<T,I>::ref(
      size_t i, size_t j) const
  {
    TMVAssert(i<colsize());
    TMVAssert(j<rowsize());
    T* mi;
    if (isrm()) mi = ptr() + int(i)*itssi + int(j);
    else if (iscm()) mi = ptr() + i + int(j)*itssj;
    else mi = ptr() + int(i)*itssi + int(j)*itssj;
#ifdef TMVFLDEBUG
    TMVAssert(mi >= first);
    TMVAssert(mi < last);
#endif
    return REF(mi,ct());
  }

  template <class T> void GenMatrix<T>::NewDivider() const
  {
    switch (this->itsdt) {
      case tmv::LU : this->itsdiv.reset(new LUDiv<T>(*this,
			   this->divinplace)); 
		     break;
      case tmv::QR : this->itsdiv.reset(new QRDiv<T>(*this,
			   this->divinplace)); 
		     break;
      case tmv::QRP : this->itsdiv.reset(new QRPDiv<T>(*this,
			    this->divinplace)); 
		      break;
      case tmv::SV : this->itsdiv.reset(new SVDiv<T>(*this,
			   this->divinplace,true,true)); 
		     break;
      case tmv::SVS : this->itsdiv.reset(new SVDiv<T>(*this,
			    this->divinplace,false,false)); 
		      break;
      case tmv::SVU : this->itsdiv.reset(new SVDiv<T>(*this,
			    this->divinplace,true,false)); 
		      break;
      case tmv::SVV : this->itsdiv.reset(new SVDiv<T>(*this,
			    this->divinplace,false,true)); 
		      break;
      default : cerr<<"dt = "<<Text(this->itsdt)<<endl; TMVAssert(FALSE);
    }
  }

#define NEW_SIZE(cs,rs) GenMatrix<T>(ColMajor,NonConj,((cs)*(rs))), \
  itsm(new T[ls()]), itscs(cs), itsrs(rs) \
  DEFFIRSTLAST(itsm.get(),itsm.get()+ls())

  template <class T, StorageType S, IndexStyle I> Matrix<T,S,I>::Matrix(
      const vector<vector<T> >& vv) :
    NEW_SIZE(vv.size(),(vv.size()>0?vv[0].size():0))
  {
    T* vi=itsm.get();
    if (S == RowMajor) {
      for(size_t i=0;i<colsize();++i) {
	TMVAssert(vv[i].size() == rowsize());
	typename vector<T>::const_iterator vvi = vv[i].begin();
	for(size_t j=0;j<rowsize();++j,++vi,++vvi) 
	  *vi = *vvi;
      }
    } else { 
      for(size_t j=0;j<rowsize();++j) {
	for(size_t i=0;i<colsize();++i,++vi) {
	  TMVAssert(vv[i].size() == rowsize());
	  *vi = vv[i][j];
	}
      }
    }
  }

#undef NEW_SIZE

  //
  // OK? (SubMatrix, SubVector)
  //

  template <class T> bool GenMatrix<T>::OKSubMatrix(
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
    return ok;
  }
  
  template <class T> bool GenMatrix<T>::OKSubVector(
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
    return ok;
  }

  template <class T> bool ConstMatrixView<T,FortranStyle>::OKSubMatrix(
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
    return ok;
  }

  template <class T> bool ConstMatrixView<T,FortranStyle>::OKSubVector(
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
    return ok;
  }

  //
  // Norms
  //

  template <class T> RealType(T) GenMatrix<T>::NormSq() const
  {
    if (CanLinearize()) return ConstLinearView().NormSq();
    else {
      RealType(T) sum(0);
      if (isrm()) 
	for(size_t i=0;i<colsize();++i) 
	  sum += row(i).NormSq();
      else 
	for(size_t j=0;j<rowsize();++j) 
	  sum += col(j).NormSq();
      return sum;
    }
  }

  template <class T> inline RealType(T) NonLapMaxAbsElement(const GenMatrix<T>& m)
  {
    if (m.CanLinearize()) return m.ConstLinearView().MaxAbsElement();
    else {
      RealType(T) max(0);
      if (m.iscm())
	for(size_t j=0;j<m.rowsize();++j) {
	  RealType(T) temp = m.col(j).NormInf();
	  if (temp > max) max = temp;
	}
      else 
	for(size_t i=0;i<m.colsize();++i) {
	  RealType(T) temp = m.row(i).NormInf();
	  if (temp > max) max = temp;
	}
      return max;
    }
  }

  template <class T> inline RealType(T) NonLapNorm1(const GenMatrix<T>& m)
  { 
    RealType(T) max(0);
    for(size_t j=0;j<m.rowsize();++j) {
      RealType(T) temp = m.col(j).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  }

  template <class T> inline RealType(T) NonLapNormF(const GenMatrix<T>& m)
  { return tmv::SQRT(m.NormSq()); }

  template <class T> inline RealType(T) NonLapNormInf(const GenMatrix<T>& m)
  { return NonLapNorm1(m.Transpose()); }

#ifdef XLAP
  template <class T> inline RealType(T) LapNorm(const char c, const GenMatrix<T>& m)
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
  template <> inline double LapNorm(const char c, const GenMatrix<double>& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
#ifndef LAPNOWORK
    auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAMEX(dlange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
	LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
  template <> inline double LapNorm(const char c, const GenMatrix<complex<double> >& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
#ifndef LAPNOWORK
    auto_array<double> work(c == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAMEX(zlange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
	LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
#ifndef NOFLOAT
  template <> inline float LapNorm(const char c, const GenMatrix<float>& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
#ifndef LAPNOWORK
    auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
    double norm = LAPNAMEX(slange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
	LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
  template <> inline float LapNorm(const char c, const GenMatrix<complex<float> >& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
#ifndef LAPNOWORK
    auto_array<float> work(c == 'I' ? new float[M] : 0);
#endif
    double norm = LAPNAMEX(clange) (LAPCM LAPV(cc),LAPV(M),LAPV(N),
	LAPP(m.cptr()),LAPV(lda) LAPWK(work.get()));
    return norm;
  }
#endif
#endif // XLAP

  template <class T> RealType(T) GenMatrix<T>::MaxAbsElement() const
  {
#ifdef XLAP
    if (isrm() && stepi() > 0) return LapNorm('M',Transpose());
    else if (iscm() && stepj() > 0) return LapNorm('M',*this);
    else
#endif
      return NonLapMaxAbsElement(*this);
  }
  template <class T> RealType(T) GenMatrix<T>::Norm1() const
  {
#ifdef XLAP
    if (isrm() && stepi() > 0) return LapNorm('I',Transpose());
    else if (iscm() && stepj() > 0) return LapNorm('1',*this);
    else
#endif
      return NonLapNorm1(*this);
  }
  template <class T> RealType(T) GenMatrix<T>::NormF() const
  {
#ifdef XLAP
    if (isrm() && stepi() > 0) return LapNorm('F',Transpose());
    else if (iscm() && stepj() > 0) return LapNorm('F',*this);
    else
#endif
      return NonLapNormF(*this);
  }

  //
  // Modifying Functions
  //

  template <class T, IndexStyle I> const MatrixView<T,I>& MatrixView<T,I>::Clip(
      RealType(T) thresh) const
  { 
    if (this->CanLinearize()) LinearView().Clip(thresh);
    else {
      if (isrm())
	for(size_t i=0;i<colsize();++i) row(i).Clip(thresh);
      else 
	for(size_t j=0;j<rowsize();++j) col(j).Clip(thresh);
    }
    return *this; 
  }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::SetAllTo(T x) const
    { 
      if (this->CanLinearize()) LinearView().SetAllTo(x);
      else {
	if (isrm())
	  for(size_t i=0;i<colsize();++i) row(i).SetAllTo(x); 
	else 
	  for(size_t j=0;j<rowsize();++j) col(j).SetAllTo(x); 
      }
      return *this; 
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::TransposeSelf() const
    {
      TMVAssert(IsSquare());
      for(size_t i=1;i<colsize();++i) Swap(row(i,0,i),col(i,0,i));
      return *this; 
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::ConjugateSelf() const
    { 
      if (IsComplex(T())) {
	if (this->CanLinearize()) LinearView().ConjugateSelf();
	else {
	  if (isrm()) 
	    for(size_t i=0;i<colsize();++i) row(i).ConjugateSelf();
	  else 
	    for(size_t j=0;j<rowsize();++j) col(j).ConjugateSelf();
	}
      }
      return *this; 
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::SetToIdentity(T x) const 
    {
      TMVAssert(IsSquare());
      Zero();
      diag().SetAllTo(x);
      return *this;
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::PermuteRows(
	const size_t* p, size_t i1, size_t i2) const
    { 
      TMVAssert(i2<=colsize());
      TMVAssert(i1<=i2);
      // This idea of doing the permutations a block at a time
      // is cribbed from the LAPack code.  It does speed things up 
      // quite a bit for large matrices.  On my machine where BLOCKSIZE=64
      // is optimal for most routines, blocks of 32 were optimal here,
      // so I use BLOCKSIZE/2 in general.
      const size_t Nx = rowsize()/PERM_BLOCKSIZE*PERM_BLOCKSIZE;
      if (Nx != 0) {
	for(size_t j=0;j<Nx;) {
	  size_t j2 = j+PERM_BLOCKSIZE;
	  const size_t* pi = p+i1;
	  for(size_t i=i1;i<i2;++i,++pi) {
            TMVAssert(*pi < colsize());
	    Cols(j,j2).SwapRows(i,*pi);
	  }
	  j = j2;
	}
      }
      if (Nx != rowsize()) {
	const size_t* pi = p+i1;
	for(size_t i=i1;i<i2;++i,++pi) {
	  TMVAssert(*pi < colsize());
	  Cols(Nx,rowsize()).SwapRows(i,*pi);
	}
      }
      return *this;
    }

  template <class T, IndexStyle I> 
    const MatrixView<T,I>& MatrixView<T,I>::ReversePermuteRows(
	const size_t* p, size_t i1, size_t i2) const
    { 
      TMVAssert(i2<=colsize());
      TMVAssert(i1<=i2);
      const size_t Nx = rowsize()/PERM_BLOCKSIZE*PERM_BLOCKSIZE;
      if (Nx != 0) {
	for(size_t j=0;j<Nx;) {
	  size_t j2 = j+PERM_BLOCKSIZE;
	  const size_t* pi = p+i2;
	  for(size_t i=i2;i>i1;) {
	    --i; --pi;
	    TMVAssert(*pi < colsize());
	    Cols(j,j2).SwapRows(i,*pi);
	  }
	  j = j2;
	}
      }
      if (Nx != rowsize()) {
	const size_t* pi = p+i2;
	for(size_t i=i2;i>i1;) {
	  --i; --pi;
	  TMVAssert(*pi < colsize());
	  Cols(Nx,rowsize()).SwapRows(i,*pi);
	}
      }
      return *this;
    }

  //
  // Copy Matrices
  //

  template <bool c1, class T> inline void NonLapCopy(
      const GenMatrix<T>& m1, const MatrixView<T>& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.ct()==NonConj);
    TMVAssert(!(m1.isrm() && m2.isrm()));
    TMVAssert(c1 == m1.isconj());
    if (m1.iscm() && m2.iscm()) 
      if (c1)
	for(size_t j=0;j<m2.rowsize();++j) 
	  DoCopySameType<c1>(m1.col(j),m2.col(j));
      else
	for(size_t j=0;j<m2.rowsize();++j) 
	  memmove(m2.col(j).ptr(),m1.col(j).cptr(),m2.colsize()*sizeof(T));
    else if (m2.colsize() > m2.rowsize())
      for(size_t j=0;j<m2.rowsize();++j) 
	DoCopySameType<c1>(m1.col(j),m2.col(j));
    else 
      for(size_t i=0;i<m2.colsize();++i) 
	DoCopySameType<c1>(m1.row(i),m2.row(i));
  }
#ifdef ELAP
  template <class T> inline void LapCopy(
      const GenMatrix<T>& m1, const MatrixView<T>& m2)
  { NonLapCopy<false>(m1,m2); }
  template <> inline void LapCopy(
      const GenMatrix<double>& m1, const MatrixView<double>& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.iscm());
    TMVAssert(m2.iscm());
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    LAPNAMEX(dlacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
	LAPP(m2.ptr()),LAPV(ld2));
  }
  template <> inline void LapCopy(
      const GenMatrix<complex<double> >& m1,
      const MatrixView<complex<double> >& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.iscm());
    TMVAssert(m2.iscm());
    TMVAssert(m1.ct() == NonConj);
    TMVAssert(m2.ct() == NonConj);
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    LAPNAMEX(zlacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
	LAPP(m2.ptr()),LAPV(ld2));
  }
#ifndef NOFLOAT
  template <> inline void LapCopy(
      const GenMatrix<float>& m1, const MatrixView<float>& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.iscm());
    TMVAssert(m2.iscm());
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    LAPNAMEX(slacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
	LAPP(m2.ptr()),LAPV(ld2));
  }
  template <> inline void LapCopy(
      const GenMatrix<complex<float> >& m1,
      const MatrixView<complex<float> >& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.iscm());
    TMVAssert(m2.iscm());
    TMVAssert(m1.ct() == NonConj);
    TMVAssert(m2.ct() == NonConj);
    char c = 'A';
    int m = m1.colsize();
    int n = m1.rowsize();
    int ld1 = m1.stepj();
    int ld2 = m2.stepj();
    LAPNAMEX(clacpy) (LAPCM LAPV(c),LAPV(m),LAPV(n),LAPP(m1.cptr()),LAPV(ld1),
	LAPP(m2.ptr()),LAPV(ld2));
  }
#endif
#endif
  template <bool c1, class T> void DoCopySameType(const GenMatrix<T>& m1, 
      const MatrixView<T>& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.rowsize() > 0);
    TMVAssert(m2.colsize() > 0);
    TMVAssert(m2.ct() == NonConj);
    TMVAssert(!(m1.isrm() && m2.isrm()));
    TMVAssert(!m2.SameAs(m1));
    TMVAssert(c1 == m1.isconj());

#ifdef ELAP
    if (!c1 && m1.iscm() && m2.iscm() && m1.stepj()>0 && m2.stepj()>0) 
      LapCopy(m1,m2);
    else
#endif
      NonLapCopy<c1>(m1,m2);
  }

  // 
  // Swap
  //

  template <class T> void Swap(
      const MatrixView<T>& m1, const MatrixView<T>& m2)
  {
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    if (m1.stor() == m2.stor() && m1.CanLinearize() && m2.CanLinearize()) {
      TMVAssert(m1.LinearView().size() == m2.LinearView().size());
      TMVAssert(m1.stepi()==m2.stepi() && m1.stepj()==m2.stepj());
      Swap(m1.LinearView(),m2.LinearView());
    }
    else {
      if (m1.isrm() && m2.isrm())
	for(size_t i=0;i<m1.colsize();++i) Swap(m1.row(i),m2.row(i)); 
      else 
	for(size_t j=0;j<m1.rowsize();++j) Swap(m1.col(j),m2.col(j)); 
    }
  }

  //
  // m1 == m2
  //

  template <class T1, class T2> bool operator==(
      const GenMatrix<T1>& m1, const GenMatrix<T2>& m2)
  {
    if (m1.colsize() != m2.colsize()) return false;
    if (m1.rowsize() != m2.rowsize()) return false;
    if (m1.SameAs(m2)) return true;
    if (m1.stepi()==m2.stepi() && m1.stepj()==m2.stepj() &&
	m1.CanLinearize() && m2.CanLinearize())
      return m1.ConstLinearView() == m2.ConstLinearView();
    else {
      for(size_t i=0;i<m1.colsize();++i) 
	if (m1.row(i) != m2.row(i)) return false;
      return true;  
    }
  }

  //
  // I/O
  //

  template <class T> void GenMatrix<T>::Write(ostream& os) const
  {
    const T* mrowi = cptr();
    const int sj = stepj();
    os << colsize() <<"  "<<rowsize()<<endl;
    for(size_t i=colsize();i>0;--i,mrowi+=stepi()) {
      os << "( ";
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=rowsize();k>0;--k,++mij) 
	    os << ' '<<CONJ(*mij)<<' ';
	else 
	  for(size_t k=rowsize();k>0;--k,++mij) 
	    os << ' '<<*mij<<' ';
      else 
	if (isconj()) 
	  for(size_t k=rowsize();k>0;--k,mij+=sj) 
	    os << ' '<<CONJ(*mij)<<' ';
	else 
	  for(size_t k=rowsize();k>0;--k,mij+=sj) 
	    os << ' '<<*mij<<' ';
      os << " )\n";
    }
  }

  template <class T> void GenMatrix<T>::Write(ostream& os,
      RealType(T) thresh) const
  {
    const T* mrowi = cptr();
    const int sj = stepj();
    os << colsize() <<"  "<<rowsize()<<endl;
    for(size_t i=colsize();i>0;--i,mrowi+=stepi()) {
      os << "( ";
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=rowsize();k>0;--k,++mij) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
	else 
	  for(size_t k=rowsize();k>0;--k,++mij) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
      else 
	if (isconj()) 
	  for(size_t k=rowsize();k>0;--k,mij+=sj) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : *mij)<<' ';
	else 
	  for(size_t k=rowsize();k>0;--k,mij+=sj) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : *mij)<<' ';
      os << " )\n";
    }
  }

  template <class T> void MatrixReadError<T>::Write(ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for Matrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (m.get() && cs != m->colsize()) {
      os<<"Wrong column size: expected "<<m->colsize()<<", got "<<cs<<".\n";
    }
    if (m.get() && rs != m->rowsize()) {
      os<<"Wrong row size: expected "<<m->rowsize()<<", got "<<rs<<".\n";
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
      os<<"The portion of the Matrix which was successfully read is: \n";
      for(size_t ii=0;ii<i;++ii) {
	os<<"( ";
	for(size_t jj=0;jj<m->rowsize();++jj)
	  os<<' '<<(*m)(ii,jj)<<' ';
	os<<" )\n";
      }
      os<<"( ";
      for(size_t jj=0;jj<j;++jj)
	os<<' '<<(*m)(i,jj)<<' ';
      os<<" )\n";
    }
  }

  template <class T, IndexStyle I> void MatrixView<T,I>::Read(
      istream& is) const
  {
    T* mrowi = ptr();
    const int sj = stepj();
    char paren;
    for(size_t i=0;i<colsize();++i,mrowi+=stepi()) {
      is >> paren;
      if (!is || paren != '(') 
	throw MatrixReadError<T>(i,0,*this,is,'(',is?paren:'(');
      T* mij = mrowi;
      if (isrm()) 
	for(size_t k=rowsize();k>0;--k,++mij) {
	  is >> *mij;
	  if (!is) 
	    throw MatrixReadError<T>(i,rowsize()-k,*this,is);
	}
      else
	for(size_t k=rowsize();k>0;--k,mij+=sj) {
	  is >> *mij;
	  if (!is) 
	    throw MatrixReadError<T>(i,rowsize()-k,*this,is);
	}
      is >> paren;
      if ((!is && i+1<colsize())  || paren != ')') 
	throw MatrixReadError<T>(i,rowsize(),*this,is,')',is?paren:')');
    }
    if (isconj()) ConjugateSelf();
  }

  template <class T, StorageType S, IndexStyle I> istream& operator>>(
      istream& is, auto_ptr<Matrix<T,S,I> >& m)
  { 
    size_t cs,rs;
    is >> cs >> rs;
    if (!is) 
      throw MatrixReadError<T>(is);
    m.reset(new Matrix<T,S,I>(cs,rs));
    m->View().Read(is); 
    return is;
  }

  template <class T> istream& operator>>(istream& is, const MatrixView<T>& m)
  { 
    size_t cs,rs;
    is >> cs >> rs;
    if (!is) 
      throw MatrixReadError<T>(is);
    if (m.colsize() != cs || m.rowsize() != rs) 
      throw MatrixReadError<T>(m,is,cs,rs);
    TMVAssert(m.colsize() == cs);
    TMVAssert(m.rowsize() == rs);
    m.Read(is);
    return is;
  }


#define InstFile "TMV_Matrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


