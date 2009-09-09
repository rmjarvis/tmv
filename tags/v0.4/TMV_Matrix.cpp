
#include "TMV.h"

namespace tmv {

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

  template <class T> RefType(T) MatrixView<T>::ref(size_t i, size_t j) const
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
      case tmv::LU : this->itsdiv = new LUDiv<T>(*this); break;
      case tmv::QR : this->itsdiv = new QRDiv<T>(*this); break;
      case tmv::QRP : this->itsdiv = new QRPDiv<T>(*this); break;
      case tmv::SV : this->itsdiv = new SVDiv<T>(*this); break;
      case tmv::SVF : this->itsdiv = new SVFDiv<T>(*this,true,true); break;
      case tmv::SVS : this->itsdiv = new SVFDiv<T>(*this,false,false); break;
      case tmv::SVU : this->itsdiv = new SVFDiv<T>(*this,true,false); break;
      case tmv::SVV : this->itsdiv = new SVFDiv<T>(*this,false,true); break;
      default : TMVAssert(false);
    }
  }

#define NEW_SIZE(cs,rs) GenMatrix<T>(ColMajor,NonConj), \
  itslen((cs)*(rs)), itsm(new T[itslen]), itscs(cs), itsrs(rs) \
  DEFFIRSTLAST(itsm,itsm+itslen)

  template <class T> Matrix<T,ColMajor>::Matrix(const vector<vector<T> >& vv) :
    NEW_SIZE(vv.size(),(vv.size()>0?vv[0].size():0))
  {
    T* vi=itsm;
    for(size_t j=0;j<rowsize();++j) {
      for(size_t i=0;i<colsize();++i,++vi) {
	TMVAssert(vv[i].size() == rowsize());
	*vi = vv[i][j];
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
      cout<<"first col element ("<<i1<<") must be in 0 -- "<<colsize()<<endl;
    }
    if (i2-istep < 0 || i2-istep >= int(colsize())) {
      ok = false;
      cout<<"last col element ("<<i2-istep<<") must be in 0 -- ";
      cout<<colsize()<<endl;
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
      cout<<"first row element ("<<j1<<") must be in 0 -- "<<rowsize()<<endl;
    }
    if (j2-jstep < 0 || j2-jstep >= int(rowsize())) {
      ok = false;
      cout<<"last row element ("<<j2-jstep<<") must be in 0 -- ";
      cout<<rowsize()<<endl;
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
      cout<<"i ("<<i<<") must be in 0 -- "<<colsize()<<endl;
    }
    if (j >= rowsize()) {
      ok = false;
      cout<<"j ("<<j<<") must be in 0 -- "<<rowsize()<<endl;
    }
    int i2 = int(i)+istep*int(size-1);
    int j2 = int(j)+jstep*int(size-1);
    if (i2 < 0 || i2 >= int(colsize())) {
      ok = false;
      cout<<"last element's i ("<<i2<<") must be in 0 -- "<<colsize()<<endl;
    }
    if (j2 < 0 || j2 >= int(rowsize())) {
      ok = false;
      cout<<"last element's j ("<<j2<<") must be in 0 -- "<<rowsize()<<endl;
    }
    return ok;
  }

  //
  // Norms
  //

  template <class T> RealType(T) GenMatrix<T>::NormSq() const
  {
    RealType(T) sum(0);
    if (isrm()) 
      for(size_t i=0;i<colsize();++i) 
	sum += row(i).NormSq();
    else 
      for(size_t j=0;j<rowsize();++j) 
	sum += col(j).NormSq();
    return sum;
  }

  template <class T> RealType(T) NonLapMaxAbsElement(const GenMatrix<T>& m)
  {
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

  template <class T> RealType(T) NonLapNorm1(const GenMatrix<T>& m)
  { 
    RealType(T) max(0);
    for(size_t j=0;j<m.rowsize();++j) {
      RealType(T) temp = m.col(j).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  }

  template <class T> RealType(T) NonLapNormF(const GenMatrix<T>& m)
  { return tmv::SQRT(m.NormSq()); }

  template <class T> RealType(T) NonLapNormInf(const GenMatrix<T>& m)
  { return NonLapNorm1(m.Transpose()); }

#ifdef LAP
  template <class T> RealType(T) LapNorm(const char c, const GenMatrix<T>& m)
  { 
    switch(c) {
      case 'M' : return NonLapMaxAbsElement(m);
      case '1' : return NonLapNorm1(m);
      case 'F' : return NonLapNormF(m);
      case 'I' : return NonLapNormInf(m);
      default : TMVAssert(false); 
    }
    return RealType(T)(0);
  }
  template <> double LapNorm(const char c, const GenMatrix<double>& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    double* work = (c == 'I' ? new double[M] : 0);
    double norm = dlange(&cc,&M,&N,const_cast<double*>(m.cptr()),&lda,work);
    if (work) delete[] work;
    return norm;
  }
  template <> double LapNorm(const char c, const GenMatrix<complex<double> >& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    double* work = (c == 'I' ? new double[M] : 0);
    double norm = zlange(&cc,&M,&N,LAP_Complex(m.cptr()),&lda,work);
    if (work) delete[] work;
    return norm;
  }
#ifndef NOFLOAT
  template <> float LapNorm(const char c, const GenMatrix<float>& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    float* work = (c == 'I' ? new float[M] : 0);
    float norm = slange(&cc,&M,&N,const_cast<float*>(m.cptr()),&lda,work);
    if (work) delete[] work;
    return norm;
  }
  template <> float LapNorm(const char c, const GenMatrix<complex<float> >& m)
  { 
    TMVAssert(m.iscm());
    char cc = c;
    int M = m.colsize();
    int N = m.rowsize();
    int lda = m.stepj();
    float* work = (c == 'I' ? new float[M] : 0);
    float norm = clange(&cc,&M,&N,LAP_Complex(m.cptr()),&lda,work);
    if (work) delete[] work;
    return norm;
  }
#endif
#endif
  template <class T> RealType(T) GenMatrix<T>::MaxAbsElement() const
  {
#ifdef LAP
    if (isrm()) return LapNorm('M',QuickTranspose());
    else if (iscm()) return LapNorm('M',*this);
    else
#endif
      return NonLapMaxAbsElement(*this);
  }
  template <class T> RealType(T) GenMatrix<T>::Norm1() const
  {
#ifdef LAP
    if (isrm()) return LapNorm('I',QuickTranspose());
    else if (iscm()) return LapNorm('1',*this);
    else
#endif
      return NonLapNorm1(*this);
  }
  template <class T> RealType(T) GenMatrix<T>::NormF() const
  {
#ifdef LAP
    if (isrm()) return LapNorm('F',QuickTranspose());
    else if (iscm()) return LapNorm('F',*this);
    else
#endif
      return NonLapNormF(*this);
  }

  //
  // Modifying Functions
  //

  template <class T> 
    const MatrixView<T>& MatrixView<T>::SetAllTo(T x) const
    { 
      if (isrm())
	for(size_t i=0;i<colsize();++i) row(i).SetAllTo(x); 
      else 
	for(size_t j=0;j<rowsize();++j) col(j).SetAllTo(x); 
      return *this; 
    }

  template <class T> 
    const MatrixView<T>& MatrixView<T>::TransposeSelf() const
    {
      TMVAssert(IsSquare());
      for(size_t i=1;i<colsize();++i) Swap(row(i,0,i),col(i,0,i));
      return *this; 
    }

  template <class T> 
    const MatrixView<T>& MatrixView<T>::ConjugateSelf() const
    { 
      if (IsComplex(T())) {
	if (isrm()) 
	  for(size_t i=0;i<colsize();++i) row(i).ConjugateSelf();
	else 
	  for(size_t j=0;j<rowsize();++j) col(j).ConjugateSelf();
      }
      return *this; 
    }

  template <class T> 
    const MatrixView<T>& MatrixView<T>::SetToIdentity(T x) const 
    {
      TMVAssert(IsSquare());
      Zero();
      diag().SetAllTo(x);
      return *this;
    }

  template <class T> 
    const MatrixView<T>& MatrixView<T>::SwapRows(
	size_t i1, size_t i2) const
    { 
      TMVAssert(i1 < colsize() && i2 < colsize());
      if (i1!=i2) Swap(row(i1),row(i2));
      return *this;
    }

  template <class T> 
    const MatrixView<T>& MatrixView<T>::SwapCols(
	size_t j1, size_t j2) const
    { 
      TMVAssert(j1 < rowsize() && j2 < rowsize());
      if (j1!=j2) Swap(col(j1),col(j2));
      return *this;
    }


  // 
  // Swap
  //

  template <class T> inline void Swap(
      const MatrixView<T>& m1, const MatrixView<T>& m2)
  {
    TMVAssert(m1.colsize() == m2.colsize());
    TMVAssert(m1.rowsize() == m2.rowsize());
    if (m1.isrm() && m2.isrm())
      for(size_t i=0;i<m1.colsize();++i) Swap(m1.row(i),m2.row(i)); 
    else 
      for(size_t j=0;j<m1.rowsize();++j) Swap(m1.col(j),m2.col(j)); 
  }

  //
  // m1 == m2
  //

  template <class T> bool operator==(
      const GenMatrix<T>& m1, const GenMatrix<T>& m2)
  {
    if (m1.colsize() != m2.colsize()) return false;
    if (m1.rowsize() != m2.rowsize()) return false;
    if (m1.SameAs(m2)) return true;
    for(size_t i=0;i<m1.colsize();++i) if (m1.row(i) != m2.row(i)) return false;
    return true;  
  }

  //
  // I/O
  //

  template <class T> void GenMatrix<T>::Write(ostream& os) const
  {
    const T* mptr = cptr();
    os << colsize() <<"  "<<rowsize()<<endl;
    for(size_t i=0;i<colsize();++i,mptr+=stepi()) {
      os << "( ";
      if (isrm()) {
	if (isconj()) {
	  CVIt<T,Unit,Conj> it(mptr,1);
	  const CVIt<T,Unit,Conj>  _end = it+rowsize();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	} else {
	  CVIt<T,Unit,NonConj> it(mptr,1);
	  const CVIt<T,Unit,NonConj>  _end = it+rowsize();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	}
      } else {
	if (isconj()) {
	  CVIt<T,Step,Conj> it(mptr,stepj());
	  const CVIt<T,Step,Conj>  _end = it+rowsize();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	} else {
	  CVIt<T,Step,NonConj> it(mptr,stepj());
	  const CVIt<T,Step,NonConj>  _end = it+rowsize();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	}
      }
      os << " )\n";
    }
  }

  template <class T> void MatrixView<T>::Read(istream& is) const
  {
    T* mptr = ptr();
    char paren;
    for(size_t i=0;i<colsize();++i,mptr+=stepi()) {
      is >> paren;
      if (!is || paren != '(') tmv_error("reading ( in Matrix::Read");
      if (isrm()) {
	if (isconj()) {
	  VIt<T,Unit,Conj> it(mptr,1 FIRSTLAST);
	  const VIt<T,Unit,Conj>  _end = it+rowsize();
	  for(; it!=_end; ++it) is >> *it;
	} else {
	  VIt<T,Unit,NonConj> it(mptr,1 FIRSTLAST);
	  const VIt<T,Unit,NonConj>  _end = it+rowsize();
	  for(; it!=_end; ++it) is >> *it;
	}
      } else {
	if (isconj()) {
	  VIt<T,Step,Conj> it(mptr,stepj() FIRSTLAST);
	  const VIt<T,Step,Conj>  _end = it+rowsize();
	  for(; it!=_end; ++it) is >> *it;
	} else {
	  VIt<T,Step,NonConj> it(mptr,stepj() FIRSTLAST);
	  const VIt<T,Step,NonConj>  _end = it+rowsize();
	  for(; it!=_end; ++it) is >> *it;
	}
      }
      if (!is) tmv_error("reading values in Matrix::Read");
      is >> paren;
      if ((!is && i+1<colsize())  || paren != ')') 
	tmv_error("reading ) in Matrix::Read");
    }
  }

  template <class T, StorageType S> istream& operator>>(istream& is, 
      Matrix<T,S>* m)
  { 
    size_t cs,rs;
    is >> cs >> rs;
    if (!is) tmv_error("reading size in Matrix::Read");
    m = new Matrix<T,S>(cs,rs);
    m->QuickView().Read(is); 
    return is;
  }

  template <class T> istream& operator>>(istream& is,
      const MatrixView<T>& m)
  { 
    size_t cs,rs;
    is >> cs >> rs;
    if (!is) tmv_error("reading size in Matrix::Read");
    if (m.colsize() != cs) tmv_error("colsize wrong on Matrix::Read");
    if (m.rowsize() != rs) tmv_error("rowsize wrong on Matrix::Read");
    m.Read(is);
    return is;
  }


#define InstFile "TMV_Matrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


