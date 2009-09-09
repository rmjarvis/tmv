
#include "TMV.h"
#include "TMV_Tri.h"

namespace tmv {

  //
  // Access
  //

  template <class T> T GenUpperTriMatrix<T>::cref(size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    if (i==j && isunit()) return T(1);
    else if (i>j) return T(0);
    else {
      const T* mi;
      if (isrm()) mi = cptr() + int(i)*stepi() + j;
      else if (iscm()) mi = cptr() + i + int(j)*stepj();
      else mi = cptr() + int(i)*stepi() + int(j)*stepj();
      return (isconj() ? CONJ(*mi) : *mi);
    }
  }

  template <class T, IndexStyle I> RefType(T) UpperTriMatrixView<T,I>::ref(
      size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    TMVAssert(this->okij(i,j));
    T* mi;
    if (isrm()) mi = ptr() + int(i)*stepi() + j;
    else if (iscm()) mi = ptr() + i + int(j)*stepj();
    else mi = ptr() + int(i)*stepi() + int(j)*stepj();
#ifdef TMVFLDEBUG
    TMVAssert(mi>=first);
    TMVAssert(mi<last);
#endif
    return REF(mi,ct());
  }

  template <class T> void GenUpperTriMatrix<T>::NewDivider() const
  {
    this->itsdiv.reset(new UpperTriDiv<T>(*this)); 
  }

  template <class T> void GenUpperTriMatrix<T>::Inverse(
      const UpperTriMatrixView<T>& minv) const
  {
    this->SetDiv();
    const UpperTriDiv<T>* tdiv =
      dynamic_cast<const UpperTriDiv<T>*>(this->GetDiv());
    TMVAssert(tdiv);
    tdiv->Inverse(minv);
  }

  template <class T> template <class T1> void GenUpperTriMatrix<T>::LDivEq(
      const UpperTriMatrixView<T1>& m) const
  {
    this->SetDiv();
    const UpperTriDiv<T>* tdiv =
      dynamic_cast<const UpperTriDiv<T>*>(this->GetDiv());
    TMVAssert(tdiv);
    tdiv->LDivEq(m);
  }

  template <class T> template <class T1> void GenUpperTriMatrix<T>::RDivEq(
      const UpperTriMatrixView<T1>& m) const
  {
    this->SetDiv();
    const UpperTriDiv<T>* tdiv =
      dynamic_cast<const UpperTriDiv<T>*>(this->GetDiv());
    TMVAssert(tdiv);
    tdiv->RDivEq(m);
  }

  template <class T> template <class T1, class T0> 
    void GenUpperTriMatrix<T>::LDiv(
      const GenUpperTriMatrix<T1>& m1,
      const UpperTriMatrixView<T0>& m0) const
    {
      this->SetDiv();
      const UpperTriDiv<T>* tdiv =
	dynamic_cast<const UpperTriDiv<T>*>(this->GetDiv());
      TMVAssert(tdiv);
      tdiv->LDiv(m1,m0);
    }

  template <class T> template <class T1, class T0>
    void GenUpperTriMatrix<T>::RDiv(
      const GenUpperTriMatrix<T1>& m1,
      const UpperTriMatrixView<T0>& m0) const
    {
      this->SetDiv();
      const UpperTriDiv<T>* tdiv =
	dynamic_cast<const UpperTriDiv<T>*>(this->GetDiv());
      TMVAssert(tdiv);
      tdiv->RDiv(m1,m0);
    }


  template <class T> T GenLowerTriMatrix<T>::cref(size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    if (i==j && isunit()) return T(1);
    else if (i<j) return T(0);
    else {
      const T* mi;
      if (isrm()) mi = cptr() + int(i)*stepi() + j;
      else if (iscm()) mi = cptr() + i + int(j)*stepj();
      else mi = cptr() + int(i)*stepi() + int(j)*stepj();
      return (isconj() ? CONJ(*mi) : *mi);
    }
  }

  template <class T, IndexStyle I> RefType(T) LowerTriMatrixView<T,I>::ref(
      size_t i, size_t j) const
  {
    TMVAssert(i < size());
    TMVAssert(j < size());
    TMVAssert(this->okij(i,j));
    T* mi;
    if (isrm()) mi = ptr() + int(i)*stepi() + j;
    else if (iscm()) mi = ptr() + i + int(j)*stepj();
    else mi = ptr() + int(i)*stepi() + int(j)*stepj();
#ifdef TMVFLDEBUG
    TMVAssert(mi>=first);
    TMVAssert(mi<last);
#endif
    return REF(mi,ct());
  }

  template <class T> void GenLowerTriMatrix<T>::NewDivider() const
  {
    this->itsdiv.reset(new LowerTriDiv<T>(*this)); 
  }

  template <class T> void GenLowerTriMatrix<T>::Inverse(
      const LowerTriMatrixView<T>& minv) const
  {
    this->SetDiv();
    const LowerTriDiv<T>* tdiv =
      dynamic_cast<const LowerTriDiv<T>*>(this->GetDiv());
    TMVAssert(tdiv);
    tdiv->Inverse(minv);
  }

  template <class T> template <class T1> void GenLowerTriMatrix<T>::LDivEq(
      const LowerTriMatrixView<T1>& m) const
  {
    this->SetDiv();
    const LowerTriDiv<T>* tdiv =
      dynamic_cast<const LowerTriDiv<T>*>(this->GetDiv());
    TMVAssert(tdiv);
    tdiv->LDivEq(m);
  }

  template <class T> template <class T1> void GenLowerTriMatrix<T>::RDivEq(
      const LowerTriMatrixView<T1>& m) const
  {
    this->SetDiv();
    const LowerTriDiv<T>* tdiv =
      dynamic_cast<const LowerTriDiv<T>*>(this->GetDiv());
    TMVAssert(tdiv);
    tdiv->RDivEq(m);
  }

  template <class T> template <class T1, class T0> 
    void GenLowerTriMatrix<T>::LDiv(
      const GenLowerTriMatrix<T1>& m1,
      const LowerTriMatrixView<T0>& m0) const
    {
      this->SetDiv();
      const LowerTriDiv<T>* tdiv =
	dynamic_cast<const LowerTriDiv<T>*>(this->GetDiv());
      TMVAssert(tdiv);
      tdiv->LDiv(m1,m0);
    }

  template <class T> template <class T1, class T0>
    void GenLowerTriMatrix<T>::RDiv(
      const GenLowerTriMatrix<T1>& m1,
      const LowerTriMatrixView<T0>& m0) const
    {
      this->SetDiv();
      const LowerTriDiv<T>* tdiv =
	dynamic_cast<const LowerTriDiv<T>*>(this->GetDiv());
      TMVAssert(tdiv);
      tdiv->RDiv(m1,m0);
    }

  //
  // OK? (SubMatrix, etc.)
  //

  template <class T> bool GenUpperTriMatrix<T>::OKSubMatrix(
      int i1, int i2, int j1, int j2, int istep, int jstep) const
  { 
    if (i1==i2 || j1==j2) return true; // no elements, so whatever...
    bool ok = true;
    if (istep == 0) {
      ok = false;
      cout<<"istep ("<<istep<<") can not be 0\n";
    }
    if (i1 < 0 || i1 >= int(size())) {
      ok = false;
      cout<<"first col element ("<<i1<<") must be in 0 -- "<<size()-1<<endl;
    }
    if (i2-istep < 0 || i2-istep >= int(size())) {
      ok = false;
      cout<<"last col element ("<<i2-istep<<") must be in 0 -- "<<size()-1<<endl;
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
    if (!this->okij(i1,j1)) {
      ok = false;
      cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i1,j2-jstep)) {
      ok = false;
      cout<<"Upper right corner ("<<i1<<','<<j2-jstep<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i2-istep,j1)) {
      ok = false;
      cout<<"Lower left corner ("<<i2-istep<<','<<j1<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i2-istep,j2-jstep)) {
      ok = false;
      cout<<"Lower right corner ("<<i2-istep<<','<<j2-jstep<<") must be in Upper Triangle\n";
    }
    return ok;
  }

  template <class T> bool GenUpperTriMatrix<T>::OKSubVector(
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
    if (!this->okij(i,j)) {
      ok = false;
      cout<<"First element ("<<i<<','<<j<<") must be in Triangle\n";
    }
    if (!this->okij(i2,j2)) {
      ok = false;
      cout<<"Last element ("<<i2<<','<<j2<<") must be in Triangle\n";
    }
    return ok;
  }

  template <class T> bool GenUpperTriMatrix<T>::OKSubTriMatrix(
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

  template <class T> bool ConstUpperTriMatrixView<T,FortranStyle>::OKSubMatrix(
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
      cout<<"last col element ("<<i2<<") must be in 1 -- "<<this->size()<<endl;
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
      cout<<"last row element ("<<j2<<") must be in 1 -- "<<this->size()<<endl;
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
      cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i1-1,j2-1)) {
      ok = false;
      cout<<"Upper right corner ("<<i1<<','<<j2<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i2-1,j1-1)) {
      ok = false;
      cout<<"Lower left corner ("<<i2<<','<<j1<<") must be in Upper Triangle\n";
    }
    if (!this->okij(i2-1,j2-1)) {
      ok = false;
      cout<<"Lower right corner ("<<i2<<','<<j2<<") must be in Upper Triangle\n";
    }
    return ok;
  }

  template <class T> bool ConstUpperTriMatrixView<T,FortranStyle>::OKSubVector(
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
    if (i < 1 || j >= this->size()) {
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
    if (!this->okij(i-1,j-1)) {
      ok = false;
      cout<<"First element ("<<i<<','<<j<<") must be in Triangle\n";
    }
    if (!this->okij(i2-1,j2-1)) {
      ok = false;
      cout<<"Last element ("<<i2<<','<<j2<<") must be in Triangle\n";
    }
    return ok;
  }

  template <class T> bool ConstUpperTriMatrixView<T,FortranStyle>::OKSubTriMatrix(
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
    if (i2<1 || i2 > int(this->size())) {
      ok = false;
      cout<<"last diag element ("<<i2<<") must be in 1 -- "<<this->size()<<endl;
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
  // Norms
  //

  template <class T> RealType(T) GenUpperTriMatrix<T>::NormSq() const
  {
    RealType(T) sum(0);
    if (isrm()) 
      if (isunit())
	for(size_t i=0;i<size();++i) 
	  sum += row(i,i+1,size()).NormSq();
      else
	for(size_t i=0;i<size();++i) 
	  sum += row(i,i,size()).NormSq();
    else
      if (isunit())
	for(size_t j=0;j<size();++j) 
	  sum += col(j,0,j).NormSq();
      else
	for(size_t j=0;j<size();++j) 
	  sum += col(j,0,j+1).NormSq();
    if (isunit()) sum += size();
    return sum;
  }

  template <class T> inline RealType(T) NonLapNormF(
      const GenUpperTriMatrix<T>& m)
  { return SQRT(m.NormSq()); }

  template <class T> inline RealType(T) NonLapMaxAbsElement(
      const GenUpperTriMatrix<T>& m)
  {
    RealType(T) max(0);
    if (m.isrm())
      for(size_t i=0;i<m.size();++i) {
	RealType(T) temp;
	if (m.isunit())
	  temp = m.row(i,i+1,m.size()).NormInf();
	else 
	  temp = m.row(i,i,m.size()).NormInf();
	if (temp > max) max = temp;
      }
    else
      for(size_t j=0;j<m.size();++j) {
	RealType(T) temp;
	if (m.isunit())
	  temp = m.col(j,0,j).NormInf();
	else 
	  temp = m.col(j,0,j+1).NormInf();
	if (temp > max) max = temp;
      }
    if (m.isunit() && max < RealType(T)(1)) max = RealType(T)(1);
    return max;
  }

  template <class T> inline RealType(T) NonLapNorm1(
      const GenUpperTriMatrix<T>& m)
  { 
    RealType(T) max(0);
    for(size_t j=0;j<m.size();++j) {
      RealType(T) temp;
      if (m.isunit()) {
	temp = m.col(j,0,j).Norm1();
	temp += RealType(T)(1);
      } else temp = m.col(j,0,j+1).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  } 

  template <class T> inline RealType(T) NonLapNormInf(
      const GenUpperTriMatrix<T>& m)
  { 
    RealType(T) max(0);
    for(size_t j=0;j<m.size();++j) {
      RealType(T) temp;
      if (m.isunit()) {
	temp = m.row(j,j+1,m.size()).Norm1();
	temp += RealType(T)(1);
      } else temp = m.row(j,j,m.size()).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  }
#ifdef XLAP
  template <class T> inline RealType(T) LapNorm(
      const char c, const GenUpperTriMatrix<T>& m)
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
      const char c, const GenUpperTriMatrix<double>& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    if (m.isrm()) {
      if (c == '1') cc = 'I';
      else if (c == 'I') cc = '1';
    }
    int M = m.size();
    int N = M;
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<double> work(cc == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAMEX(dlantr) (LAPCM LAPV(cc),
	m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
	LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
	LAP1 LAP1 LAP1);
    return norm;
  }
  template <> inline double LapNorm(
      const char c, const GenUpperTriMatrix<complex<double> >& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    if (m.isrm()) {
      if (c == '1') cc = 'I';
      else if (c == 'I') cc = '1';
    }
    int M = m.size();
    int N = M;
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<double> work(cc == 'I' ? new double[M] : 0);
#endif
    double norm = LAPNAMEX(zlantr) (LAPCM LAPV(cc),
	m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
	LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
	LAP1 LAP1 LAP1);
    return norm;
  }
#ifndef NOFLOAT
  template <> inline float LapNorm(
      const char c, const GenUpperTriMatrix<float>& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    if (m.isrm()) {
      if (c == '1') cc = 'I';
      else if (c == 'I') cc = '1';
    }
    int M = m.size();
    int N = M;
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<float> work(cc == 'I' ? new float[M] : 0);
#endif
    float norm = LAPNAMEX(slantr) (LAPCM LAPV(cc),
	m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
	LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
	LAP1 LAP1 LAP1);
    return norm;
  }
  template <> inline float LapNorm(
      const char c, const GenUpperTriMatrix<complex<float> >& m)
  {
    TMVAssert(m.iscm() || m.isrm());
    char cc = c;
    if (m.isrm()) {
      if (c == '1') cc = 'I';
      else if (c == 'I') cc = '1';
    }
    int M = m.size();
    int N = M;
    int lda = m.iscm() ? m.stepj() : m.stepi();
#ifndef LAPNOWORK
    auto_array<float> work(cc == 'I' ? new float[M] : 0);
#endif
    float norm = LAPNAMEX(clantr) (LAPCM LAPV(cc),
	m.iscm() ? LAPCH_UP : LAPCH_LO, m.isunit() ? LAPCH_U : LAPCH_NU,
	LAPV(M),LAPV(N),LAPP(m.cptr()),LAPV(lda) LAPWK(work.get())
	LAP1 LAP1 LAP1);
    return norm;
  }
#endif
#endif // XLAP

  template <class T> RealType(T) GenUpperTriMatrix<T>::MaxAbsElement() const
  {
#ifdef XLAP
    return LapNorm('M',*this);
#else
    return NonLapMaxAbsElement(*this);
#endif
  }
  template <class T> RealType(T) GenUpperTriMatrix<T>::Norm1() const
  {
#ifdef XLAP
    return LapNorm('1',*this);
#else
    return NonLapNorm1(*this);
#endif
  }
  template <class T> RealType(T) GenUpperTriMatrix<T>::NormInf() const
  {
#ifdef XLAP
    return LapNorm('I',*this);
#else
    return NonLapNormInf(*this);
#endif
  }
  template <class T> RealType(T) GenUpperTriMatrix<T>::NormF() const
  {
#ifdef XLAP
    return LapNorm('F',*this);
#else
    return NonLapNormF(*this);
#endif
  }

  //
  // Modifying Functions
  //

  template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::SetAllTo(T x) const
    { 
      const size_t N = size();

      if (isrm())
	if (isunit())
	  for(size_t i=0;i<N;++i) row(i,i+1,N).SetAllTo(x); 
	else
	  for(size_t i=0;i<N;++i) row(i,i,N).SetAllTo(x); 
      else 
	if (isunit())
	  for(size_t j=0;j<N;++j) col(j,0,j).SetAllTo(x); 
	else
	  for(size_t j=0;j<N;++j) col(j,0,j+1).SetAllTo(x); 
      return *this; 
    }

  template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::Clip(
	RealType(T) thresh) const
    { 
      const size_t N = size();

      if (isrm())
	if (isunit())
	  for(size_t i=0;i<N;++i) row(i,i+1,N).Clip(thresh);
	else
	  for(size_t i=0;i<N;++i) row(i,i,N).Clip(thresh);
      else 
	if (isunit())
	  for(size_t j=0;j<N;++j) col(j,0,j).Clip(thresh);
	else
	  for(size_t j=0;j<N;++j) col(j,0,j+1).Clip(thresh);
      return *this; 
    }

  template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::ConjugateSelf() const
    { 
      const size_t N = size();

      if (IsComplex(T())) {
	if (isrm())
	  if (isunit())
	    for(size_t i=0;i<N;++i) row(i,i+1,N).ConjugateSelf();
	  else
	    for(size_t i=0;i<N;++i) row(i,i,N).ConjugateSelf();
	else
	  if (isunit())
	    for(size_t j=0;j<N;++j) col(j,0,j).ConjugateSelf();
	  else
	    for(size_t j=0;j<N;++j) col(j,0,j+1).ConjugateSelf();
      }
      return *this; 
    }

  template <class T, IndexStyle I> 
    const UpperTriMatrixView<T,I>& UpperTriMatrixView<T,I>::SetToIdentity(
	T x) const 
    {
      TMVAssert(!isunit() || x==T(1));
      Zero();
      if (!isunit()) diag().SetAllTo(x);
      return *this;
    }

  //
  // Swap
  //

  template <class T> void Swap(
      const UpperTriMatrixView<T>& m1, const UpperTriMatrixView<T>& m2)
  {
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.dt() == m2.dt());
    const size_t N = m1.size();

    if (m1.isrm() && m2.isrm())
      if (m1.isunit()) 
	for(size_t i=0;i<N;++i) 
	  Swap(m1.row(i,i+1,N), m2.row(i,i+1,N));
      else
	for(size_t i=0;i<N;++i) 
	  Swap(m1.row(i,i,N), m2.row(i,i,N));
    else
      if (m1.isunit()) 
	for(size_t j=0;j<N;++j) 
	  Swap(m1.col(j,0,j), m2.col(j,0,j));
      else
	for(size_t j=0;j<N;++j) 
	  Swap(m1.col(j,0,j+1), m2.col(j,0,j+1));
  }

  //
  // m1 == m2
  //

  template <class T1, class T2> bool operator==(
      const GenUpperTriMatrix<T1>& m1, const GenUpperTriMatrix<T2>& m2)
  {
    if (m1.size() != m2.size()) return false;
    if (m1.SameAs(m2)) return true;
    for(size_t j=0;j<m1.size();++j) {
      if (m1.col(j,0,j) != m2.col(j,0,j)) return false;
    }

    if (m1.isunit() && !m2.isunit()) {
      for(size_t i=0;i<m1.size();++i) if (m2(i,i) != T2(1)) return false;
    }
    else if (m2.isunit() && !m1.isunit()) {
      for(size_t i=0;i<m1.size();++i) if (m1(i,i) != T1(1)) return false;
    }
    else if (!m1.isunit() && !m2.isunit()) {
      if (m1.diag() != m2.diag()) return false;
    }

    return true;  
  }

  //
  // I/O
  //

  template <class T> void GenUpperTriMatrix<T>::Write(ostream& os) const
  {
    const T* mrowi = cptr();
    const int sj = stepj();
    const int ds = stepi()+sj;
    size_t len = size();
    if (isunit()) {
      mrowi += sj;
      --len;
    }
    os << size() <<' '<<size() << endl;
    for(size_t i=0;i<size();++i,--len,mrowi+=ds) {
      os << "( ";
      for(size_t j=0;j<i;j++) os << ' '<<T(0)<<' ';
      if (isunit()) os << ' '<<T(1)<<' ';
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,++mij) os << ' '<<CONJ(*mij)<<' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) os << ' '<<*mij<<' ';
      else 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj) os << ' '<<CONJ(*mij)<<' ';
	else 
	  for(size_t k=len;k>0;--k,mij+=sj) os << ' '<<*mij<<' ';
      os << " )\n";
    }
  }

  template <class T> void GenUpperTriMatrix<T>::Write(ostream& os,
      RealType(T) thresh) const
  {
    const T* mrowi = cptr();
    const int sj = stepj();
    const int ds = stepi()+sj;
    size_t len = size();
    if (isunit()) {
      mrowi += sj;
      --len;
    }
    os << size() <<' '<<size() << endl;
    for(size_t i=0;i<size();++i,--len,mrowi+=ds) {
      os << "( ";
      for(size_t j=0;j<i;j++) os << ' '<<T(0)<<' ';
      if (isunit()) os << ' '<<T(1)<<' ';
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : *mij)<<' ';
      else 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
	else 
	  for(size_t k=len;k>0;--k,mij+=sj) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : *mij)<<' ';
      os << " )\n";
    }
  }

  template <class T> void GenUpperTriMatrix<T>::WriteCompact(ostream& os) const
  {
    const T* mrowi = cptr();
    const int sj = stepj();
    const int ds = stepi()+sj;
    size_t len = size();
    if (isunit()) {
      mrowi += sj;
      --len;
    }
    os << "U " << size() << endl;
    for(size_t i=0;i<size();++i,--len,mrowi+=ds) {
      os << "( ";
      if (isunit()) os << ' '<<T(1)<<' ';
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,++mij) os << ' '<<CONJ(*mij)<<' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) os << ' '<<*mij<<' ';
      else 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj) os << ' '<<CONJ(*mij)<<' ';
	else 
	  for(size_t k=len;k>0;--k,mij+=sj) os << ' '<<*mij<<' ';
      os << " )\n";
    }
  }

  template <class T> void GenUpperTriMatrix<T>::WriteCompact(ostream& os,
      RealType(T) thresh) const
  {
    const T* mrowi = cptr();
    const int sj = stepj();
    const int ds = stepi()+sj;
    size_t len = size();
    if (isunit()) {
      mrowi += sj;
      --len;
    }
    os << "U " << size() << endl;
    for(size_t i=0;i<size();++i,--len,mrowi+=ds) {
      os << "( ";
      if (isunit()) os << ' '<<T(1)<<' ';
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : *mij)<<' ';
      else 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
	else 
	  for(size_t k=len;k>0;--k,mij+=sj) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : *mij)<<' ';
      os << " )\n";
    }
  }

  template <class T> void UpperTriMatrixReadError<T>::Write(
      ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for UpperTriMatrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (unitgot != T(1)) {
      os<<"Wrong format: expected 1 on the diagonal, got "<<unitgot<<".\n";
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
      os<<"The portion of the UpperTriMatrix which was successfully read is: \n";
      ConstUpperTriMatrixView<T> mm = m->View();
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
  
  template <class T> void LowerTriMatrixReadError<T>::Write(
      ostream& os) const throw()
  {
    os<<"TMV Read Error: Reading istream input for LowerTriMatrix\n";
    if (exp != got) {
      os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
    }
    if (unitgot != T(1)) {
      os<<"Wrong format: expected 1 on the diagonal, got "<<unitgot<<".\n";
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
      os<<"The portion of the LowerTriMatrix which was successfully read is: \n";
      ConstLowerTriMatrixView<T> mm = m->View();
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

  template <class T, IndexStyle I> void UpperTriMatrixView<T,I>::Read(
      istream& is) const
  {
    T* mrowi = ptr();
    const int sj = stepj();
    const int ds = stepi()+sj;
    size_t len = size();
    if (isunit()) {
      mrowi += sj;
      --len;
    }
    char paren;
    for(size_t i=0;i<size();++i,--len,mrowi+=ds) {
      is >> paren;
      if (!is || paren != '(') 
	throw UpperTriMatrixReadError<T>(i,0,*this,is,'(',is?paren:'(');
      if (isunit()) {
	T unit;
	is >> unit;
	if (!is || unit != T(1))
	  throw UpperTriMatrixReadError<T>(i,i,*this,is,is?unit:T(1));
      }
      T* mij = mrowi;
      if (isrm()) 
	for(size_t k=len;k>0;--k,++mij) {
	  is >> *mij;
	  if (!is) 
	    throw UpperTriMatrixReadError<T>(i,size()-k,*this,is);
	}
      else 
	for(size_t k=len;k>0;--k,mij+=sj) {
	  is >> *mij;
	  if (!is) 
	    throw UpperTriMatrixReadError<T>(i,size()-k,*this,is);
	}
      is >> paren;
      if ((!is && i+1<size())  || paren != ')') 
	throw UpperTriMatrixReadError<T>(i,size(),*this,is,')',is?paren:')');
    }
    if (isconj()) ConjugateSelf();
  }

  template <class T> void GenLowerTriMatrix<T>::Write(ostream& os) const
  {
    const T* mrowi = cptr();
    const int sj = stepj();
    size_t len = isunit() ? 0 : 1;
    os << size()<<' '<<size() << endl;
    for(size_t i=0;i<size();++i,++len,mrowi+=stepi()) {
      os << "( ";
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,++mij) os << ' '<<CONJ(*mij)<<' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) os << ' '<<*mij<<' ';
      else
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj) os << ' '<<CONJ(*mij)<<' ';
	else 
	  for(size_t k=len;k>0;--k,mij+=sj) os << ' '<<*mij<<' ';
      if (isunit()) os << ' '<<T(1)<<' ';
      for(size_t j=isunit()?len+1:len;j<size();j++) 
	os << ' '<<T(0)<<' ';
      os << " )\n";
    }
  }

  template <class T> void GenLowerTriMatrix<T>::Write(ostream& os,
      RealType(T) thresh) const
  {
    const T* mrowi = cptr();
    const int sj = stepj();
    size_t len = isunit() ? 0 : 1;
    os << size()<<' '<<size() << endl;
    for(size_t i=0;i<size();++i,++len,mrowi+=stepi()) {
      os << "( ";
      const T* mij = mrowi;
      if (isrm()) 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : *mij)<<' ';
      else
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
	else 
	  for(size_t k=len;k>0;--k,mij+=sj) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : *mij)<<' ';
      if (isunit()) os << ' '<<T(1)<<' ';
      for(size_t j=isunit()?len+1:len;j<size();j++) 
	os << ' '<<T(0)<<' ';
      os << " )\n";
    }
  }

  template <class T> void GenLowerTriMatrix<T>::WriteCompact(ostream& os) const
  {
    const T* mrowi = cptr();
    const int sj = stepj();
    size_t len = isunit() ? 0 : 1;
    os << "L " << size() << endl;
    for(size_t i=0;i<size();++i,++len,mrowi+=stepi()) {
      os << "( ";
      const T* mij = mrowi;
      if (isrm())
	if (isconj())
	  for(size_t k=len;k>0;--k,++mij) os << ' '<<CONJ(*mij)<<' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) os << ' '<<*mij<<' ';
      else 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj) os << ' '<<CONJ(*mij)<<' ';
	else
	  for(size_t k=len;k>0;--k,mij+=sj) os << ' '<<*mij<<' ';
      if (isunit()) os << ' '<<T(1)<<' ';
      os << " )\n";
    }
  }

  template <class T> void GenLowerTriMatrix<T>::WriteCompact(ostream& os,
      RealType(T) thresh) const
  {
    const T* mrowi = cptr();
    const int sj = stepj();
    size_t len = isunit() ? 0 : 1;
    os << "L " << size() << endl;
    for(size_t i=0;i<size();++i,++len,mrowi+=stepi()) {
      os << "( ";
      const T* mij = mrowi;
      if (isrm())
	if (isconj())
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
	else 
	  for(size_t k=len;k>0;--k,++mij) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : *mij)<<' ';
      else 
	if (isconj()) 
	  for(size_t k=len;k>0;--k,mij+=sj) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : CONJ(*mij))<<' ';
	else
	  for(size_t k=len;k>0;--k,mij+=sj) 
	    os << ' '<<(abs(*mij) < thresh ? T(0) : *mij)<<' ';
      if (isunit()) os << ' '<<T(1)<<' ';
      os << " )\n";
    }
  }

  template <class T, IndexStyle I> void LowerTriMatrixView<T,I>::Read(
      istream& is) const
  {
    T* mrowi = ptr();
    const int sj = stepj();
    size_t len = isunit() ? 0 : 1;
    char paren;
    for(size_t i=0;i<size();++i,++len,mrowi+=stepi()) {
      is >> paren;
      if (!is || paren != '(') 
	throw LowerTriMatrixReadError<T>(i,0,*this,is,'(',is?paren:'(');
      T* mij = mrowi;
      if (isrm()) 
	for(size_t k=len;k>0;--k,++mij) {
	  is >> *mij;
	  if (!is) 
	    throw LowerTriMatrixReadError<T>(i,len-k,*this,is);
	}
      else 
	for(size_t k=len;k>0;--k,mij+=sj) {
	  is >> *mij;
	  if (!is) 
	    throw LowerTriMatrixReadError<T>(i,len-k,*this,is);
	}
      if (isunit()) {
	T unit;
	is >> unit;
	if (!is || unit != T(1))
	  throw LowerTriMatrixReadError<T>(i,i,*this,is,is?unit:T(1));
      }
      is >> paren;
      if ((!is && i+1<size())  || paren != ')') 
	throw LowerTriMatrixReadError<T>(i,size(),*this,is,')',is?paren:')');
    }
    if (isconj()) ConjugateSelf();
  }

  template <class T, DiagType D, StorageType S, IndexStyle I> 
    istream& operator>>(istream& is, auto_ptr<UpperTriMatrix<T,D,S,I> >& m)
    { 
      char ul;
      is >> ul;
      if (!is)
	throw UpperTriMatrixReadError<T>(is);
      if (ul != 'U')
	throw UpperTriMatrixReadError<T>(is,'U',ul);
      size_t size;
      is >> size;
      if (!is) 
	throw UpperTriMatrixReadError<T>(is);
      m.reset(new UpperTriMatrix<T,D,S,I>(size));
      m->View().Read(is); 
      return is;
    }

  template <class T, DiagType D, StorageType S, IndexStyle I>
    istream& operator>>(istream& is, auto_ptr<LowerTriMatrix<T,D,S,I> >& m)
    { 
      char ul;
      is >> ul;
      if (!is)
	throw LowerTriMatrixReadError<T>(is);
      if (ul != 'L')
	throw LowerTriMatrixReadError<T>(is,'L',ul);
      size_t size;
      is >> size;
      if (!is) 
	throw LowerTriMatrixReadError<T>(is);
      m.reset(new LowerTriMatrix<T,D,S,I>(size));
      m->View().Read(is); 
      return is;
    }

  template <class T> istream& operator>>(
      istream& is, const UpperTriMatrixView<T>& m)
  { 
    char ul;
    is >> ul;
    if (!is)
      throw UpperTriMatrixReadError<T>(is);
    if (ul != 'U')
      throw UpperTriMatrixReadError<T>(is,'U',ul);
    size_t s;
    is >> s;
    if (!is)
      throw UpperTriMatrixReadError<T>(is);
    if (m.size() != s)
      throw UpperTriMatrixReadError<T>(m,is,s);
    TMVAssert(m.size() == s);
    m.Read(is);
    return is;
  }

  template <class T> istream& operator>>(
      istream& is, const LowerTriMatrixView<T>& m)
  { 
    char ul;
    is >> ul;
    if (!is)
      throw LowerTriMatrixReadError<T>(is);
    if (ul != 'L')
      throw LowerTriMatrixReadError<T>(is,'L',ul);
    size_t s;
    is >> s;
    if (!is)
      throw LowerTriMatrixReadError<T>(is);
    if (m.size() != s)
      throw LowerTriMatrixReadError<T>(m,is,s);
    TMVAssert(m.size() == s);
    m.Read(is);
    return is;
  }


#define InstFile "TMV_TriMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


