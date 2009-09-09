
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

  template <class T> RefType(T) UpperTriMatrixView<T>::ref(
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
    switch (this->itsdt) {
      case tmv::LU : case tmv::QR :
	this->itsdiv = new UpperTriLUDiv<T>(*this);
	break;
      case tmv::SV :
	this->itsdiv = new UpperTriSVDiv<T>(*this);
	break;
      default : TMVAssert(false);
    }
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

  template <class T> RefType(T) LowerTriMatrixView<T>::ref(
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
    switch (this->itsdt) {
      case tmv::LU : case tmv::QR :
	this->itsdiv = new LowerTriLUDiv<T>(*this);
	break;
      case tmv::SV :
	this->itsdiv = new LowerTriSVDiv<T>(*this);
	break;
      default : TMVAssert(false);
    }
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
    if (i1 < 0 || i1 >= int(colsize())) {
      ok = false;
      cout<<"first col element ("<<i1<<") must be in 0 -- "<<size()<<endl;
    }
    if (i2-istep < 0 || i2-istep >= int(colsize())) {
      ok = false;
      cout<<"last col element ("<<i2-istep<<") must be in 0 -- "<<size()<<endl;
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
      cout<<"first row element ("<<j1<<") must be in 0 -- "<<size()<<endl;
    }
    if (j2-jstep < 0 || j2-jstep >= int(rowsize())) {
      ok = false;
      cout<<"last row element ("<<j2-jstep<<") must be in 0 -- "<<size()<<endl;
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
      cout<<"Upper left corner ("<<i1<<','<<j1<<") must be in Triangle\n";
    }
    if (!okij(i1,j2-jstep)) {
      ok = false;
      cout<<"Upper right corner ("<<i1<<','<<j2-jstep;
      cout<<") must be in Triangle\n";
    }
    if (!okij(i2-istep,j1)) {
      ok = false;
      cout<<"Lower left corner ("<<i2-istep<<','<<j1;
      cout<<") must be in  Triangle\n";
    }
    if (!okij(i2-istep,j2-jstep)) {
      ok = false;
      cout<<"Lower right corner ("<<i2-istep<<','<<j2-jstep;
      cout<<") must be in  Triangle\n";
    }
    return ok;
  }

  template <class T> bool GenUpperTriMatrix<T>::OKSubVector(
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
    if (!okij(i,j)) {
      ok = false;
      cout<<"First element ("<<i<<','<<j<<") must be in Triangle\n";
    }
    if (!okij(i2,j2)) {
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
      cout<<"first diag element ("<<i1<<") must be in 0 -- "<<size()<<endl;
    }
    if (i2-istep<0 || i2-istep >= int(size())) {
      ok = false;
      cout<<"last diag element ("<<i2-istep<<") must be in 0 -- "<<size()<<endl;
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

  //
  // Norms
  //

  template <class T> RealType(T) GenUpperTriMatrix<T>::NormSq() const
  {
    RealType(T) sum(0);
    if (isrm()) 
      if (isunit())
	for(size_t i=0;i<colsize();++i) 
	  sum += row(i,i+1,rowsize()).NormSq();
      else
	for(size_t i=0;i<colsize();++i) 
	  sum += row(i,i,rowsize()).NormSq();
    else
      if (isunit())
	for(size_t j=0;j<rowsize();++j) 
	  sum += col(j,0,j).NormSq();
      else
	for(size_t j=0;j<rowsize();++j) 
	  sum += col(j,0,j+1).NormSq();
    if (isunit()) sum += size();
    return sum;
  }

  template <class T> RealType(T) GenUpperTriMatrix<T>::MaxAbsElement() const
  {
    RealType(T) max(0);
    if (isrm())
      for(size_t i=0;i<colsize();++i) {
	RealType(T) temp;
	if (isunit())
	  temp = row(i,i+1,rowsize()).NormInf();
	else 
	  temp = row(i,i,rowsize()).NormInf();
	if (temp > max) max = temp;
      }
    else
      for(size_t j=0;j<rowsize();++j) {
	RealType(T) temp;
	if (isunit())
	  temp = col(j,0,j).NormInf();
	else 
	  temp = col(j,0,j+1).NormInf();
	if (temp > max) max = temp;
      }
    if (isunit() && max < RealType(T)(1)) max = RealType(T)(1);
    return max;
  }

  template <class T> RealType(T) GenUpperTriMatrix<T>::Norm1() const
  { 
    RealType(T) max(0);
    for(size_t j=0;j<rowsize();++j) {
      RealType(T) temp;
      if (isunit()) {
	temp = col(j,0,j).Norm1();
	if (isunit()) temp += RealType(T)(1);
      } else temp = col(j,0,j+1).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  } 

  template <class T> RealType(T) GenUpperTriMatrix<T>::NormInf() const
  { 
    RealType(T) max(0);
    for(size_t j=0;j<rowsize();++j) {
      RealType(T) temp;
      if (isunit()) {
	temp = row(j,j+1,rowsize()).Norm1();
	if (isunit()) temp += RealType(T)(1);
      } else temp = row(j,j,rowsize()).Norm1();
      if (temp > max) max = temp;
    }
    return max;
  }

  //
  // Modifying Functions
  //

  template <class T> 
    const UpperTriMatrixView<T>& UpperTriMatrixView<T>::SetAllTo(T x) const
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

  template <class T> 
    const UpperTriMatrixView<T>& UpperTriMatrixView<T>::ConjugateSelf() const
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

  template <class T> 
    const UpperTriMatrixView<T>& UpperTriMatrixView<T>::SetToIdentity(
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

  template <class T> bool operator==(
      const GenUpperTriMatrix<T>& m1, const GenUpperTriMatrix<T>& m2)
  {
    if (m1.size() != m2.size()) return false;
    if (m1.SameAs(m2)) return true;
    for(size_t j=0;j<m1.size();++j) 
      if (m1.col(j,0,j) != m2.col(j,0,j)) return false;
    if (m1.isunit() && !m2.isunit())
      for(size_t i=0;i<m1.size();++i) if (m2(i,i) != T(1)) return false;
    else if (m2.isunit() && !m1.isunit())
      for(size_t i=0;i<m1.size();++i) if (m1(i,i) != T(1)) return false;
    else if (!m1.isunit() && !m2.isunit()) 
      if (m1.diag() != m2.diag()) return false;

    return true;  
  }

  //
  // I/O
  //

  // MJ: Add Aptr stuff
  template <class T> void GenUpperTriMatrix<T>::WriteCompact(
      ostream& os) const
  {
    os << "UT " << size() << endl;
    for(size_t i=0;i<colsize();++i) {
      os << "( ";
      size_t j1 = (isunit() ? i+1 : i);
      size_t j2 = rowsize();
      if (isunit()) os << ' '<<T(1)<<' ';
      if (isrm()) {
	if (isconj()) {
	  const CVIt<T,Unit,Conj> _end = row(i,j1,j2).end();
	  CVIt<T,Unit,Conj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	} else {
	  const CVIt<T,Unit,NonConj> _end = row(i,j1,j2).end();
	  CVIt<T,Unit,NonConj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	}
      } else {
	if (isconj()) {
	  const CVIt<T,Step,Conj> _end = row(i,j1,j2).end();
	  CVIt<T,Step,Conj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	} else {
	  const CVIt<T,Step,NonConj> _end = row(i,j1,j2).end();
	  CVIt<T,Step,NonConj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	}
      }
      os << " )\n";
    }
  }

  template <class T> void GenLowerTriMatrix<T>::WriteCompact(
      ostream& os) const
  {
    os << "LT " << size() << endl;
    for(size_t i=0;i<colsize();++i) {
      os << "( ";
      size_t j2 = (isunit() ? i : i+1);
      if (isrm()) {
	if (isconj()) {
	  const CVIt<T,Unit,Conj> _end = row(i,0,j2).end();
	  CVIt<T,Unit,Conj> it = row(i,0,j2).begin();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	} else {
	  const CVIt<T,Unit,NonConj> _end = row(i,0,j2).end();
	  CVIt<T,Unit,NonConj> it = row(i,0,j2).begin();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	}
      } else {
	if (isconj()) {
	  const CVIt<T,Step,Conj> _end = row(i,0,j2).end();
	  CVIt<T,Step,Conj> it = row(i,0,j2).begin();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	} else {
	  const CVIt<T,Step,NonConj> _end = row(i,0,j2).end();
	  CVIt<T,Step,NonConj> it = row(i,0,j2).begin();
	  for(; it!=_end; ++it) os << ' '<<*it<<' ';
	}
      }
      if (isunit()) os << ' '<<T(1)<<' ';
      os << " )\n";
    }
  }

  template <class T> inline bool GetUnit(istream& is, DiagType dt)
  {
    TMVAssert(dt == UnitDiag);
    T unit;
    is >> unit;
    return (unit == T(1));
  }

  template <class T> void UpperTriMatrixView<T>::Read(istream& is) const
  {
    char paren;
    for(size_t i=0;i<size();++i) {
      is >> paren;
      if (!is || paren != '(') tmv_error("reading ( in TriMatrix::Read");
      size_t j1 = (isunit() ? i+1 : i);
      if (isunit())
	if (!GetUnit<T>(is,dt())) 
	  tmv_error("reading diagonal 1/0 in TriMatrix::Read");
      if (isrm()) {
	if (isconj()) {
	  const VIt<T,Unit,Conj> _end = row(i,j1,size()).end();
	  VIt<T,Unit,Conj> it = row(i,j1,size()).begin();
	  for(; it!=_end; ++it) is >> *it;
	} else {
	  const VIt<T,Unit,NonConj> _end = row(i,j1,size()).end();
	  VIt<T,Unit,NonConj> it = row(i,j1,size()).begin();
	  for(; it!=_end; ++it) is >> *it;
	}
      } else {
	if (isconj()) {
	  const VIt<T,Step,Conj> _end = row(i,j1,size()).end();
	  VIt<T,Step,Conj> it = row(i,j1,size()).begin();
	  for(; it!=_end; ++it) is >> *it;
	} else {
	  const VIt<T,Step,NonConj> _end = row(i,j1,size()).end();
	  VIt<T,Step,NonConj> it = row(i,j1,size()).begin();
	  for(; it!=_end; ++it) is >> *it;
	}
      }
      if (!is) tmv_error("reading values in TriMatrix::Read");
      is >> paren;
      if ((!is && i+1<size())  || paren != ')') 
	tmv_error("reading ) in TriMatrix::Read");
    }
  }

  template <class T> void LowerTriMatrixView<T>::Read(istream& is) const
  {
    char paren;
    for(size_t i=0;i<size();++i) {
      is >> paren;
      if (!is || paren != '(') tmv_error("reading ( in TriMatrix::Read");
      size_t j2 = (isunit() ? i : i+1);
      if (isrm()) {
	if (isconj()) {
	  const VIt<T,Unit,Conj> _end = row(i,0,j2).end();
	  VIt<T,Unit,Conj> it = row(i,0,j2).begin();
	  for(; it!=_end; ++it) is >> *it;
	} else {
	  const VIt<T,Unit,NonConj> _end = row(i,0,j2).end();
	  VIt<T,Unit,NonConj> it = row(i,0,j2).begin();
	  for(; it!=_end; ++it) is >> *it;
	}
      } else {
	if (isconj()) {
	  const VIt<T,Step,Conj> _end = row(i,0,j2).end();
	  VIt<T,Step,Conj> it = row(i,0,j2).begin();
	  for(; it!=_end; ++it) is >> *it;
	} else {
	  const VIt<T,Step,NonConj> _end = row(i,0,j2).end();
	  VIt<T,Step,NonConj> it = row(i,0,j2).begin();
	  for(; it!=_end; ++it) is >> *it;
	}
      }
      if (isunit())
	if (!GetUnit<T>(is,dt())) 
	  tmv_error("reading diagonal 1/0 in TriMatrix::Read");
      if (!is) tmv_error("reading values in TriMatrix::Read");
      is >> paren;
      if ((!is && i+1<size())  || paren != ')') 
	tmv_error("reading ) in TriMatrix::Read");
    }
  }

  template <class T, DiagType D, StorageType S> istream& operator>>(
      istream& is, UpperTriMatrix<T,D,S>* m)
  { 
    size_t size;
    is >> size;
    if (!is) tmv_error("reading size in TriMatrix::Read");
    m = new UpperTriMatrix<T,D,S>(size);
    m->QuickView().Read(is); 
    return is;
  }

  template <class T, DiagType D, StorageType S> istream& operator>>(
      istream& is, LowerTriMatrix<T,D,S>* m)
  { 
    size_t size;
    is >> size;
    if (!is) tmv_error("reading size in TriMatrix::Read");
    m = new LowerTriMatrix<T,D,S>(size);
    m->QuickView().Read(is); 
    return is;
  }

  template <class T> istream& operator>>(
      istream& is, const UpperTriMatrixView<T>& m)
  { 
    size_t size;
    is >> size;
    if (!is) tmv_error("reading size in TriMatrix::Read");
    TMVAssert(m.size() == size);
    m.Read(is);
    return is;
  }

  template <class T> istream& operator>>(
      istream& is, const LowerTriMatrixView<T>& m)
  { 
    size_t size;
    is >> size;
    if (!is) tmv_error("reading size in TriMatrix::Read");
    TMVAssert(m.size() == size);
    m.Read(is);
    return is;
  }


#define InstFile "TMV_TriMatrix.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


