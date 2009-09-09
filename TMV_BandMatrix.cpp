
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

  template <class T> RefType(T) BandMatrixView<T>::ref(size_t i, size_t j) const
  {
    TMVAssert(this->okij(i,j));
    T* mi;
    if (this->isrm()) mi = ptr() + int(i)*stepi() + j;
    else if (this->iscm()) mi = ptr() + i + int(j)*stepj();
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
      case tmv::LU : this->itsdiv = new BandLUDiv<T>(*this); break;
      case tmv::QR : this->itsdiv = new BandQRDiv<T>(*this); break;
      case tmv::SV : this->itsdiv = new BandSVDiv<T>(*this); break;
      default : TMVAssert(false);
    }
  }

  template <StorageType S> size_t BandStorageLength(
      size_t cs, size_t rs, int lo, int hi)
  { TMVAssert(false); }

  template <> size_t BandStorageLength<RowMajor>(
      size_t cs, size_t rs, int lo, int hi)
  {
    if (cs == rs) return (cs-1)*(lo+hi+1)+1;
    else if (cs > rs) return min(rs+lo,cs)*(lo+hi)+rs;
    else return (cs-1)*(lo+hi)+min(cs+hi,rs);
  }

  template <> size_t BandStorageLength<ColMajor>(
      size_t cs, size_t rs, int lo, int hi)
  {
    if (cs == rs) return (cs-1)*(lo+hi+1)+1;
    else if (cs > rs) return (rs-1)*(lo+hi)+min(rs+lo,cs);
    else return min(cs+hi,rs)*(lo+hi)+cs;
  }

  template <> size_t BandStorageLength<DiagMajor>(
      size_t cs, size_t rs, int lo, int hi)
  {
    if (cs == rs) return (cs-1)*(lo+hi+1)+1;
    else if (cs > rs) return rs*(lo+hi+1);
    else return min(cs+hi,rs)*(lo+hi)+cs;
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
    if (i1 < 0 || i1 >= int(this->colsize())) {
      ok = false;
      cout<<"first col element ("<<i1<<") must be in 0 -- "<<this->colsize()<<endl;
    }
    if (i2-istep < 0 || i2-istep >= int(this->colsize())) {
      ok = false;
      cout<<"last col element ("<<i2-istep<<") must be in 0 -- "<<this->colsize()<<endl;
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
    if (j1 < 0 || j1 >= int(this->rowsize())) {
      ok = false;
      cout<<"first row element ("<<j1<<") must be in 0 -- "<<this->rowsize()<<endl;
    }
    if (j2-jstep < 0 || j2-jstep >= int(this->rowsize())) {
      ok = false;
      cout<<"last row element ("<<j2-jstep<<") must be in 0 -- "<<this->rowsize()<<endl;
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
    if (i >= this->colsize()) {
      ok = false;
      cout<<"i ("<<i<<") must be in 0 -- "<<this->colsize()<<endl;
    }
    if (j >= this->rowsize()) {
      ok = false;
      cout<<"j ("<<j<<") must be in 0 -- "<<this->rowsize()<<endl;
    }
    int i2 = int(i)+istep*int(size-1);
    int j2 = int(j)+jstep*int(size-1);
    if (i2 < 0 || i2 >= int(this->colsize())) {
      ok = false;
      cout<<"last element's i ("<<i2<<") must be in 0 -- "<<this->colsize()<<endl;
    }
    if (j2 < 0 || j2 >= int(this->rowsize())) {
      ok = false;
      cout<<"last element's j ("<<j2<<") must be in 0 -- "<<this->rowsize()<<endl;
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
    if (i1 < 0 || i1 >= int(this->colsize())) {
      ok = false;
      cout<<"first col element ("<<i1<<") must be in 0 -- "<<this->colsize()<<endl;
    }
    if (i2-istep < 0 || i2-istep >= int(this->colsize())) {
      ok = false;
      cout<<"last col element ("<<i2-istep<<") must be in 0 -- "<<this->colsize()<<endl;
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
    if (j1 < 0 || j1 >= int(this->rowsize())) {
      ok = false;
      cout<<"first row element ("<<j1<<") must be in 0 -- "<<this->rowsize()<<endl;
    }
    if (j2-jstep < 0 || j2-jstep >= int(this->rowsize())) {
      ok = false;
      cout<<"last row element ("<<j2-jstep<<") must be in 0 -- "<<this->rowsize()<<endl;
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

  template <class T> RealType(T) GenBandMatrix<T>::MaxAbsElement() const
  { 
    RealType(T) max = 0;
    if (isrm()) {
      size_t j1=0;
      size_t j2=nhi()+1;
      size_t k=nlo();
      for(size_t i=0;i<colsize();++i) {
	RealType(T) temp = row(i,j1,j2).MaxAbsElement();
	if (temp > max) max = temp;
	if (k>0) --k; else ++j1;
	if (j2<rowsize()) ++j2;
	else if (j1==rowsize()) break;
      }
    } else if (iscm()) {
      size_t i1=0;
      size_t i2=nlo()+1;
      size_t k=nhi();
      for(size_t j=0;j<rowsize();++j) {
	RealType(T) temp = col(j,i1,i2).MaxAbsElement();
	if (temp > max) max = temp;
	if (k>0) --k; else ++i1;
	if (i2<colsize()) ++i2;
	else if (i1==colsize()) break;
      }
    } else {
      for(int i=-nlo();i<=nhi();++i) {
	RealType(T) temp = diag(i).MaxAbsElement();
	if (temp > max) max = temp;
      }
    }
    return max;
  }

  // 1-Norm = max_j (sum_i |a_ij|)
  template <class T> RealType(T) GenBandMatrix<T>::Norm1() const
  { 
    RealType(T) max = 0;
    size_t i1=0;
    size_t i2=nlo()+1;
    size_t k=nhi();
    for(size_t j=0;j<rowsize();++j) {
      RealType(T) temp = col(j,i1,i2).Norm1();
      if (temp > max) max = temp;
      if (k>0) --k; else ++i1;
      if (i2<colsize()) ++i2;
      else if (i1==colsize()) break;
    }
    return max;
  }

  //
  // Modifying Functions
  //

  template <class T> const BandMatrixView<T>& BandMatrixView<T>::SetAllTo(
      T x) const
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

  template <class T> 
    const BandMatrixView<T>& BandMatrixView<T>::ConjugateSelf() const
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

  template <class T> inline BandMatrix<T,DiagMajor> UpperBiDiagMatrix(
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

  template <class T> inline BandMatrix<T,DiagMajor> LowerBiDiagMatrix(
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

  template <class T> inline BandMatrix<T,DiagMajor> TriDiagMatrix(
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

  template <class T> bool operator==(const GenBandMatrix<T>& m1,
      const GenBandMatrix<T>& m2)
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
 
  // MJ: Add Aptr stuff
  template <class T> void GenBandMatrix<T>::WriteCompact(ostream& os) const
  {
    os << "B "<<this->colsize()<<' '<<this->rowsize()<<' '<<nlo()<<' '<<nhi()<<endl;
    size_t j1=0;
    size_t j2=nhi()+1;
    for(size_t i=0;i<this->colsize();++i) {
      os << "( ";
      if (isrm()) {
	if (isconj()) {
	  const CVIt<T,Unit,Conj> _end = row(i,j1,j2).end();
	  CVIt<T,Unit,Conj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) os << ' ' << *it << ' ';
	} else {
	  const CVIt<T,Unit,NonConj> _end = row(i,j1,j2).end();
	  CVIt<T,Unit,NonConj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) os << ' ' << *it << ' ';
	}
      } else {
	if (isconj()) {
	  const CVIt<T,Step,Conj> _end = row(i,j1,j2).end();
	  CVIt<T,Step,Conj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) os << ' ' << *it << ' ';
	} else {
	  const CVIt<T,Step,NonConj> _end = row(i,j1,j2).end();
	  CVIt<T,Step,NonConj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) os << ' ' << *it << ' ';
	}
      }
      os << " )\n";
      if (int(i)>=nlo() && j1<this->rowsize()) ++j2;
      if (j2<this->rowsize()) ++j2;
    }
  }

  template <class T> void BandMatrixView<T>::Read(istream& is) const
  {
    char paren;
    size_t j1=0;
    size_t j2=nhi()+1;
    for(size_t i=0;i<colsize();++i) {
      is >> paren;
      if (!is || paren != '(') tmv_error("reading ( in BandMatrix::Read");
      if (isrm()) {
	if (isconj()) {
	  const VIt<T,Unit,Conj> _end = row(i,j1,j2).end();
	  VIt<T,Unit,Conj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) is >> *it;
	} else {
	  const VIt<T,Unit,NonConj> _end = row(i,j1,j2).end();
	  VIt<T,Unit,NonConj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) is >> *it;
	}
      } else {
	if (isconj()) {
	  const VIt<T,Step,Conj> _end = row(i,j1,j2).end();
	  VIt<T,Step,Conj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) is >> *it;
	} else {
	  const VIt<T,Step,NonConj> _end = row(i,j1,j2).end();
	  VIt<T,Step,NonConj> it = row(i,j1,j2).begin();
	  for(; it!=_end; ++it) is >> *it;
	}
      }
      if (!is) tmv_error("reading values in BandMatrix::Read");
      is >> paren;
      if (!is || paren != ')') tmv_error("reading ) in BandMatrix::Read");
      if (int(i)>=nlo() && j1<rowsize()) ++j1;
      if (j2<rowsize()) ++j2;
    }
  }

  template <class T, StorageType S> istream& operator>>(istream& is,
      BandMatrix<T,S>* m)
  { 
    size_t cs,rs;
    int lo,hi;
    is >> cs >> rs >> lo >> hi;
    if (!is) tmv_error("reading size in BandMatrix::Read");
    m = new BandMatrix<T,S>(cs,rs,lo,hi);
    m->QuickView().Read(is); 
    return is;
  }

  template <class T> istream& operator>>(istream& is,
      const BandMatrixView<T>& m)
  { 
    size_t cs,rs;
    int lo,hi;
    is >> cs >> rs >> lo >> hi;
    if (!is) tmv_error("reading size in BandMatrix::Read");
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


