
#include "TMV.h"

namespace tmv {

  Permutation& Permutation::InvertStorage()
  {
    if (!isident) {
      TMVAssert(IsValid());
      MakeCycles();
      for(size_t i=0;i<cycles.size();++i) {
	vector<size_t> ci = cycles[i];
	TMVAssert(ci.size() > 1);
	size_t last = ci.back();
	for(size_t j=1;j<ci.size();++j) 
	  itsp[ci[j]] = ci[j-1];
	itsp[ci.front()] = last;
      }
      isinv = !isinv;
      changed = true;
      TMVAssert(IsValid());
    }
    return *this;
  }

  bool Permutation::IsValid() const 
  {
    vector<bool> used(size(),false);
    for(size_t i=0;i<size();++i) used[itsp[i]] = true;
    for(size_t i=0;i<size();++i) if (!used[i]) {
      cout<<*this<<endl;
      cout<<"isident = "<<isident<<endl;
      cout<<"isinv = "<<isinv<<endl;
      cout<<i<<" is not used.\n";
      return false;
    }
    if (isident) for(size_t i=0;i<size();++i) if(itsp[i]!=i) {
      cout<<*this<<endl;
      cout<<"isident = "<<isident<<endl;
      cout<<"isinv = "<<isinv<<endl;
      cout<<"itsp["<<i<<"] = "<<itsp[i]<<", but isident, so should = "<<i<<endl;
      return false;
    }
    if (negdet != DetermineNegDet()) {
      cout<<*this<<endl;
      cout<<"isident = "<<isident<<endl;
      cout<<"isinv = "<<isinv<<endl;
      cout<<"negdet = "<<negdet<<endl;
      cout<<"actual negdet = "<<DetermineNegDet()<<endl;
      return false;
    }
    return true;
  }

  void Permutation::MakeCycles() const
  {
    if (changed) {
      cycles.clear();
      if (!isident) {
	vector<bool> done(size(),false);
	for(size_t i=0; i<size(); ++i) if (!done[i]) {
	  if (itsp[i] != i) {
	    cycles.push_back(vector<size_t>());
	    vector<size_t>& cycle = cycles.back();
	    cycle.push_back(i);
	    size_t k = itsp[i];
	    while (k != i) {
	      done[k] = true;
	      cycle.push_back(k);
	      k = itsp[k];
	    }
	  }
	}
      }
      changed = false;
    }
  }

  int Permutation::DetermineNegDet() const
  {
    MakeCycles();
    bool isneg = false;
    for(size_t i=0; i<cycles.size(); ++i) 
      if (cycles[i].size() %2 == 0) isneg = !isneg;
    return isneg;
  }

  size_t gcd(size_t a, size_t b)
  {
    if (b > a) swap(a,b);
    while (b > 0) {
      size_t c = a%b;
      a = b;
      b = c;
    }
    return a;
  }

  Permutation& Permutation::operator^=(int n)
    // p = p^n
  {
    TMVAssert(IsValid());
    if (n == 0) SetToIdentity();
    else if (n == 1) {}
    else if (n == -1) InvertSelf();
    else if (!isident) {
      MakeCycles();
      for(size_t i=0;i<cycles.size();++i) {
        vector<size_t> ci = cycles[i];
	TMVAssert(ci.size() > 1);
	int nn = n % ci.size();
	if (nn == 0) {
	  for(size_t j = 0;j<ci.size();++j)
	    itsp[ci[j]] = ci[j];
	} else if (nn > 1) {
	  size_t g = gcd(ci.size(),nn);
	  for(size_t jstart = 0;jstart<g;++jstart) {
	    size_t j = jstart;
	    size_t k = jstart + nn;
	    TMVAssert(k < ci.size());
	    size_t first = ci[jstart];
	    while (k != jstart) {
	      itsp[ci[j]] = ci[k];
	      size_t k2 = (k + nn)%ci.size();
	      j = k;
	      k = k2;
	    }
	    itsp[ci[j]] = first;
	  }
	}
      }
      if (n%2 == 0) negdet = false;
      changed = true;
      // else det stays the same
      TMVAssert(IsValid());
    }
    return *this;
  }

  void Permutation::SetToPP(const Permutation& p1, const Permutation& p2) 
    // p = p1 * p2
    // This assumes all Permutations are row-based.
  {
    TMVAssert(p1.IsValid());
    TMVAssert(p2.IsValid());
    isident = p1.isident && p2.isident;
    if (isident) SetToIdentity();
    else {
      isinv = false;
      TMVAssert(size() == p1.size());
      for(size_t i=0;i<size();++i) itsp[i] = p2[p1[i]];
      negdet = (p1.negdet != p2.negdet);
      isident = false;
      changed = true;
    }
    TMVAssert(IsValid());
  }

  void Permutation::RMultEq(Permutation& p1) const
    // p1 = p1 * p
    // This assumes all Permutations are row-based.
  {
    TMVAssert(IsValid());
    TMVAssert(p1.IsValid());
    if (!isident) {
      if (this == &p1) {
	p1 ^= 2;
      } else {
	for(size_t i=0;i<size();++i) p1.itsp[i] = itsp[p1[i]];
	if (negdet) p1.negdet = !p1.negdet;
	p1.isident = false;
	p1.changed = true;
      }
      TMVAssert(p1.IsValid());
    }
  }

  void Permutation::LMultEq(Permutation& p2) const
    // p2 = p * p2
    // This assumes all Permutations are row-based.
  {
    TMVAssert(IsValid());
    TMVAssert(p2.IsValid());
    if (!isident) {
      if (this == &p2) {
	p2 ^= 2;
      } else {
	MakeCycles();
	for(size_t i=0;i<cycles.size();++i) {
	  vector<size_t> ci = cycles[i];
	  TMVAssert(ci.size() > 1);
	  size_t first = p2[ci.front()];
	  for(size_t j=1;j<ci.size();++j) 
	    p2.itsp[ci[j-1]] = p2[ci[j]];
	  p2.itsp[ci.back()] = first;
	}
	if (negdet) p2.negdet = !p2.negdet;
	p2.isident = false;
	p2.changed = true;
      }
      TMVAssert(p2.IsValid());
    }
  }

  void Permutation::SetToPtP(const Permutation& p1, const Permutation& p2) 
    // p = p1t * p2
  {
    TMVAssert(p1.IsValid());
    TMVAssert(p2.IsValid());
    isinv = false;
    TMVAssert(size() == p1.size());
    isident = p1.isident && p2.isident;
    if (&p1 == &p2) isident = true;
    if (isident) SetToIdentity();
    else {
      for(size_t i=0;i<p1.size();++i) itsp[p1[i]] = p2[i];
      negdet = (p1.negdet != p2.negdet);
      isident = false;
      changed = true;
    }
    TMVAssert(IsValid());
  }

  void Permutation::LDivEq(Permutation& p2) const
    // p2 = pt * p2
  {
    TMVAssert(IsValid());
    TMVAssert(p2.IsValid());
    TMVAssert(size() == p2.size());
    if (this == &p2) { 
      p2.SetToIdentity();
    } else if (!isident) {
      MakeCycles();
      for(size_t i=0;i<cycles.size();++i) {
	const vector<size_t>& ci = cycles[i];
	TMVAssert(ci.size() > 1);
	size_t last = p2[ci.back()];
	for(size_t j=ci.size()-1;j>0;j--) 
	  p2.itsp[ci[j]] = p2[ci[j-1]];
	p2.itsp[ci.front()] = last;
      }
      if (negdet) p2.negdet = !p2.negdet;
      p2.isident = false;
      p2.changed = true;
      TMVAssert(p2.IsValid());
    }
  }

  void Permutation::SetToPPt(const Permutation& p1, const Permutation& p2) 
    // p = p1 * p2t
  { 
    TMVAssert(p1.IsValid());
    TMVAssert(p2.IsValid());
    isinv = false;
    TMVAssert(size() == p1.size());
    isident = p1.isident && p2.isident;
    if (&p1 == &p2) isident = true;
    if (isident) SetToIdentity();
    else {
      Permutation p2t = p2.Inverse();
      SetToPP(p1,p2t);
    }
    TMVAssert(IsValid());
  }

  void Permutation::RDivEq(Permutation& p2) const
    // p2 = p2 * pt
  { 
    TMVAssert(IsValid());
    if (!isident) {
      Permutation pt = Inverse();
      pt.RMultEq(p2);
    }
  }

  template <class T> void Permutation::LMultEq(const VectorView<T>& v) const
    // v = p * v
  { 
    TMVAssert(IsValid());
    TMVAssert(size() == v.size());
    if (!isident) {
      MakeCycles();
      for(size_t i=0;i<cycles.size();++i) {
	vector<size_t> ci = cycles[i];
	TMVAssert(ci.size() > 1);
	T first = v[ci.front()];
	for(size_t j=1;j<ci.size();++j) 
	  v[ci[j-1]] = v[ci[j]];
	v[ci.back()] = first;
      }
    }
  }

  template <class T> void Permutation::LDivEq(const VectorView<T>& v) const
    // v = pt * v
  { 
    TMVAssert(IsValid());
    TMVAssert(size() == v.size());
    if (!isident) {
      MakeCycles();
      for(size_t i=0;i<cycles.size();++i) {
	vector<size_t> ci = cycles[i];
	TMVAssert(ci.size() > 1);
	T last = v[ci.back()];
	for(size_t j=ci.size()-1;j>0;j--) 
	  v[ci[j]] = v[ci[j-1]];
	v[ci.front()] = last;
      }
    }
  }

  template <class T> void Permutation::LMult(
      const GenVector<T>& v1, const VectorView<T>& v0) const
    // v0 = p * v1
  {
    TMVAssert(IsValid());
    TMVAssert(size() == v1.size());
    TMVAssert(v0.size() == v1.size());
    if (isident) v0 = v1;
    else for(size_t i=0;i<size();++i) v0(i) = v1(itsp[i]); 
  }

  template <class T> void Permutation::LDiv(
      const GenVector<T>& v1, const VectorView<T>& v0) const
    // v0 = pt * v1
  {
    TMVAssert(IsValid());
    TMVAssert(size() == v1.size());
    TMVAssert(v0.size() == v1.size());
    if (isident) v0 = v1;
    else for(size_t i=0;i<size();++i) v0(itsp[i]) = v1(i); 
  }

  template <class T> void Permutation::LMultEq(const MatrixView<T>& m) const
    // m = p * m
  { 
    TMVAssert(IsValid());
    TMVAssert(size() == m.colsize());
    if (!isident) {
      MakeCycles();
      for(size_t i=0;i<cycles.size();++i) {
	vector<size_t> ci = cycles[i];
	TMVAssert(ci.size() > 1);
	Vector<T> first = m.row(ci.front());
	for(size_t j=1;j<ci.size();++j) 
	  m.row(ci[j-1]) = m.row(ci[j]);
	m.row(ci.back()) = first;
      }
    }
  }

  template <class T> void Permutation::LDivEq(const MatrixView<T>& m) const
    // m = pt * m
  { 
    TMVAssert(IsValid());
    TMVAssert(size() == m.colsize());
    if (!isident) {
      MakeCycles();
      for(size_t i=0;i<cycles.size();++i) {
	vector<size_t> ci = cycles[i];
	TMVAssert(ci.size() > 1);
	Vector<T> last = m.row(ci.back());
	for(size_t j=ci.size()-1;j>0;j--) 
	  m.row(ci[j]) = m.row(ci[j-1]);
	m.row(ci.front()) = last;
      }
    }
  }

  template <class T> void Permutation::LMult(
      const GenMatrix<T>& m1, const MatrixView<T>& m0) const
    // m0 = p * m1
  {
    TMVAssert(IsValid());
    TMVAssert(size() == m1.colsize());
    TMVAssert(m0.colsize() == m1.colsize());
    TMVAssert(m0.rowsize() == m1.rowsize());
    if (isident) m0 = m1;
    else for(size_t i=0;i<size();++i) m0.row(i) = m1.row(itsp[i]); 
  }

  template <class T> void Permutation::LDiv(
      const GenMatrix<T>& m1, const MatrixView<T>& m0) const
    // m0 = pt * m1
  {
    TMVAssert(IsValid());
    TMVAssert(size() == m1.colsize());
    TMVAssert(m0.colsize() == m1.colsize());
    TMVAssert(m0.rowsize() == m1.rowsize());
    if (isident) m0 = m1;
    else for(size_t i=0;i<size();++i) m0.row(itsp[i]) = m1.row(i); 
  }

  void Permutation::Write(ostream& fout) const
  {
    fout << size() << " " << isinv << " (";
    for(size_t i=0;i<size();++i) {
      fout << " "<< itsp[i] << " ";
    }
    fout << " )";
  }

  void Permutation::WriteCycles(ostream& os) const
  {
    MakeCycles();
    os << "[ ";
    for(size_t i=0; i<cycles.size(); ++i) {
      vector<size_t> ci = cycles[i];
      os << "(";
      for(size_t j=0; j<ci.size(); ++j) {
	if (j>0) os <<",";
	os << ci[j];
      }
      os << ") ";
    }
    os <<"]";
  }

  void Permutation::WriteMatrix(ostream& os) const
  { 
    os << size() <<"  "<<size()<<endl;
    for(size_t i=0; i<size(); ++i) {
      os << "( ";
      for(size_t j=0; j<size(); ++j) os <<" "<<(*this)(i,j)<<" ";
      os << " )\n";
    }
  }

  void Permutation::Read(istream& fin)
  {
    fin >> isinv;
    if (!fin) tmv_error("reading isinv in Permutation::Read");
    char paren;
    fin >> paren;
    if (!fin || paren != '(') tmv_error("reading ( in Permutation::Read");
    for(size_t i=0;i<size();++i) {
      fin >> itsp[i];
      if (!fin) tmv_error("reading value in Permutation::Read");
    }
    fin >> paren;
    if (!fin || paren != ')') tmv_error("reading ) in Permutation::Read");
    if (!IsValid()) tmv_error("Permutation read is invalid");
  }

  bool operator==(const Permutation& p1, const Permutation& p2) 
  {
    if (p1.size() != p2.size()) return false;
    for(size_t i=0;i<p1.size();++i) if (p1[i]!=p2[i]) return false;
    return true; 
  }

  template <class T> bool operator==(
      const GenMatrix<T>& m, const Permutation& p) 
  {
    if (m.colsize() != p.size()) return false;
    if (m.rowsize() != p.size()) return false;
    for(size_t i=0;i<m.colsize();++i) for(size_t j=0;j<m.rowsize();++j) 
      if (p(i,j)==0) {
	if (REAL(m(i,j)) != RealType(T)(0)) return false;
	if (IMAG(m(i,j)) != RealType(T)(0)) return false;
      } else {
	if (REAL(m(i,j)) != RealType(T)(1)) return false;
	if (IMAG(m(i,j)) != RealType(T)(0)) return false;
      }
    return true; 
  }

  istream& operator>>(istream& fin, Permutation& p)
  { 
    size_t n;
    fin >> n;
    if (!fin) tmv_error("reading size in Permutation::Read");
    TMVAssert(p.size() == n);
    p.Read(fin);
    return fin;
  }

  istream& operator>>(istream& fin, Permutation* p)
  { 
    size_t n;
    fin >> n;
    if (!fin) tmv_error("reading size in Permutation::Read");
    p = new Permutation(n);
    p->Read(fin);
    return fin;
  }

  void MultPP(const Permutation& p1, const bool inv1, const Permutation& p2,
      const bool inv2, Permutation& p0)
  {
    TMVAssert(p0.size() == p1.size());
    TMVAssert(p0.size() == p2.size());

    if (&p0 == &p1) {
      if (inv1) p0.InvertSelf();
      if (!p0.IsInverse()) {
	if ((p2.IsInverse() == inv2)) p2.RMultEq(p0);
	// p1 = p1 * p2
	else p2.RDivEq(p0);
	// p1 = p1 * p2t
      } else {
	if ((p2.IsInverse() == inv2)) p2.LDivEq(p0);
	// p1t = p1t * p2
	// ==> p1 = p2t * p1
	else p2.LMultEq(p0);
	// p1t = p1t * p2t
	// ==> p1 = p2 * p1
      }
    } else if (&p0 == &p2) {
      if (inv2) p0.InvertSelf();
      if (!p0.IsInverse()) {
	if ((p1.IsInverse() == inv1)) p1.LMultEq(p0);
	// p2 = p1 * p2
	else p1.LDivEq(p0);
	// p2 = p1t * p2
      } else {
	if ((p1.IsInverse() == inv1)) p1.RDivEq(p0);
	// p2t = p1 * p2t
	// ==> p2 = p2 * p1t
	else p1.RMultEq(p0);
	// p2t = p1t * p2t
	// ==> p2 = p2 * p1
      }
    } else {
      if (p1.IsInverse() == inv1) {
	if (p2.IsInverse() == inv2) p0.SetToPP(p1,p2);
	// p = p1 * p2
	else p0.SetToPPt(p1,p2);
	// p = p1 * p2t
      } else {
	if (p2.IsInverse() == inv2) p0.SetToPtP(p1,p2);
	// p = p1t * p2
	else {
	  // p = p1t * p2t
	  // ==> pt = p2 * p1
	  p0.SetToPP(p2,p1);
	  p0.InvertSelf();
	}
      }
    }
  }

  template <class T> void MultPV(const Permutation& p1, const bool inv1,
      const GenVector<T>& v2, const VectorView<T>& v0) 
  {
    TMVAssert(v0.size() == v2.size()); 
    if (v2.SameAs(v0)) {
      if (inv1 == p1.IsInverse()) p1.LMultEq(v0);
      // v0 = p * v0
      else p1.LDivEq(v0);
      // v0 = pt * v0
    } else if (v2.SameStorageAs(v0)) {
      Vector<T> temp(v2);
      if (inv1 == p1.IsInverse()) p1.LMult(temp,v0);
      // v0 = p * v2
      else p1.LDiv(temp,v0);
      // v0 = pt * v2
    } else {
      if (inv1 == p1.IsInverse()) p1.LMult(v2,v0);
      // v0 = p * v2
      else p1.LDiv(v2,v0);
      // v0 = pt * v2
    }
  };

  //
  // Permutations of Matrices
  //

  template <class T> void MultPM(const Permutation& p1, const bool inv1,
      const GenMatrix<T>& m2, const MatrixView<T>& m0)
  {
    TMVAssert(m0.colsize() == m2.colsize()); 
    TMVAssert(m0.rowsize() == m2.rowsize()); 
    if (m2.SameAs(m0)) {
      if (inv1 == p1.IsInverse()) p1.LMultEq(m0);
      // m0 = p * m0
      else p1.LDivEq(m0);
      // m0 = pt * m0
    } else if (m2.SameStorageAs(m0)) {
      if (m0.isrm()) {
	Matrix<T,RowMajor> temp(m2);
	if (inv1 == p1.IsInverse()) p1.LMult(temp.QuickView(),m0);
	else p1.LDiv(temp.QuickView(),m0);
      } else {
	Matrix<T,ColMajor> temp(m2);
	if (inv1 == p1.IsInverse()) p1.LMult(temp.QuickView(),m0);
	else p1.LDiv(temp.QuickView(),m0);
      }
    } else {
      if (inv1 == p1.IsInverse()) p1.LMult(m2,m0);
      else p1.LDiv(m2,m0);
    }
  };

#define InstFile "TMV_Permutation.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace mv


