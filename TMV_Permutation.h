//---------------------------------------------------------------------------
//
// This file defines the TMV Permutation class.
//
// The Permutation class and all associated functions are contained
// in the namespace tmv, so to use the Permutation class, you need
// to write tmv::Permutation, or put using tmv::Permutation somewhere
// in your code.
//
// A Permutation is an operator that rearranges the elements of 
// a Vector.  We treat the Permutation as a Matrix P whose elements
// are all zeros except for one 1 per row/column.
//
// Note that Transpose(P) x P = I, so P is an orthogonal Matrix.
// Thus, Inverse() and Transpose() are the same thing.
// Also, InvertSelf() and TransposeSelf() both turn the Permutation
// into its inverse.
//
// Constructors:
//
//    Permutation(size_t n)
//        Make an identity permutation of size n
//        This is equivalent to the nxn identity matrix
//
//    Permutation(const vector<size_t>& pp, bool isinv, int det)
//    Permutation(const valarray<size_t>& pp, bool isinv, int det)
//    Permutation(size_t n, const size_t* pp, bool isinv, int det)
//        Make a Permutation P from vector pp.
//        (For the third version, the size of pp is needed)
//        If isinv = false, the unit elements of P are (i,pp[i])
//          This is what might be considered 'normal' form.
//          v' = Pv means v'[i] = v[p[i]]
//        If isinv = true, the unit elements of P are (pp[j],j)
//          This is the inverse of the previous form.
//        If det is given, it is taken to be the determinant of P. 
//        If it is omitted, the determinant will be determined for you.
//
// Special Constructors:
//
//    SortPermutation(Vector<T> v, AD, COMP)
//        Makes a permutations such that P * v is sorted.
//        Technically, this is written as a template, so a vector, deque, etc.
//          are all valid types for v.  The requirements are just that 
//          v.begin() and v.end() are valid iterators.
//        However, only Vector and vector work as given.  For other 
//          "vector-like" types write (for example): SortPermutation<float>(v)
//        AD = ASCEND or DESCEND (ASCEND=default)
//        COMP = REAL_COMP, ABS_COMP, IMAG_COMP, or ARG_COMP (REAL_COMP=default)
//        Note: the resulting permutation, p, is such that:
//          p[0] = the index in v that belongs first
//          p[1] = the index in v that belongs second, etc.
//        
//    ReversePermutation(size_t n)
//        Makes a permutation that reverses the order of elements in 
//        a vector of length n
//
// Access Functions
//
//    size_t size() const
//    size_t colsize() const
//    size_t rowsize() const
//        Return the size of the Permutation.  
//        Since P is a square matrix, all three of these return the same value.
//
//    size_t operator[](size_t i) const
//        If !IsInverse(), returns the j for which P(i,j) = 1
//        If IsInverse(), returns the j for which P(j,i) = 1
//
//    size_t operator()(size_t i, size_t j) const
//        Returns the (i,j) element of P (always either 1 or 0).
//
//    bool IsInverse() const 
//        Returns whether the data storage is the inverse of P
//
//
// Modifying Functions
//
//    Permutation& SetToIdentity()
//        Set to Identity Permutation
//
//    Permutation& InvertSelf()
//    Permutation& TransposeSelf()
//        Invert the Permutation.
//        Note, this does the quick method of switching the nominal
//        storage to the inverse of the current storage method
//        and leaving the p[i] values the same.
//        This has implications for the SwapRows and SwapCols 
//        functions below.
//
//    Permutation& SwapRows(size_t i1, size_t i2)
//    Permutation& SwapCols(size_t i1, size_t i2)
//        Swap p[i] and p[j].  
//        From an identity matrix, you can do either of these commands.
//        However, if the internal storage is of the normal form
//        described above, then it is inefficient to swap to columns.
//        Likewise, if P is storing the inverse permutation,
//        then it is inefficient to swap two rows.
//        Thus, starting from an identity matrix (the usual time
//        you would want to use these functions), you can use either
//        one, but you cannot mix them.
//        If you must mix them, you need to use the function
//        InvertStorage() below.
//
//    void Swap(Permutation& p1, Permutation& p2)
//        Swap the values of two Permutations
//
//    Permutation& InvertStorage()
//        Invert the internal storage of P, but leave the effect of 
//        P as a matrix unchanged.
//
//
// Functions of Permutations:
//        (These are all both member functions and functions of a Permutation,
//         so Det(p) and p.Det() for example are equivalent.)
//
//    Det(p)
//        Returns the determinant of the Matrix version of P.
//
//    Trace(p)
//        Returns the trace of the Matrix version of P.
//        = sum_i ( p_ii )
//
//    Inverse(p)
//    Transpose(p)
//        Returns the inverse Permutation of p 
//        (in the same storage format as p)
//
//
// I/O: 
//
//    os << p 
//        Writes p to ostream os in the Vector format:
//          size inv ( p(0) p(1) p(2) ... p(n-1) )
//
//    is >> p
//        Reads p from istream is in the same format
//
//    WriteCycles(os)
//        Write the cycle notation for p
//
//    WriteMatrix(os)
//        Write p as a Matrix
//
// Operators:
//       Here we use p for a Permutation, m for a Matrix and v for a Vector
//
//    p = p2
//
//    p * p2
//    p / p2  // = p2^-1 * p
//    p % p2  // = p * p2^-1
//
//    p ^= n
//    p ^ n
//
//    v * p
//    p * v
//    v * p
//    v / p // = p^-1 * v
//    v % p // = v * p^-1
//
//    p * m
//    m * p
//    m / p // = p^-1 * m
//    m % p // = m * p^-1
//
//    p == p
//    p != p
//
//       These all act in the way you would expect for the Matrix
//       formulation of p.
//       The operators ^ and ^= are powers.  So p^4 = p*p*p*p
//       Also, we omit operators such as v *= p, since usually you
//       will want v = p*v, but v *= p is normally v = v*p,
//       so we make all x = p*x or x = x*p operations efficiently
//       operate in place without temporary storage.
//       (See Stroustrup 22.4.7 for the basics of this technique.)
//
//


#ifndef TMV_Permutation_H
#define TMV_Permutation_H

#include "TMV_Base.h"
#include "TMV_Matrix.h"
#include <algorithm>

namespace tmv {

  class Permutation;

  class PermutationComposite
  {
    public:
      PermutationComposite() {}
      virtual ~PermutationComposite() {};
      virtual size_t size() const = 0;
      virtual void AssignTo(Permutation&) const = 0;
  };

  class Permutation 
  {

    public:

      //
      // Constructors
      //

      explicit Permutation(size_t n) :
	itsp(new size_t[n]), itsn(n), negdet(false), 
	isident(true), isinv(false), changed(true) 
	{ SetToIdentity(); }

      Permutation(const valarray<size_t>& pp, bool isinv, int dd=0) :
	itsp(new size_t[pp.size()]), itsn(pp.size()), negdet(dd<0), 
	isident(false), isinv(isinv), changed(true)
      { 
	for(size_t i=0;i<pp.size();++i) itsp[i] = pp[i]; 
	if (!dd) negdet = DetermineNegDet();
      }

      Permutation(const valarray<int>& pp, bool isinv, int dd=0) :
	itsp(new size_t[pp.size()]), itsn(pp.size()), negdet(dd<0), 
	isident(false), isinv(isinv), changed(true)
      { 
	for(size_t i=0;i<pp.size();++i) itsp[i] = pp[i]; 
	if (!dd) negdet = DetermineNegDet();
      }

      Permutation(size_t n, const size_t* pp, bool isinv, int dd=0) :
	itsp(new size_t[n]), itsn(n), negdet(dd<0), 
	isident(false), isinv(isinv), changed(true)
      { 
	for(size_t i=0;i<n;++i) itsp[i] = pp[i]; 
	if (!dd) negdet = DetermineNegDet();
      }

      Permutation(size_t n, const int* pp, bool isinv, int dd=0) :
	itsp(new size_t[n]), itsn(n), negdet(dd<0), 
	isident(false), isinv(isinv), changed(true)
      { 
	for(size_t i=0;i<n;++i) itsp[i] = pp[i]; 
	if (!dd) negdet = DetermineNegDet();
      }

      Permutation(const vector<size_t>& pp, bool isinv, int dd=0) :
	itsp(new size_t[pp.size()]), itsn(pp.size()), negdet(dd<0), 
	isident(false), isinv(isinv), changed(true)
      { 
	for(size_t i=0;i<pp.size();++i) itsp[i] = pp[i]; 
	if (!dd) negdet = DetermineNegDet();
      }

      Permutation(const vector<int>& pp, bool isinv, int dd=0) :
	itsp(new size_t[pp.size()]), itsn(pp.size()), negdet(dd<0), 
	isident(false), isinv(isinv), changed(true)
      { 
	for(size_t i=0;i<pp.size();++i) itsp[i] = pp[i]; 
	if (!dd) negdet = DetermineNegDet();
      }

      Permutation(const Permutation& rhs) :
	itsp(new size_t[rhs.size()]), itsn(rhs.size()), negdet(rhs.negdet), 
	isident(rhs.isident), isinv(rhs.isinv), changed(true) 
      { for(size_t i=0;i<rhs.size();++i) itsp[i] = rhs.itsp[i]; }

      Permutation(const PermutationComposite& pcomp) :
	itsp(new size_t[pcomp.size()]), itsn(pcomp.size()), negdet(false), 
	isident(false), isinv(false), changed(true)
      { pcomp.AssignTo(*this); }

      ~Permutation() { TMVAssert(itsp); delete[] itsp; }

      //
      // Op=
      //

      inline Permutation& operator=(const Permutation& rhs)
      {
	for(size_t i=0;i<rhs.size();++i) itsp[i] = rhs.itsp[i]; 
	itsn = rhs.itsn;
	isident = rhs.isident;
	isinv = rhs.isinv;
	negdet = rhs.negdet;
	changed = true;
	return *this;
      }

      inline Permutation& operator=(const PermutationComposite& pcomp)
      { pcomp.AssignTo(*this); return *this; }


      //
      // Access 
      //

      inline size_t size() const { return itsn; }

      inline size_t colsize() const { return itsn; }

      inline size_t rowsize() const { return itsn; }

      inline size_t operator[](size_t i) const 
      { TMVAssert(i<size()); return itsp[i]; }

      inline size_t operator()(size_t i,size_t j) const
      { 
	TMVAssert(i<size()&&j<size());
	return (isinv ? i == itsp[j] : itsp[i]==j) ? 1 : 0; 
      }

      inline bool IsInverse() const { return isinv; }

      template <class T, StorageType S> operator Matrix<T,S>() const 
      {
	Matrix<T,S> temp(size(),size(),T(0));
	if (isinv) for(size_t i=0;i<size();++i) temp(itsp[i],i) = 1;
	else  for(size_t i=0;i<size();++i) temp(i,itsp[i]) = 1;
	return temp;
      }

      //
      // Modifying Functions
      //

      inline Permutation& SetToIdentity()
      { 
	for(size_t i=0;i<size();++i) itsp[i] = i; 
	negdet = false;
	isident = true;
	isinv = false;
	changed = true;
	return *this;
      }

      inline Permutation& SwapRows(size_t i1, size_t i2)
      { 
	TMVAssert(isident || !isinv);
	TMVAssert(i1<size());
	TMVAssert(i2<size());
	swap(itsp[i1],itsp[i2]); 
	negdet = !negdet; 
	isident = false;
	isinv = false;
	changed = true; 
	return *this; 
      }

      inline Permutation& SwapCols(size_t i1, size_t i2)
      { 
	TMVAssert(isident || isinv);
	TMVAssert(i1<size());
	TMVAssert(i2<size());
	swap(itsp[i1],itsp[i2]); 
	negdet = !negdet; 
	isident = false;
	isinv = true;
	changed = true; 
	return *this; 
      }

      inline Permutation& SwapWith(Permutation& p2)
      { 
	TMVAssert(p2.size() == size());
	for(size_t i=0;i<size();++i) swap(itsp[i],p2.itsp[i]);
	swap(negdet,p2.negdet);
	swap(isident,p2.isident);
	swap(isinv,p2.isinv);
	changed = true;
	p2.changed = true;
	return *this;
      }

      inline Permutation& InvertSelf()
      { isinv = !isinv; return *this; }

      inline Permutation& TransposeSelf()
      { return InvertSelf(); }

      Permutation& InvertStorage();
      // Defined in TMV_Permutation.cpp


      //
      // Functions of Permutations
      //

      inline Permutation Inverse() const
      { 
	Permutation temp(size(),negdet);
	for(size_t i=0;i<size();++i) temp.itsp[itsp[i]] = i;
	temp.isinv = isinv;
	temp.isident = isident;
	temp.changed = true;
	return temp;
      }

      inline Permutation Transpose() const
      { return Inverse(); }

      inline int Det() const
      { return negdet ? -1 : 1; }

      inline int Trace() const
      { 
	int trace = 0; 
	for(size_t i=0;i<size();++i) if (itsp[i] == i) ++trace;
	return trace;
      }


      //
      // Operators
      //

      Permutation& operator^=(int n);
      // Defined in TMV_Permutation.cpp


      //
      // I/O
      //

      void Read(istream& fin);
      // Defined in TMV_Permutation.cpp

      void Write(ostream& fout) const;
      // Defined in TMV_Permutation.cpp

      void WriteCycles(ostream& fout) const;
      // Defined in TMV_Permutation.cpp

      void WriteMatrix(ostream& fout) const;
      // Defined in TMV_Permutation.cpp


      //
      // Check validity
      //

      bool IsValid() const;
      // Defined in TMV_Permutation.cpp

      // 
      // Helpers for arithmetic functions
      //
      template <class T> void LMultEq(const VectorView<T>& v) const;
      template <class T> void LMult(
	  const GenVector<T>& v1, const VectorView<T>& v0) const;
      template <class T> void LDivEq(const VectorView<T>& v) const;
      template <class T> void LDiv(
	  const GenVector<T>& v1, const VectorView<T>& v0) const;

      template <class T> void LMultEq(const MatrixView<T>& m) const;
      template <class T> void LMult(
	  const GenMatrix<T>& m1, const MatrixView<T>& m0) const;
      template <class T> void LDivEq(const MatrixView<T>& m) const;
      template <class T> void LDiv(
	  const GenMatrix<T>& m1, const MatrixView<T>& m0) const;

      void RMultEq(Permutation& p) const;
      void LMultEq(Permutation& p) const;
      void RDivEq(Permutation& p) const;
      void LDivEq(Permutation& p) const;

      void SetToPP(const Permutation& p1, const Permutation& p2);
      void SetToPPt(const Permutation& p1, const Permutation& p2);
      void SetToPtP(const Permutation& p1, const Permutation& p2);

      size_t* ptr() { return itsp; }
      const size_t* cptr() const { return itsp; }

    protected:

      size_t*const itsp;
      size_t itsn;
      bool negdet; // det = 1 or -1

      bool isident;
      bool isinv;
      mutable bool changed;
      mutable vector<vector<size_t> > cycles;

      Permutation(size_t n, bool _negdet) :
	itsp(new size_t[n]), itsn(n), negdet(_negdet),
	isident(true), isinv(false), changed(true)  {}

      // 
      // Auxiliary Functions
      //
      int DetermineNegDet() const;
      void MakeCycles() const;

  }; // Permutation


  // 
  // SortPermutation
  //

  enum ADType { ASCEND, DESCEND };
  enum COMPType { REAL_COMP, ABS_COMP, IMAG_COMP, ARG_COMP };

  template <class T, class VT> class VTIndex 
  {

    private :

      RealType(T) itsvalue;
      size_t itsi;

    public :

      VTIndex() : itsvalue(RealType(T)(0)), itsi(0) {}

      VTIndex(const VT& v, size_t i, ADType ad, COMPType comp) : itsi(i)
      {
	bool neg = ad==DESCEND;
	switch(comp) {
	  case REAL_COMP : itsvalue = neg ? -REAL(v[i]) : REAL(v[i]); break;
	  case ABS_COMP : itsvalue = neg ? -abs(v[i]) : abs(v[i]); break;
	  case IMAG_COMP : itsvalue = neg ? -IMAG(v[i]) : IMAG(v[i]); break;
	  case ARG_COMP : itsvalue = neg ? -ARG(v[i]) : ARG(v[i]); break;
	  default : TMVAssert(false);
	}
      }

      // Use default copy, op=, destructor

      size_t GetI() const { return itsi; }
      bool operator<(const VTIndex& rhs) const 
      { return itsvalue < rhs.itsvalue; }

  };

  template <class T, class VT> inline Permutation SortPermutation(
      const VT& v, ADType ad=ASCEND, COMPType comp=REAL_COMP)
  { 
    vector<VTIndex<T,VT> > index(v.size());
    for(size_t i=0;i<index.size();++i) index[i] = VTIndex<T,VT>(v,i,ad,comp);
    std::sort(index.begin(),index.end());

    vector<size_t> parray(index.size());
    for(size_t i=0;i<index.size();++i) parray[i] = index[i].GetI();
    return Permutation(parray,false);
  }

  template <class T> inline Permutation SortPermutation(
      const GenVector<T>& v, ADType ad=ASCEND, COMPType comp=REAL_COMP)
  { return SortPermutation<T,GenVector<T> >(v,ad,comp); }

  template <class T> inline Permutation SortPermutation(
      const vector<T>& v, ADType ad=ASCEND, COMPType comp=REAL_COMP)
  { return SortPermutation<T, vector<T> >(v,ad,comp); }
  
  //
  // ReversePermutation
  // 
  
  inline Permutation ReversePermutation(const size_t n)
  {
    vector<size_t> v(n);
    for(size_t i=0,j=n-1;i<n;++i,--j) v[i] = j;
    return Permutation(v,false);
  }

  //
  // Swap
  //

  inline void Swap(Permutation& p1, Permutation& p2)
  { 
    TMVAssert(p1.size()==p2.size());
    p1.SwapWith(p2); 
  }

  //
  // Functions of Permutations
  //

  inline Permutation Inverse(const Permutation& p)
  { return p.Inverse(); }

  inline Permutation Transpose(const Permutation& p)
  { return p.Inverse(); }

  inline int Det(const Permutation& p)
  { return p.Det(); }

  inline int Trace(const Permutation& p)
  { return p.Trace(); }


  //
  // Permutation ==, != Permutation
  //

  bool operator==(const Permutation& p1, const Permutation& p2);

  inline bool operator!=(const Permutation& p1, const Permutation& p2) 
  { return !(p1 == p2); }

  template <class T> bool operator==(
      const GenMatrix<T>& m, const Permutation& p);

  template <class T> inline bool operator==(
      const Permutation& p, const GenMatrix<T>& m)
  { return (m == p); }

  template <class T> inline bool operator!=(
      const GenMatrix<T>& m, const Permutation& p) 
  { return !(m==p); }

  template <class T> inline bool operator!=(
      const Permutation& p, const GenMatrix<T>& m)
  { return !(m == p); }


  //
  // I/O
  //

  inline ostream& operator<<(ostream& fout, const Permutation& p)
  { p.Write(fout); return fout;}

  istream& operator>>(istream& fin, Permutation& p);

  istream& operator>>(istream& fin, Permutation* p);

} // namespace tmv

#endif
