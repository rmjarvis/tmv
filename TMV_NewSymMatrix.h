//---------------------------------------------------------------------------
//
// This file defines the TMV SymMatrix and HermMatrix classes.
//
// Constructors:
//
//    SymMatrix is used for symmetric matrices, and HermMatrix is used
//    for Hermitian matrices.  For real matrices, these are the same thing:
//    A = A.Transpose().
//    But for complex, they are different:
//    A_sym = A_sym.Transpose()
//    A_herm = A_herm.Adjoint()
//
//    For these notes, I will always write SymMatrix, but (except where
//    otherwise indicated) everything applies the same for Sym and Herm.
//
//    Caveat: Complex Hermitian matrices are such that A = At, which 
//    implies taht their diagonal elements are real.  Most routines
//    involving HermMatrixes assume the reality of the diagonal.
//    However, it is possible to assign a non-real value to a diagonal
//    element.  If the user does this, incorrect answers may result.
//
//    In addition to the type template parameter (T), SymMatrixes have two
//    additional template parameters:
//        UpLoType uplo = Upper || Lower
//        StorageType stor = RowMajor || ColMajor
//
//        They both have default values, so you can omit both, or
//        just stor.  The default values are: {Upper, RowMajor}
//
//        The first, uplo, refers to which triangular half stores the actual 
//        data (the other triangle is unreferenced).
//
//        The storage follows the same meaning as for regular Matrices.
//
//    SymMatrix<T,uplo,stor>(size_t n)
//        Makes a Symangular Matrix with column size = row size = n
//        with _uninitialized_ values.
//
//    SymMatrix<T,uplo,stor>(const Matrix<T>& m)
//    SymMatrix<T,uplo,stor>(const SymMatrix<T>& m)
//        Makes a SymMatrix which copies the corresponding elements of m.
//
//    ConstSymMatrixView<T> SymMatrixViewOf(const Matrix<T>& m, uplo)
//        Makes a constant SymMatrix view of the corresponding part of m.
//        While this view cannot be modified, changing the original m
//        will cause corresponding changes in this view of m.
//
//    SymMatrixView<T> SymMatrixViewOf(Matrix<T>& m, uplo)
//        Makes a modifiable SymMatrix view of the corresponding part of m.
//        Modifying this matrix will change the corresponding elements in
//        the original Matrix.
//
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//        Return the dimensions of the SymMatrix
//
//    T& operator()(size_t i, size_t j)
//    T operator()(size_t i, size_t j) const
//        Return the (i,j) element of the SymMatrix
//
//    Vector& row(size_t i, size_t j1, size_t j2)
//        Return a portion of the ith row 
//        This range must be a valid range for the requested row
//        and be entirely in the upper or lower triangle.
//
//    Vector& col(size_t j, size_t i1, size_t i2)
//        Return a portion of the jth column
//        This range must be a valid range for the requested column
//        and be entirely in the upper or lower triangle.
//
//    Vector& diag()
//        Return the main diagonal
//
//    Vector& diag(int i, size_t j1, size_t j2)
//    Vector& diag(int i)
//        Return the super- or sub-diagonal i
//        If i > 0 return the super diagonal starting at m_0i
//        If i < 0 return the sub diagonal starting at m_|i|0
//        If j1,j2 are given, it returns the diagonal SubVector 
//        either from m_j1,i+j1 to m_j2,i+j2 (for i>0) 
//        or from m_|i|+j1,j1 to m_|i|+j2,j2 (for i<0)
//
// Modifying Functions:
//
//    Zero()
//    SetAllTo(T x) 
//        For HermMatrix, x must be real.
//    ConjugateSelf()
//    SetToIdentity(x = 1)
//    void Swap(SymMatrix& m1, SymMatrix& m2)
//        The SymMatrices must be the same size.
//
// Views of a SymMatrix:
//
//    SubMatrix(int i1, int i2, int j1, int j2, int istep=1, int jstep=1)
//        This member function will return a submatrix using rows i1 to i2
//        and columns j1 to j2 which refers
//        to the same physical elements as the original.
//        The submatrix must be completely contained within the upper 
//        or lower triangle.
//
//    SubVector(int i, int j, int istep, int jstep, int size)
//        Returns a SubVector which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.  The subvector must be completely with the upper or
//        lower triangle.
//
//    SubSymMatrix(int i1, int i2, int istep)
//        Returns the SymMatrix which runs from i1 to i2 along the diagonal
//        (not including i2) with an optional step, and includes the 
//        off diagonal in the same rows/cols.
//
//        For example, with an UpperSymMatrix of size 10, the x's below
//        are the original data, the O's are the SubSymMatrix returned
//        with the command SubSymMatrix(3,11,2), and the #'s are the 
//        SubSymMatrix returned with SubSymMatrix(0,3)
//
//        ###xxxxxxx
//         ##xxxxxxx
//          #xxxxxxx
//           OxOxOxO
//            xxxxxx
//             OxOxO
//              xxxx
//               OxO
//                xx
//                 O
//
//    Transpose(m)
//    Adjoint(m)
//    Conjugate(m)
//
//
// Functions of Matrixs:
//
//    Det(m)
//    Trace(m)
//    Norm(m) or NormF(m)
//    NormSq(m)
//    Norm1(m) 
//    Norm2(m) 
//    NormInf(m) 
//
//    Inverse(m)
//    InverseATA(m)
//
//    m.NewTranspose()
//    m.NewConjugate()
//    m.NewAdjoint()
//    m.NewInverse()
//    m.NewCopy()
//
//
// I/O: 
//
//    os << m 
//        Writes m to ostream os in the usual Matrix format
//
//    m.WriteCompact(os)
//        Writes m to ostream os in the following compact format:
//          size 
//          ( m(0,0) m(0,1) ... m(0,size) )
//          ( m(1,1) .. m(1,size) )
//          ...
//          ( m(size,size) )
//
//    is >> m
//        Reads m from istream is in the compact format
//        m must already be the correct size for this to work.
//
//    is >> mptr
//        If you do not know the size of the SymMatrix to be read, you can
//        use this form where mptr is a pointer to an undefined SymMatrix.
//        Then *mptr will be created with the right size using new,
//        so you should subsequently delete it.
//        (Note: if the DiagType for the SymMatrix is UnitDiag, then
//        all of the diagonals read in must be = 1.)
//
//
// Division Control Functions:
//
//    For SymMatrixes, LU actually does an LDLT decomposition.
//        ie. L(Lower Triangle) * Diagonal * L.Transpose()
//        For non-positive definite matrices, D may have 2x2 blocks along
//        the diagonal, which can then be further decomposed via normal LU
//        decomposition.
//        
//    We add the option of CH for SymMatrixes - Cholskey decomposition.
//        This is only appropriate if you know that your SymMatrix is 
//        positive definite (ie. all eigenvalues are positive).
//        This is guaranteed, for example, if all diagonal elements are 
//        positive and all are greater than the largest off-diagonal element
//        (in absolute value).
//        In this case, the SymMatrix can be decomposed into L*L.Transpose()
//        (for SymMatrixes) or L*L.Adjoint() (for HermMatrixes).
//
//    QR and SV are also still options, but they are not fast for SymMatrixes.  
//        For these, we simply copy the SymMatrix to a normal Matrix
//        and decompose that.
//


#ifndef TMV_SymMatrix_H
#define TMV_SymMatrix_H

#include "TMV_BaseMatrix.h"

namespace tmv {

  template <class T> class GenSymMatrix;
  template <class T> class ConstSymMatrixView;
  template <class T> class SymMatrixView;
  template <class T, DiagType D=NonUnitDiag, StorageType S=RowMajor> 
    class SymMatrix;
  template <class T> class SymMatrixComposite; 
  template <class T> class SymDivider;
  template <class T> class MetaSymDivider;

  template <class T> class SymLUDiv;
  template <class T> class SymCHDiv;
  template <class T> class HermLUDiv;
  template <class T> class HermCHDiv;
  template <class T> class SymQRDiv;
  template <class T> class SymSVDiv;

  template <class T1, class T2> void Copy(
      const GenSymMatrix<T1>& m1, const SymMatrixView<T2>& m2);

  template <class T> class GenSymMatrix : 
    public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      GenSymMatrix(SymType h, UpLoType u, StorageType s, ConjItType c) : 
	itssym(h), itsuplo(u), itsstor(s), itsct(c) {}
      GenSymMatrix(SymType h, UpLoType u, StorageType s, ConjItType c,
	  const MetaDivider<T>* div) : 
	BaseMatrix<T>(div), itssym(h), itsuplo(u), itsstor(s), itsct(c) {}
      GenSymMatrix(const GenSymMatrix<T>& rhs) : 
	BaseMatrix<T>(rhs), itssym(rhs.itssym), itsuplo(rhs.itsuplo), 
	itsstor(rhs.itsstor), itsct(rhs.itsct) {}
      ~GenSymMatrix() {}

      //
      // Access Functions
      //

      inline size_t colsize() const { return size(); }
      inline size_t rowsize() const { return size(); }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert( i<=j1 || j2<=i );
	if ((i<=j1 && uplo()==Upper) || (j2<=i && uplo()==Lower))
	  return ConstVectorView<T>(cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
	      itsct); 
	else
	  return ConstVectorView<T>(cptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
	      isherm() ? ConjOf(T,itsct) : itsct);
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert( i2<=j || j<=i1 );
	if ((i2<=j && uplo()==Upper) || (j<=i1 && uplo()==Lower))
	  return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
	      itsct); 
	else 
	  return ConstVectorView<T>(cptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
	      isherm() ? ConjOf(T,itsct) : itsct);
      }

      inline ConstVectorView<T> diag() const
      { return ConstVectorView<T>(cptr(),size(),stepi()+stepj(),itsct); }

      inline ConstVectorView<T> diag(int i=0) const
      {
	TMVAssert(i<=int(size())); 
	if (i>=0)
	  if (uplo()==Upper)
	    return ConstVectorView<T>(cptr()+i*stepj(),size()-i,
		stepi()+stepj(),itsct); 
	  else
	    return ConstVectorView<T>(cptr()+i*stepi(),size()-i,
		stepi()+stepj(), isherm() ? ConjOf(T,itsct) : itsct);
	else
	  if (uplo()==Upper)
	    return ConstVectorView<T>(cptr()-i*stepj(),size()+i,
		stepi()+stepj(), isherm() ? ConjOf(T,itsct) : itsct);
	  else
	    return ConstVectorView<T>(cptr()-i*stepi(),size()+i,
		stepi()+stepj(),itsct); 
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(isunit() ? i>0 : i>=0);
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()-abs(i));
	if (i>=0)
	  if (uplo()==Upper)
	    return ConstVectorView<T>(cptr()+i*stepj(),j2-j1,
		stepi()+stepj(),itsct); 
	  else
	    return ConstVectorView<T>(cptr()+i*stepi(),j2-j1,
		stepi()+stepj(), isherm() ? ConjOf(T,itsct) : itsct);
	else
	  if (uplo()==Upper)
	    return ConstVectorView<T>(cptr()-i*stepj(),j2-j1,
		stepi()+stepj(), isherm() ? ConjOf(T,itsct) : itsct);
	  else
	    return ConstVectorView<T>(cptr()-i*stepi(),j2-j1,
		stepi()+stepj(),itsct); 
      }

      template <class T2> inline bool SameStorageAs(
	  const BaseMatrix<T2>& m2) const
      { return false; }

      inline bool SameStorageAs(const GenMatrix<T>& m2) const
      { return (cptr()==m2.cptr()); }

      inline bool SameStorageAs(const GenSymMatrix<T>& m2) const
      { return (cptr()==m2.cptr()); }

      template <class T2> inline bool SameAs(const BaseMatrix<T2>& m2) const
      { return false; }

      inline bool SameAs(const GenSymMatrix<T>& m2) const
      { 
	return (this == &m2 || (cptr()==m2.cptr() && size()==m2.size() && 
	      st() == m2.st() && uplo() == m2.uplo() && itsct == m2.itsct &&
	      stepi()==m2.stepi() && stepj()==m2.stepj()));
      }

      template <class T2> inline void DoCopyToMatrix(
	  const MatrixView<T2>& m2) const
      {
	TMVAssert(m2.colsize() == size());
	TMVAssert(m2.rowsize() == size());
	m2.diag() = diag();
	UpperTriMatrixViewOf(m2).OffDiag() = UpperOffDiag();
	LowerTriMatrixViewOf(m2).OffDiag() = LowerOffDiag();
      }
      inline void CopyToMatrix(const MatrixView<RealType(T)>& m2) const
      { TMVAssert(IsReal(T())); DoCopyToMatrix(m2); }
      inline void CopyToMatrix(const MatrixView<ComplexType(T)>& m2) const
      { DoCopyToMatrix(m2); }

      //
      // SubMatrix
      //

      bool OKSubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const;

      inline ConstMatrixView<T> SubMatrix(int i1, int i2,
	  int j1, int j2) const
      {
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,1,1));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(), stepj(), stor(), itsct);
      }

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	const StorageType newstor =
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	    isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(), 
	    newstor, itsct);
      }

      bool OKSubVector(size_t i, size_t j, int istep, int jstep, 
	  size_t size) const;

      inline ConstVectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(),itsct);
      }

      bool OKSubSymMatrix(int i1, int i2, int istep) const;

      inline ConstSymMatrixView<T> SubSymMatrix(int i1, int i2) const
      {
	TMVAssert(OKSubSymMatrix(i1,i2,1));
	return ConstSymMatrixView<T>(cptr()+i1*
	    stepi()+stepj(),i2-i1,
	    stepi(),stepj(),st(),uplo(),stor(),itsct);
      }

      inline ConstSymMatrixView<T> SubSymMatrix(int i1, int i2,
	  int istep) const
      {
	TMVAssert(OKSubSymMatrix(i1,i2,istep));
	const int diagstep = stepi()+stepj();
	return ConstSymMatrixView<T>(cptr()+i1*diagstep,(i2-i1)/istep,
	    istep*stepi(),istep*stepj(),istep*diagstep,st(),uplo(),
	    istep==1 ? stor() : NoMajor,itsct);
      }

      inline ConstUpperTriMatrixView<T> UpperOffDiag() const
      {
	TMVAssert(size() > 0);
	if (uplo() == Upper)
	  return ConstUpperTriMatrixView<T>(ptr()+stepj(),size()-1,
	      stepi(),stepj(),NonUnitDiag,stor(),itsct FIRSTLAST);
	else
	  return ConstUpperTriMatrixView<T>(ptr()+stepi(),size()-1,
	      stepj(),stepi(),NonUnitDiag,TransOf(stor()),
	      isherm()?ConjOf(T,itsct):itsct FIRSTLAST);
      }

      inline ConstLowerTriMatrixView<T> LowerOffDiag() const
      {
	TMVAssert(size() > 0);
	if (uplo() == Lower)
	  return ConstLowerTriMatrixView<T>(cptr()+stepi(),size()-1,
	      stepi(),stepj(),NonUnitDiag,stor(),itsct);
	else
	  return ConstLowerTriMatrixView<T>(cptr()+stepj(),size()-1,
	      stepj(),stepi(),NonUnitDiag,TransOf(stor()),
	      isherm()?ConjOf(T,itsct):itsct);
      }

      inline ConstSymMatrixView<T> View() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepi(),stepj(),
	    st(),uplo(),stor(),itsct,
	    MakeMetaDivider(false,false,this));
      }

      inline ConstSymMatrixView<T> Transpose() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepj(),stepi(),
	    st(),TransOf(uplo()),TransOf(stor()),itsct,
	    MakeMetaDivider(true,false,this));
      }

      inline ConstSymMatrixView<T> Conjugate() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepi(),stepj(),
	    st(),uplo(),stor(),ConjOf(T,itsct),
	    MakeMetaDivider(false,true,this));
      }

      inline ConstSymMatrixView<T> Adjoint() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepj(),stepi(),
	    st(),TransOf(uplo()),TransOf(stor()),ConjOf(T,itsct),
	    MakeMetaDivider(true,true,this));
      }

      inline ConstSymMatrixView<T> QuickView() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepi(),stepj(),
	    st(),uplo(),stor(),itsct);
      }

      inline ConstSymMatrixView<T> QuickTranspose() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepj(),stepi(),
	    st(),TransOf(uplo()),TransOf(stor()),itsct);
      }

      inline ConstSymMatrixView<T> QuickConjugate() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepi(),stepj(),
	    st(),uplo(),stor(),ConjOf(T,itsct));
      }

      inline ConstSymMatrixView<T> QuickAdjoint() const
      { 
	return ConstSymMatrixView<T>(cptr(),size(),stepj(),stepi(),
	    st(),TransOf(uplo()),TransOf(stor()),ConjOf(T,itsct));
      }

      inline ConstSymMatrixView<RealType(T)> Real() const
      {
	return ConstSymMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    Sym, uplo(), IsReal(T()) ? stor() : NoMajor,NonConj);
      }

      inline ConstSymMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(!isherm());
	// The imaginary part of a Hermitian matrix is anti-symmetric
	return ConstSymMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,
	    size(),2*stepi(),2*stepj(), Sym,uplo(),NoMajor,NonConj);
      }

      //
      // Functions of Matrix
      //

      inline T Trace() const
      { return diag().SumElements(); }

      inline RealType(T) Norm() const 
      { return SQRT(NormSq()); }

      inline RealType(T) NormF() const 
      { return SQRT(NormSq()); }

      // NormF()^2
      RealType(T) NormSq() const;

      // 1-Norm = max_j (sum_i |a_ij|)
      RealType(T) Norm1() const;

      // 2-Norm defined in BaseMatrix.h

      // inf-Norm = max_i (sum_j |a_ij|)
      RealType(T) NormInf() const
      { return Norm1(); }

      // = max_i,j (|a_ij|)
      RealType(T) MaxAbsElement() const;

      inline BaseMatrix<T>* NewTranspose() const 
      { return new ConstSymMatrixView<T>(Transpose()); }

      inline BaseMatrix<T>* NewConjugate() const
      { return new ConstSymMatrixView<T>(Conjugate()); }

      inline BaseMatrix<T>* NewAdjoint() const
      { return new ConstSymMatrixView<T>(Adjoint()); }

      inline BaseMatrix<T>* NewView() const
      { return new ConstSymMatrixView<T>(View()); }

      inline BaseMatrix<T>* NewCopy() const
      { 
	if (uplo() == Upper)
	  if (isrm()) return new SymMatrix<T,Upper,RowMajor>(*this); 
	  else return new SymMatrix<T,Upper,ColMajor>(*this); 
	else
	  if (isrm()) return new SymMatrix<T,Lower,RowMajor>(*this); 
	  else return new SymMatrix<T,Lower,ColMajor>(*this); 
      }

      // 
      // Division Control
      //

      inline void DivideUsing(DivType dt) const
      { BaseMatrix<T>::DivideUsing(dt); }

      //
      // I/O
      //

      void WriteCompact(ostream& os) const;

      //
      // Arithmetic Helpers
      //

      using BaseMatrix<T>::LDivEq;
      using BaseMatrix<T>::RDivEq;
      using BaseMatrix<T>::LDiv;
      using BaseMatrix<T>::RDiv;

      virtual const T* cptr() const = 0;
      virtual int stepi() const = 0;
      virtual int stepj() const = 0;
      virtual size_t size() const = 0;
      virtual inline SymType st() const { return itssym; }
      virtual inline UpLoType uplo() const { return itsuplo; }
      virtual inline StorageType stor() const { return itsstor; }
      virtual inline ConjItType ct() const { return itsct; }
      virtual inline bool isrm() const { return itsstor == RowMajor; }
      virtual inline bool iscm() const { return itsstor == ColMajor; }
      virtual inline bool isconj() const
      { 
	TMVAssert(IsComplex(T()) || itsct==NonConj);
	return IsComplex(T()) && itsct == Conj; 
      }
      virtual inline bool isherm() const 
      { 
	TMVAssert(IsComplex(T()) || itssym==Sym);
	return IsComplex(T()) && itssym == Herm; 
      }

    protected :

      virtual T cref(size_t i, size_t j) const;

      void NewDivider() const;

    private :

      SymType itssym;
      UpLoType itsuplo;
      StorageType itsstor;
      ConjItType itsct;

      void operator=(const GenSymMatrix<T>&) { TMVAssert(false); }

  }; // GenSymMatrix

  template <class T> class ConstSymMatrixView : 
    public GenSymMatrix<T>
  {
    public :

      ConstSymMatrixView(const ConstSymMatrixView<T>& rhs) :
	GenSymMatrix<T>(rhs), itsm(rhs.itsm),
	itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj) {}

      ConstSymMatrixView(const GenSymMatrix<T>& rhs) :
	GenSymMatrix<T>(rhs), itsm(rhs.cptr()),
	itss(rhs.size()), itssi(rhs.stepi()), itssj(rhs.stepj()) {}

      ConstSymMatrixView(const T* _m, size_t _s, int _si, int _sj,
	  SymType insym, UpLoType inuplo, StorageType instor,
	  ConjItType inct) : 
	GenSymMatrix<T>(insym,inuplo,instor,inct), itsm(_m), itss(_s),
	itssi(_si), itssj(_sj)
      { 
	TMVAssert(instor==RowMajor ? _sj == 1 : instor==ColMajor ?
	    _si==1 : true);
      }

      ConstSymMatrixView(const T* _m, size_t _s, int _si, int _sj, 
	  SymType insym, UpLoType inuplo, StorageType instor, ConjItType inct,
	  const MetaDivider<T>* mdiv) : 
	GenSymMatrix<T>(insym,inuplo,instor,inct,mdiv), 
	itsm(_m), itss(_s), itssi(_si), itssj(_sj)
      { 
	TMVAssert(instor==RowMajor ? _sj == 1 : instor==ColMajor ?
	    _si==1 : true);
      }

      ConstSymMatrixView(const GenMatrix<T>& rhs, 
	  SymType insym, UpLoType inuplo) : 
	GenSymMatrix<T>(insym,inuplo,rhs.stor(),rhs.ct()),
	itsm(rhs.cptr()), itss(rhs.colsize()), 
	itssi(rhs.stepi()), itssj(rhs.stepj())
      { TMVAssert(rhs.IsSquare()); }

      ~ConstSymMatrixView() {}

      inline size_t size() const { return itss; }

      inline const T* cptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }

    protected :

      const T*const itsm;
      const size_t itss;
      const int itssi;
      const int itssj;

    private :

      inline void operator=(const ConstSymMatrixView<T>&) 
      { TMVAssert(false); }

  }; // ConstSymMatrixView

  template <class T> class SymMatrixView : 
    public GenSymMatrix<T>
  {

    public:

      //
      // Constructors
      //

      SymMatrixView(const SymMatrixView<T>& rhs) : 
	GenSymMatrix<T>(rhs), itsm(rhs.itsm), itss(rhs.itss), 
	itssi(rhs.itssi), itssj(rhs.itssj)
	  DEFFIRSTLAST(rhs.first,rhs.last) {}

      SymMatrixView(const SymMatrixView<T>& rhs, DiagType d) : 
	GenSymMatrix<T>(d,rhs.stor(),rhs.ct()), itsm(rhs.itsm), 
	itss(rhs.itss), itssi(rhs.itssi), itssj(rhs.itssj)
	  DEFFIRSTLAST(rhs.first,rhs.last) 
      { TMVAssert(d!=NonUnitDiag || !rhs.isunit()); }

      SymMatrixView(T* _m, size_t _s, int _si, int _sj, 
	  DiagType indt, StorageType instor, ConjItType inct 
	  PARAMFIRSTLAST(T) ) :
	GenSymMatrix<T>(indt,instor,inct), itsm(_m), itss(_s),
	itssi(_si), itssj(_sj) DEFFIRSTLAST(_first,_last)
      {
	TMVAssert(instor==RowMajor ? _sj==1 : instor==ColMajor ?
	  _si==1 : true); 
      }

      SymMatrixView(T* _m, size_t _s, int _si, int _sj,
	  DiagType indt, StorageType instor, ConjItType inct,
	  const MetaDivider<T>* mdiv PARAMFIRSTLAST(T) ) :
	GenSymMatrix<T>(indt,instor,inct,mdiv), itsm(_m), itss(_s),
	itssi(_si), itssj(_sj) DEFFIRSTLAST(_first,_last)
      {
	TMVAssert(instor==RowMajor ? _sj==1 : instor==ColMajor ?
	  _si==1 : true); 
      }

      SymMatrixView(const MatrixView<T>& rhs, DiagType dt) : 
	GenSymMatrix<T>(dt,rhs.stor(),rhs.ct()),
	itsm(rhs.ptr()), itss(rhs.colsize()), 
	itssi(rhs.stepi()), itssj(rhs.stepj())
	  DEFFIRSTLAST(rhs.first,rhs.last) 
      { TMVAssert(rhs.IsSquare()); }

      ~SymMatrixView() {} 

      //
      // Op=
      //

      template <class T2> inline void DoCopy(
	  const GenSymMatrix<T2>& m2) const
      { 
	diag() = m2.diag();
	UpperOffDiag() = m2.UpperOffDiag();
	LowerOffDiag() = m2.LowerOffDiag();
      }

      inline const SymMatrixView<T>& operator=(
	  const SymMatrixView<T>& m2) const
      { if (!SameAs(m2)) DoCopy(m2); return *this; }

      inline const SymMatrixView<T>& operator=(
	  const GenSymMatrix<T>& m2) const
      { if (!SameAs(m2)) DoCopy(m2); return *this; }

      template <class T2> inline const SymMatrixView<T>& operator=(
	  const GenSymMatrix<T2>& m2) const
      { DoCopy(m2); return *this; }

      inline const SymMatrixView<T>& operator=(T x) const 
      { return SetToIdentity(x); }

      inline const SymMatrixView<T>& operator=(
	  const SymMatrixComposite<T>& mcomp) const
      { 
	TMVAssert(size() == mcomp.size());
	mcomp.AssignTo(QuickView());
	return *this;
      }

      //
      // Access
      //

      inline size_t size() const { return itss; }

      inline VectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert( i<=j1 || j2<=i );
	if ((i<=j1 && uplo()==Upper) || (j2<=i && uplo()==Lower))
	  return VectorView<T>(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),
	      ct() FIRSTLAST); 
	else
	  return VectorView<T>(ptr()+i*stepj()+j1*stepi(),j2-j1,stepi(),
	      isherm() ? ConjOf(T,ct()) : ct() FIRSTLAST);
      }

      inline VectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert( i2<=j || j<=i1 );
	if ((i2<=j && uplo()==Upper) || (j<=i1 && uplo()==Lower))
	  return VectorView<T>(ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),
	      ct() FIRSTLAST); 
	else 
	  return VectorView<T>(ptr()+i1*stepj()+j*stepi(),i2-i1,stepj(),
	      isherm() ? ConjOf(T,ct()) : ct() FIRSTLAST);
      }

      inline VectorView<T> diag() const
      { return VectorView<T>(ptr(),size(),stepi()+stepj(),ct() FIRSTLAST); }

      inline VectorView<T> diag(int i=0) const
      {
	TMVAssert(i<=int(size())); 
	if (i>=0)
	  if (uplo()==Upper)
	    return VectorView<T>(ptr()+i*stepj(),size()-i,
		stepi()+stepj(),ct() FIRSTLAST); 
	  else
	    return VectorView<T>(ptr()+i*stepi(),size()-i,
		stepi()+stepj(), isherm() ? ConjOf(T,ct()) : ct() FIRSTLAST);
	else
	  if (uplo()==Upper)
	    return VectorView<T>(ptr()-i*stepj(),size()+i,
		stepi()+stepj(), isherm() ? ConjOf(T,ct()) : ct() FIRSTLAST);
	  else
	    return VectorView<T>(ptr()-i*stepi(),size()+i,
		stepi()+stepj(),ct() FIRSTLAST); 
      }

      inline VectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(isunit() ? i>0 : i>=0);
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()-abs(i));
	if (i>=0)
	  if (uplo()==Upper)
	    return VectorView<T>(ptr()+i*stepj(),j2-j1,
		stepi()+stepj(),ct() FIRSTLAST); 
	  else
	    return VectorView<T>(ptr()+i*stepi(),j2-j1,
		stepi()+stepj(), isherm() ? ConjOf(T,ct()) : ct() FIRSTLAST);
	else
	  if (uplo()==Upper)
	    return VectorView<T>(ptr()-i*stepj(),j2-j1,
		stepi()+stepj(), isherm() ? ConjOf(T,ct()) : ct() FIRSTLAST);
	  else
	    return VectorView<T>(ptr()-i*stepi(),j2-j1,
		stepi()+stepj(),ct() FIRSTLAST); 
      }

      inline RefType(T) operator()(size_t i,size_t j) const 
      { return ref(i,j); }

      //
      // Modifying Functions
      //

      inline const SymMatrixView<T>& Zero() const 
      { return SetAllTo(T(0)); }
      const SymMatrixView<T>& SetAllTo(T x) const;
      const SymMatrixView<T>& ConjugateSelf() const;
      const SymMatrixView<T>& SetToIdentity(T x=T(1)) const;

      //
      // SubMatrix
      //

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) const
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1, j2-j1, stepi(),stepj(),stor(),ct() FIRSTLAST );
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	const StorageType newstor =
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	  isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
	    newstor,ct() FIRSTLAST );
      }

      inline VectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return VectorView<T>(ptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep*stepj(),ct() FIRSTLAST );
      }

      inline SymMatrixView<T> SubSymMatrix(int i1, int i2) const
      {
	TMVAssert(this->OKSubSymMatrix(i1,i2,1));
	return SymMatrixView<T>(ptr()+i1*(stepi()+stepj()),i2-i1,
	    stepi(),stepj(),sym(),uplo(),stor(),ct() FIRSTLAST);
      }

      inline SymMatrixView<T> SubSymMatrix(int i1, int i2,
	  int istep) const
      {
	TMVAssert(this->OKSubSymMatrix(i1,i2,istep));
	const int diagstep = stepi()+stepj();
	return SymMatrixView<T>(ptr()+i1*(stepi()+stepj()),(i2-i1)/istep,
	    istep*stepi(),istep*stepj(),sym(),uplo(),
	    istep==1 ? stor() : NoMajor,ct() FIRSTLAST);
      }

      inline UpperTriMatrixView<T> UpperOffDiag() const
      {
	TMVAssert(size() > 0);
	if (uplo() == Upper)
	  return UpperTriMatrixView<T>(ptr()+stepj(),size()-1,
	      stepi(),stepj(),NonUnitDiag,stor(),ct() FIRSTLAST);
	else
	  return UpperTriMatrixView<T>(ptr()+stepi(),size()-1,
	      stepj(),stepi(),NonUnitDiag,TransOf(stor()),
	      isherm()?ConjOf(T,ct()):ct() FIRSTLAST);
      }

      inline LowerTriMatrixView<T> LowerOffDiag() const
      {
	TMVAssert(size() > 0);
	if (uplo() == Lower)
	  return LowerTriMatrixView<T>(ptr()+stepi(),size()-1,
	      stepi(),stepj(),NonUnitDiag,stor(),ct() FIRSTLAST);
	else
	  return LowerTriMatrixView<T>(ptr()+stepj(),size()-1,
	      stepj(),stepi(),NonUnitDiag,TransOf(stor()),
	      isherm()?ConjOf(T,ct()):ct() FIRSTLAST);
      }

      inline SymMatrixView<T> View() const
      { return *this; }

      inline SymMatrixView<T> Transpose() const
      {
	return SymMatrixView<T>(ptr(),size(),stepj(),stepi(),
	    sym(),TransOf(uplo()),TransOf(stor()),ct(),
	    MakeMetaDivider(true,false,this) FIRSTLAST );
      }

      inline SymMatrixView<T> Conjugate() const
      {
	return SymMatrixView<T>(ptr(),size(),stepi(),stepj(),
	    sym(),uplo(),stor(),ConjOf(T,ct()),
	    MakeMetaDivider(false,true,this) FIRSTLAST ); 
      }

      inline SymMatrixView<T> Adjoint() const
      {
	return SymMatrixView<T>(ptr(),size(),stepj(),stepi(),
	    sym(),TransOf(uplo()),TransOf(stor()),ConjOf(T,ct()),
	    MakeMetaDivider(true,true,this) FIRSTLAST );
      }

      inline SymMatrixView<T> QuickView() const
      { 
	return SymMatrixView<T>(ptr(),size(),stepi(),stepj(),
	    sym(),uplo(),stor(),ct() FIRSTLAST);
      }

      inline SymMatrixView<T> QuickTranspose() const
      {
	return SymMatrixView<T>(ptr(),size(),stepj(),stepi(),
	    sym(),TransOf(uplo()),TransOf(stor()),ct() FIRSTLAST);
      }

      inline SymMatrixView<T> QuickConjugate() const
      {
	return SymMatrixView<T>(ptr(),size(),stepi(),stepj(),
	    sym(),uplo(),stor(),ConjOf(T,ct()) FIRSTLAST);
      }

      inline SymMatrixView<T> QuickAdjoint() const
      {
	return SymMatrixView<T>(ptr(),size(),stepj(),stepi(),
	    sym(),TransOf(uplo()),TransOf(stor()),ConjOf(T,ct()) FIRSTLAST);
      }

      inline SymMatrixView<RealType(T)> Real() const
      {
	return SymMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr()),size(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    sym(),uplo(), IsReal(T()) ? stor() : NoMajor, NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)
	    ,reinterpret_cast<const RealType(T)*>(last)
#endif
	    );
      }

      inline SymMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(!isherm());
	return SymMatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,
	    size(),2*stepi(),2*stepj(),sym(),uplo(),NoMajor, NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      //
      // I/O
      //

      void Read(istream& is) const;

      inline const T* cptr() const { return itsm; }
      inline T* ptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      using GenSymMatrix<T>::sym;
      using GenSymMatrix<T>::uplo;
      using GenSymMatrix<T>::stor;
      using GenSymMatrix<T>::ct;
      using GenSymMatrix<T>::isrm;
      using GenSymMatrix<T>::iscm;
      using GenSymMatrix<T>::isherm;

      // This makes it easier for some things to compile
      // The ones that work are overridden
      template <class T2> void operator=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }
      template <class T2> void operator+=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }
      template <class T2> void operator-=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }
      template <class T2> void operator*=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }
      template <class T2> void operator/=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }
      template <class T2> void operator%=(
	  const BaseMatrix<T2>&) const { TMVAssert(false); }

    protected :

      T*const itsm;
      const size_t itss;
      const int itssi;
      const int itssj;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
    protected:
#endif

      RefType(T) ref(size_t i, size_t j) const;

  }; // SymMatrixView

  template <class T, UpLoType U, StorageType S> class SymMatrix : 
    public GenSymMatrix<T>
  {

    public:

      //
      // Constructors
      //

#define NEW_SIZE(s) GenSymMatrix<T>(Sym,U,S,NonConj), \
      itslen((s)*(s)), itsm(new T[itslen]), itss(s)  \
      DEFFIRSTLAST(itsm,itsm+itslen)

      explicit UpperTriMatrix(size_t _size) : NEW_SIZE(_size) 
      { TMVAssert(S==RowMajor); }

      UpperTriMatrix(const UpperTriMatrix<T,D,S>& rhs) : NEW_SIZE(rhs.size())
      {
	TMVAssert(S==RowMajor); 
	memmove(itsm,rhs.itsm,itslen*sizeof(T));
      }

      template <DiagType D2> UpperTriMatrix(
	  const UpperTriMatrix<T,D2,S>& rhs) : NEW_SIZE(rhs.size())
      {
	TMVAssert(S==RowMajor); 
	memmove(itsm,rhs.cptr(),itslen*sizeof(T));
	if (D==NonUnitDiag && D2==UnitDiag) diag().SetAllTo(T(1));
      }

      template <DiagType D2, StorageType S2> UpperTriMatrix(
	  const UpperTriMatrix<T,D2,S2>& rhs) : NEW_SIZE(rhs.size())
      {
	TMVAssert(S==RowMajor); 
	if (D==D2) Copy(rhs,QuickView());
	else {
	  Copy(rhs.OffDiag(),OffDiag());
	  if (D==NonUnitDiag && D2==UnitDiag) diag().SetAllTo(T(1));
	}
      }

      template <StorageType S2> UpperTriMatrix(const Matrix<T,S2>& rhs) :
	NEW_SIZE(rhs.rowsize())
      {
	TMVAssert(S==RowMajor); 
	if (S==S2 && rhs.IsSquare()) 
	  memmove(itsm,rhs.cptr(),itslen*sizeof(T));
	else 
	  Copy(UpperTriMatrixViewOf(rhs,D),QuickView());
      }

      template <class T2> UpperTriMatrix(const GenMatrix<T2>& rhs) :
	NEW_SIZE(rhs.rowsize())
      { 
	TMVAssert(S==RowMajor); 
	Copy(UpperTriMatrixViewOf(rhs,D),QuickView()); 
      }

      template <class T2> UpperTriMatrix(const GenUpperTriMatrix<T2>& rhs) :
	NEW_SIZE(rhs.size())
      { 
	TMVAssert(S==RowMajor); 
	if (D==rhs.dt()) Copy(rhs,QuickView());
	else {
	  Copy(rhs.OffDiag(),OffDiag());
	  if (D==NonUnitDiag && rhs.dt()==UnitDiag) diag().SetAllTo(T(1));
	}
      }

      UpperTriMatrix(const UpperTriMatrixComposite<T>& mcomp) :
	NEW_SIZE(mcomp.size())
      {
	TMVAssert(S==RowMajor); 
	TMVAssert(!(mcomp.dt()==NonUnitDiag && D==UnitDiag));
	if (mcomp.dt()==UnitDiag && D==NonUnitDiag) {
	  mcomp.AssignTo(UpperTriMatrixViewOf(*this,UnitDiag));
	  diag().SetAllTo(T(1));
	} else 
	  mcomp.AssignTo(QuickView());
      }

#undef NEW_SIZE

      ~UpperTriMatrix() { TMVAssert(itsm); delete[] itsm; }


      //
      // Op=
      //

      inline UpperTriMatrix<T,D,S>& operator=(const UpperTriMatrix<T,D,S>& m2)
      { 
	if (&m2 != this) memmove(itsm,m2.itsm,size()*size()*sizeof(T));
	return *this;
      }

      template <class T2> inline UpperTriMatrix<T,D,S>& operator=(
	  const GenUpperTriMatrix<T2>& m2)
      { 
	TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	if (m2.dt()==UnitDiag && D==NonUnitDiag) {
	  Copy(m2.OffDiag(),OffDiag());
	  diag().SetAllTo(T(1));
	} else 
	  Copy(m2,QuickView());
	return *this;
      }

      inline UpperTriMatrix<T,D,S>& operator=(T x) 
      { 
	TMVAssert(!this->isunit() || x==T(1));
	return SetToIdentity(x); 
      }

      inline UpperTriMatrix<T,D,S>& operator=(
	  const UpperTriMatrixComposite<T>& mcomp)
      { 
	TMVAssert(size() == mcomp.size());
	TMVAssert(!(mcomp.dt()==NonUnitDiag && D==UnitDiag));
	if (mcomp.dt()==UnitDiag && D==NonUnitDiag) {
	  mcomp.AssignTo(UpperTriMatrixViewOf(*this,UnitDiag));
	  diag().SetAllTo(T(1));
	} else 
	  mcomp.AssignTo(QuickView());
	return *this;
      }

      //
      // Access
      //

      inline size_t size() const { return itss; }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert(j1==j2 || okij(i,j1));
	return ConstVectorView<T>(cptr()+i*stepi()+j1,j2-j1,1,NonConj);
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert(i1==i2 || okij(i2-1,j));
	return ConstVectorView<T>(cptr()+i1*stepi()+j,i2-i1,stepi(),NonConj);
      }

      inline ConstVectorView<T> diag(int i=0) const
      {
	TMVAssert(i<=int(size())); 
	TMVAssert(isunit() ? i>0 : i>=0);
	return ConstVectorView<T>(cptr()+i,size()-i,NonConj);
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(isunit() ? i>0 : i>=0);
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()+i);
	return ConstVectorView<T>(cptr()+i+j1*(stepi()+1),j2-j1,NonConj);
      }

      inline VectorView<T> row(size_t i, size_t j1, size_t j2)
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert(j1==j2 || okij(i,j1));
	return VectorView<T>(ptr()+i*stepi()+j1,j2-j1,1,NonConj FIRSTLAST);
      }

      inline VectorView<T> col(size_t j, size_t i1, size_t i2)
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert(i1==i2 || okij(i2-1,j));
	return VectorView<T>(ptr()+i1*stepi()+j,i2-i1,stepi(),NonConj
	    FIRSTLAST);
      }

      inline VectorView<T> diag(int i=0)
      {
	TMVAssert(i<=int(size())); 
	TMVAssert(isunit() ? i>0 : i>=0);
	return VectorView<T>(ptr()+i,size()-i,NonConj FIRSTLAST);
      }

      inline VectorView<T> diag(int i, size_t j1, size_t j2)
      {
	TMVAssert(isunit() ? i>0 : i>=0);
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()+i);
	return VectorView<T>(ptr()+i+j1*(stepi()+1),j2-j1,NonConj FIRSTLAST);
      }

      inline T operator()(size_t i, size_t j) const
      { return cref(i,j); }

      inline T& operator()(size_t i, size_t j) 
      { return ref(i,j); }

      //
      // Modifying Functions
      //

      inline UpperTriMatrix<T,D,S>& Zero() { return SetAllTo(0); }
      inline UpperTriMatrix<T,D,S>& SetAllTo(T x) 
      { QuickView().SetAllTo(x); return *this; }
      inline UpperTriMatrix<T,D,S>& ConjugateSelf() 
      { QuickView().ConjugateSelf(); return *this; }
      inline UpperTriMatrix<T,D,S>& SetToIdentity(T x=T(1)) 
      { QuickView().SetToIdentity(x);  return *this; }

      //
      // SubMatrix
      //

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) const
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1,
	    i2-i1, j2-j1,stepi(),1,RowMajor,NonConj);
      }

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1,
	    (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep,
	    jstep==1 ? RowMajor : NoMajor,NonConj);
      }

      inline ConstVectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i*stepi()+j*stepj(),size,
	    istep*stepi()+jstep,NonConj);
      }

      inline ConstUpperTriMatrixView<T> SubTriMatrix(int i1, int i2) const
      {
	TMVAssert(this->OKSubTriMatrix(i1,i2,1));
	return ConstUpperTriMatrixView<T>(cptr()+i1*(stepi()+1),i2-i1,
	    stepi(),1,D,RowMajor,NonConj);
      }

      inline ConstUpperTriMatrixView<T> SubTriMatrix(int i1, int i2,
	  int istep) const
      {
	TMVAssert(this->OKSubTriMatrix(i1,i2,istep));
	return ConstUpperTriMatrixView<T>(cptr()+i1*(stepi()+1),(i2-i1)/istep,
	    istep*stepi(),istep,D,
	    istep==1 ? RowMajor : NoMajor, NonConj);
      }

      inline ConstUpperTriMatrixView<T> OffDiag() const
      {
	TMVAssert(size() > 0);
	return ConstUpperTriMatrixView<T>(cptr()+1,size()-1,
	    stepi(),1,NonUnitDiag,RowMajor,NonConj);
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) 
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T>(ptr()+i1*stepi()+j1,
	    i2-i1, j2-j1, stepi(),1,RowMajor,NonConj FIRSTLAST );
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep)
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T>(ptr()+i1*stepi()+j1,
	    (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep,
	    jstep==1?RowMajor:NoMajor,NonConj FIRSTLAST );
      }

      inline VectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size)
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return VectorView<T>(ptr()+i*stepi()+j,size,
	    istep*stepi()+jstep,NonConj FIRSTLAST );
      }

      inline UpperTriMatrixView<T> SubTriMatrix(int i1, int i2)
      {
	TMVAssert(this->OKSubTriMatrix(i1,i2,1));
	return UpperTriMatrixView<T>(ptr()+i1*(stepi()+1),i2-i1,
	    stepi(),1,dt(),RowMajor,NonConj FIRSTLAST);
      }

      inline UpperTriMatrixView<T> SubTriMatrix(int i1, int i2, int istep) 
      {
	TMVAssert(this->OKSubTriMatrix(i1,i2,istep));
	return UpperTriMatrixView<T>(ptr()+i1*(stepi()+1),(i2-i1)/istep,
	    istep*stepi(),istep,dt(),
	    istep==1 ? RowMajor : NoMajor,NonConj FIRSTLAST);
      }

      inline UpperTriMatrixView<T> OffDiag()
      {
	TMVAssert(size() > 0);
	return UpperTriMatrixView<T>(ptr()+1,size()-1,
	    stepi(),1,NonUnitDiag,RowMajor,NonConj FIRSTLAST);
      }

      inline ConstUpperTriMatrixView<T> View() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    stepi(),1,D,RowMajor,NonConj,
	    new MetaDivider<T>(false,false,this));
      }

      inline ConstLowerTriMatrixView<T> Transpose() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    1,stepi(),D,ColMajor,NonConj,
	    new MetaDivider<T>(true,false,this));
      }

      inline ConstUpperTriMatrixView<T> Conjugate() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    stepi(),1,D,RowMajor,ConjOf(T,NonConj),
	    new MetaUpperTriDivider<T>(false,true,this));
      }

      inline ConstLowerTriMatrixView<T> Adjoint() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    1,stepi(),D,ColMajor,ConjOf(T,NonConj),
	    new MetaLowerTriDivider<T>(true,true,this));
      }

      inline ConstUpperTriMatrixView<T> QuickView() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    stepi(),1,D,RowMajor,NonConj);
      }

      inline ConstLowerTriMatrixView<T> QuickTranspose() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    1,stepi(),D,ColMajor,NonConj);
      }

      inline ConstUpperTriMatrixView<T> QuickConjugate() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    stepi(),1,D,RowMajor,ConjOf(T,NonConj));
      }

      inline ConstLowerTriMatrixView<T> QuickAdjoint() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    1,stepi(),D,ColMajor,ConjOf(T,NonConj));
      }

      inline UpperTriMatrixView<T> View() 
      { 
	return UpperTriMatrixView<T>(ptr(),size(),
	    stepi(),1,D,RowMajor,NonConj,
	    new MetaUpperTriDivider<T>(false,false,this) FIRSTLAST);
      }

      inline LowerTriMatrixView<T> Transpose() 
      { 
	return LowerTriMatrixView<T>(ptr(),size(),
	    1,stepi(),D,ColMajor,NonConj,
	    new MetaLowerTriDivider<T>(true,false,this) FIRSTLAST);
      }

      inline UpperTriMatrixView<T> Conjugate() 
      { 
	return UpperTriMatrixView<T>(ptr(),size(),
	    stepi(),1,D,RowMajor,ConjOf(T,NonConj),
	    new MetaUpperTriDivider<T>(false,true,this) FIRSTLAST);
      }

      inline LowerTriMatrixView<T> Adjoint() 
      { 
	return LowerTriMatrixView<T>(ptr(),size(),
	    1,stepi(),D,ColMajor,ConjOf(T,NonConj),
	    new MetaLowerTriDivider<T>(true,true,this) FIRSTLAST);
      }

      inline UpperTriMatrixView<T> QuickView() 
      { 
	return UpperTriMatrixView<T>(ptr(),size(),
	    stepi(),1,D,RowMajor,NonConj FIRSTLAST);
      }

      inline LowerTriMatrixView<T> QuickTranspose() 
      { 
	return LowerTriMatrixView<T>(ptr(),size(),
	    1,stepi(),D,ColMajor,NonConj FIRSTLAST);
      }

      inline UpperTriMatrixView<T> QuickConjugate() 
      { 
	return UpperTriMatrixView<T>(ptr(),size(),
	    stepi(),1,D,RowMajor,ConjOf(T,NonConj) FIRSTLAST);
      }

      inline LowerTriMatrixView<T> QuickAdjoint() 
      { 
	return LowerTriMatrixView<T>(ptr(),size(),
	    1,stepi(),D,ColMajor,ConjOf(T,NonConj) FIRSTLAST);
      }


      inline const T* cptr() const { return itsm; }
      inline T* ptr() { return itsm; }
      inline int stepi() const { return itss; }
      inline int stepj() const { return 1; }
      inline DiagType dt() const { return D; }
      inline StorageType stor() const { return RowMajor; }
      inline ConjItType ct() const { return NonConj; }
      inline bool isrm() const { return true; }
      inline bool iscm() const { return false; }
      inline bool isunit() const { return D == UnitDiag; }
      inline bool isconj() const { return false; }

      // This makes it easier for some things to compile
      // The ones that work are overridden
      template <class T2> void operator=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }
      template <class T2> void operator+=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }
      template <class T2> void operator-=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }
      template <class T2> void operator*=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }
      template <class T2> void operator/=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }
      template <class T2> void operator%=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }

    protected :

      T*const itsm;
      const size_t itss;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
#endif

      inline bool okij(size_t i, size_t j) const
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	if (isunit()) return i<j; else return i<=j;
      }

      inline T& ref(size_t i, size_t j)
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	TMVAssert(okij(i,j));
	return *(ptr() + int(i)*stepi() + j);
      }

      inline T cref(size_t i, size_t j) const 
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	if (i==j && isunit()) return T(1);
	else if (i>j) return T(0);
	else return *(cptr() + int(i)*stepi() + j);
      }

  }; // UpperTriMatrix - RowMajor

  template <class T, DiagType D> class UpperTriMatrix<T,D,ColMajor> : 
    public GenUpperTriMatrix<T>
  {

    public:

      //
      // Constructors
      //

#define NEW_SIZE(s) GenUpperTriMatrix<T>(D,ColMajor,NonConj), \
      itsm(new T[(s)*(s)]), itss(s) DEFFIRSTLAST(itsm,itsm+(s)*(s))

      explicit UpperTriMatrix(size_t _size) : NEW_SIZE(_size) {}

      UpperTriMatrix(const UpperTriMatrix<T,D,ColMajor>& rhs) : 
	NEW_SIZE(rhs.size())
      {
	memmove(itsm,rhs.itsm,size()*size()*sizeof(T));
      }

      template <DiagType D2> UpperTriMatrix(
	  const UpperTriMatrix<T,D2,ColMajor>& rhs) : NEW_SIZE(rhs.size())
      {
	memmove(itsm,rhs.cptr(),size()*size()*sizeof(T));
	if (D==NonUnitDiag && D2==UnitDiag) diag().SetAllTo(T(1));
      }

      template <DiagType D2, StorageType S2> UpperTriMatrix(
	  const UpperTriMatrix<T,D2,S2>& rhs) : NEW_SIZE(rhs.size())
      {
	if (D==D2) Copy(rhs,QuickView());
	else {
	  Copy(rhs.OffDiag(),OffDiag());
	  if (D==NonUnitDiag && D2==UnitDiag) diag().SetAllTo(T(1));
	}
      }

      template <StorageType S2> UpperTriMatrix(const Matrix<T,S2>& rhs) :
	NEW_SIZE(rhs.rowsize())
      {
	if (ColMajor==S2 && rhs.IsSquare()) 
	  memmove(itsm,rhs.cptr(),size()*size()*sizeof(T));
	else 
	  Copy(UpperTriMatrixViewOf(rhs,D),QuickView());
      }

      template <class T2> UpperTriMatrix(const GenMatrix<T2>& rhs) :
	NEW_SIZE(rhs.rowsize())
      { 
	Copy(UpperTriMatrixViewOf(rhs,D),QuickView()); 
      }

      template <class T2> UpperTriMatrix(const GenUpperTriMatrix<T2>& rhs) :
	NEW_SIZE(rhs.size())
      { 
	if (D==rhs.dt()) Copy(rhs,QuickView());
	else {
	  Copy(rhs.OffDiag(),OffDiag());
	  if (D==NonUnitDiag && rhs.dt()==UnitDiag) diag().SetAllTo(T(1));
	}
      }

      UpperTriMatrix(const UpperTriMatrixComposite<T>& mcomp) :
	NEW_SIZE(mcomp.size())
      {
	TMVAssert(!(mcomp.dt()==NonUnitDiag && D==UnitDiag));
	if (mcomp.dt()==UnitDiag && D==NonUnitDiag) {
	  mcomp.AssignTo(UpperTriMatrixViewOf(*this,UnitDiag));
	  diag().SetAllTo(T(1));
	} else 
	  mcomp.AssignTo(QuickView());
      }

#undef NEW_SIZE

      ~UpperTriMatrix() { TMVAssert(itsm); delete[] itsm; }


      //
      // Op=
      //

      inline UpperTriMatrix<T,D,ColMajor>& operator=(
	  const UpperTriMatrix<T,D,ColMajor>& m2)
      { 
	if (&m2 != this) memmove(itsm,m2.itsm,size()*size()*sizeof(T));
	return *this;
      }

      template <class T2> inline UpperTriMatrix<T,D,ColMajor>& operator=(
	  const GenUpperTriMatrix<T2>& m2)
      { 
	TMVAssert(!(m2.dt()==NonUnitDiag && D==UnitDiag));
	if (m2.dt()==UnitDiag && D==NonUnitDiag) {
	  Copy(m2.OffDiag(),OffDiag());
	  diag().SetAllTo(T(1));
	} else 
	  Copy(m2,QuickView());
	return *this;
      }

      inline UpperTriMatrix<T,D,ColMajor>& operator=(T x) 
      { 
	TMVAssert(!this->isunit() || x==T(1));
	return SetToIdentity(x); 
      }

      inline UpperTriMatrix<T,D,ColMajor>& operator=(
	  const UpperTriMatrixComposite<T>& mcomp)
      { 
	TMVAssert(size() == mcomp.size());
	TMVAssert(!(mcomp.dt()==NonUnitDiag && D==UnitDiag));
	if (mcomp.dt()==UnitDiag && D==NonUnitDiag) {
	  mcomp.AssignTo(UpperTriMatrixViewOf(*this,UnitDiag));
	  diag().SetAllTo(T(1));
	} else 
	  mcomp.AssignTo(QuickView());
	return *this;
      }

      //
      // Access
      //

      inline size_t size() const { return itss; }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert(j1==j2 || okij(i,j1));
	return ConstVectorView<T>(cptr()+i+j1*stepj(),j2-j1,stepj(),NonConj);
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert(i1==i2 || okij(i2-1,j));
	return ConstVectorView<T>(cptr()+i1+j*stepj(),i2-i1,1,NonConj);
      }

      inline ConstVectorView<T> diag(int i=0) const
      {
	TMVAssert(i<=int(size())); 
	TMVAssert(isunit() ? i>0 : i>=0);
	return ConstVectorView<T>(cptr()+i*stepj(),size()-i,(1+stepj()),
	    NonConj);
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(isunit() ? i>0 : i>=0);
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()+i);
	return ConstVectorView<T>(cptr()+(i+j1)*stepj()+j1,j2-j1,(1+stepj()),
	    NonConj);
      }

      inline VectorView<T> row(size_t i, size_t j1, size_t j2)
      { 
	TMVAssert(i<size());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=size());
	TMVAssert(j1==j2 || okij(i,j1));
	return VectorView<T>(ptr()+i+j1*stepj(),j2-j1,stepj(),NonConj
	    FIRSTLAST);
      }

      inline VectorView<T> col(size_t j, size_t i1, size_t i2)
      {
	TMVAssert(j<size());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=size());
	TMVAssert(i1==i2 || okij(i2-1,j));
	return VectorView<T>(ptr()+i1+j*stepj(),i2-i1,1,NonConj FIRSTLAST);
      }

      inline VectorView<T> diag(int i=0)
      {
	TMVAssert(i<=int(size())); 
	TMVAssert(isunit() ? i>0 : i>=0);
	return VectorView<T>(ptr()+i*stepj(),size()-i,(1+stepj()),NonConj
	    FIRSTLAST);
      }

      inline VectorView<T> diag(int i, size_t j1, size_t j2)
      {
	TMVAssert(isunit() ? i>0 : i>=0);
	TMVAssert(j1 <= j2);
	TMVAssert(j2 <= size()+i);
	return VectorView<T>(ptr()+(i+j1)*stepj()+j1,j2-j1,(1+stepj()),
	    NonConj FIRSTLAST);
      }

      inline T operator()(size_t i, size_t j) const
      { return cref(i,j); }

      inline T& operator()(size_t i, size_t j) 
      { return ref(i,j); }

      //
      // Modifying Functions
      //

      inline UpperTriMatrix<T,D,ColMajor>& Zero() { return SetAllTo(0); }
      inline UpperTriMatrix<T,D,ColMajor>& SetAllTo(T x) 
      { QuickView().SetAllTo(x); return *this; }
      inline UpperTriMatrix<T,D,ColMajor>& ConjugateSelf() 
      { QuickView().ConjugateSelf(); return *this; }
      inline UpperTriMatrix<T,D,ColMajor>& SetToIdentity(T x=T(1)) 
      { QuickView().SetToIdentity(x);  return *this; }

      //
      // SubMatrix
      //

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) const
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return ConstMatrixView<T>(cptr()+i1+j1*stepj(),
	    i2-i1, j2-j1,1,stepj(),ColMajor,NonConj);
      }

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, istep, jstep*stepj(),
	    istep==1 ? ColMajor : NoMajor,NonConj);
      }

      inline ConstVectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i+j*stepj(),size,
	    istep+jstep*stepj(),NonConj);
      }

      inline ConstUpperTriMatrixView<T> SubTriMatrix(int i1, int i2) const
      {
	TMVAssert(this->OKSubTriMatrix(i1,i2,1));
	return ConstUpperTriMatrixView<T>(cptr()+i1*(1+stepj()),i2-i1,
	    1,stepj(),D,ColMajor,NonConj);
      }

      inline ConstUpperTriMatrixView<T> SubTriMatrix(int i1, int i2,
	  int istep) const
      {
	TMVAssert(this->OKSubTriMatrix(i1,i2,istep));
	return ConstUpperTriMatrixView<T>(cptr()+i1*(1+stepj()),(i2-i1)/istep,
	    istep,istep*stepj(),D,istep==1 ? ColMajor : NoMajor, NonConj);
      }

      inline ConstUpperTriMatrixView<T> OffDiag() const
      {
	TMVAssert(size() > 0);
	return ConstUpperTriMatrixView<T>(cptr()+stepj(),size()-1,
	    1,stepj(),NonUnitDiag,ColMajor,NonConj);
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) 
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T>(ptr()+i1+j1*stepj(),
	    i2-i1, j2-j1, 1,stepj(),ColMajor,NonConj FIRSTLAST );
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep)
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T>(ptr()+i1+j1*stepj(),
	    (i2-i1)/istep, (j2-j1)/jstep, istep, jstep*stepj(),
	    istep==1?ColMajor:NoMajor,NonConj FIRSTLAST );
      }

      inline VectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size)
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return VectorView<T>(ptr()+i+j*stepj(),size,
	    istep+jstep*stepj(),NonConj FIRSTLAST );
      }

      inline UpperTriMatrixView<T> SubTriMatrix(int i1, int i2)
      {
	TMVAssert(this->OKSubTriMatrix(i1,i2,1));
	return UpperTriMatrixView<T>(ptr()+i1*(1+stepj()),i2-i1,
	    1,stepj(),dt(),ColMajor,NonConj FIRSTLAST);
      }

      inline UpperTriMatrixView<T> SubTriMatrix(int i1, int i2, int istep) 
      {
	TMVAssert(this->OKSubTriMatrix(i1,i2,istep));
	return UpperTriMatrixView<T>(ptr()+i1*(1+stepj()),(i2-i1)/istep,
	    istep,istep*stepj(),dt(),
	    istep==1 ? ColMajor : NoMajor,NonConj FIRSTLAST);
      }

      inline UpperTriMatrixView<T> OffDiag()
      {
	TMVAssert(size() > 0);
	return UpperTriMatrixView<T>(ptr()+stepj(),size()-1,
	    1,stepj(),NonUnitDiag,ColMajor,NonConj FIRSTLAST);
      }

      inline ConstUpperTriMatrixView<T> View() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    1,stepj(),D,ColMajor,NonConj,
	    new MetaUpperTriDivider<T>(false,false,this));
      }

      inline ConstLowerTriMatrixView<T> Transpose() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    stepj(),1,D,RowMajor,NonConj,
	    new MetaLowerTriDivider<T>(true,false,this));
      }

      inline ConstUpperTriMatrixView<T> Conjugate() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    1,stepj(),D,ColMajor,ConjOf(T,NonConj),
	    new MetaUpperTriDivider<T>(false,true,this));
      }

      inline ConstLowerTriMatrixView<T> Adjoint() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    stepj(),1,D,RowMajor,ConjOf(T,NonConj),
	    new MetaLowerTriDivider<T>(true,true,this));
      }

      inline ConstUpperTriMatrixView<T> QuickView() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    1,stepj(),D,ColMajor,NonConj);
      }

      inline ConstLowerTriMatrixView<T> QuickTranspose() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    stepj(),1,D,RowMajor,NonConj);
      }

      inline ConstUpperTriMatrixView<T> QuickConjugate() const
      { 
	return ConstUpperTriMatrixView<T>(cptr(),size(),
	    1,stepj(),D,ColMajor,ConjOf(T,NonConj));
      }

      inline ConstLowerTriMatrixView<T> QuickAdjoint() const
      { 
	return ConstLowerTriMatrixView<T>(cptr(),size(),
	    stepj(),1,D,RowMajor,ConjOf(T,NonConj));
      }

      inline UpperTriMatrixView<T> View() 
      { 
	return UpperTriMatrixView<T>(ptr(),size(),
	    1,stepj(),D,ColMajor,NonConj,
	    new MetaUpperTriDivider<T>(false,false,this) FIRSTLAST);
      }

      inline LowerTriMatrixView<T> Transpose() 
      { 
	return LowerTriMatrixView<T>(ptr(),size(),
	    stepj(),1,D,RowMajor,NonConj,
	    new MetaLowerTriDivider<T>(true,false,this) FIRSTLAST);
      }

      inline UpperTriMatrixView<T> Conjugate() 
      { 
	return UpperTriMatrixView<T>(ptr(),size(),
	    1,stepj(),D,ColMajor,ConjOf(T,NonConj),
	    new MetaUpperTriDivider<T>(false,true,this) FIRSTLAST);
      }

      inline LowerTriMatrixView<T> Adjoint() 
      { 
	return LowerTriMatrixView<T>(ptr(),size(),
	    stepj(),1,D,RowMajor,ConjOf(T,NonConj),
	    new MetaLowerTriDivider<T>(true,true,this) FIRSTLAST);
      }

      inline UpperTriMatrixView<T> QuickView() 
      { 
	return UpperTriMatrixView<T>(ptr(),size(),
	    1,stepj(),D,ColMajor,NonConj FIRSTLAST);
      }

      inline LowerTriMatrixView<T> QuickTranspose() 
      { 
	return LowerTriMatrixView<T>(ptr(),size(),
	    stepj(),1,D,RowMajor,NonConj FIRSTLAST);
      }

      inline UpperTriMatrixView<T> QuickConjugate() 
      { 
	return UpperTriMatrixView<T>(ptr(),size(),
	    1,stepj(),D,ColMajor,ConjOf(T,NonConj) FIRSTLAST);
      }

      inline LowerTriMatrixView<T> QuickAdjoint() 
      { 
	return LowerTriMatrixView<T>(ptr(),size(),
	    stepj(),1,D,RowMajor,ConjOf(T,NonConj) FIRSTLAST);
      }

      inline const T* cptr() const { return itsm; }
      inline T* ptr() { return itsm; }
      inline int stepi() const { return 1; }
      inline int stepj() const { return itss; }
      inline DiagType dt() const { return D; }
      inline StorageType stor() const { return ColMajor; }
      inline ConjItType ct() const { return NonConj; }
      inline bool isrm() const { return false; }
      inline bool iscm() const { return true; }
      inline bool isunit() const { return D == UnitDiag; }
      inline bool isconj() const { return false; }

      // This makes it easier for some things to compile
      // The ones that work are overridden
      template <class T2> void operator=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }
      template <class T2> void operator+=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }
      template <class T2> void operator-=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }
      template <class T2> void operator*=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }
      template <class T2> void operator/=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }
      template <class T2> void operator%=(
	  const BaseMatrix<T2>&) { TMVAssert(false); }

    protected :

      T*const itsm;
      const size_t itss;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
#endif

      inline bool okij(size_t i, size_t j) const
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	if (isunit()) return i<j; else return i<=j;
      }

      inline T& ref(size_t i, size_t j)
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	TMVAssert(okij(i,j));
	return *(ptr() + i + int(j)*stepj());
      }

      inline T cref(size_t i, size_t j) const 
      {
	TMVAssert(i < size());
	TMVAssert(j < size());
	if (i==j && isunit()) return T(1);
	else if (i>j) return T(0);
	else return *(cptr() + i + int(j)*stepj());
      }

  }; // UpperTriMatrix - ColMajor

//---------------------------------------------------------------------------

  //
  // Special Creators: 
  //   UpperTriMatrixViewOf(m)
  //   LowerTriMatrixViewOf(m)
  //   UnitTriMatrixViewOf(t)
  //

  template <class T> inline ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
      const GenMatrix<T>& m, DiagType dt=NonUnitDiag)
  {
    TMVAssert(m.colsize()>=m.rowsize());
    return ConstUpperTriMatrixView<T>(m.Rows(0,m.rowsize()),dt); 
  }

  template <class T> inline UpperTriMatrixView<T> UpperTriMatrixViewOf(
      const MatrixView<T>& m, DiagType dt=NonUnitDiag)
  { 
    TMVAssert(m.colsize()>=m.rowsize());
    return UpperTriMatrixView<T>(m.Rows(0,m.rowsize()),dt); 
  }

  template <class T, StorageType S>
    inline UpperTriMatrixView<T> UpperTriMatrixViewOf(
	Matrix<T,S>& m, DiagType dt=NonUnitDiag)
    {
      TMVAssert(m.colsize()>=m.rowsize());
      return UpperTriMatrixView<T>(m.Rows(0,m.rowsize()),dt); 
    }

  template <class T> inline ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
      const GenUpperTriMatrix<T>& m, DiagType dt)
  {
    TMVAssert(m.colsize()>=m.rowsize());
    TMVAssert(!(m.dt()==UnitDiag && dt==NonUnitDiag));
    return ConstUpperTriMatrixView<T>(m,dt); 
  }

  template <class T> inline UpperTriMatrixView<T> UpperTriMatrixViewOf(
      const UpperTriMatrixView<T>& m, DiagType dt)
  { 
    TMVAssert(m.colsize()>=m.rowsize());
    TMVAssert(!(m.dt()==UnitDiag && dt==NonUnitDiag));
    return UpperTriMatrixView<T>(m,dt); 
  }

  template <class T, DiagType D, StorageType S>
    inline UpperTriMatrixView<T> UpperTriMatrixViewOf(
	UpperTriMatrix<T,D,S>& m, DiagType dt)
    {
      TMVAssert(m.colsize()>=m.rowsize());
      TMVAssert(!(D==UnitDiag && dt==NonUnitDiag));
      return UpperTriMatrixView<T>(m.QuickView(),dt); 
    }

  template <class T> inline ConstUpperTriMatrixView<T> UpperTriMatrixViewOf(
      const T* m, size_t size, DiagType dt, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstUpperTriMatrixView<T>(m,size,size,1,size+1,
	  dt,RowMajor,NonConj);
    else
      return ConstUpperTriMatrixView<T>(m,size,1,size,size+1,
	  dt,ColMajor,NonConj);
  }

  template <class T> inline UpperTriMatrixView<T> UpperTriMatrixViewOf(
      T* m, size_t size, DiagType dt, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return UpperTriMatrixView<T>(m,size,size,1,size+1,
	  dt,RowMajor,NonConj FIRSTLAST1(m,m+size*size));
    else
      return UpperTriMatrixView<T>(m,size,1,size,size+1,
	  dt,ColMajor,NonConj FIRSTLAST1(m,m+size*size));
  }

  template <class T> inline ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
      const GenMatrix<T>& m, DiagType dt=NonUnitDiag)
  {
    TMVAssert(m.colsize()<=m.rowsize());
    return ConstLowerTriMatrixView<T>(m.Cols(0,m.colsize()),dt); 
  }

  template <class T> inline LowerTriMatrixView<T> LowerTriMatrixViewOf(
      const MatrixView<T>& m, DiagType dt=NonUnitDiag)
  { 
    TMVAssert(m.colsize()<=m.rowsize());
    return LowerTriMatrixView<T>(m.Cols(0,m.colsize()),dt); 
  }

  template <class T,StorageType S> 
    inline LowerTriMatrixView<T> LowerTriMatrixViewOf(
	Matrix<T,S>& m, DiagType dt=NonUnitDiag)
    {
      TMVAssert(m.colsize()<=m.rowsize());
      return LowerTriMatrixView<T>(m.Cols(0,m.colsize()),dt); 
    }

  template <class T> inline ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
      const GenLowerTriMatrix<T>& m, DiagType dt)
  {
    TMVAssert(m.colsize()<=m.rowsize());
    TMVAssert(!(m.dt()==UnitDiag && dt==NonUnitDiag));
    return ConstLowerTriMatrixView<T>(m,dt); 
  }

  template <class T> inline LowerTriMatrixView<T> LowerTriMatrixViewOf(
      const LowerTriMatrixView<T>& m, DiagType dt)
  { 
    TMVAssert(m.colsize()<=m.rowsize());
    TMVAssert(!(m.dt()==UnitDiag && dt==NonUnitDiag));
    return LowerTriMatrixView<T>(m,dt); 
  }

  template <class T, DiagType D, StorageType S> 
    inline LowerTriMatrixView<T> LowerTriMatrixViewOf(
	LowerTriMatrix<T,D,S>& m, DiagType dt)
    {
      TMVAssert(m.colsize()<=m.rowsize());
      TMVAssert(!(D==UnitDiag && dt==NonUnitDiag));
      return LowerTriMatrixView<T>(m.QuickView(),dt); 
    }

  template <class T> inline ConstLowerTriMatrixView<T> LowerTriMatrixViewOf(
      const T* m, size_t size, DiagType dt, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstLowerTriMatrixView<T>(m,size,size,1,size+1,
	  dt,RowMajor,NonConj);
    else
      return ConstLowerTriMatrixView<T>(m,size,1,size,size+1,
	  dt,ColMajor,NonConj);
  }

  template <class T> inline LowerTriMatrixView<T> LowerTriMatrixViewOf(
      T* m, size_t size, DiagType dt, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return LowerTriMatrixView<T>(m,size,size,1,size+1,
	  dt,RowMajor,NonConj FIRSTLAST1(m,m+size*size));
    else
      return LowerTriMatrixView<T>(m,size,1,size,size+1,
	  dt,ColMajor,NonConj FIRSTLAST1(m,m+size*size));
  }

  //
  // Copy
  //

  template <class T1, class T2> inline void Copy(
      const GenUpperTriMatrix<T1>& m1, const UpperTriMatrixView<T2>& m2)
  {
    TMVAssert(m1.size() == m2.size());
    TMVAssert(m1.dt() == m2.dt());
    const size_t N = m1.size();

    if (m1.dt() == UnitDiag)
      if (m1.isrm() && m2.isrm()) 
	for(size_t i=0;i<N;++i) m2.row(i,i+1,N) = m1.row(i,i+1,N);
      else 
	for(size_t j=0;j<N;++j) m2.col(j,0,j) = m1.col(j,0,j);
    else
      if (m1.iscm() && m2.iscm()) 
	for(size_t j=0;j<N;++j) m2.col(j,0,j+1) = m1.col(j,0,j+1);
      else 
	for(size_t i=0;i<N;++i) m2.row(i,i,N) = m1.row(i,i,N);
  }

  //
  // Swap Matrices
  //

  template <class T> void Swap(
      const UpperTriMatrixView<T>& m1, const UpperTriMatrixView<T>& m2);
  template <class T, DiagType D, StorageType S> inline void Swap(
      const UpperTriMatrixView<T>& m1, UpperTriMatrix<T,D,S>& m2)
  { Swap(m1,m2.QuickView()); }
  template <class T, DiagType D, StorageType S> inline void Swap(
      UpperTriMatrix<T,D,S>& m1, const UpperTriMatrixView<T>& m2)
  { Swap(m1.QuickView(),m2); }
  template <class T, DiagType D, StorageType S1, StorageType S2> 
    inline void Swap(UpperTriMatrix<T,D,S1>& m1, UpperTriMatrix<T,D,S2>& m2)
    { Swap(m1.QuickView(),m2.QuickView()); }

  template <class T> inline void Swap(
      const LowerTriMatrixView<T>& m1, const LowerTriMatrixView<T>& m2)
  { Swap(m1.QuickTranspose(),m2.QuickTranspose()); }
  template <class T, DiagType D, StorageType S> inline void Swap(
      const LowerTriMatrixView<T>& m1, LowerTriMatrix<T,D,S>& m2)
  { Swap(m1.QuickTranspose(),m2.QuickTranspose()); }
  template <class T, DiagType D, StorageType S> inline void Swap(
      LowerTriMatrix<T,D,S>& m1, const LowerTriMatrixView<T>& m2)
  { Swap(m1.QuickTranspose(),m2.QuickTranspose()); }
  template <class T, DiagType D, StorageType S1, StorageType S2> 
    inline void Swap(LowerTriMatrix<T,D,S1>& m1, LowerTriMatrix<T,D,S2>& m2)
    { Swap(m1.QuickTranspose(),m2.QuickTranspose()); }

  //
  // Functions of Matrices:
  //

  template <class T> inline T Det(const GenUpperTriMatrix<T>& m)
  { return m.Det(); }
  template <class T> inline T Det(const GenLowerTriMatrix<T>& m)
  { return m.Det(); }

  template <class T> inline T Trace(const GenUpperTriMatrix<T>& m)
  { return m.Trace(); }
  template <class T> inline T Trace(const GenLowerTriMatrix<T>& m)
  { return m.Trace(); }

  template <class T> inline RealType(T) Norm(const GenUpperTriMatrix<T>& m)
  { return m.Norm(); }
  template <class T> inline RealType(T) Norm(const GenLowerTriMatrix<T>& m)
  { return m.Norm(); }

  template <class T> inline RealType(T) NormF(const GenUpperTriMatrix<T>& m)
  { return m.NormF(); }
  template <class T> inline RealType(T) NormF(const GenLowerTriMatrix<T>& m)
  { return m.NormF(); }

  template <class T> inline RealType(T) Norm2(const GenUpperTriMatrix<T>& m)
  { return m.Norm2(); }
  template <class T> inline RealType(T) Norm2(const GenLowerTriMatrix<T>& m)
  { return m.Norm2(); }

  template <class T> inline RealType(T) NormInf(const GenUpperTriMatrix<T>& m)
  { return m.NormInf(); }
  template <class T> inline RealType(T) NormInf(const GenLowerTriMatrix<T>& m)
  { return m.NormInf(); }

  template <class T> inline UpperTriMatrix<T,NonUnitDiag,ColMajor> Inverse(
      const GenUpperTriMatrix<T>& m)
  { return m.TInverse(); }
  template <class T> inline LowerTriMatrix<T,NonUnitDiag,ColMajor> Inverse(
      const GenLowerTriMatrix<T>& m)
  { return m.TInverse(); }

  template <class T> inline Matrix<T,ColMajor> InverseATA(
      const GenUpperTriMatrix<T>& m)
  { return m.InverseATA(); }
  template <class T> inline Matrix<T,ColMajor> InverseATA(
      const GenLowerTriMatrix<T>& m)
  { return m.InverseATA(); }

  template <class T> inline ConstLowerTriMatrixView<T> Transpose(
      const GenUpperTriMatrix<T>& m)
  { return m.Transpose(); }
  template <class T> inline ConstUpperTriMatrixView<T> Transpose(
      const GenLowerTriMatrix<T>& m)
  { return m.Transpose(); }

  template <class T> inline LowerTriMatrixView<T> Transpose(
      const UpperTriMatrixView<T>& m)
  { return m.Transpose(); }
  template <class T> inline UpperTriMatrixView<T> Transpose(
      const LowerTriMatrixView<T>& m)
  { return m.Transpose(); }

  template <class T, DiagType D, StorageType S> 
    inline LowerTriMatrixView<T> Transpose(UpperTriMatrix<T,D,S>& m)
    { return m.Transpose(); }
  template <class T, DiagType D, StorageType S> 
    inline UpperTriMatrixView<T> Transpose(LowerTriMatrix<T,D,S>& m)
    { return m.Transpose(); }

  template <class T> inline ConstUpperTriMatrixView<T> Conjugate(
      const GenUpperTriMatrix<T>& m)
  { return m.Conjugate(); }
  template <class T> inline ConstLowerTriMatrixView<T> Conjugate(
      const GenLowerTriMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T> inline UpperTriMatrixView<T> Conjugate(
      const UpperTriMatrixView<T>& m)
  { return m.Conjugate(); }
  template <class T> inline LowerTriMatrixView<T> Conjugate(
      const LowerTriMatrixView<T>& m)
  { return m.Conjugate(); }

  template <class T, DiagType D, StorageType S> 
    inline UpperTriMatrixView<T> Conjugate(UpperTriMatrix<T,D,S>& m)
  { return m.Conjugate(); }
  template <class T, DiagType D, StorageType S> 
    inline LowerTriMatrixView<T> Conjugate(LowerTriMatrix<T,D,S>& m)
  { return m.Conjugate(); }

  template <class T> inline ConstLowerTriMatrixView<T> Adjoint(
      const GenUpperTriMatrix<T>& m)
  { return m.Adjoint(); }
  template <class T> inline ConstUpperTriMatrixView<T> Adjoint(
      const GenLowerTriMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T> inline LowerTriMatrixView<T> Adjoint(
      const UpperTriMatrixView<T>& m)
  { return m.Adjoint(); }
  template <class T> inline UpperTriMatrixView<T> Adjoint(
      const LowerTriMatrixView<T>& m)
  { return m.Adjoint(); }

  template <class T, DiagType D, StorageType S> 
    inline LowerTriMatrixView<T> Adjoint(UpperTriMatrix<T,D,S>& m)
  { return m.Adjoint(); }
  template <class T, DiagType D, StorageType S> 
    inline UpperTriMatrixView<T> Adjoint(LowerTriMatrix<T,D,S>& m)
  { return m.Adjoint(); }


  //
  // TriMatrix ==, != TriMatrix
  //

  template <class T> bool operator==(
      const GenUpperTriMatrix<T>& m1, const GenUpperTriMatrix<T>& m2);
  template <class T> inline bool operator==(
      const GenLowerTriMatrix<T>& m1, const GenLowerTriMatrix<T>& m2)
  { return m1.QuickTranspose() == m2.QuickTranspose(); }
  template <class T> bool operator!=(
      const GenUpperTriMatrix<T>& m1, const GenUpperTriMatrix<T>& m2)
  { return !(m1 == m2); }
  template <class T> inline bool operator!=(
      const GenLowerTriMatrix<T>& m1, const GenLowerTriMatrix<T>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T, DiagType D, StorageType S> istream& operator>>(
      istream& is, UpperTriMatrix<T,D,S>* m);
  template <class T, DiagType D, StorageType S> istream& operator>>(
      istream& is, LowerTriMatrix<T,D,S>* m);

  template <class T> istream& operator>>(
      istream& is, const UpperTriMatrixView<T>& m);
  template <class T> istream& operator>>(
      istream& is, const LowerTriMatrixView<T>& m);

  template <class T, DiagType D, StorageType S> istream& operator>>(
      istream& is, UpperTriMatrix<T,D,S>& m)
  { return is>>m.QuickView(); }
  template <class T, DiagType D, StorageType S> istream& operator>>(
      istream& is, LowerTriMatrix<T,D,S>& m)
  { return is>>m.QuickView(); }

  template <class T, DiagType D, StorageType S> inline std::string Type(
      const UpperTriMatrix<T,D,S>& m)
  { 
    return std::string("UpperTriMatrix<")+Type(T())+","
      + Text(D)+","+Text(S)+">"; 
  }
  template <class T, DiagType D, StorageType S> inline std::string Type(
      const LowerTriMatrix<T,D,S>& m)
  { 
    return std::string("LowerTriMatrix<")+Type(T())+","
      + Text(D)+","+Text(S)+">"; 
  }

  template <class T> inline std::string Type(const GenUpperTriMatrix<T>& m)
  { 
    return std::string("GenUpperTriMatrix<")+Type(T())+","
      + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+">"; 
  }
  template <class T> inline std::string Type(const GenLowerTriMatrix<T>& m)
  { 
    return std::string("GenLowerTriMatrix<")+Type(T())+","
      + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+">"; 
  }

  template <class T> inline std::string Type(
      const ConstUpperTriMatrixView<T>& m)
  {
    return std::string("ConstUpperTriMatrixView<")+Type(T())+","
      + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+">"; 
  }
  template <class T> inline std::string Type(
      const ConstLowerTriMatrixView<T>& m)
  {
    return std::string("ConstLowerTriMatrixView<")+Type(T())+","
      + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+">"; 
  }

  template <class T> inline std::string Type(const UpperTriMatrixView<T>& m)
  {
    return std::string("UpperTriMatrixView<")+Type(T())+","
      + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+">"; 
  }
  template <class T> inline std::string Type(const LowerTriMatrixView<T>& m)
  {
    return std::string("LowerTriMatrixView<")+Type(T())+","
      + Text(m.dt())+","+Text(m.ct())+","+Text(m.stor())+">"; 
  }

}; // namespace tmv

#endif
