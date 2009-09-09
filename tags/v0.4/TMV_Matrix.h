//---------------------------------------------------------------------------
//
// This file defines the TMV Matrix class.
//
// The Matrix class and all associated functions are contained
// in the namespace tmv.  Alse, the Matrix class is a template, 
// so for a Matrix of long doubles, one would write 
// tmv::Matrix<long double>.  
//
// There is also an optional second template parameter which is 
// either RowMajor or ColMajor, referring to the storage pattern of 
// the matrix elements.
//
// Constructors:
//
//    Matrix<T>(size_t colsize, size_t rowsize)
//        Makes a Matrix with column size = colsize and row size = rowsize
//        with _uninitialized_ values
//
//    Matrix<T>(size_t colsize, size_t rowsize, T x)
//        Makes a Matrix of size n with all values = x
//
//    Matrix<T>(const vector<vector<T> >& m)
//        Makes a Matrix with a_ij = m[i][j]
//
//    Matrix<T,stor>(size_t colsize, size_t rowsize, const valarray<T>& m)
//    Matrix<T,stor>(size_t colsize, size_t rowsize, const T* m)
//    Matrix<T,stor>(size_t colsize, size_t rowsize, const vector<T>& m)
//        Make a Matrix which copies the elements of m.
//        If stor is tmv::RowMajor then the elements are taken in row major
//        order (m00,m01,..m0n,m10,m11...).  If stor is tmv::ColMajor
//        then the elements are taken in column major order.
//        If stor is omitted, then tmv::RowMajor is assumed.
//        (stor is also an optional last parameter on the other above 
//        constructors as well.)
//
//
// Special Constructors
//
//    Eye<T,stor>(size_t n)
//        Returns the nxn identity matrix.
//        stor = RowMajor by default.
//
//    Eye<T,stor>(size_t n, T x)
//        Returns x times the nxn identity matrix
//        stor = RowMajor by default.
//
//    Diag<T,stor>(const Vector& v)  
//        Returns a diagonal Matrix with Vector v along the diagonal
//        stor = RowMajor by default.
//
//    RowVector(const Vector& v)
//        Return a 1xn (RowMajor) Matrix whose only row is v
//
//    ColVector(const Vector& v)
//        Return an nx1 (ColMajor) Matrix whose only column is v
//
//    OuterProduct<T,stor>(const Vector& v1, const Vector& v2)
//        Returns the outer product of v1 and v2.
//        ie. m_ij = v1[i] v2[j]
//        stor = RowMajor by default.
//
//    MatrixView RowVectorViewOf(Vector& v)
//    MatrixView RowVectorViewOf(const VectorView& v)
//    ConstMatrixView RowVectorViewOf(const Vector& v)
//    ConstMatrixView RowVectorViewOf(const ConstVectorView& v)
//        Returns a 1xn MatrixView with v in the (only) row.
//        Unlike creating a Matrix with RowVector(v), this uses the same
//        data storage as the actual Vector v.
//        For a const Vector or a ConstVectorView, this returns a 
//        ConstMatrixView.
//
//    MatrixView ColVectorViewOf(Vector& v)
//    MatrixView ColVectorViewOf(const VectorView& v)
//    ConstMatrixView ColVectorViewOf(const Vector& v)
//    ConstMatrixView ColVectorViewOf(const ConstVectorView& v)
//        Returns an nx1 MatrixView with v in the (only) column.
//
//    MatrixView MatrixViewOf(T* m, size_t colsize, size_t rowsize,
//            StorageType stor)
//    ConstMatrixView MatrixViewOf(const T* m, size_t colsize, size_t rowsize,
//            StorageType stor)
//        Returns a MatrixView of the elements in m, using the actual
//        elements m for the storage.  This is essentially the same as the 
//        constructor with (const T*m), except that the data isn't duplicated.
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//        Return the dimensions of the Matrix
//
//    T& operator()(size_t i, size_t j)
//    T operator()(size_t i, size_t j) const
//        Return the (i,j) element of the Matrix
//
//    VectorView& row(size_t i, size_t j1, size_t j2)
//    ConstVectorView& row(size_t i) const
//        Return the ith row of the Matrix as a Vector
//        If j1,j2 are given, it returns the SubVector from j1 to j2 
//        (not including j2) within the row.
//
//    VectorView& col(size_t j, size_t i1, size_t i2)
//    ConstVectorView& col(size_t j) const
//        Return the jth column of the Matrix as a Vector
//        If i1,i2 are given, it returns the SubVector from i1 to i2 
//        (not including i2) within the column.
//
//    VectorView& diag()
//    ConstVectorView& diag() const
//        Return the diagonal of the Matrix as a Vector
//
//    VectorView& diag(int i, size_t j1, size_t j2)
//    ConstVectorView& diag(i, j1, j2) const
//        Return the super- or sub-diagonal i
//        If i > 0 return the super diagonal starting at m_0i
//        If i < 0 return the sub diagonal starting at m_|i|0
//        If j1,j2 are given, it returns the diagonal SubVector 
//        either from m_j1,i+j1 to m_j2,i+j2 (for i>0) 
//        or from m_|i|+j1,j1 to m_|i|+j2,j2 (for i<0)
//
// Modifying Functions
//
//    Matrix& Zero()
//        Sets all elements to 0
//
//    Matrix& SetAllTo(T x)
//        Sets all elements to x
//
//    Matrix<T>& TransposeSelf() 
//        Transposes the elements of a square Matrix or SubMatrix
//
//    Matrix& ConjugateSelf()
//        Sets all elements to its conjugate
//
//    Matrix& SetToIdentity(x = 1)
//        Set to Identity Matrix, or 
//        with a parameter, set to x times Identity Matrix
//
//    Matrix& SwapRows(i1, i2)
//        Swap two rows
//
//    Matrix& SwapCols(j1, j2)
//        Swap two columns
//
//    void Swap(Matrix& m1, Matrix& m2)
//        Swap the values of two Matrices
//        The Matrices must be the same size
//
// MatrixViews:
//
//    As with VectorViews, we have both constant and mutable views of Matrices.
//    A ConstMatrixView object can only view the underlying Matrix object
//    whereas a MatrixView object can modify it as well.
//    For the below routines, a ConstMatrixView is returned from 
//    ConstMatrixViews and const Matrix objects.
//    A MatrixView is returned from MatrixViews and (non-const) Matrix objects.
//
//    MatrixView SubMatrix(int i1, int i2, int j1, int j2,
//            int istep=1, int jstep=1)
//        This member function will return a submatrix using rows i1 to i2
//        and columns j1 to j2 (not inclusive of i2,j2) which refers
//        to the same physical elements as the original.
//        Thus, to make the matrix:
//           ( 0 0 1 0 )
//           ( 0 0 0 1 )
//           ( 2 2 0 0 )
//           ( 2 2 0 0 )
//        you could write:
//           Matrix<int> A(4,4,0);
//           A.SubMatrix(0,1,2,3) = Eye<int>(2);
//           A.SubMatrix(2,3,0,1).SetAllTo(2);
//        The substep values allow you to space the elements of 
//        the submatrix at steps larger than 1.
//        eg. To make an 8x8 checkerboard of 1's and 0's, you could write:
//           Matrix<int> A(8,8,0);
//           A.SubMatrix(0,3,1,4,2,2) = 1;
//           A.SubMatrix(1,4,0,3,2,2) = 1;
//
//    VectorView SubVector(int i, int j, int istep, int jstep, int size)
//        Returns a SubVector which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.
//        For example, the cross-diagonal from the lower left to the upper
//        right of a 6x6 matrix could be accessed using:
//        m.SubVector(5,0,-1,1,6)
//
//    MatrixView ColPair(int j1, int j2)
//        This returns an Mx2 submatrix which consists of the 
//        columns j1 and j2.  This is useful for multiplying two 
//        (not necessarily adjacent) columns of a matrix by a 2x2 matrix.
//
//    MatrixView RowPair(int i1, int i2)
//        Same as above, but two rows.
//
//    MatrixView Cols(int j1, int j2)
//        This is shorthand for SubMatrix(0,m.colsize(),j1,j2)
//
//    MatrixView Rows(int i1, int i2)
//        This is shorthand for SubMatrix(i1,i2,0,m.rowsize())
//
//    MatrixView Real()
//    MatrixView Imag()
//        Returns a view to the real/imag elements of a complex Matrix.
//
//    MatrixView View()
//        Returns a view of a matrix.
//
//    MatrixView QuickView()
//        Same thing, but doesn't make a MetaDivider.  This is quicker
//        if you are not going to divide by the View, but the MetaDivider
//        can speed up routines that divide by Views, Transposes, etc.
//        There is a "Quick" version of each of the below: 
//        QuickConjugate, QuickTranspose, and QuickAdjoint which have
//        the same tradeoffs.
//
//    MatrixView Conjugate()
//        Returns a view to the conjugate of a Matrix.
//
//    MatrixView Transpose()
//        Returns a view to the transpose of a Matrix.
//
//    MatrixView Adjoint()
//        Returns a view to the adjoint (conjugate transpose) of a Matrix.
//        Note: Some people define the adjoint as the cofactor matrix.
//              This is not the same as our definition of the Adjoint.
//
// Functions of Matrixs:
//        (These are all both member functions and functions of a Matrix,
//         so Norm(m) and m.Norm() for example are equivalent.)
//
//    Det(m)
//        Returns the determinant of a Matrix.
//        Note: If the matrix is not square, the determinant is not
//              well defined.  The returned value is such that
//              conj(det) * det = Det(Adjoint(m) * m)
//              So for real nonsquare matrices, the sign is arbitrary,
//              and for complex nonsquare matrices, it is multiplied
//              by an arbitrary phase.
//
//    Trace(m)
//        Returns the trace of a Matrix.
//        = sum_i ( a_ii )
//
//    Norm(m) or NormF(m)
//        Return the Frobenius norm of a Matrix.
//        = sqrt( sum_ij |a_ij|^2 )
//
//    NormSq(m)
//        Returns the square of Norm().
//
//    Norm1(m) 
//        Returns the 1-norm of a Matrix.
//        = max_j (sum_i |a_ij|)
//
//    Norm2(m) 
//        Returns the 2-norm of a Matrix.
//        = sqrt( Max Eigenvalue of (A.Adjoint * A) )
//        = Maximum singular value
//        Note: This norm is costly to calculate if one is not 
//              otherwise doing a singular value decomposition
//              of the Matrix.
//
//    NormInf(m) 
//        Returns the infinity-norm of a Matrix.
//        = max_i (sum_j |a_ij|)
//
//    Inverse(m)
//        Returns the inverse of m if it exists.
//        If m is singular and square, and LU is set for dividing
//          (LU is default for square matrices)
//          then a run-time error will occur.
//        If m is singular or not square and SV is set 
//          then the returned matrix is the pseudo-inverse which satisfies:
//          MXM = M
//          XMX = X
//          (MX)T = MX
//          (XM)T = XM
//          (If QR or QRP is set, all but the last of these will be true.)
//
//    InverseATA(m)
//        Returns the Inverse of (A.Adjoint * A) for matrix m = A
//        If Ax=b is solved for x, then (AT A)^-1 is the 
//        covariance matrix of the least-squares solution x
//        This is only efficient to calculate if m is using
//        SVD for its divisions and SVD has been set (explicitly through
//        SetDiv() or implicitly through matrix division or Inverse()).
//
//    m.NewTranspose()
//    m.NewConjugate()
//    m.NewAdjoint()
//    m.NewView()
//    m.NewInverse()
//    m.NewCopy()
//        These all return pointers to new BaseMatrix's equal to the 
//        Transpose, Conjugate, Adjoint, View, Inverse, or itself respectively.
//        Inverse and Copy allocate new memory, the others are view.
//
// Operators:
//        Here we use m for a Matrix, v for a Vector and x for a Scalar.
//
//        You can also mix real and complex Vectors of the same
//        underlying type.  eg. Matrix<double> and Matrix<complex<double> >
//
//    -m
//
//    m = m
//
//    m += x
//    m += m
//    m + x
//    x + m
//    m + m
//
//    m -= x
//    m -= m
//    m - x
//    x - m
//    m - m
//
//    m *= x
//    m *= m
//    m * x
//    x * m
//    m * v
//    v * m
//    m * m
//
//    m /= x
//    v /= m
//    v %= m
//    m /= m
//    m %= m
//    m / x
//    x / m
//    v / m
//    v % m
//    m / m
//    m % m
//
//    m == m
//    m != m
//
//       Most of these behave in the logical way for dealing with matrices.
//       Three comments about behavior that might not be obvious:
//
//       1) Sometimes x should be thought of as x * I
//          eg. m+x means m+xI,  x/m means xI/m
//
//       2) Vectors are either row or column Vectors as appropriate.
//          eg. For m*v, v is a column Vector, but for v*m, v is a row Vector
//
//       3) Division by a matrix can be from either the left or the right.
//          ie. v/m can mean either m^-1 v or v m^-1
//          The first case is the solution of mx=v, the second is of xm = v.
//          Since the first case is the way one generally poses a problem
//          for solving a set of equations, we take v/m to be left-division.
//          If you want right-division (v m^-1), then we supply the % operator
//          to do so.  
//          ie. v%m means v m^-1
//
//
// I/O: 
//
//    os << m 
//        Writes m to ostream os in the following format:
//          colsize rowsize
//          m(0,0) m(0,1) m(0,2) ... m(0,rowsize)
//          m(1,0) m(1,1) m(1,2) ... m(1,rowsize)
//          ...
//          m(colsize,0) ... 
//
//    is >> m
//        Reads m from istream is in the same format
//        m must already be the correct size for this to work.
//
//    is >> mptr
//        If you do not know the size of the Matrix to be read, you can
//        use this form where mptr is a pointer to an undefined Matrix.
//        Then *mptr will be created with the right size using new,
//        so you should subsequently delete it.
//
//
// Division Control Functions:
//
//    There are a number of algorithms available for dividing
//    matrices.  We provide functions to allow you to 
//    change the algorithm used by the code on the fly.
//    In particular, you can write:
//    m.DivideUsing(dt)
//    where dt is LU, QR, QRP, SV, SVF, SVS, SVU or SVV
//    (ie. anything but CH)
//    Note: the last 3 cannot actually perform divisions.  They are 
//    available for when you want an SV decomposition, but don't need to 
//    store both of the U,V matrices.  SVS only stores S, 
//    SVU stores S,U, and SVV stores S,V.
//
//    The default method is LU (LU Decomposition) for square
//    matrices and QR (QR Decomposition) for non-square.
//   
//    See the header comment to TMV_Divider.h for more info about
//    the different algorithms.
//


#ifndef TMV_Matrix_H
#define TMV_Matrix_H

#include "TMV_Vector.h"

namespace tmv {

  template <class T> class GenMatrix;
  template <class T> class ConstMatrixView;
  template <class T> class MatrixView;
  template <class T, StorageType S=RowMajor> class Matrix;
  template <class T> class MatrixComposite; // Defined in TMV_MatrixArith.h

  inline StorageType TransOf(StorageType s)
  { return s==RowMajor ? ColMajor : s==ColMajor ? RowMajor : s; }

  template <class M> inline StorageType BaseStorOf(const M& m)
  {
    return (m.stor()==RowMajor || m.stor()==ColMajor) ? m.stor() : 
      m.stepi() > m.stepj() ? RowMajor : ColMajor;
  }

  template <class T> class LUDiv;
  template <class T> class QRDiv;
  template <class T> class QRPDiv;
  template <class T> class SVDiv;
  template <class T> class SVFDiv;

}

#include "TMV_BaseMatrix.h"

namespace tmv {

  // Defined in TMV_Matrix.cpp
  template <class T1, class T2> inline void Copy(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m2);

  template <class T> class GenMatrix : 
    public BaseMatrix<T>
  {

    public:

      //
      // Constructors
      //

      GenMatrix(StorageType s, ConjItType c) : itsstor(s), itsct(c) {}
      GenMatrix(StorageType s, ConjItType c, const MetaDivider<T>* div) : 
	BaseMatrix<T>(div), itsstor(s), itsct(c) {}
      GenMatrix(const GenMatrix<T>& rhs) : 
	BaseMatrix<T>(rhs), itsstor(rhs.itsstor), itsct(rhs.itsct) {}
      ~GenMatrix() {}

      //
      // Access Functions
      //

      inline ConstVectorView<T> row(size_t i) const 
      { 
	TMVAssert(i<colsize());
	return ConstVectorView<T>(cptr()+i*stepi(),rowsize(),stepj(),itsct); 
      }

      inline ConstVectorView<T> col(size_t j) const
      {
	TMVAssert(j<rowsize());
        return ConstVectorView<T>(cptr()+j*stepj(),colsize(),stepi(),itsct); 
      }

      inline ConstVectorView<T> diag() const
      {
	return ConstVectorView<T>(cptr(),min(colsize(),rowsize()),
	      stepi()+stepj(),itsct);
      }

      inline ConstVectorView<T> diag(int i) const
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	const size_t diagstep = stepi() + stepj();
	if (i >= 0) {
	  const size_t diagsize = min(colsize(),rowsize()-i);
	  return ConstVectorView<T>(cptr()+i*stepj(),diagsize,diagstep,itsct);
	} else {
	  const size_t diagsize = min(colsize()+i,rowsize());
	  return ConstVectorView<T>(cptr()-i*stepi(),diagsize,diagstep,itsct);
	}
      }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return ConstVectorView<T>(cptr()+i*stepi()+j1*stepj(),
	    j2-j1,stepj(),itsct); 
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return ConstVectorView<T>(cptr()+i1*stepi()+j*stepj(),
	    i2-i1,stepi(),itsct); 
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	const size_t diagstep = stepi() + stepj();
	if (i >= 0) {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize(),rowsize()-i));
	  return ConstVectorView<T>(cptr()+i*stepj()+j1*diagstep, 
	      j2-j1, diagstep,itsct);
	} else {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize()+i,rowsize()));
	  return ConstVectorView<T>(cptr()-i*stepi()+j1*diagstep, 
	      j2-j1, diagstep,itsct);
	}
      }

      inline ConstVectorView<T> operator[](size_t i) const
      { return row(i); }

      template <class T2> inline bool SameStorageAs(
	  const BaseMatrix<T2>& m2) const
      { return false; }

      inline bool SameStorageAs(const GenMatrix<T>& m2) const
      { return (cptr()==m2.cptr()); }

      template <class T2> inline bool SameAs(
	  const BaseMatrix<T2>& m2) const
      { return false; }

      inline bool SameAs(const GenMatrix<T>& m2) const
      { 
	return (this==&m2 || (cptr()==m2.cptr() && 
	    rowsize()==m2.rowsize() && colsize()==m2.colsize() &&
	    stepi()==m2.stepi() && stepj()==m2.stepj() && itsct==m2.itsct));
      }

      inline void CopyToMatrix(const MatrixView<RealType(T)>& m) const
      { 
	TMVAssert(IsReal(T()));
	TMVAssert(m.colsize() == colsize());
	TMVAssert(m.rowsize() == rowsize());
	Copy(*this,m);
      }
      inline void CopyToMatrix(const MatrixView<ComplexType(T)>& m) const
      {
	TMVAssert(m.colsize() == colsize());
	TMVAssert(m.rowsize() == rowsize());
	Copy(*this,m);
      }

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
	    i2-i1,j2-j1,stepi(),stepj(),stor(),itsct);
      }

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	StorageType newstor = 
	  iscm() ? (istep == 1 ? ColMajor : NoMajor) :
	  isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1*stepj(),
	    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
	    newstor,itsct);
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

      inline ConstMatrixView<T> ColPair(size_t j1, size_t j2) const
      {
	StorageType newstor = 
	  iscm() ? ColMajor : 
	  isrm() ? (j2==j1+1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(j1<rowsize() && j2<rowsize());
	return ConstMatrixView<T>(cptr()+j1*stepj(),colsize(),2,
	    stepi(),(int(j2)-int(j1))*stepj(),newstor,itsct);
      }

      inline ConstMatrixView<T> RowPair(size_t i1, size_t i2) const
      {
	StorageType newstor = 
	  isrm() ? RowMajor : 
	  iscm() ? (i2==i1+1 ? ColMajor : NoMajor) : NoMajor;
	TMVAssert(i1<colsize() && i2<colsize());
	return ConstMatrixView<T>(cptr()+i1*stepi(),2,rowsize(),
	    (int(i2)-int(i1))*stepi(),stepj(),newstor,itsct);
      }

      inline ConstMatrixView<T> Cols(size_t j1, size_t j2) const
      {
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return ConstMatrixView<T>(cptr()+j1*stepj(),colsize(),j2-j1,
	    stepi(),stepj(),stor(),itsct);
      }

      inline ConstMatrixView<T> Rows(size_t i1, size_t i2) const
      {
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return ConstMatrixView<T>(cptr()+i1*stepi(),i2-i1,rowsize(),
	    stepi(),stepj(),stor(),itsct);
      }

      inline ConstMatrixView<RealType(T)> Real() const
      {
	return ConstMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),colsize(),rowsize(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    IsReal(T()) ? stor() : NoMajor,NonConj);
      }

      inline ConstMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(itsct==NonConj);
	return ConstMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,
	    colsize(),rowsize(),2*stepi(),2*stepj(),NoMajor,NonConj);
      }

      //
      // Views
      //

      inline ConstMatrixView<T> View() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    stepi(),stepj(),stor(),itsct,
	    MakeMetaDivider(false,false,this));
      }

      inline ConstMatrixView<T> Transpose() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    stepj(),stepi(),TransOf(stor()),itsct,
	    MakeMetaDivider(true,false,this));
      }
  
      inline ConstMatrixView<T> Conjugate() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    stepi(),stepj(),stor(),ConjOf(T,itsct),
	    MakeMetaDivider(false,true,this));
      }

      inline ConstMatrixView<T> Adjoint() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    stepj(),stepi(),TransOf(stor()),ConjOf(T,itsct),
	    MakeMetaDivider(true,true,this));
      }

      inline ConstMatrixView<T> QuickView() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    stepi(),stepj(),stor(),itsct);
      }

      inline ConstMatrixView<T> QuickTranspose() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    stepj(),stepi(),TransOf(stor()),itsct);
      }
  
      inline ConstMatrixView<T> QuickConjugate() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    stepi(),stepj(),stor(),ConjOf(T,itsct));
      }

      inline ConstMatrixView<T> QuickAdjoint() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    stepj(),stepi(),TransOf(stor()),ConjOf(T,itsct));
      }

      inline ConstVectorView<T> LinearView() const
      {
	TMVAssert(stepi()==stepj()*int(rowsize()) || 
	    stepj()==stepi()*int(colsize()));
	return ConstVectorView<T>(cptr(),colsize()*rowsize(),
	    std::min(stepi(),stepj()),itsct);
      }

      //
      // Functions of Matrix
      //

      inline T Trace() const
      { return diag().SumElements(); }

      inline RealType(T) Norm() const 
      { return NormF(); }

      RealType(T) NormF() const;

      // NormF()^2
      RealType(T) NormSq() const;

      // 1-Norm = max_j (sum_i |a_ij|)
      RealType(T) Norm1() const;

      // 2-Norm defined in BaseMatrix.h

      // inf-Norm = max_i (sum_j |a_ij|)
      inline RealType(T) NormInf() const
      { return Transpose().Norm1(); }

      // = max_i,j (|a_ij|)
      RealType(T) MaxAbsElement() const;

      inline BaseMatrix<T>* NewTranspose() const 
      { return new ConstMatrixView<T>(Transpose()); }

      inline BaseMatrix<T>* NewConjugate() const
      { return new ConstMatrixView<T>(Conjugate()); }

      inline BaseMatrix<T>* NewAdjoint() const
      { return new ConstMatrixView<T>(Adjoint()); }

      inline BaseMatrix<T>* NewView() const
      { return new ConstMatrixView<T>(View()); }

      inline BaseMatrix<T>* NewCopy() const
      { 
	if (isrm()) return new Matrix<T,RowMajor>(*this); 
	else return new Matrix<T,ColMajor>(*this);
      }

      // 
      // Division Control
      //

      inline void DivideUsing(DivType dt) const
      {
	TMVAssert(dt == LU || dt == QR || dt == QRP || SV_Type(dt));
	BaseMatrix<T>::DivideUsing(dt); 
      }

      //
      // I/O
      //

      void Write(ostream& os) const;

      using BaseMatrix<T>::colsize;
      using BaseMatrix<T>::rowsize;
      using BaseMatrix<T>::IsSquare;
      virtual const T* cptr() const = 0;
      virtual int stepi() const = 0;
      virtual int stepj() const = 0;
      virtual inline StorageType stor() const { return itsstor; }
      virtual inline bool isrm() const { return itsstor == RowMajor; }
      virtual inline bool iscm() const { return itsstor == ColMajor; }
      virtual inline ConjItType ct() const { return itsct; }
      inline bool isconj() const 
      { 
	TMVAssert(IsComplex(T()) || itsct==NonConj);
	return IsComplex(T()) && itsct==Conj; 
      }

    protected :

      virtual T cref(size_t i, size_t j) const;

      void NewDivider() const;

    private :

      const StorageType itsstor;
      const ConjItType itsct;

      void operator=(const GenMatrix<T>&) { TMVAssert(false); }

  }; // GenMatrix

  template <class T> class ConstMatrixView : 
    public GenMatrix<T>
  {
    public :

      ConstMatrixView(const ConstMatrixView<T>& rhs) :
	GenMatrix<T>(rhs),
	itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs), 
	itssi(rhs.itssi), itssj(rhs.itssj) {}

      ConstMatrixView(const GenMatrix<T>& rhs) :
	GenMatrix<T>(rhs),
	itsm(rhs.cptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
	itssi(rhs.stepi()), itssj(rhs.stepj()) {}

      ConstMatrixView(const T* _m, size_t _cs, size_t _rs, int _si, int _sj, 
	  StorageType instor, ConjItType inct) : GenMatrix<T>(instor,inct),
	itsm(_m), itscs(_cs), itsrs(_rs), itssi(_si), itssj(_sj) 
      { 
	TMVAssert(instor==RowMajor ? _sj==1 : instor==ColMajor ? _si==1 : true);
      }

      ConstMatrixView(const T* _m, size_t _cs, size_t _rs, int _si, int _sj, 
	  StorageType instor, ConjItType inct, const MetaDivider<T>* mdiv) :
	GenMatrix<T>(instor,inct,mdiv),
	itsm(_m), itscs(_cs), itsrs(_rs), itssi(_si), itssj(_sj)
      { 
	TMVAssert(instor==RowMajor ? _sj==1 : instor==ColMajor ? _si==1 : true);
      }

      ~ConstMatrixView() {}

      inline size_t colsize() const { return itscs; }
      inline size_t rowsize() const { return itsrs; }
      inline const T* cptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }

    protected :

      const T*const itsm;
      const size_t itscs;
      const size_t itsrs;
      const int itssi;
      const int itssj;

    private :

      inline void operator=(const ConstMatrixView<T>&) { TMVAssert(false); }

  }; // ConstMatrixView

  template <class T> class MatrixView : 
    public GenMatrix<T>
  {

    public:

      //
      // Constructors
      //

      MatrixView(const MatrixView<T>& rhs) : 
	GenMatrix<T>(rhs), itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs),
	itssi(rhs.itssi), itssj(rhs.itssj) DEFFIRSTLAST(rhs.first,rhs.last) {}

      MatrixView(T* _m, size_t _cs, size_t _rs, int _si, int _sj,
	  StorageType instor, ConjItType inct PARAMFIRSTLAST(T) ) :
	GenMatrix<T>(instor,inct), itsm(_m), itscs(_cs), itsrs(_rs), 
	itssi(_si), itssj(_sj) DEFFIRSTLAST(_first,_last)
      {
	TMVAssert(instor==RowMajor ? _sj==1 : instor==ColMajor ? _si==1 : true);
      }

      MatrixView(T* _m, size_t _cs, size_t _rs, int _si, int _sj, 
	  StorageType instor, ConjItType inct, const MetaDivider<T>* mdiv 
	  PARAMFIRSTLAST(T) ) :
	GenMatrix<T>(instor,inct,mdiv), itsm(_m), itscs(_cs), itsrs(_rs),
	itssi(_si), itssj(_sj) DEFFIRSTLAST(_first,_last) 
      {
	TMVAssert(instor==RowMajor ? _sj==1 : instor==ColMajor ? _si==1 : true);
      }

      ~MatrixView() {} 

      //
      // Op=
      //

      inline const MatrixView<T>& operator=(const MatrixView<T>& m2) const
      { if (!SameAs(m2)) Copy(m2,*this); return *this; }

      inline const MatrixView<T>& operator=(const GenMatrix<T>& m2) const
      { if (!SameAs(m2)) Copy(m2,*this); return *this; }

      template <class T2> inline const MatrixView<T>& operator=(
	  const GenMatrix<T2>& m2) const
      { Copy(m2,*this); return *this; }

      template <class T2> inline const MatrixView<T>& operator=(
	  const BaseMatrix<T2>& m2) const
      { m2.CopyToMatrix(*this); return *this; }

      inline const MatrixView<T>& operator=(T x) const 
      { return SetToIdentity(x); }

      inline const MatrixView<T>& operator=(
	  const MatrixComposite<T>& mcomp) const
      { 
	TMVAssert(colsize() == mcomp.colsize());
	TMVAssert(rowsize() == mcomp.rowsize());
	mcomp.AssignTo(*this);
	return *this;
      }

      //
      // Access
      //
 
      inline VectorView<T> row(size_t i) const
      {
	TMVAssert(i<colsize());
	return VectorView<T>(ptr()+i*stepi(),rowsize(),stepj(),ct() 
	    FIRSTLAST ); 
      }

      inline VectorView<T> col(size_t j) const
      {
	TMVAssert(j<rowsize());
	return VectorView<T>(ptr()+j*stepj(),colsize(),stepi(),ct() 
	    FIRSTLAST ); 
      }

      inline VectorView<T> diag() const
      {
	return VectorView<T>(ptr(),min(colsize(),rowsize()),
	      stepi()+stepj(),ct() FIRSTLAST);
      }

      inline VectorView<T> diag(int i) const
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	const size_t diagstep = stepi() + stepj();
	if (i >= 0) {
	  const size_t diagsize = min(colsize(),rowsize()-i);
	  return VectorView<T>(ptr()+i*stepj(),diagsize,diagstep,ct() 
	      FIRSTLAST );
	} else {
	  const size_t diagsize = min(colsize()+i,rowsize());
	  return VectorView<T>(ptr()-i*stepi(),diagsize,diagstep,ct()
	      FIRSTLAST );
	}
      }

      inline VectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return VectorView<T>(ptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),ct()
	    FIRSTLAST ); 
      }

      inline VectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return VectorView<T>(ptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),ct()
	    FIRSTLAST ); 
      }

      inline VectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	const size_t diagstep = stepi() + stepj();
	if (i >= 0) {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize(),rowsize()-i));
	  return VectorView<T>(ptr()+i*stepj()+j1*diagstep,j2-j1,diagstep,ct() 
	      FIRSTLAST );
	} else {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize()+i,rowsize()));
	  return VectorView<T>(ptr()-i*stepi()+j1*diagstep,j2-j1,diagstep,ct()
	      FIRSTLAST );
	}
      }

      inline RefType(T) operator()(size_t i,size_t j) const 
      { return ref(i,j); }
      inline VectorView<T> operator[](size_t i) const 
      { return row(i); }

      //
      // Modifying Functions
      //

      inline const MatrixView<T>& Zero() const { return SetAllTo(T(0)); }
      const MatrixView<T>& SetAllTo(T x) const;
      const MatrixView<T>& TransposeSelf() const;
      const MatrixView<T>& ConjugateSelf() const;
      const MatrixView<T>& SetToIdentity(T x=T(1)) const;
      const MatrixView<T>& SwapRows(size_t i1, size_t i2) const;
      const MatrixView<T>& SwapCols(size_t j1, size_t j2) const;

      //
      // SubMatrix
      //

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) const
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T>(ptr()+i1*stepi()+j1*stepj(),
	    i2-i1,j2-j1,stepi(),stepj(),stor(),ct() FIRSTLAST );
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	StorageType newstor = 
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

      inline MatrixView<T> ColPair(size_t j1, size_t j2) const
      {
	StorageType newstor = 
	  iscm() ? ColMajor : 
	  isrm() ? (j2==j1+1 ? RowMajor : NoMajor) : NoMajor;
	TMVAssert(j1<rowsize() && j2<rowsize());
	return MatrixView<T>(ptr()+j1*stepj(),colsize(),2,
	    stepi(),(int(j2)-int(j1))*stepj(),newstor,ct() FIRSTLAST );
      }

      inline MatrixView<T> RowPair(size_t i1, size_t i2) const
      {
	StorageType newstor = 
	  isrm() ? RowMajor : 
	  iscm() ? (i2==i1+1 ? ColMajor : NoMajor) : NoMajor;
	TMVAssert(i1<colsize() && i2<colsize());
	return MatrixView<T>(ptr()+i1*stepi(),2,rowsize(),
	    (int(i2)-int(i1))*stepi(),stepj(),newstor,ct() FIRSTLAST );
      }

      inline MatrixView<T> Cols(size_t j1, size_t j2) const
      {
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return MatrixView<T>(ptr()+j1*stepj(),colsize(),j2-j1,
	    stepi(),stepj(),stor(),ct() FIRSTLAST);
      }

      inline MatrixView<T> Rows(size_t i1, size_t i2) const
      {
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return MatrixView<T>(ptr()+i1*stepi(),i2-i1,rowsize(),
	    stepi(),stepj(),stor(),ct() FIRSTLAST);
      }

      inline MatrixView<RealType(T)> Real() const
      {
	return MatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr()),colsize(),rowsize(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? stepj() : 2*stepj(),
	    IsReal(T()) ? stor() : NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)
	    ,reinterpret_cast<const RealType(T)*>(last)
#endif
	    );
      }

      inline MatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	TMVAssert(ct()==NonConj);
	return MatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,
	    colsize(),rowsize(),2*stepi(),2*stepj(), NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      //
      // Views
      //

      inline MatrixView<T> View() const
      { return *this; }

      inline MatrixView<T> Transpose() const
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    stepj(),stepi(),TransOf(stor()),ct(),
	    MakeMetaDivider(true,false,this) FIRSTLAST);
      }
  
      inline MatrixView<T> Conjugate() const
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    stepi(),stepj(),stor(),ConjOf(T,ct()),
	    MakeMetaDivider(false,true,this) FIRSTLAST);
      }

      inline MatrixView<T> Adjoint() const
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    stepj(),stepi(),TransOf(stor()),ConjOf(T,ct()),
	    MakeMetaDivider(true,true,this) FIRSTLAST);
      }

      inline MatrixView<T> QuickView() const
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    stepi(),stepj(),stor(),ct() FIRSTLAST);
      }

      inline MatrixView<T> QuickTranspose() const
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    stepj(),stepi(),TransOf(stor()),ct() FIRSTLAST);
      }
  
      inline MatrixView<T> QuickConjugate() const
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    stepi(),stepj(),stor(),ConjOf(T,ct()) FIRSTLAST);
      }

      inline MatrixView<T> QuickAdjoint() const
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    stepj(),stepi(),TransOf(stor()),ConjOf(T,ct()) FIRSTLAST);
      }

      inline VectorView<T> LinearView() const
      {
	TMVAssert(stepi()==stepj()*int(rowsize()) || 
	    stepj()==stepi()*int(colsize()));
	return VectorView<T>(ptr(),colsize()*rowsize(),
	    std::min(stepi(),stepj()),ct() FIRSTLAST);
      }

      
      //
      // I/O
      //

      void Read(istream& is) const;

      inline size_t colsize() const { return itscs; }
      inline size_t rowsize() const { return itsrs; }
      using BaseMatrix<T>::IsSquare;
      inline const T* cptr() const { return itsm; }
      inline T* ptr() const { return itsm; }
      inline int stepi() const { return itssi; }
      inline int stepj() const { return itssj; }
      using GenMatrix<T>::stor;
      using GenMatrix<T>::isrm;
      using GenMatrix<T>::iscm;
      using GenMatrix<T>::ct;
      using GenMatrix<T>::isconj;

    protected:

      T*const itsm;
      const size_t itscs;
      const size_t itsrs;
      const int itssi;
      const int itssj;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
    protected:
#endif

      RefType(T) ref(size_t i, size_t j) const;

  }; // MatrixView

  // The first definition is for the RowMajor version
  template <class T, StorageType S> class Matrix : 
    public GenMatrix<T> 
  {

    public:

      //
      // Constructors
      //

#define NEW_SIZE(cs,rs) GenMatrix<T>(RowMajor,NonConj), \
      itslen((cs)*(rs)), itsm(new T[itslen]), itscs(cs), itsrs(rs) \
      DEFFIRSTLAST(itsm,itsm+itslen)

      Matrix(size_t _colsize, size_t _rowsize) :
	NEW_SIZE(_colsize,_rowsize)
      {
	TMVAssert(S==RowMajor); 
#ifdef TMVDEBUG
	SetAllTo(T(888));
#endif
      }

      Matrix(size_t _colsize, size_t _rowsize, T x) :
	NEW_SIZE(_colsize,_rowsize) 
      { 
	TMVAssert(S==RowMajor);
	std::fill(itsm,itsm+itslen,x);
      }

      Matrix(size_t _colsize, size_t _rowsize, const valarray<T>& vv) : 
	NEW_SIZE(_colsize,_rowsize) 
      { 
	TMVAssert(S==RowMajor);
	TMVAssert(vv.size() == itslen);
	T* vi=itsm;
	for(size_t i=0;i<itslen;++i) *vi=vv[i];
      }

      Matrix(size_t _colsize,size_t _rowsize, const T* vv) : 
	NEW_SIZE(_colsize,_rowsize) 
      { 
	TMVAssert(S==RowMajor);
	memmove(itsm,vv,itslen*sizeof(T));
      }

      Matrix(size_t _colsize, size_t _rowsize, const vector<T>& vv) : 
	NEW_SIZE(_colsize,_rowsize) 
      { 
	TMVAssert(S==RowMajor);
	TMVAssert(vv.size() == itslen);
	std::copy(vv.begin(),vv.end(),itsm);
      }

      explicit Matrix(const vector<vector<T> >& vv) :
	NEW_SIZE(vv.size(),(vv.size()>0?vv[0].size():0)) 
	{ 
	  TMVAssert(S==RowMajor);
	  T* vi=itsm;
	  for(size_t i=0;i<colsize();++i,vi+=stepi()) {
	    TMVAssert(vv[i].size() == rowsize());
	    std::copy(vv[i].begin(),vv[i].end(),vi);
	  }
	}

      Matrix(const Matrix<T,S>& rhs) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize())
      { memmove(itsm,rhs.itsm,itslen*sizeof(T)); }

      template <class T2> Matrix(const GenMatrix<T2>& rhs) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize()) 
      { Copy(rhs,QuickView()); }

      template <class T2> Matrix(const BaseMatrix<T2>& rhs) : 
	  NEW_SIZE(rhs.colsize(),rhs.rowsize()) 
      {
	TMVAssert(S==RowMajor);
	rhs.CopyToMatrix(QuickView()); 
      }

      Matrix(const MatrixComposite<T>& mcomp) :
	NEW_SIZE(mcomp.colsize(),mcomp.rowsize())
      {
	TMVAssert(S==RowMajor);
	mcomp.AssignTo(QuickView()); 
      }

#undef NEW_SIZE

      ~Matrix() { TMVAssert(itsm); delete[] itsm; }

      //
      // Op=
      //

      inline Matrix<T,S>& operator=(const Matrix<T,S>& m2)
      { 
	TMVAssert(m2.colsize() == colsize() && m2.rowsize()==rowsize());
	if (&m2 != this) memmove(itsm,m2.itsm,itslen*sizeof(T));
	return *this; 
      }

      template <class T2> inline Matrix<T,S>& operator=(
	  const GenMatrix<T2>& m2)
      { Copy(m2,QuickView()); return *this; }

      template <class T2> inline Matrix<T,S>& operator=(
	  const BaseMatrix<T2>& m2)
      { m2.CopyToMatrix(QuickView()); return *this; }

      inline Matrix<T,S>& operator=(T x) { return SetToIdentity(x); }

      inline Matrix<T,S>& operator=(const MatrixComposite<T>& mcomp)
      { mcomp.AssignTo(QuickView()); return *this; }


      //
      // Access
      //

      inline ConstVectorView<T> row(size_t i) const 
      { 
	TMVAssert(i<colsize());
	return ConstVectorView<T>(cptr()+i*stepi(),rowsize(),1,NonConj); 
      }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return ConstVectorView<T>(cptr()+i*stepi()+j1,j2-j1,1,NonConj); 
      }

      inline ConstVectorView<T> operator[](size_t i) const
      { return row(i); }

      inline ConstVectorView<T> col(size_t j) const
      {
	TMVAssert(j<rowsize());
        return ConstVectorView<T>(cptr()+j,colsize(),stepi(),NonConj); 
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return ConstVectorView<T>(cptr()+i1*stepi()+j,i2-i1,stepi(),NonConj); 
      }

      inline ConstVectorView<T> diag() const
      {
	return ConstVectorView<T>(cptr(),min(colsize(),rowsize()),
	      stepi()+1,NonConj);
      }

      inline ConstVectorView<T> diag(int i) const
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	const size_t diagstep = stepi() + 1;
	if (i >= 0) {
	  const size_t diagsize = min(colsize(),rowsize()-i);
	  return ConstVectorView<T>(cptr()+i,diagsize,diagstep,NonConj);
	} else {
	  const size_t diagsize = min(colsize()+i,rowsize());
	  return ConstVectorView<T>(cptr()-i*stepi(),diagsize,diagstep,NonConj);
	}
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	const size_t diagstep = stepi() + 1;
	if (i >= 0) {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize(),rowsize()-i));
	  return ConstVectorView<T>(cptr()+i+j1*diagstep,j2-j1,diagstep,
	      NonConj);
	} else {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize()+i,rowsize()));
	  return ConstVectorView<T>(cptr()-i*stepi()+j1*diagstep,j2-j1,
	      diagstep,NonConj);
	}
      }

      inline T operator()(size_t i,size_t j) const
      { return cref(i,j); }

      inline VectorView<T> row(size_t i)
      { 
	TMVAssert(i<colsize());
	return VectorView<T>(ptr()+i*stepi(),rowsize(),1,NonConj FIRSTLAST); 
      }

      inline VectorView<T> row(size_t i, size_t j1, size_t j2)
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return VectorView<T>(ptr()+i*stepi()+j1,j2-j1,1,NonConj FIRSTLAST); 
      }

      inline VectorView<T> operator[](size_t i)
      { return row(i); }

      inline VectorView<T> col(size_t j)
      {
	TMVAssert(j<rowsize());
        return VectorView<T>(ptr()+j,colsize(),stepi(),NonConj FIRSTLAST); 
      }

      inline VectorView<T> col(size_t j, size_t i1, size_t i2)
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return VectorView<T>(ptr()+i1*stepi()+j,i2-i1,stepi(),NonConj 
	    FIRSTLAST); 
      }

      inline VectorView<T> diag()
      {
	return VectorView<T>(ptr(),min(colsize(),rowsize()),
	      stepi()+1,NonConj FIRSTLAST);
      }

      inline VectorView<T> diag(int i)
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	const size_t diagstep = stepi() + 1;
	if (i >= 0) {
	  const size_t diagsize = min(colsize(),rowsize()-i);
	  return VectorView<T>(ptr()+i,diagsize,diagstep,NonConj FIRSTLAST);
	} else {
	  const size_t diagsize = min(colsize()+i,rowsize());
	  return VectorView<T>(ptr()-i*stepi(),diagsize,diagstep,NonConj 
	      FIRSTLAST);
	}
      }

      inline VectorView<T> diag(int i, size_t j1, size_t j2) 
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	const size_t diagstep = stepi() + 1;
	if (i >= 0) {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize(),rowsize()-i));
	  return VectorView<T>(ptr()+i+j1*diagstep,j2-j1,diagstep,NonConj 
	      FIRSTLAST);
	} else {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize()+i,rowsize()));
	  return VectorView<T>(ptr()-i*stepi()+j1*diagstep,j2-j1,diagstep,
	      NonConj FIRSTLAST);
	}
      }

      inline T& operator()(size_t i,size_t j) 
      { return ref(i,j); }

      //
      // Modifying Functions
      //

      inline Matrix<T,S>& Zero() { return SetAllTo(0); }
      inline Matrix<T,S>& SetAllTo(T x) 
      { std::fill(itsm,itsm+itslen,x); return *this; }
      inline Matrix<T,S>& TransposeSelf() 
      { QuickView().TransposeSelf(); return *this; }
      inline Matrix<T,S>& ConjugateSelf() 
      { LinearView().ConjugateSelf(); return *this; }
      inline Matrix<T,S>& SetToIdentity(T x=T(1)) 
      { 
	TMVAssert(this->IsSquare());
	Zero(); diag().SetAllTo(x);
	return *this; 
      }
      inline Matrix<T,S>& SwapRows(size_t i1, size_t i2)
      { 
	TMVAssert(i1<colsize() && i2<colsize());
	if (i1!=i2) Swap(row(i1),row(i2));
	return *this; 
      }
      inline Matrix<T,S>& SwapCols(size_t j1, size_t j2)
      { 
	TMVAssert(j1<rowsize() && j2<rowsize());
	if (j1!=j2) Swap(col(j1),col(j2));
	return *this; 
      }

      //
      // SubMatrix
      //

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, 
	  int j1, int j2) const
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1,
	    i2-i1,j2-j1,stepi(),1,stor(),NonConj);
      }

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	StorageType newstor = jstep == 1 ? RowMajor : NoMajor;
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1*stepi()+j1,
	    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep, newstor,NonConj);
      }

      inline ConstVectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i*stepi()+j,size,
	    istep*stepi()+jstep,NonConj);
      }

      inline ConstMatrixView<T> ColPair(size_t j1, size_t j2) const
      {
	StorageType newstor = j2==j1+1 ? RowMajor : NoMajor;
	TMVAssert(j1<rowsize() && j2<rowsize());
	return ConstMatrixView<T>(cptr()+j1,colsize(),2,
	    stepi(),(int(j2)-int(j1)), newstor,NonConj);
      }

      inline ConstMatrixView<T> RowPair(size_t i1, size_t i2) const
      {
	TMVAssert(i1<colsize() && i2<colsize());
	return ConstMatrixView<T>(cptr()+i1*stepi(),2,rowsize(),
	    (int(i2)-int(i1))*stepi(),1, RowMajor,NonConj);
      }

      inline ConstMatrixView<T> Cols(size_t j1, size_t j2) const
      {
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return ConstMatrixView<T>(cptr()+j1,colsize(),j2-j1,stepi(),1,RowMajor,
	    NonConj);
      }

      inline ConstMatrixView<T> Rows(size_t i1, size_t i2) const
      {
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return ConstMatrixView<T>(cptr()+i1*stepi(),i2-i1,rowsize(),
	    stepi(),1,RowMajor,NonConj);
      }

      inline ConstMatrixView<RealType(T)> Real() const
      {
	return ConstMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),colsize(),rowsize(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? 1 : 2,
	    IsReal(T()) ? RowMajor : NoMajor,NonConj);
      }

      inline ConstMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	return ConstMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,
	    colsize(),rowsize(),2*stepi(),2,NoMajor,NonConj);
      }
      inline MatrixView<T> SubMatrix(int i1, int i2, 
	  int j1, int j2)
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T>(ptr()+i1*stepi()+j1,
	    i2-i1,j2-j1,stepi(),1,stor(),NonConj FIRSTLAST);
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) 
      {
	StorageType newstor = jstep == 1 ? RowMajor : NoMajor;
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T>(ptr()+i1*stepi()+j1,
	    (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep, newstor,NonConj 
	    FIRSTLAST);
      }

      inline VectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) 
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return VectorView<T>(ptr()+i*stepi()+j,size,
	    istep*stepi()+jstep,NonConj FIRSTLAST);
      }

      inline MatrixView<T> ColPair(size_t j1, size_t j2) 
      {
	StorageType newstor = j2==j1+1 ? RowMajor : NoMajor;
	TMVAssert(j1<rowsize() && j2<rowsize());
	return MatrixView<T>(ptr()+j1,colsize(),2,
	    stepi(),(int(j2)-int(j1)), newstor,NonConj FIRSTLAST);
      }

      inline MatrixView<T> RowPair(size_t i1, size_t i2) 
      {
	TMVAssert(i1<colsize() && i2<colsize());
	return MatrixView<T>(ptr()+i1*stepi(),2,rowsize(),
	    (int(i2)-int(i1))*stepi(),1, RowMajor,NonConj FIRSTLAST);
      }

      inline MatrixView<T> Cols(size_t j1, size_t j2) 
      {
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return MatrixView<T>(ptr()+j1,colsize(),j2-j1,
	    stepi(),1,RowMajor,NonConj FIRSTLAST);
      }

      inline MatrixView<T> Rows(size_t i1, size_t i2) 
      {
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return MatrixView<T>(ptr()+i1*stepi(),i2-i1,rowsize(),
	    stepi(),1,RowMajor,NonConj FIRSTLAST);
      }

      inline MatrixView<RealType(T)> Real() 
      {
	return MatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr()),colsize(),rowsize(),
	    IsReal(T()) ? stepi() : 2*stepi(),
	    IsReal(T()) ? 1 : 2,
	    IsReal(T()) ? RowMajor : NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      inline MatrixView<RealType(T)> Imag() 
      {
	TMVAssert(IsComplex(T()));
	return MatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,
	    colsize(),rowsize(),2*stepi(),2,NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      //
      // Views
      //

      inline ConstMatrixView<T> View() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    stepi(),1,RowMajor,NonConj,
	    new MetaDivider<T>(false,false,this));
      }

      inline ConstMatrixView<T> Transpose() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    1,stepi(),ColMajor,NonConj,
	    new MetaDivider<T>(true,false,this));
      }
  
      inline ConstMatrixView<T> Conjugate() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    stepi(),1,RowMajor,ConjOf(T,NonConj),
	    new MetaDivider<T>(false,true,this));
      }

      inline ConstMatrixView<T> Adjoint() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    1,stepi(),ColMajor,ConjOf(T,NonConj),
	    new MetaDivider<T>(true,true,this));
      }

      inline MatrixView<T> View()
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    stepi(),1,RowMajor,NonConj,
	    new MetaDivider<T>(false,false,this) FIRSTLAST);
      }

      inline MatrixView<T> Transpose()
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    1,stepi(),ColMajor,NonConj,
	    new MetaDivider<T>(true,false,this) FIRSTLAST);
      }
  
      inline MatrixView<T> Conjugate()
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    stepi(),1,RowMajor,ConjOf(T,NonConj),
	    new MetaDivider<T>(false,true,this) FIRSTLAST);
      }

      inline MatrixView<T> Adjoint()
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    1,stepi(),ColMajor,ConjOf(T,NonConj),
	    new MetaDivider<T>(true,true,this) FIRSTLAST);
      }

      inline ConstMatrixView<T> QuickView() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    stepi(),1,RowMajor,NonConj);
      }

      inline ConstMatrixView<T> QuickTranspose() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    1,stepi(),ColMajor,NonConj);
      }
  
      inline ConstMatrixView<T> QuickConjugate() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    stepi(),1,RowMajor,ConjOf(T,NonConj));
      }

      inline ConstMatrixView<T> QuickAdjoint() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    1,stepi(),ColMajor,ConjOf(T,NonConj));
      }

      inline MatrixView<T> QuickView()
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    stepi(),1,RowMajor,NonConj FIRSTLAST);
      }

      inline MatrixView<T> QuickTranspose()
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    1,stepi(),ColMajor,NonConj FIRSTLAST);
      }
  
      inline MatrixView<T> QuickConjugate()
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    stepi(),1,RowMajor,ConjOf(T,NonConj) FIRSTLAST);
      }

      inline MatrixView<T> QuickAdjoint()
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    1,stepi(),ColMajor,ConjOf(T,NonConj) FIRSTLAST);
      }

      inline ConstVectorView<T> LinearView() const
      { return ConstVectorView<T>(cptr(),itslen,1,NonConj); }

      inline VectorView<T> LinearView()
      { return VectorView<T>(ptr(),itslen,1,NonConj FIRSTLAST); }

      inline size_t colsize() const { return itscs; }
      inline size_t rowsize() const { return itsrs; }
      using BaseMatrix<T>::IsSquare;
      inline const T* cptr() const { return itsm; }
      inline T* ptr() { return itsm; }
      inline int stepi() const { return itsrs; }
      inline int stepj() const { return 1; }
      inline StorageType stor() const { return RowMajor; }
      inline bool isrm() const { return true; }
      inline bool iscm() const { return false; }
      inline ConjItType ct() const { return NonConj; }
      inline bool isconj() const { return false; }

    protected :

      const size_t itslen;
      T*const itsm;
      const size_t itscs;
      const size_t itsrs;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
    protected:
#endif

      inline T cref(size_t i, size_t j) const
      { 
	TMVAssert(i<colsize());
	TMVAssert(j<rowsize());
	return *(cptr() + int(i)*itsrs + int(j));
      }

      inline T& ref(size_t i, size_t j)
      { 
	TMVAssert(i<colsize());
	TMVAssert(j<rowsize());
	T*const mi = ptr() + int(i)*itsrs + int(j);
#ifdef TMVFLDEBUG
	TMVAssert(mi >= first);
	TMVAssert(mi < last);
#endif
	return *mi;
      }

  }; // RowMajor Matrix

  template <class T> class Matrix<T,ColMajor> : 
    public GenMatrix<T> 
  {

    public:

      //
      // Constructors
      //

#define NEW_SIZE(cs,rs) GenMatrix<T>(ColMajor,NonConj), \
      itslen((cs)*(rs)), itsm(new T[itslen]), itscs(cs), itsrs(rs) \
      DEFFIRSTLAST(itsm,itsm+itslen)

      Matrix(size_t _colsize, size_t _rowsize) :
	NEW_SIZE(_colsize,_rowsize)
      {
#ifdef TMVDEBUG
	SetAllTo(T(888));
#endif
      }

      Matrix(size_t _colsize, size_t _rowsize, T x) :
	NEW_SIZE(_colsize,_rowsize) 
      { std::fill(itsm,itsm+itslen,x); }

      Matrix(size_t _colsize, size_t _rowsize, const valarray<T>& vv) : 
	NEW_SIZE(_colsize,_rowsize) 
      { 
	TMVAssert(vv.size() == itslen); 
	T* vi=itsm;
	for(size_t i=0;i<vv.size();++i) *vi=vv[i];
      }

      Matrix(size_t _colsize,size_t _rowsize, const T* vv) : 
	NEW_SIZE(_colsize,_rowsize) 
      { memmove(itsm,vv,itslen*sizeof(T)); }

      Matrix(size_t _colsize, size_t _rowsize, const vector<T>& vv) : 
	NEW_SIZE(_colsize,_rowsize) 
      { 
	TMVAssert(vv.size() == itslen); 
	std::copy(vv.begin(),vv.end(),itsm);
      }

      explicit Matrix(const vector<vector<T> >& vv);

      Matrix(const Matrix<T,ColMajor>& rhs) :
	NEW_SIZE(rhs.colsize(),rhs.rowsize())
      { memmove(itsm,rhs.itsm,colsize()*rowsize()*sizeof(T)); }

      template <class T2> Matrix(
	  const GenMatrix<T2>& rhs) :
	  NEW_SIZE(rhs.colsize(),rhs.rowsize()) 
      { Copy(rhs,QuickView()); }

      template <class T2> Matrix(const BaseMatrix<T2>& rhs) : 
	  NEW_SIZE(rhs.colsize(),rhs.rowsize()) 
      { rhs.CopyToMatrix(QuickView()); }

      Matrix(const MatrixComposite<T>& mcomp) :
	NEW_SIZE(mcomp.colsize(),mcomp.rowsize())
      { mcomp.AssignTo(QuickView()); }
#undef NEW_SIZE

      ~Matrix() { TMVAssert(itsm); delete[] itsm; }

      //
      // Op=
      //

      inline Matrix<T,ColMajor>& operator=(const Matrix<T,ColMajor>& m2)
      { 
	TMVAssert(m2.colsize() == colsize() && m2.rowsize()==rowsize());
	if (&m2 != this) memmove(itsm,m2.itsm,itslen*sizeof(T));
	return *this; 
      }

      template <class T2> inline Matrix<T,ColMajor>& operator=(
	  const GenMatrix<T2>& m2)
      { Copy(m2,QuickView()); return *this; }

      template <class T2> inline Matrix<T,ColMajor>& operator=(
	  const BaseMatrix<T2>& m2)
      { m2.CopyToMatrix(QuickView()); return *this; }

      inline Matrix<T,ColMajor>& operator=(T x) { return SetToIdentity(x); }

      inline Matrix<T,ColMajor>& operator=(const MatrixComposite<T>& mcomp)
      { mcomp.AssignTo(QuickView()); return *this; }


      //
      // Access
      //

      inline ConstVectorView<T> row(size_t i) const 
      { 
	TMVAssert(i<colsize());
	return ConstVectorView<T>(cptr()+i,rowsize(),colsize(),NonConj); 
      }

      inline ConstVectorView<T> row(size_t i, size_t j1, size_t j2) const 
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return ConstVectorView<T>(cptr()+i+j1*colsize(),j2-j1,colsize(),
	    NonConj); 
      }

      inline ConstVectorView<T> operator[](size_t i) const
      { return row(i); }

      inline ConstVectorView<T> col(size_t j) const
      {
	TMVAssert(j<rowsize());
        return ConstVectorView<T>(cptr()+j*colsize(),colsize(),1,NonConj); 
      }

      inline ConstVectorView<T> col(size_t j, size_t i1, size_t i2) const
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return ConstVectorView<T>(cptr()+i1+j*colsize(),i2-i1,1,NonConj); 
      }

      inline ConstVectorView<T> diag() const
      {
	return ConstVectorView<T>(cptr(),min(colsize(),rowsize()),
	      1+stepj(),NonConj);
      }

      inline ConstVectorView<T> diag(int i) const
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	size_t diagstep = 1 + colsize();
	if (i >= 0) {
	  size_t diagsize = min(colsize(),rowsize()-i);
	  return ConstVectorView<T>(cptr()+i*colsize(),diagsize,diagstep,
	      NonConj);
	} else {
	  i = -i;
	  size_t diagsize = min(colsize()-i,rowsize());
	  return ConstVectorView<T>(cptr()+i,diagsize,diagstep,NonConj);
	}
      }

      inline ConstVectorView<T> diag(int i, size_t j1, size_t j2) const
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	size_t diagstep = 1 + colsize();
	if (i >= 0) {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize(),rowsize()-i));
	  return ConstVectorView<T>(cptr()+i*colsize()+j1*diagstep, 
	      j2-j1, diagstep,NonConj);
	} else {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize()+i,rowsize()));
	  return ConstVectorView<T>(cptr()-i+j1*diagstep,j2-j1,diagstep,
	      NonConj);
	}
      }

      inline T operator()(size_t i, size_t j) const
      { return cref(i,j); }

      inline VectorView<T> row(size_t i) 
      { 
	TMVAssert(i<colsize());
	return VectorView<T>(ptr()+i,rowsize(),colsize(),NonConj FIRSTLAST); 
      }

      inline VectorView<T> row(size_t i, size_t j1, size_t j2) 
      { 
	TMVAssert(i<colsize());
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return VectorView<T>(ptr()+i+j1*colsize(),j2-j1,colsize(),NonConj 
	    FIRSTLAST); 
      }

      inline VectorView<T> operator[](size_t i)
      { return row(i); }

      inline VectorView<T> col(size_t j) 
      {
	TMVAssert(j<rowsize());
        return VectorView<T>(ptr()+j*colsize(),colsize(),1,NonConj FIRSTLAST); 
      }

      inline VectorView<T> col(size_t j, size_t i1, size_t i2) 
      {
	TMVAssert(j<rowsize());
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return VectorView<T>(ptr()+i1+j*colsize(),i2-i1,1,NonConj FIRSTLAST); 
      }

      inline VectorView<T> diag() 
      {
	return VectorView<T>(ptr(),min(colsize(),rowsize()),
	      1+stepj(),NonConj FIRSTLAST);
      }

      inline VectorView<T> diag(int i) 
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	size_t diagstep = 1 + colsize();
	if (i >= 0) {
	  size_t diagsize = min(colsize(),rowsize()-i);
	  return VectorView<T>(ptr()+i*colsize(),diagsize,diagstep,NonConj 
	      FIRSTLAST);
	} else {
	  size_t diagsize = min(colsize()+i,rowsize());
	  return VectorView<T>(ptr()-i,diagsize,diagstep,NonConj FIRSTLAST);
	}
      }

      inline VectorView<T> diag(int i, size_t j1, size_t j2) 
      {
	TMVAssert(i<=int(rowsize()));
	TMVAssert(-i<=int(colsize()));
	size_t diagstep = 1 + colsize();
	if (i >= 0) {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize(),rowsize()-i));
	  return VectorView<T>(ptr()+i*colsize()+j1*diagstep,j2-j1,diagstep,
	      NonConj FIRSTLAST);
	} else {
	  TMVAssert(j1 <= j2);
	  TMVAssert(j2 <= min(colsize()+i,rowsize()));
	  return VectorView<T>(ptr()-i+j1*diagstep,j2-j1,diagstep,NonConj 
	      FIRSTLAST);
	}
      }

      inline T& operator()(size_t i, size_t j) 
      { return ref(i,j); }

      //
      // Modifying Functions
      //

      inline Matrix<T,ColMajor>& Zero() { return SetAllTo(0); }
      inline Matrix<T,ColMajor>& SetAllTo(T x) 
      { std::fill(itsm,itsm+itslen,x); return *this; }
      inline Matrix<T,ColMajor>& TransposeSelf() 
      { QuickView().TransposeSelf(); return *this; }
      inline Matrix<T,ColMajor>& ConjugateSelf() 
      { LinearView().ConjugateSelf(); return *this; }
      inline Matrix<T,ColMajor>& SetToIdentity(T x=T(1)) 
      { 
	TMVAssert(this->IsSquare());
	Zero(); diag().SetAllTo(x);
	return *this; 
      }
      inline Matrix<T,ColMajor>& SwapRows(size_t i1, size_t i2)
      { 
	TMVAssert(i1<colsize() && i2<colsize());
	if (i1!=i2) Swap(row(i1),row(i2));
	return *this; 
      }
      inline Matrix<T,ColMajor>& SwapCols(size_t j1, size_t j2)
      { 
	TMVAssert(j1<rowsize() && j2<rowsize());
	if (j1!=j2) Swap(col(j1),col(j2));
	return *this; 
      }

      //
      // SubMatrix
      //

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) const
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return ConstMatrixView<T>(cptr()+i1+j1*stepj(),
	    i2-i1,j2-j1,1,stepj(),ColMajor,NonConj);
      }

      inline ConstMatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) const
      {
	StorageType newstor = istep == 1 ? ColMajor : NoMajor;
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return ConstMatrixView<T>(cptr()+i1+j1*stepj(),
	    (i2-i1)/istep,(j2-j1)/jstep,istep,jstep*stepj(), newstor,NonConj);
      }

      inline ConstVectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) const
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return ConstVectorView<T>(cptr()+i+j*stepj(),size,
	    istep+jstep*stepj(),NonConj);
      }

      inline ConstMatrixView<T> ColPair(size_t j1, size_t j2) const
      {
	TMVAssert(j1<rowsize() && j2<rowsize());
	return ConstMatrixView<T>(cptr()+j1*stepj(),colsize(),2,
	    1,(int(j2)-int(j1))*stepj(),ColMajor,NonConj);
      }

      inline ConstMatrixView<T> RowPair(size_t i1, size_t i2) const
      {
	StorageType newstor = i2==i1+1 ? ColMajor : NoMajor;
	TMVAssert(i1<colsize() && i2<colsize());
	return ConstMatrixView<T>(cptr()+i1,2,rowsize(),
	    (int(i2)-int(i1)),stepj(),newstor,NonConj);
      }

      inline ConstMatrixView<T> Cols(size_t j1, size_t j2) const
      {
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return ConstMatrixView<T>(cptr()+j1*stepj(),colsize(),j2-j1,
	    1,stepj(),ColMajor,NonConj);
      }

      inline ConstMatrixView<T> Rows(size_t i1, size_t i2) const
      {
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return ConstMatrixView<T>(cptr()+i1,i2-i1,rowsize(),
	    1,stepj(),ColMajor,NonConj);
      }

      inline ConstMatrixView<RealType(T)> Real() const
      {
	return ConstMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr()),colsize(),rowsize(),
	    IsReal(T()) ? 1 : 2,
	    IsReal(T()) ? stepj() : 2*stepj(),
	    IsReal(T()) ? ColMajor : NoMajor,NonConj);
      }

      inline ConstMatrixView<RealType(T)> Imag() const
      {
	TMVAssert(IsComplex(T()));
	return ConstMatrixView<RealType(T)>(
	    reinterpret_cast<const RealType(T)*>(cptr())+1,
	    colsize(),rowsize(),2,2*stepj(),NoMajor,NonConj);
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2) 
      {
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,1,1));
	return MatrixView<T>(ptr()+i1+j1*stepj(),
	    i2-i1,j2-j1,1,stepj(),ColMajor,NonConj FIRSTLAST);
      }

      inline MatrixView<T> SubMatrix(int i1, int i2, int j1, int j2,
	  int istep, int jstep) 
      {
	StorageType newstor = istep == 1 ? ColMajor : NoMajor;
	TMVAssert(this->OKSubMatrix(i1,i2,j1,j2,istep,jstep));
	return MatrixView<T>(ptr()+i1+j1*stepj(),
	    (i2-i1)/istep,(j2-j1)/jstep,istep,jstep*stepj(), newstor,
	    NonConj FIRSTLAST);
      }

      inline VectorView<T> SubVector(size_t i, size_t j,
	  int istep, int jstep, size_t size) 
      {
	TMVAssert(this->OKSubVector(i,j,istep,jstep,size));
	return VectorView<T>(ptr()+i+j*stepj(),size,
	    istep+jstep*stepj(),NonConj FIRSTLAST);
      }

      inline MatrixView<T> ColPair(size_t j1, size_t j2) 
      {
	TMVAssert(j1<rowsize() && j2<rowsize());
	return MatrixView<T>(ptr()+j1*stepj(),colsize(),2,
	    1,(int(j2)-int(j1))*stepj(), ColMajor,NonConj FIRSTLAST);
      }

      inline MatrixView<T> RowPair(size_t i1, size_t i2) 
      {
	StorageType newstor = i2==i1+1 ? ColMajor : NoMajor;
	TMVAssert(i1<colsize() && i2<colsize());
	return MatrixView<T>(ptr()+i1,2,rowsize(),
	    (int(i2)-int(i1)),stepj(), newstor,NonConj FIRSTLAST);
      }

      inline MatrixView<T> Cols(size_t j1, size_t j2) 
      {
	TMVAssert(j1<=j2);
	TMVAssert(j2<=rowsize());
	return MatrixView<T>(ptr()+j1*stepj(),colsize(),j2-j1,
	    1,stepj(),ColMajor,NonConj FIRSTLAST);
      }

      inline MatrixView<T> Rows(size_t i1, size_t i2) 
      {
	TMVAssert(i1<=i2);
	TMVAssert(i2<=colsize());
	return MatrixView<T>(ptr()+i1,i2-i1,rowsize(),
	    1,stepj(),ColMajor,NonConj FIRSTLAST);
      }

      inline MatrixView<RealType(T)> Real() 
      {
	return MatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr()),colsize(),rowsize(),
	    IsReal(T()) ? 1 : 2,
	    IsReal(T()) ? stepj() : 2*stepj(),
	    IsReal(T()) ? ColMajor : NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      inline MatrixView<RealType(T)> Imag() 
      {
	TMVAssert(IsComplex(T()));
	return MatrixView<RealType(T)>(
	    reinterpret_cast<RealType(T)*>(ptr())+1,
	    colsize(),rowsize(),2,2*stepj(),NoMajor,NonConj
#ifdef TMVFLDEBUG
	    ,reinterpret_cast<const RealType(T)*>(first)+1
	    ,reinterpret_cast<const RealType(T)*>(last)+1
#endif
	    );
      }

      //
      // Views
      //
      
      inline ConstMatrixView<T> View() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    1,colsize(),ColMajor,NonConj,
	    new MetaDivider<T>(false,false,this));
      }

      inline ConstMatrixView<T> Transpose() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    colsize(),1,RowMajor,NonConj,
	    new MetaDivider<T>(true,false,this));
      }
  
      inline ConstMatrixView<T> Conjugate() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    1,colsize(),ColMajor,ConjOf(T,NonConj),
	    new MetaDivider<T>(false,true,this));
      }

      inline ConstMatrixView<T> Adjoint() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    colsize(),1,RowMajor,ConjOf(T,NonConj),
	    new MetaDivider<T>(true,true,this));
      }

      inline MatrixView<T> View()
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    1,colsize(),ColMajor,NonConj,
	    new MetaDivider<T>(false,false,this) FIRSTLAST);
      }

      inline MatrixView<T> Transpose()
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    colsize(),1,RowMajor,NonConj,
	    new MetaDivider<T>(true,false,this) FIRSTLAST);
      }
  
      inline MatrixView<T> Conjugate()
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    1,colsize(),ColMajor,ConjOf(T,NonConj),
	    new MetaDivider<T>(false,true,this) FIRSTLAST);
      }

      inline MatrixView<T> Adjoint()
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    colsize(),1,RowMajor,ConjOf(T,NonConj),
	    new MetaDivider<T>(true,true,this) FIRSTLAST);
      }

      inline ConstMatrixView<T> QuickView() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    1,colsize(),ColMajor,NonConj);
      }

      inline ConstMatrixView<T> QuickTranspose() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    colsize(),1,RowMajor,NonConj);
      }
  
      inline ConstMatrixView<T> QuickConjugate() const
      { 
	return ConstMatrixView<T>(cptr(),colsize(),rowsize(),
	    1,colsize(),ColMajor,ConjOf(T,NonConj));
      }

      inline ConstMatrixView<T> QuickAdjoint() const
      { 
	return ConstMatrixView<T>(cptr(),rowsize(),colsize(),
	    colsize(),1,RowMajor,ConjOf(T,NonConj));
      }

      inline MatrixView<T> QuickView()
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    1,colsize(),ColMajor,NonConj FIRSTLAST);
      }

      inline MatrixView<T> QuickTranspose()
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    colsize(),1,RowMajor,NonConj FIRSTLAST);
      }
  
      inline MatrixView<T> QuickConjugate()
      { 
	return MatrixView<T>(ptr(),colsize(),rowsize(),
	    1,colsize(),ColMajor,ConjOf(T,NonConj) FIRSTLAST);
      }

      inline MatrixView<T> QuickAdjoint()
      { 
	return MatrixView<T>(ptr(),rowsize(),colsize(),
	    colsize(),1,RowMajor,ConjOf(T,NonConj) FIRSTLAST);
      }

      inline ConstVectorView<T> LinearView() const
      { return ConstVectorView<T>(cptr(),itslen,1,NonConj); }

      inline VectorView<T> LinearView()
      { return VectorView<T>(ptr(),itslen,1,NonConj FIRSTLAST); }

      inline size_t colsize() const { return itscs; }
      inline size_t rowsize() const { return itsrs; }
      using BaseMatrix<T>::IsSquare;
      inline const T* cptr() const { return itsm; }
      inline T* ptr() { return itsm; }
      inline int stepi() const { return 1; }
      inline int stepj() const { return itscs; }
      inline StorageType stor() const { return ColMajor; }
      inline bool isrm() const { return false; }
      inline bool iscm() const { return true; }
      inline ConjItType ct() const { return NonConj; }
      inline bool isconj() const { return false; }

    protected :

      const size_t itslen;
      T*const itsm;
      const size_t itscs;
      const size_t itsrs;

#ifdef TMVFLDEBUG
    public:
      const T*const first;
      const T*const last;
    protected:
#endif

      inline T cref(size_t i, size_t j) const
      { 
	TMVAssert(i<colsize());
	TMVAssert(j<rowsize());
	return *(cptr() + i + int(j)*itscs);
      }

      inline T& ref(size_t i, size_t j)
      { 
	TMVAssert(i<colsize());
	TMVAssert(j<rowsize());
	T*const mi = ptr() + i + int(j)*itscs;
#ifdef TMVFLDEBUG
	TMVAssert(mi >= first);
	TMVAssert(mi < last);
#endif
	return *mi;
      }
  };  // ColMajor Matrix

//---------------------------------------------------------------------------

  //
  // Special Creators: 
  //   Eye(n) = Identity
  //   Eye(n,x) = x * Identity
  //   Diag(v) = Diagonal Matrix with vector v along the diagonal
  //   OuterProduct(v1,v2) = Outer Product of two Vectors
  //   RowVector(v) = 1xn Matrix with v in only row
  //   ColVector(v) = nx1 Matrix with v in only col
  //   RowVectorViewOf(v) = 1xn Matrix with v in only row - Same Storage
  //   ColVectorViewOf(v) = nx1 Matrix with v in only col - Same Storage
  //   MatrixViewOf(m,colsize,rowsize,storage) = MatrixView of m 
  //

  template <class T, StorageType S> inline Matrix<T,S> Eye(
      size_t n, T x=1)
  {
    Matrix<T,S> temp(n,n,T(0));
    temp.diag().SetAllTo(x);
    return temp;
  }
  template <class T> inline Matrix<T,RowMajor> Eye(size_t n, T x=1)
  { return Eye<T,RowMajor>(n,x); }

  template <class T, StorageType S> inline Matrix<T,S> Diag(
      const GenVector<T>& v)
  { 
    Matrix<T,S> temp(v.size(),v.size(),T(0));
    temp.diag() = v;
    return temp;
  }
  template <class T> inline Matrix<T,RowMajor> Diag(const GenVector<T>& v)
  { return Diag<T,RowMajor>(v); }

  template <class T, StorageType S> inline Matrix<T,S> OuterProduct(
      const GenVector<T>& v1, const GenVector<T>& v2)
  {
    Matrix<T,S> temp(v1.size(),v2.size(),T(0));
    Rank1Update(T(1),v1,v2,temp.QuickView());
    return temp;
  }
  template <class T> inline Matrix<T,RowMajor> OuterProduct(
      const GenVector<T>& v1, const GenVector<T>& v2)
  { return OuterProduct<T,RowMajor>(v1,v2); }

  template <class T> inline Matrix<T,RowMajor> RowVector(const GenVector<T>& v)
  { 
    Matrix<T,RowMajor> temp(1,v.size());
    temp.row(0) = v;
    return temp;
  }

  template <class T> inline Matrix<T,ColMajor> ColVector(const GenVector<T>& v)
  { 
    Matrix<T,ColMajor> temp(v.size(),1);
    temp.col(0) = v;
    return temp;
  }

  template <class T> inline ConstMatrixView<T> RowVectorViewOf(
      const GenVector<T>& v)
  {
    if (v.step() == 1)
      return ConstMatrixView<T>(v.cptr(),1,v.size(),v.size(),1,RowMajor,v.ct());
    else
      return ConstMatrixView<T>(v.cptr(),1,v.size(),1,v.step(),ColMajor,v.ct());
  }

  template <class T> inline ConstMatrixView<T> ColVectorViewOf(
      const GenVector<T>& v)
  { 
    if (v.step() == 1)
      return ConstMatrixView<T>(v.cptr(),v.size(),1,1,v.size(),ColMajor,v.ct());
    else
      return ConstMatrixView<T>(v.cptr(),v.size(),1,v.step(),1,RowMajor,v.ct());
  }

  template <class T> inline MatrixView<T> RowVectorViewOf(
      const VectorView<T>& v)
  {
    if (v.step() == 1)
      return MatrixView<T>(v.ptr(),1,v.size(),v.size(),1,RowMajor,v.ct()
	  FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    else
      return MatrixView<T>(v.ptr(),1,v.size(),1,v.step(),ColMajor,v.ct()
	  FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
  }

  template <class T> inline MatrixView<T> ColVectorViewOf(
      const VectorView<T>& v)
  { 
    if (v.step() == 1)
      return MatrixView<T>(v.ptr(),v.size(),1,1,v.size(),ColMajor,v.ct()
	  FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    else
      return MatrixView<T>(v.ptr(),v.size(),1,v.step(),1,RowMajor,v.ct()
	  FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
  }

  template <class T> inline MatrixView<T> RowVectorViewOf(Vector<T>& v)
  {
    return MatrixView<T>(v.ptr(),1,v.size(),v.size(),1,RowMajor,NonConj
	FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
  }

  template <class T> inline MatrixView<T> ColVectorViewOf(Vector<T>& v)
  { 
    return MatrixView<T>(v.ptr(),v.size(),1,1,v.size(),ColMajor,NonConj
	FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
  }

  template <class T> inline MatrixView<T> MatrixViewOf(
      T* m, size_t colsize, size_t rowsize, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return MatrixView<T>(m,colsize,rowsize,rowsize,1,RowMajor,NonConj
	  FIRSTLAST1(m,m+colsize*rowsize));
    else 
      return MatrixView<T>(m,colsize,rowsize,1,colsize,ColMajor,NonConj
	  FIRSTLAST1(m,m+colsize*rowsize));
  }

  template <class T> inline ConstMatrixView<T> MatrixViewOf(
      const T* m, size_t colsize, size_t rowsize, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstMatrixView<T>(m,colsize,rowsize,rowsize,1,RowMajor,NonConj);
    else 
      return ConstMatrixView<T>(m,colsize,rowsize,1,colsize,ColMajor,NonConj);
  }

  //
  // Copy Matrices
  //

  template <class T1, class T2> inline void NonLapCopy(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.ct()==NonConj);
    if (m1.isrm() && m2.isrm()) {
      const size_t size = m1.rowsize();
      if (m1.isconj())
	for(size_t i=0;i<m2.colsize();++i) 
	  DoCopy(CVIt<T1,Unit,Conj>(m1.row(i).begin()),
	      VIt<T2,Unit,NonConj>(m2.row(i).begin()),size);
      else
	for(size_t i=0;i<m2.colsize();++i) 
	  DoCopy(CVIt<T1,Unit,NonConj>(m1.row(i).begin()),
	      VIt<T2,Unit,NonConj>(m2.row(i).begin()),size);
    }
    else if (m1.iscm() && m2.iscm()) {
      const size_t size = m1.colsize();
      if (m1.isconj())
	for(size_t j=0;j<m2.rowsize();++j) 
	  DoCopy(CVIt<T1,Unit,Conj>(m1.col(j).begin()),
	      VIt<T2,Unit,NonConj>(m2.col(j).begin()),size);
      else
	for(size_t j=0;j<m2.rowsize();++j) 
	  DoCopy(CVIt<T1,Unit,NonConj>(m1.col(j).begin()),
	      VIt<T2,Unit,NonConj>(m2.col(j).begin()),size);
    }
    else if (m2.colsize() <= m2.rowsize())
      for(size_t i=0;i<m2.colsize();++i) m2.row(i) = m1.row(i);
    else
      for(size_t j=0;j<m2.rowsize();++j) m2.col(j) = m1.col(j);
  }
#ifdef LAP
  template <class T> inline void LapCopy(
      const GenMatrix<T>& m1, const MatrixView<T>& m2)
  { NonLapCopy(m1,m2); }
  template <> inline void LapCopy(
      const GenMatrix<double>& m1, const MatrixView<double>& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert((m1.iscm()&&m2.iscm()) || (m1.isrm()&&m2.isrm()));
    char c = 'A';
    int m = m1.isrm() ? m1.rowsize() : m1.colsize();
    int n = m1.isrm() ? m1.colsize() : m1.rowsize();
    int ld1 = m1.isrm() ? m1.stepi() : m1.stepj();
    int ld2 = m1.isrm() ? m2.stepi() : m2.stepj();
    dlacpy(&c,&m,&n,const_cast<double*>(m1.cptr()),&ld1,
	const_cast<double*>(m2.ptr()),&ld2);
  }
  template <> inline void LapCopy(
      const GenMatrix<complex<double> >& m1,
      const MatrixView<complex<double> >& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.iscm());
    TMVAssert(m2.iscm());
    char c = 'A';
    int m = m1.isrm() ? m1.rowsize() : m1.colsize();
    int n = m1.isrm() ? m1.colsize() : m1.rowsize();
    int ld1 = m1.isrm() ? m1.stepi() : m1.stepj();
    int ld2 = m1.isrm() ? m2.stepi() : m2.stepj();
    zlacpy(&c,&m,&n,LAP_Complex(m1.cptr()),&ld1,LAP_Complex(m2.ptr()),&ld2);
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
    int m = m1.isrm() ? m1.rowsize() : m1.colsize();
    int n = m1.isrm() ? m1.colsize() : m1.rowsize();
    int ld1 = m1.isrm() ? m1.stepi() : m1.stepj();
    int ld2 = m1.isrm() ? m2.stepi() : m2.stepj();
    slacpy(&c,&m,&n,const_cast<float*>(m1.cptr()),&ld1,
	const_cast<float*>(m2.ptr()),&ld2);
  }
  template <> inline void LapCopy(
      const GenMatrix<complex<float> >& m1,
      const MatrixView<complex<float> >& m2)
  { 
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m1.iscm());
    TMVAssert(m2.iscm());
    char c = 'A';
    int m = m1.isrm() ? m1.rowsize() : m1.colsize();
    int n = m1.isrm() ? m1.colsize() : m1.rowsize();
    int ld1 = m1.isrm() ? m1.stepi() : m1.stepj();
    int ld2 = m1.isrm() ? m2.stepi() : m2.stepj();
    clacpy(&c,&m,&n,LAP_Complex(m1.cptr()),&ld1,LAP_Complex(m2.ptr()),&ld2);
  }
#endif
#endif
  template <class T> inline void CallCopy(const GenMatrix<T>& m1, 
      const MatrixView<T>& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.ct()==NonConj);
#ifdef LAP
    if (!m1.isconj()) 
      if (m1.isrm() && m2.isrm())
	LapCopy(m1.QuickTranspose(),m2.QuickTranspose());
      else if (m1.iscm() && m2.iscm()) LapCopy(m1,m2);
      else NonLapCopy(m1,m2);
    else
#endif
      NonLapCopy(m1,m2);
  }
  template <class T1, class T2> inline void CallCopy(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.ct()==NonConj);
    NonLapCopy(m1,m2);
  }
  template <class T> inline void CallCopy(
      const GenMatrix<T>& m1, const MatrixView<complex<T> >& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    TMVAssert(m2.ct()==NonConj);
    m2.Imag().Zero();
    CallCopy(m1,m2.Real());
  }
  template <class T> inline void CallCopy(
      const GenMatrix<complex<T> >& m1, const MatrixView<T>& m2)
  { TMVAssert(false); }

  template <class T1, class T2> inline void Copy(
      const GenMatrix<T1>& m1, const MatrixView<T2>& m2)
  {
    TMVAssert(m2.rowsize() == m1.rowsize());
    TMVAssert(m2.colsize() == m1.colsize());
    if (m1.stepi()==m2.stepi() && m1.stepj()==m2.stepj() && 
	((m2.stepi()==int(1) && m2.stepj()==int(m2.colsize())) ||
	 (m2.stepj()==int(1) && m2.stepi()==int(m2.rowsize())) ||
	 (m2.stepi()<m2.stepj()&&m2.stepj()==m2.stepi()*int(m2.colsize())) || 
	 (m2.stepj()<m2.stepi()&&m2.stepi()==m2.stepj()*int(m2.rowsize())))) {
      m2.LinearView() = m1.LinearView();
    } else  {
      if (m2.isconj()) CallCopy(m1.Conjugate(),m2.Conjugate());
      else CallCopy(m1,m2);
    }
  }

  //
  // Swap Matrices
  //

  template <class T> void Swap(const MatrixView<T>& m1,
      const MatrixView<T>& m2);
  template <class T, StorageType S> inline void Swap(
      const MatrixView<T>& m1, Matrix<T,S>& m2)
  { Swap(m1,m2.QuickView()); }
  template <class T, StorageType S> inline void Swap(
      Matrix<T,S>& m1, const MatrixView<T>& m2)
  { Swap(m1.QuickView(),m2); }
  template <class T, StorageType S1, StorageType S2> inline void Swap(
      Matrix<T,S1>& m1, Matrix<T,S2>& m2)
  { Swap(m1.QuickView(),m2.QuickView()); }

  //
  // Functions of Matrices:
  //

  template <class T> inline T Det(const GenMatrix<T>& m)
  { return m.Det(); }

  template <class T> inline T Trace(const GenMatrix<T>& m)
  { return m.Trace(); }

  template <class T> inline RealType(T) Norm(const GenMatrix<T>& m)
  { return m.Norm(); }

  template <class T> inline RealType(T) NormSq(const GenMatrix<T>& m)
  { return m.NormSq(); }

  template <class T> inline RealType(T) NormF(const GenMatrix<T>& m)
  { return m.NormF(); }

  template <class T> RealType(T) Norm1(const GenMatrix<T>& m)
  { return m.Norm1(); }

  template <class T> inline RealType(T) Norm2(const GenMatrix<T>& m)
  { return m.Norm2(); }

  template <class T> inline RealType(T) NormInf(const GenMatrix<T>& m)
  { return m.NormInf(); }

  template <class T> RealType(T) MaxAbsElement(const GenMatrix<T>& m)
  { return m.MaxAbsElement(); }

  template <class T> inline Matrix<T,ColMajor> Inverse(const GenMatrix<T>& m)
  { return m.Inverse(); }

  template <class T> inline Matrix<T,ColMajor> InverseATA(const GenMatrix<T>& m)
  { return m.InverseATA(); }

  template <class T> inline ConstMatrixView<T> Transpose(const GenMatrix<T>& m)
  { return m.Transpose(); }

  template <class T> inline MatrixView<T> Transpose(const MatrixView<T>& m)
  { return m.Transpose(); }

  template <class T, StorageType S> inline MatrixView<T> Transpose(
      Matrix<T,S>& m)
  { return m.Transpose(); }

  template <class T> inline ConstMatrixView<T> Conjugate(const GenMatrix<T>& m)
  { return m.Conjugate(); }

  template <class T> inline MatrixView<T> Conjugate(const MatrixView<T>& m)
  { return m.Conjugate(); }

  template <class T, StorageType S> inline MatrixView<T> Conjugate(
      Matrix<T,S>& m)
  { return m.Conjugate(); }

  template <class T> inline ConstMatrixView<T> Adjoint(const GenMatrix<T>& m)
  { return m.Adjoint(); }

  template <class T> inline MatrixView<T> Adjoint(const MatrixView<T>& m)
  { return m.Adjoint(); }

  template <class T, StorageType S> inline MatrixView<T> Adjoint(
      Matrix<T,S>& m)
  { return m.Adjoint(); }

  //
  // Matrix ==, != Matrix
  //

  template <class T> bool operator==(const GenMatrix<T>& m1,
      const GenMatrix<T>& m2);
  template <class T> inline bool operator!=(const GenMatrix<T>& m1,
      const GenMatrix<T>& m2)
  { return !(m1 == m2); }


  //
  // I/O
  //
 
  template <class T, StorageType S> istream& operator>>(istream& is, 
      Matrix<T,S>* m);

  template <class T> istream& operator>>(istream& is, const MatrixView<T>& m);

  template <class T, StorageType S> inline istream& operator>>(istream& is, 
      Matrix<T,S>& m)
  { return is>>m.QuickView(); }

  inline std::string Text(StorageType s)
  { 
    return s == RowMajor ? "RowMajor" :
      s == ColMajor ? "ColMajor" :
      s == DiagMajor ? "DiagMajor" :
      s == NoMajor ? "NoMajor" : "Unknown";
  }

  template <class T, StorageType S> inline std::string Type(
      const Matrix<T,S>& m)
  { 
    return std::string("Matrix<")+Type(T())+","+Text(S)+">"; 
  }
  template <class T> inline std::string Type(const GenMatrix<T>& m)
  {
    return std::string("GenMatrix<")+Type(T())+","+Text(m.ct())+","
      + Text(m.stor())+">"; 
  }
  template <class T> inline std::string Type(const ConstMatrixView<T>& m)
  { 
    return std::string("ConstMatrixView<")+Type(T())+","+Text(m.ct())+","
      + Text(m.stor())+">"; 
  }
  template <class T> inline std::string Type(const MatrixView<T>& m)
  { 
    return std::string("MatrixView<")+Type(T())+","+Text(m.ct())+","
      + Text(m.stor())+">"; 
  }

}; // namespace tmv

#endif
