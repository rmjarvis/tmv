///////////////////////////////////////////////////////////////////////////////
// vim:et:ts=2:sw=2:ci:cino=f0,g0,t0,+0:
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2009                                                 //
//                                                                           //
// The project is hosted at http://sourceforge.net/projects/tmv-cpp/         //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis@users.sourceforge.net                                         //
//                                                                           //
// This program is free software; you can redistribute it and/or             //
// modify it under the terms of the GNU General Public License               //
// as published by the Free Software Foundation; either version 2            //
// of the License, or (at your option) any later version.                    //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program in the file LICENSE.                              //
//                                                                           //
// If not, write to:                                                         //
// The Free Software Foundation, Inc.                                        //
// 51 Franklin Street, Fifth Floor,                                          //
// Boston, MA  02110-1301, USA.                                              //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------------
//
// This file defines the TMV Matrix class.
//
// A Matrix is a mathematical matrix whose sizes are not necessarily
// known at compile time.  (c.f. SmallMatrix in TMV_SmallMatrix.h.)
// Note, that this means that multiplication, for example, is real
// matrix multiplication, not multiplication of the corresponding elements.
// The other arithmetic operators are similarly correct for a mathematical
// matrix, but not necessarily for a simple 2-d array of numbers.
//
// The Matrix class and all associated functions are contained
// in the namespace tmv.  Alse, the Matrix class is a template, 
// so for a matrix of doubles, one would write 
// tmv::Matrix<double>.  
//
// An optional second template parameter, which can be 
// either RowMajor or ColMajor, sets which order to use for the
// storage of the matrix elements.  If omitted, ColMajor is assumed.
//
// Finally, an optional third template parameter, either CStyle or
// FortranStyle, indicate whether to use C-style or Fortran-style indexing.  
// With C-style (the default), the upper left corner of an MxN matrix is 
// m(0,0), and the lower right is m(M-1,N-1).  
// With Fortran-style, these are m(1,1) and m(M,N) respectively.  
// Also, when a function takes an index range, i1,i2, 
// then with C-style, this means elements from i1...i2-1 inclusive. 
// With Fortran-style, this means i1..i2 inclusive.
//
// If the third template parameter is omitted, CStyle is assumed.
// Note that if you want to specify Fortran-style, then you must also 
// specify the storage pattern, since that parameter precedes it.  For example:
// Matrix<T,RowMajor,FortranStyle> m(10,10);
//
// The return type of many of the below methods are given generically
// as value_type, real_type, iterator, etc.
// All of these types are typedefs in Matrix.  
// So, you could write, for example:
//
// tmv::Matrix<double> m(30,20);
// typename tmv::Matrix<double>::row_type row10 = m.row(10);
// typename tmv::Matrix<double>::rows_type lastrows = m.Rows(20,30);
//
//
// Constructors:
//
//    Matrix<T,stor,I>(size_t colsize, size_t rowsize)
//        Makes a matrix with column size = colsize and row size = rowsize
//        with _uninitialized_ values
//
//    Matrix<T,stor,I>(size_t colsize, size_t rowsize, T x)
//        Makes a matrix of size n with all values = x
//
//    Matrix<T,stor,I>(const vector<vector<T> >& m)
//        Makes a matrix with a_ij = m[i][j]
//
//    Matrix<T,stor,I>(size_t colsize, size_t rowsize, const T* vv)
//    Matrix<T,stor,I>(size_t colsize, size_t rowsize, 
//            const std::vector<T>& vv)
//        Make a matrix which copies the elements of vv.
//        If stor is tmv::RowMajor then the elements are taken in row major
//        order (m00,m01,..m0n,m10,m11...).  If stor is tmv::ColMajor
//        then the elements are taken in column major order.
//        If stor is omitted, then tmv::ColMajor is assumed.
//        (stor is also an optional last parameter on the other above 
//        constructors as well.)
//
//
// Special Creators:
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
//        Return the dimensions of the matrix
//
//    value_type operator()(int i, int j) const
//    value_type cref(int i, int j) const
//        Return the (i,j) element of the matrix
//        The first one respects the index-style of the underlying matrix.
//        The second, cref, always uses CStyle indexing and does not 
//        do any checking of the valididty of i,j.
//
//    reference operator()(int i, int j)
//    reference ref(int i, int j)
//        Return a reference to the (i,j) element of the matrix
//        The first one respects the index-style of the underlying matrix.
//        The second, ref, always uses CStyle indexing and does not 
//        do any checking of the valididty of i,j.
//        
//    row_type operator[](int i)
//    const_row_type operator[](int i) const
//        Return the ith row of the matrix as a vector.
//        This allows for m[i][j] style access, which shouldn't be
//        any slower than m(i,j) style access.
//
//    row_type row(int i)
//    const_row_type row(int i) const
//    row_range_type row(int i, int j1, int j2)
//    const_row_range_type row(int i, int j1, int j2) const
//        Return the ith row of the matrix as a vector
//        If j1,j2 are given, it returns the SubVector from j1 to j2 
//        (not including j2) within the row.
//
//    col_type col(int j)
//    const_col_type col(int j) const
//    col_range_type col(int j, int i1, int i2)
//    const_col_range_type col(int j, int j1, int j2) const
//        Return the jth column of the matrix as a Vector
//        If i1,i2 are given, it returns the SubVector from i1 to i2 
//        (not including i2) within the column.
//
//    diag_type diag()
//    const_diag_type diag() const
//        Return the diagonal of the matrix as a Vector
//
//    diag_range_type diag(int i)
//    diag_range_type diag(int i, int j1, int j2)
//    const_diag_range_type diag(i) const
//    const_diag_range_type diag(i, j1, j2) const
//        Return the super- or sub-diagonal i
//        If i > 0 return the super diagonal starting at m_0i
//        If i < 0 return the sub diagonal starting at m_|i|0
//        If j1,j2 are given, it returns the diagonal SubVector 
//        either from m_j1,i+j1 to m_j2,i+j2 (for i>0) 
//        or from m_|i|+j1,j1 to m_|i|+j2,j2 (for i<0)
//
// Functions of Matrices:
//
//    Most of these are both member functions and functions of a matrix,
//    so Norm(m) and m.Norm(), for example, are equivalent.
//
//    value_type m.Det() const    or Det(m)
//        Returns the determinant of a matrix.
//        Note: If the matrix is not square, the determinant is not
//              well defined.  The returned value is such that
//              conj(det) * det = Det(Adjoint(m) * m)
//              So for real nonsquare matrices, the sign is arbitrary,
//              and for complex nonsquare matrices, it is multiplied
//              by an arbitrary phase.
//
//    real_type m.LogDet(value_type* sign=0) const   or LogDet(m,sign)
//        Returns the logarithm of the absolute value of the determinant.
//        For many large matrices, the determinant yields to overflow.
//        Hence, this function is provided, which stably calculates the
//        natural logarithm of the absolute value of the determinant.
//        The optional sign argument returns the sign of the determinant
//        if value_type is real, or the factor exp(it) by which exp(logdet) 
//        would have to be multiplied to get the actual determinant.
//
//    value_type m.Trace() const    or Trace(m)
//        Returns the trace of a matrix.
//        = sum_i ( a_ii )
//
//    real_type m.Norm() const    or Norm(m)
//    real_type m.NormF() const    or NormF(m)
//        Return the Frobenius norm of a matrix.
//        = sqrt( sum_ij |a_ij|^2 )
//
//    real_type m.NormSq() const    or NormSq()
//    real_type m.NormSq(real_type scale) const
//        Returns the square of Norm().
//        In the method version, you can provide an optional scale, in
//        which case the output is equal to NormSq(scale*m).
//
//    real_type m.Norm1() const    or Norm1(m)
//        Returns the 1-norm of a matrix.
//        = max_j (sum_i |a_ij|)
//
//    real_type m.Norm2() const    or Norm2(m)
//        Returns the 2-norm of a matrix.
//        = sqrt( Max Eigenvalue of (A.Adjoint * A) )
//        = Maximum singular value
//        Note: This norm is costly to calculate if one is not 
//              otherwise doing a singular value decomposition
//              of the matrix.
//
//    real_type m.NormInf() const    or NormInf(m)
//        Returns the infinity-norm of a matrix.
//        = max_i (sum_j |a_ij|)
//
//    value_type SumElements() const    or SumElements(m) 
//        Returns the sum of all elements in the matrix.
//
//    real_type SumAbsElements() const    or SumAbsElements(m) 
//        Returns the sum of absolute values of elements in the matrix.
//
//    value_type MaxElement() const    or MaxElement(m) 
//        Returns the maximum value of any element in the matrix.
//        As "max" doesn't make sense for complex values, for these
//        we use just the real components.
//
//    value_type MinElement() const    or MinElement(m) 
//        Returns the minimum value of any element in the matrix.
//
//    real_type MaxAbsElement() const    or MaxAbsElement(m) 
//        The same as MaxElement, except absolute values are used.
//
//    real_type MinAbsElement() const    or MinAbsElement(m) 
//        The same as MinElement, except absolute values are used.
//
//    void m.Inverse(minv) const
//        Sets minv to the inverse of m if it exists.
//        If m is singular and square, and LU is set for dividing
//          (LU is default for square matrices)
//          then a run-time error will occur.
//        If m is singular or not square and SV is set 
//          then the returned matrix is the pseudo-inverse which satisfies:
//          MXM = M
//          XMX = X
//          (MX)T = MX
//          (XM)T = XM
//        [If dividing using QR or QRP, all but the last of these will 
//         be true.]
//
//    void m.InverseATA(invata) const
//        Sets invata to the Inverse of (A.Adjoint * A) for matrix m = A
//        If Ax=b is solved for x, then (AT A)^-1 is the 
//        covariance matrix of the least-squares solution x
//
//    inverse_type m.Inverse() const    or Inverse(m)
//        Returns an auxiliary object that delays the calculation of the 
//        inverse until there is appropriate storage for it.
//        m.Inverse(minv) is equivalent to minv = m.Inverse().
//
//
// Modifying Functions
//
//    type& Zero()
//        Sets all elements to 0
//
//    type& SetAllTo(value_type x)
//        Sets all elements to x
//
//    type& AddToAll(value_type x)
//        Add x to each element of m
//
//    type& Clip(real_type thresh)
//        Set to 0 all elements whose abolute value is < thresh
//
//    type& ConjugateSelf()
//        Sets all elements to its conjugate
//
//    type& TransposeSelf() 
//        Transposes the elements of a square matrix
//
//    type& SetToIdentity(value_type x = 1)
//        Set to the identity matrix, or 
//        with a parameter, set to x times identity matrix
//
//    type& SwapRows(int i1, int i2)
//        Swap two rows
//
//    type& SwapCols(int j1, int j2)
//        Swap two columns
//
//    type& PermuteRows(const int* p)
//    type& PermuteRows(const int* p, int i1, int i2)
//        Perform a series of row swaps (0,p[0]), (1,p[1])...
//        In the second case, only do (i1,p[i1)...(i2-1,p[i2-1])
//    type& ReversePermuteRows(const int* p)
//    type& ReversePermuteRows(const int* p, int i1, int i2)
//        The same, but perform the swaps in reverse order
//
//    type& PermuteCols(const int* p)
//    type& PermuteCols(const int* p, int j1, int j2)
//        Perform a series of column swaps (0,p[0]), (1,p[1])...
//        In the second case, only do (j1,p[j1)...(j2-1,p[j2-1])
//    type& ReversePermuteCols(const int* p)
//    type& ReversePermuteCols(const int* p, int j1, int j2)
//        The same, but perform the swaps in reverse order
//
//    Swap(Matrix& m1, Matrix& m2)
//        Swap the values of two matrices
//        The matrices must be the same size
//
// MatrixView:
//
//    A MatrixView object refers to the elements of a regular matrix
//    so that altering the elements in the view alters the corresponding
//    elements in the original object.  
//
//    As for VectorView, we have both mutable and non-mutable views,
//    called MatrixView and ConstMatrixView respectively.
//
//    Both of these take several template arguments:
//    T = the underlying data type.
//    Si = (int) the step size along a column, if known. 
//    Sj = (int) the step size along a row, if known. 
//    C = (bool) whether the view is the conjugate of the data.
//    I = (CStyle or FortranStyle) which indexing style to use.
// 
//    There are default arguments for all but T:
//    Si = Sj = UNKNOWN, C = false, I = CStyle
//    so you can omit these, although you will generally get better results
//    by specifying the step sizes, since these are usually known.
//
//    ConstMatrixView<T,Si,Sj,C,I>(const T* p, size_t m, size_t n, 
//            int si=Si, int sj=Sj)
//    MatrixView<T,Si,Sj,C,I>(T* p, size_t m, size_t n, 
//            int si=Si, int sj=Sj)
//        Make a matrix view starting at memory location p, 
//        with m rows and n cols, stepping over the data with step sizes 
//        si and sj along the columns and rows respectively.
//        Note, that si,sj are generally omitted if you are specifying Si,Sj
//        as a template argument.  But if either S == UNKNOWN, then you would
//        need to specify that s as an argument.
//
//    There are also copy constructors that change Si,Sj (e.g. from UNKNOWN
//    to some compile-time-known value), or change MatrixView to 
//    ConstMatrixView.  e.g.:
//
//    ConstMatrixView<T,Si,Sj,C,I>(const MatrixView<T,Si,Sj,C,I>& m);
//    ConstMatrixView<T,Si,Sj,C,I>(const MatrixView<T,Si,UNKNOWN,C,I>& m);
//    MatrixView<T,Si,Sj,C,I>(const MatrixView<T,UNKNOWN,Sj,C,I>& m);
//
// Views:
//
//    All of these methods return some kind of view into the matrix
//    data.  The return types are either ConstMatrixView or MatrixView
//    with the appropriate values for the template parameters Si, Sj, and C.
//    (Except for SubVector and LinearView, which return either 
//    ConstVectorView or VectorView.)
//    The template parameter I is preserved from the original matrix.
//
//    All of these have a const and non-const version.  For the const
//    versions, prepend const_ to the return type name.  
//    e.g. const_submatrix_type.  
//
//    submatrix_type SubMatrix(int i1, int i2, int j1, int j2,
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
//           A.SubMatrix(0,2,2,4).SetToIdentity();
//           A.SubMatrix(2,4,0,2).SetAllTo(2);
//        The substep values allow you to space the elements of 
//        the submatrix at steps larger than 1.
//        eg. To make an 8x8 checkerboard of 1's and 0's, you could write:
//           Matrix<int> A(8,8,0);
//           A.SubMatrix(0,8,1,9,2,2) = 1;
//           A.SubMatrix(1,9,0,8,2,2) = 1;
//        Note that in this case the i2,j2 values need to be the 
//        "one past the end" value given the step size, so 9 here when
//        starting at 1.
//
//    subvector_type SubVector(int i, int j, int istep, int jstep, int size)
//        Returns a subvector which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.
//        For example, the cross-diagonal from the lower left to the upper
//        right of a 6x6 matrix could be accessed using:
//        m.SubVector(5,0,-1,1,6)
//
//    colpair_type ColPair(int j1, int j2)
//        This returns an Mx2 submatrix which consists of the 
//        columns j1 and j2.  This is useful for multiplying two 
//        (not necessarily adjacent) columns of a matrix by a 2x2 matrix.
//
//    rowpair_type RowPair(int i1, int i2)
//        Same as above, but two rows.
//
//    cols_type Cols(int j1, int j2)
//        This is shorthand for SubMatrix(0,m.colsize(),j1,j2)
//
//    rows_type Rows(int i1, int i2)
//        This is shorthand for SubMatrix(i1,i2,0,m.rowsize())
//
//    realview_type Real()
//    imagview_type Imag()
//        Returns a view to the real/imag elements of a complex matrix.
//
//    view_type View()
//        Returns a view of a matrix.
//
//    conjugate_type Conjugate()
//        Returns a view to the conjugate of a matrix.
//
//    transpose_type Transpose()
//        Returns a view to the transpose of a matrix.
//
//    adjoint_type Adjoint()
//        Returns a view to the adjoint (conjugate transpose) of a matrix.
//        Note: Some people define the adjoint as the cofactor matrix.
//              This is not the same as our definition of the Adjoint.
//
//    uppertri_type UpperTri()
//        Returns a view to the upper triangle portion of the matrix.
//
//    unit_uppertri_type UnitUpperTri()
//        Returns a view to the upper triangle portion of the matrix
//        where the diagonal elements are all taken to be equal to 1.
//
//    lowertri_type LowerTri()
//        Returns a view to the lower triangle portion of the matrix.
//
//    unit_lowertri_type UnitLowerTri()
//        Returns a view to the lower triangle portion of the matrix
//        where the diagonal elements are all taken to be equal to 1.
//
//    linearview_type LinearView()
//        Returns a VectorView with all the elements of the matrix.
//        This is mostly used internally for things like MaxElement
//        and ConjugateSelf, where the matrix structure is irrelevant,
//        and we just want to do something to all the elements.
//        The correlary function CanLinearize() returns whether this is 
//        allowed.
//
//    nonconj_type NonConj()
//        Returns a view of the underlying memory elements, removing
//        any mconj that might be set in the current view.
//        This is sometimes useful when an operation has the same
//        effect regardless of mconj, so it is better to just ignore
//        any mconj value.
//
//
//
// Operators:
//    Here we use m for a matrix, v for a vector and x for a scalar.
//
//    You can also mix real and complex Vectors of the same
//    underlying type.  eg. Matrix<double> and Matrix<complex<double> >
//
//    -m
//
//    m = m
//
//    m += m
//    m + m
//    m += x
//    m + x
//    x + m
//
//    m -= m
//    m - m
//    m -= x
//    m - x
//    x - m
//
//    m *= x
//    m *= m
//    m * x
//    x * m
//    m * v
//    v * m
//    v *= m
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
//    Most of these behave in the logical way for dealing with matrices.
//    Two comments about behavior that might not be obvious:
//
//    1) Vectors are either row or column Vectors as appropriate.
//       eg. For m*v, v is a column Vector, but for v*m, v is a row Vector
//
//    2) Sometimes x should be thought of a x*I, where I is an appropriately
//       sized identity matrix.  For example m-1 means m-I,
//       and m += 3 means m += 3*I (or equivalently, m.diag().AddToAll(3)).
//
//    3) Division by a matrix can be from either the left or the right.
//       ie. v/m can mean either m^-1 v or v m^-1
//       The first case is the solution of mx=v, the second is of xm = v.
//       Since the first case is the way one generally poses a problem
//       for solving a set of equations, we take v/m to be left-division.
//       If you want right-division (v m^-1), then we supply the % operator
//       to do so.  
//       ie. v%m means v m^-1
//       If you want to be more explicit, you can write:
//       v/m as m.Inverse() * v and v%m as v * m.Inverse().
//       In all cases, the actual calculation is delayed until there is
//       storage to put it.  (Unless you string too many calculations 
//       together, in which case it will use a temporary.)
//
//
// I/O: 
//
//    void Write(std::ostream& os) const    or os << m
//        Writes m to ostream os in the following format:
//          colsize rowsize
//          ( m(0,0) m(0,1) m(0,2) ... m(0,rowsize) )
//          ( m(1,0) m(1,1) m(1,2) ... m(1,rowsize) )
//          ...
//          ( m(colsize,0)  ...  m(colsize,rowsize) )
//
//
//    m.Write(ostream& os, real_type thresh)
//        Write m to os as above, but if |m(i,j)| < thresh, write 0 instead
//
//    void Read(std::istream& is)    or is >> m
//        Reads the matrix from istream is in the same format
//        Note: the matrix must already be the correct size
//
//    std::auto_ptr<tmv::Matrix<T> > mptr;
//    is >> mptr
//        If you do not know the size of the matrix to be read in, you can
//        use this form which will allocate the matrix to be the correct
//        size according to the input data.
//
//
// Division Control Functions:
//
//    There are a number of algorithms available for dividing
//    matrices.  We provide functions to allow you to 
//    change the algorithm used by the code on the fly.
//    In particular, you can write:
//
//    m.DivideUsing(dt)
//    where dt is LU, QR, QRP, or SV
//
//    Each of these also has an in-place version whcih overwrites the
//    current matrix memory with the decomposition needed for 
//    doing the division.  Obviously, if you try to use the matrix
//    after doing SetDiv (implicitly or explicitly), the results will
//    be invalid.
//
//    The default method is LU (LU Decomposition) for square
//    matrices and QR (QR Decomposition) for non-square.
//   
//    See the header comment to TMV_Divider.h for more info about
//    the different algorithms.
//
//    There are also shorthands for accessing the decomposition.
//    If dt = LU, then LUD() returns the LUDiv<T> class, performing
//    the decomposition if it hasn't be done yet.
//
//    Likewise:
//    QRD(), QRPD(), SVD() return the corresponding Divider classes for
//    those algorithms.
//


#ifndef TMV_Matrix_H
#define TMV_Matrix_H

#include "tmv/TMV_BaseMatrix_Rec.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_BaseMatrix_Tri.h"
//#include "tmv/TMV_DivHelper.h"
#include <vector>

namespace tmv {

  template <class T, StorageType S, IndexStyle I>
  struct Traits<Matrix<T,S,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef Matrix<T,S,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef type copy_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { mshape = Rec };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (S == RowMajor) };
    enum { mcolmajor = (S == ColMajor) };
    enum { mstor = S };
    enum { mstepi = (S==ColMajor ? 1 : UNKNOWN) };
    enum { mstepj = (S==RowMajor ? 1 : UNKNOWN) };
    enum { mdiagstep = UNKNOWN };
    enum { mconj = false };
    enum { mcanlin = true };
    enum { twoSi = misreal ? int(mstepi) : int(IntTraits<mstepi>::twoS) };
    enum { twoSj = misreal ? int(mstepj) : int(IntTraits<mstepj>::twoS) };
    enum { notC = miscomplex };

    typedef ConstVectorView<T,mstepi,false,I> const_col_type;
    typedef ConstVectorView<T,mstepi,false,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,false,I> const_row_type;
    typedef ConstVectorView<T,mstepj,false,I> const_row_range_type;
    typedef ConstVectorView<T,mdiagstep,false,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,false,I> const_diag_range_type;

    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;
    typedef ConstSmallMatrixView<T,UNKNOWN,2,mstepi,UNKNOWN,false,I> const_colpair_type;
    typedef ConstSmallMatrixView<T,2,UNKNOWN,UNKNOWN,mstepj,false,I> const_rowpair_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_cols_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_rows_type;

    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_view_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,CStyle> const_cview_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,FortranStyle> const_fview_type;
    typedef ConstMatrixView<T> const_xview_type;
    typedef ConstMatrixView<T,1,mstepj,false,I> const_cmview_type;
    typedef ConstMatrixView<T,mstepi,1,false,I> const_rmview_type;
    typedef ConstMatrixView<T,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstMatrixView<T,mstepj,mstepi,false,I> const_transpose_type;
    typedef ConstMatrixView<T,mstepj,mstepi,notC,I> const_adjoint_type;
    typedef ConstUpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> const_uppertri_type;
    typedef ConstUpperTriMatrixView<T,UnitDiag,mstepi,mstepj,false,I> const_unit_uppertri_type;
    typedef ConstLowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> const_lowertri_type;
    typedef ConstLowerTriMatrixView<T,UnitDiag,mstepi,mstepj,false,I> const_unit_lowertri_type;
    typedef ConstMatrixView<real_type,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstVectorView<T,1,false,I> const_linearview_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_nonconj_type;
    typedef MatrixView<T,mstepi,mstepj,false,I> nonconst_type;

    typedef T& reference;

    typedef VectorView<T,mstepi,false,I> col_type;
    typedef VectorView<T,mstepi,false,I> col_range_type;
    typedef VectorView<T,mstepj,false,I> row_type;
    typedef VectorView<T,mstepj,false,I> row_range_type;
    typedef VectorView<T,mdiagstep,false,I> diag_type;
    typedef VectorView<T,mdiagstep,false,I> diag_range_type;

    typedef MatrixView<T,mstepi,mstepj,false,I> submatrix_type;
    typedef MatrixView<T,UNKNOWN,UNKNOWN,false,I> submatrix_step_type;
    typedef VectorView<T,UNKNOWN,false,I> subvector_type;
    typedef SmallMatrixView<T,UNKNOWN,2,mstepi,UNKNOWN,false,I> colpair_type;
    typedef SmallMatrixView<T,2,UNKNOWN,UNKNOWN,mstepj,false,I> rowpair_type;
    typedef MatrixView<T,mstepi,mstepj,false,I> cols_type;
    typedef MatrixView<T,mstepi,mstepj,false,I> rows_type;

    typedef MatrixView<T,mstepi,mstepj,false,I> view_type;
    typedef MatrixView<T,mstepi,mstepj,false,CStyle> cview_type;
    typedef MatrixView<T,mstepi,mstepj,false,FortranStyle> fview_type;
    typedef MatrixView<T> xview_type;
    typedef MatrixView<T,1,mstepj,false,I> cmview_type;
    typedef MatrixView<T,mstepi,1,false,I> rmview_type;
    typedef MatrixView<T,mstepi,mstepj,notC,I> conjugate_type;
    typedef MatrixView<T,mstepj,mstepi,false,I> transpose_type;
    typedef MatrixView<T,mstepj,mstepi,notC,I> adjoint_type;
    typedef UpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> uppertri_type;
    typedef UpperTriMatrixView<T,UnitDiag,mstepi,mstepj,false,I> unit_uppertri_type;
    typedef LowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,false,I> lowertri_type;
    typedef LowerTriMatrixView<T,UnitDiag,mstepi,mstepj,false,I> unit_lowertri_type;
    typedef MatrixView<real_type,twoSi,twoSj,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef VectorView<T,1,false,I> linearview_type;
    typedef MatrixView<T,mstepi,mstepj,false,I> nonconj_type;
  };

#ifdef XTEST
#ifdef TMV_DEBUG
#define XTEST_DEBUG
#endif
#endif

  template <class T, StorageType S, IndexStyle I>
  class Matrix :
    public BaseMatrix_Rec_Mutable<Matrix<T,S,I> >
    //public DivHelper<Matrix<T,S,I> >
  {
  public:

    typedef Matrix<T,S,I> type;
    typedef BaseMatrix_Rec_Mutable<type> base_mut;
    enum { misreal = Traits<type>::misreal };
    enum { miscomplex = Traits<type>::miscomplex };
    enum { mshape = Traits<type>::mshape };

#if 0
    typedef DivHelper<type> base_div;

    // 
    // Division Control
    //

    using DivHelper<T>::DivideInPlace;
    using DivHelper<T>::SaveDiv;
    using DivHelper<T>::SetDiv;
    using DivHelper<T>::UnSetDiv;
    using DivHelper<T>::ReSetDiv;
    using DivHelper<T>::DivIsSet;
    using DivHelper<T>::CheckDecomp;

    inline void DivideUsing(DivType dt) const
    {
      TMVAssert(dt == LU || dt == QR || dt == QRP || dt == SV);
      DivHelper<T>::DivideUsing(dt); 
    }

    inline const LUDiv<T>& LUD() const
    {
      DivideUsing(LU);
      SetDiv();
      TMVAssert(GetDiv());
      TMVAssert(dynamic_cast<const LUDiv<T>*>(GetDiv()));
      return *dynamic_cast<const LUDiv<T>*>(GetDiv());
    }

    inline const QRDiv<T>& QRD() const
    {
      DivideUsing(QR);
      SetDiv();
      TMVAssert(GetDiv());
      TMVAssert(dynamic_cast<const QRDiv<T>*>(GetDiv()));
      return *dynamic_cast<const QRDiv<T>*>(GetDiv());
    }

    inline const QRPDiv<T>& QRPD() const
    {
      DivideUsing(QRP);
      SetDiv();
      TMVAssert(GetDiv());
      TMVAssert(dynamic_cast<const QRPDiv<T>*>(GetDiv()));
      return *dynamic_cast<const QRPDiv<T>*>(GetDiv());
    }

    inline const SVDiv<T>& SVD() const
    {
      DivideUsing(SV);
      SetDiv();
      TMVAssert(GetDiv());
      TMVAssert(dynamic_cast<const SVDiv<T>*>(GetDiv()));
      return *dynamic_cast<const SVDiv<T>*>(GetDiv());
    }

    template <class T1> 
    inline void LDivEq(const VectorView<T1>& v) const 
    { DivHelper<T>::LDivEq(v); }

    template <class T1> 
    inline void LDivEq(const MatrixView<T1>& m) const 
    { DivHelper<T>::LDivEq(m); }

    template <class T1> 
    inline void RDivEq(const VectorView<T1>& v) const 
    { DivHelper<T>::RDivEq(v); }

    template <class T1> 
    inline void RDivEq(const MatrixView<T1>& m) const 
    { DivHelper<T>::RDivEq(m); }

    template <class T1, class T0> 
    inline void LDiv(const GenVector<T1>& v1, const VectorView<T0>& v0) const
    { DivHelper<T>::LDiv(v1,v0); }

    template <class T1, class T0> 
    inline void LDiv(const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
    { DivHelper<T>::LDiv(m1,m0); }

    template <class T1, class T0> 
    inline void RDiv(const GenVector<T1>& v1, const VectorView<T0>& v0) const
    { DivHelper<T>::RDiv(v1,v0); }

    template <class T1, class T0> 
    inline void RDiv(const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
    { DivHelper<T>::RDiv(m1,m0); }
    inline const BaseMatrix<T>& GetMatrix() const { return *this; }
#endif


    //
    // Constructors
    //

    inline Matrix(size_t cs, size_t rs) :
      itscs(cs), itsrs(rs),
      linsize(cs*rs), itsm(new T[linsize])
    {
      TMVStaticAssert(S==RowMajor || S==ColMajor);
#ifdef TMV_DEBUG
      this->SetAllTo(T(888));
#endif
    }

    inline Matrix(size_t cs, size_t rs, T x) :
      itscs(cs), itsrs(rs),
      linsize(cs*rs), itsm(new T[linsize])
    {
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      this->SetAllTo(x);
    }

    inline Matrix(size_t cs, size_t rs, const T* vv) :
      itscs(cs), itsrs(rs),
      linsize(cs*rs), itsm(new T[linsize])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      std::copy(vv,vv+linsize,itsm.get());
    }

    inline Matrix(size_t cs, size_t rs, const std::vector<T>& vv) : 
      itscs(cs), itsrs(rs),
      linsize(cs*rs), itsm(new T[linsize])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVAssert(vv.size() == linsize);
      std::copy(vv.begin(),vv.end(),itsm.get());
    }

    inline Matrix(const std::vector<std::vector<T> >& vv) :
      itscs(vv.size())
    {
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      if (itscs == 0) {
        itsrs = 0;
        linsize = 0;
        itsm.reset(new T[0]);
      } else {
        itsrs = vv[0].size();
        linsize = itscs * itsrs;
        itsm.reset(new T[linsize]);
        for(int i=0;i<itscs;++i) {
          TMVAssert(vv[i].size() == itsrs);
          std::copy(vv[i].begin(),vv[i].end(),this->get_row(i).begin());
        }
      }
    }

    inline Matrix(const type& m2) :
      itscs(m2.itscs), itsrs(m2.itsrs),
      linsize(m2.linsize), itsm(new T[linsize])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      m2.AssignTo(*this);
    }

    template <class M2>
    inline Matrix(const BaseMatrix<M2>& m2) :
      itscs(m2.colsize()), itsrs(m2.rowsize()),
      linsize(itscs * itsrs), itsm(new T[linsize])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVStaticAssert(M2::misreal || miscomplex);
      TMVStaticAssert((ShapeTraits2<M2::mshape,mshape>::assignable));
      m2.AssignTo(*this);
    }

    template <class M2>
    inline Matrix(const BaseMatrix_Rec<M2>& m2) :
      itscs(m2.colsize()), itsrs(m2.rowsize()),
      linsize(m2.ls()), itsm(new T[linsize])
    {
#ifdef XTEST_DEBUG
      this->SetAllTo(T(888));
#endif
      TMVStaticAssert(S==RowMajor || S==ColMajor);
      TMVStaticAssert(M2::misreal || miscomplex);
      m2.AssignTo(*this);
    }

    inline ~Matrix() 
    {
#ifdef TMV_DEBUG
      this->SetAllTo(T(999));
#endif
    }

    //
    // Op=
    //

    inline type& operator=(const type& m2)
    {
      if (&m2 != this)
        base_mut::operator=(m2);
      return *this; 
    }

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    {
      base_mut::operator=(m2);
      return *this; 
    }

    inline type& operator=(T x)
    {
      base_mut::operator=(x);
      return *this; 
    }


    //
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm.get(); }
    inline T* ptr() { return itsm.get(); }

    inline T cref(int i, int j) const
    { return itsm[S==RowMajor ? i*stepi()+j : i+j*stepj()]; }

    inline T& ref(int i, int j)
    { return itsm[S==RowMajor ? i*stepi()+j : i+j*stepj()]; }

    inline size_t ls() const { return linsize; }
    inline size_t colsize() const { return itscs; }
    inline size_t rowsize() const { return itsrs; }
    inline int stepi() const { return S == RowMajor ? itsrs : 1; }
    inline int stepj() const { return S == RowMajor ? 1 : itscs; }
    inline bool isconj() const { return false; }
    inline bool isrm() const { return S == RowMajor; }
    inline bool iscm() const { return S == ColMajor; }
    inline StorageType stor() const { return S; }

  private:

    const size_t itscs;
    const size_t itsrs;
    const size_t linsize;
    auto_array<T> itsm;

  }; // Matrix

  template <class T, StorageType S>
  class MatrixF :
    public Matrix<T,S,FortranStyle>
  {
  public:
    typedef MatrixF<T,S> type;
    typedef Matrix<T,S,FortranStyle> mtype;

    inline MatrixF(size_t cs, size_t rs) : mtype(cs,rs) {}
    inline MatrixF(size_t cs, size_t rs, T x) : mtype(cs,rs,x) {}
    inline MatrixF(size_t cs, size_t rs, const T* vv) : mtype(cs,rs,vv) {}
    inline MatrixF(size_t cs, size_t rs, const std::vector<T>& vv) : 
      mtype(cs,rs,vv) {}
    inline MatrixF(const std::vector<std::vector<T> >& vv) : mtype(vv) {}
    inline MatrixF(const type& m2) : mtype(m2) {}
    template <class M2>
    inline MatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
    inline ~MatrixF() {}

    inline type& operator=(const type& m2)
    { mtype::operator=(m2); return *this; }
    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    { mtype::operator=(m2); return *this; }
    inline type& operator=(T x)
    { mtype::operator=(x); return *this; }

  }; // MatrixF

  template <class T, int Si, int Sj, bool C, IndexStyle I>
  struct Traits<ConstMatrixView<T,Si,Sj,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef ConstMatrixView<T,Si,Sj,C,I> type;
    typedef const type& calc_type;
    typedef const type& eval_type;
    typedef Matrix<T,Sj==1?RowMajor:ColMajor,I> copy_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { mshape = Rec };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (Sj == 1) };
    enum { mcolmajor = (Si == 1) };
    enum { mstor = (mrowmajor ? RowMajor : ColMajor) };
    enum { mstepi = Si };
    enum { mstepj = Sj };
    enum { mdiagstep = IntTraits2<Si,Sj>::sum };
    enum { mconj = C };
    enum { mcanlin = false };
    enum { twoSi = misreal ? Si : IntTraits<Si>::twoS };
    enum { twoSj = misreal ? Sj : IntTraits<Sj>::twoS };
    enum { notC = !C && miscomplex };

    typedef ConstVectorView<T,mstepi,C,I> const_col_type;
    typedef ConstVectorView<T,mstepi,C,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,C,I> const_row_type;
    typedef ConstVectorView<T,mstepj,C,I> const_row_range_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_range_type;

    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;
    typedef ConstSmallMatrixView<T,UNKNOWN,2,mstepi,UNKNOWN,C,I> const_colpair_type;
    typedef ConstSmallMatrixView<T,2,UNKNOWN,UNKNOWN,mstepj,C,I> const_rowpair_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_cols_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_rows_type;

    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_view_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,CStyle> const_cview_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,FortranStyle> const_fview_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C> const_xview_type;
    typedef ConstMatrixView<T,1,mstepj,C,I> const_cmview_type;
    typedef ConstMatrixView<T,mstepi,1,C,I> const_rmview_type;
    typedef ConstMatrixView<T,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstMatrixView<T,mstepj,mstepi,C,I> const_transpose_type;
    typedef ConstMatrixView<T,mstepj,mstepi,notC,I> const_adjoint_type;
    typedef ConstUpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> const_uppertri_type;
    typedef ConstUpperTriMatrixView<T,UnitDiag,mstepi,mstepj,C,I> const_unit_uppertri_type;
    typedef ConstLowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> const_lowertri_type;
    typedef ConstLowerTriMatrixView<T,UnitDiag,mstepi,mstepj,C,I> const_unit_lowertri_type;
    typedef ConstMatrixView<real_type,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstVectorView<T,1,C,I> const_linearview_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_nonconj_type;
    typedef MatrixView<T,mstepi,mstepj,C,I> nonconst_type;
  };

  template <class T, int Si, int Sj, bool C, IndexStyle I>
  class ConstMatrixView :
    public BaseMatrix_Rec<ConstMatrixView<T,Si,Sj,C,I> >
  {
  public:

    typedef ConstMatrixView<T,Si,Sj,C,I> type;
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };

    //
    // Constructors
    //

    inline ConstMatrixView(const T* m, size_t cs, size_t rs, int si, int sj) :
      itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) {}

    inline ConstMatrixView(const T* m, size_t cs, size_t rs, int si) :
      itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj) 
    { TMVStaticAssert(Sj != UNKNOWN); }

    inline ConstMatrixView(const T* m, size_t cs, size_t rs) :
      itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
    { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

    inline ConstMatrixView(const type& m2) :
      itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstMatrixView(const ConstMatrixView<T,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstMatrixView(const MatrixView<T,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstMatrixView(
        const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstMatrixView(
        const SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) :
      itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    inline ~ConstMatrixView() { 
#ifdef TMV_DEBUG
      itsm = 0; 
#endif
    }

  private :
    inline void operator=(const type& m2);
  public :

    //
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm; }

    inline T cref(int i, int j) const
    { return DoConj<C>(itsm[i*stepi()+j*stepj()]); }

    inline size_t ls() const { return itscs*itsrs; }
    inline size_t colsize() const { return itscs; }
    inline size_t rowsize() const { return itsrs; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline bool isconj() const { return C; }
    inline bool isrm() const 
    { return mrowmajor || (!mcolmajor &&  stepj() == 1); }
    inline bool iscm() const 
    { return mcolmajor || (!mrowmajor &&  stepi() == 1); }

  private :

#ifdef TMV_DEBUG
    const T* itsm;
#else
    const T*const itsm;
#endif
    const size_t itscs;
    const size_t itsrs;
    const StepInt<Si> itssi;
    const StepInt<Sj> itssj;

  }; // ConstMatrixView

  template <class T, int Si, int Sj, bool C>
  class ConstMatrixViewF :
    public ConstMatrixView<T,Si,Sj,C,FortranStyle>
  {
  public:

    typedef ConstMatrixViewF<T,Si,Sj,C> type;
    typedef ConstMatrixView<T,Si,Sj,C,FortranStyle> mtype;

    inline ConstMatrixViewF(const T* m, size_t cs, size_t rs, int si, int sj) :
      mtype(m,cs,rs,si,sj) {}
    inline ConstMatrixViewF(const T* m, size_t cs, size_t rs, int si) :
      mtype(m,cs,rs,si) {}
    inline ConstMatrixViewF(const T* m, size_t cs, size_t rs) : 
      mtype(m,cs,rs) {}
    inline ConstMatrixViewF(const type& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstMatrixViewF(const ConstMatrixView<T,Si2,Sj2,C,I2>& m2) :
      mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline ConstMatrixViewF(const MatrixView<T,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstMatrixViewF(
        const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
    inline ConstMatrixViewF(
        const SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) : mtype(m2) {}
    inline ~ConstMatrixViewF() {}

  private :
    inline void operator=(const type& m2);

  }; // ConstMatrixViewF


  template <class T, int Si, int Sj, bool C, IndexStyle I>
  struct Traits<MatrixView<T,Si,Sj,C,I> >
  {
    typedef T value_type;

    typedef typename Traits<T>::real_type real_type;
    typedef typename Traits<T>::complex_type complex_type;
    enum { misreal = Traits<T>::isreal };
    enum { miscomplex = Traits<T>::iscomplex };

    typedef MatrixView<T,Si,Sj,C,I> type;
    typedef const ConstMatrixView<T,Si,Sj,C,I> calc_type;
    typedef calc_type eval_type;
    typedef Matrix<T,Sj==1?RowMajor:ColMajor,I> copy_type;
    typedef InvalidType inverse_type;

    enum { mcolsize = UNKNOWN };
    enum { mrowsize = UNKNOWN };
    enum { mshape = Rec };
    enum { mfort = (I == FortranStyle) };
    enum { mcalc = true };
    enum { mrowmajor = (Sj == 1) };
    enum { mcolmajor = (Si == 1) };
    enum { mstor = (mrowmajor ? RowMajor : ColMajor) };
    enum { mstepi = Si };
    enum { mstepj = Sj };
    enum { mdiagstep = IntTraits2<Si,Sj>::sum };
    enum { mconj = C };
    enum { mcanlin = false };
    enum { twoSi = misreal ? Si : IntTraits<Si>::twoS };
    enum { twoSj = misreal ? Sj : IntTraits<Sj>::twoS };
    enum { notC = !C && miscomplex };

    typedef ConstVectorView<T,mstepi,C,I> const_col_type;
    typedef ConstVectorView<T,mstepi,C,I> const_col_range_type;
    typedef ConstVectorView<T,mstepj,C,I> const_row_type;
    typedef ConstVectorView<T,mstepj,C,I> const_row_range_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_type;
    typedef ConstVectorView<T,mdiagstep,C,I> const_diag_range_type;

    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_submatrix_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> const_submatrix_step_type;
    typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;
    typedef ConstSmallMatrixView<T,UNKNOWN,2,mstepi,UNKNOWN,C,I> const_colpair_type;
    typedef ConstSmallMatrixView<T,2,UNKNOWN,UNKNOWN,mstepj,C,I> const_rowpair_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_cols_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_rows_type;

    typedef ConstMatrixView<T,mstepi,mstepj,C,I> const_view_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,CStyle> const_cview_type;
    typedef ConstMatrixView<T,mstepi,mstepj,C,FortranStyle> const_fview_type;
    typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C> const_xview_type;
    typedef ConstMatrixView<T,1,mstepj,C,I> const_cmview_type;
    typedef ConstMatrixView<T,mstepi,1,C,I> const_rmview_type;
    typedef ConstMatrixView<T,mstepi,mstepj,notC,I> const_conjugate_type;
    typedef ConstMatrixView<T,mstepj,mstepi,C,I> const_transpose_type;
    typedef ConstMatrixView<T,mstepj,mstepi,notC,I> const_adjoint_type;
    typedef ConstUpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> const_uppertri_type;
    typedef ConstUpperTriMatrixView<T,UnitDiag,mstepi,mstepj,C,I> const_unit_uppertri_type;
    typedef ConstLowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> const_lowertri_type;
    typedef ConstLowerTriMatrixView<T,UnitDiag,mstepi,mstepj,C,I> const_unit_lowertri_type;
    typedef ConstMatrixView<real_type,twoSi,twoSj,false,I> const_realview_type;
    typedef const_realview_type const_imagview_type;
    typedef ConstVectorView<T,1,C,I> const_linearview_type;
    typedef ConstMatrixView<T,mstepi,mstepj,false,I> const_nonconj_type;
    typedef MatrixView<T,mstepi,mstepj,C,I> nonconst_type;

    typedef typename AuxRef<T,C>::reference reference;

    typedef VectorView<T,mstepi,C,I> col_type;
    typedef VectorView<T,mstepi,C,I> col_range_type;
    typedef VectorView<T,mstepj,C,I> row_type;
    typedef VectorView<T,mstepj,C,I> row_range_type;
    typedef VectorView<T,mdiagstep,C,I> diag_type;
    typedef VectorView<T,mdiagstep,C,I> diag_range_type;

    typedef MatrixView<T,mstepi,mstepj,C,I> submatrix_type;
    typedef MatrixView<T,UNKNOWN,UNKNOWN,C,I> submatrix_step_type;
    typedef VectorView<T,UNKNOWN,C,I> subvector_type;
    typedef SmallMatrixView<T,UNKNOWN,2,mstepi,UNKNOWN,C,I> colpair_type;
    typedef SmallMatrixView<T,2,UNKNOWN,UNKNOWN,mstepj,C,I> rowpair_type;
    typedef MatrixView<T,mstepi,mstepj,C,I> cols_type;
    typedef MatrixView<T,mstepi,mstepj,C,I> rows_type;

    typedef MatrixView<T,mstepi,mstepj,C,I> view_type;
    typedef MatrixView<T,mstepi,mstepj,C,CStyle> cview_type;
    typedef MatrixView<T,mstepi,mstepj,C,FortranStyle> fview_type;
    typedef MatrixView<T,UNKNOWN,UNKNOWN,C> xview_type;
    typedef MatrixView<T,1,mstepj,C,I> cmview_type;
    typedef MatrixView<T,mstepi,1,C,I> rmview_type;
    typedef MatrixView<T,mstepi,mstepj,notC,I> conjugate_type;
    typedef MatrixView<T,mstepj,mstepi,C,I> transpose_type;
    typedef MatrixView<T,mstepj,mstepi,notC,I> adjoint_type;
    typedef UpperTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> uppertri_type;
    typedef UpperTriMatrixView<T,UnitDiag,mstepi,mstepj,C,I> unit_uppertri_type;
    typedef LowerTriMatrixView<T,NonUnitDiag,mstepi,mstepj,C,I> lowertri_type;
    typedef LowerTriMatrixView<T,UnitDiag,mstepi,mstepj,C,I> unit_lowertri_type;
    typedef MatrixView<real_type,twoSi,twoSj,false,I> realview_type;
    typedef realview_type imagview_type;
    typedef VectorView<T,1,C,I> linearview_type;
    typedef MatrixView<T,mstepi,mstepj,false,I> nonconj_type;
  };

  template <class T, int Si, int Sj, bool C, IndexStyle I>
  class MatrixView :
    public BaseMatrix_Rec_Mutable<MatrixView<T,Si,Sj,C,I> >
  {
  public:

    typedef MatrixView<T,Si,Sj,C,I> type;
    typedef BaseMatrix_Rec_Mutable<type> base_mut;
    typedef typename base_mut::reference reference;
    enum { mrowmajor = Traits<type>::mrowmajor };
    enum { mcolmajor = Traits<type>::mcolmajor };

    //
    // Constructors
    //

    inline MatrixView(T* m, size_t cs, size_t rs, int si, int sj) :
      itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) {}

    inline MatrixView(T* m, size_t cs, size_t rs, int si) :
      itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj) 
    { TMVStaticAssert(Sj != UNKNOWN); }

    inline MatrixView(T* m, size_t cs, size_t rs) :
      itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
    { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

    inline MatrixView(const type& m2) :
      itsm(m2.itsm), itscs(m2.colsize()), itsrs(m2.rowsize()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int Si2, int Sj2, IndexStyle I2>
    inline MatrixView(MatrixView<T,Si2,Sj2,C,I2> m2) :
      itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
    inline MatrixView(SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2> m2) :
      itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
      itssi(m2.stepi()), itssj(m2.stepj()) {}

    inline ~MatrixView() { 
#ifdef TMV_DEBUG
      itsm = 0; 
#endif
    }


    //
    // Op = 
    //

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    {
      base_mut::operator=(m2); 
      return *this;
    }

    inline type& operator=(const type& m2)
    {
      base_mut::operator=(m2); 
      return *this;
    }

    inline type& operator=(const T x)
    {
      base_mut::operator=(x); 
      return *this;
    }


    //
    // Auxilliary Functions
    //

    inline const T* cptr() const { return itsm; }
    inline T* ptr() { return itsm; }

    inline T cref(int i, int j) const
    { return DoConj<C>(itsm[i*stepi()+j*stepj()]); }

    inline reference ref(int i, int j) 
    { return reference(itsm[i*stepi()+j*stepj()]); }

    inline size_t ls() const { return itscs*itsrs; }
    inline size_t colsize() const { return itscs; }
    inline size_t rowsize() const { return itsrs; }
    inline int stepi() const { return itssi; }
    inline int stepj() const { return itssj; }
    inline bool isconj() const { return C; }
    inline bool isrm() const 
    { return mrowmajor || (!mcolmajor &&  stepj() == 1); }
    inline bool iscm() const 
    { return mcolmajor || (!mrowmajor &&  stepi() == 1); }

  private :

#ifdef TMV_DEBUG
    T* itsm;
#else
    T*const itsm;
#endif
    const size_t itscs;
    const size_t itsrs;
    const StepInt<Si> itssi;
    const StepInt<Sj> itssj;

  }; // MatrixView

  template <class T, int Si, int Sj, bool C>
  class MatrixViewF :
    public MatrixView<T,Si,Sj,C,FortranStyle>
  {
  public:

    typedef MatrixViewF<T,Si,Sj,C> type;
    typedef MatrixView<T,Si,Sj,C,FortranStyle> mtype;

    inline MatrixViewF(T* m, size_t cs, size_t rs, int si, int sj) :
      mtype(m,cs,rs,si,sj) {}
    inline MatrixViewF(T* m, size_t cs, size_t rs, int si) :
      mtype(m,cs,rs,si) {}
    inline MatrixViewF(T* m, size_t cs, size_t rs) : mtype(m,cs,rs) {}
    inline MatrixViewF(const type& m2) : mtype(m2) {}
    template <int Si2, int Sj2, IndexStyle I2>
    inline MatrixViewF(MatrixView<T,Si2,Sj2,C,I2> m2) : mtype(m2) {}
    template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
    inline MatrixViewF(SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2> m2) : mtype(m2) {}
    inline ~MatrixViewF() {}

    template <class M2>
    inline type& operator=(const BaseMatrix<M2>& m2)
    { mtype::operator=(m2); return *this; }
    inline type& operator=(const type& m2)
    { mtype::operator=(m2); return *this; }
    inline type& operator=(const T x)
    { mtype::operator=(x); return *this; }

  }; // MatrixViewF

  

  //---------------------------------------------------------------------------

  //
  // Special Creators: 
  //   MatrixViewOf(m,colsize,rowsize,stor) = MatrixView of m 
  //   MatrixViewOf(m,colsize,rowsize,stepi,stepj) = MatrixView of m 
  //

  // MatrixView of raw memory:
  template <class T> 
  inline MatrixView<T> MatrixViewOf(
      T* m, size_t colsize, size_t rowsize, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor) 
      return MatrixView<T>(m,colsize,rowsize,rowsize,1);
    else 
      return MatrixView<T>(m,colsize,rowsize,1,colsize);
  }

  template <class T> 
  inline ConstMatrixView<T> MatrixViewOf(
      const T* m, size_t colsize, size_t rowsize, StorageType stor)
  {
    TMVAssert(stor == RowMajor || stor == ColMajor);
    if (stor == RowMajor)
      return ConstMatrixView<T>(m,colsize,rowsize,rowsize,1);
    else 
      return ConstMatrixView<T>(m,colsize,rowsize,1,colsize);
  }

  template <class T> 
  inline MatrixView<T> MatrixViewOf(
      T* m, size_t colsize, size_t rowsize, int stepi, int stepj)
  { return MatrixView<T>(m,colsize,rowsize,stepi,stepj); }

  template <class T> 
  inline ConstMatrixView<T> MatrixViewOf(
      const T* m, size_t colsize, size_t rowsize, int stepi, int stepj)
  { return ConstMatrixView<T>(m,colsize,rowsize,stepi,stepj); }


  //
  // Swap
  //

  template <class T, int Si, int Sj, bool C, IndexStyle I, class M>
  inline void Swap(
      BaseMatrix_Rec_Mutable<M>& m1, MatrixView<T,Si,Sj,C,I> m2)
  { DoSwap(m1,m2); }
  template <class T, int Si, int Sj, bool C, IndexStyle I, class M>
  inline void Swap(
      MatrixView<T,Si,Sj,C,I> m1, BaseMatrix_Rec_Mutable<M>& m2)
  { DoSwap(m1,m2); }
  template <class T, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
  inline void Swap(
      MatrixView<T,Si1,Sj1,C1,I1> m1, MatrixView<T,Si2,Sj2,C2,I2> m2)
  { DoSwap(m1,m2); }


  //
  // TypeText 
  //

  template <class T, StorageType S, IndexStyle I>
  inline std::string TypeText(const Matrix<T,S,I>& )
  {
    std::ostringstream s;
    s << "Matrix<"<<TypeText(T())<<","<<Text(S)<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, int Si, int Sj, bool C, IndexStyle I>
  inline std::string TypeText(const ConstMatrixView<T,Si,Sj,C,I>& m)
  {
    std::ostringstream s;
    s << "ConstMatrixView<"<<TypeText(T());
    s << ","<<IntTraits<Si>::text();
    if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
    s << ","<<IntTraits<Sj>::text();
    if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }

  template <class T, int Si, int Sj, bool C, IndexStyle I>
  inline std::string TypeText(const MatrixView<T,Si,Sj,C,I>& m)
  {
    std::ostringstream s;
    s << "MatrixView<"<<TypeText(T());
    s << ","<<IntTraits<Si>::text();
    if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
    s << ","<<IntTraits<Sj>::text();
    if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
    s << ","<<C<<","<<Text(I)<<">";
    return s.str();
  }

} // namespace tmv

#include "tmv/TMV_TriMatrix.h"

#endif
