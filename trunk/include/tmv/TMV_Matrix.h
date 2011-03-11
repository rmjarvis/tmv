///////////////////////////////////////////////////////////////////////////////
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
// typename tmv::Matrix<double>::rowrange_type lastrows = m.rowRange(20,30);
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
//    MatrixView MatrixViewOf(T* m, size_t colsize, size_t rowsize,
//            int stepi, int stepj)
//    ConstMatrixView MatrixViewOf(const T* m, size_t colsize, size_t rowsize,
//            int stepi, int stepj)
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
//    row_sub_type row(int i, int j1, int j2)
//    const_row_sub_type row(int i, int j1, int j2) const
//        Return the ith row of the matrix as a vector
//        If j1,j2 are given, it returns the SubVector from j1 to j2 
//        (not including j2) within the row.
//
//    col_type col(int j)
//    const_col_type col(int j) const
//    col_sub_type col(int j, int i1, int i2)
//    const_col_sub_type col(int j, int j1, int j2) const
//        Return the jth column of the matrix as a Vector
//        If i1,i2 are given, it returns the SubVector from i1 to i2 
//        (not including i2) within the column.
//
//    diag_type diag()
//    const_diag_type diag() const
//        Return the diagonal of the matrix as a Vector
//
//    diag_sub_type diag(int i)
//    diag_sub_type diag(int i, int j1, int j2)
//    const_diag_sub_type diag(i) const
//    const_diag_sub_type diag(i, j1, j2) const
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
//    so Norm(m) and m.norm(), for example, are equivalent.
//
//    value_type m.det() const    or Det(m)
//        Returns the determinant of a matrix.
//        Note: If the matrix is not square, the determinant is not
//              well defined.  The returned value is such that
//              conj(det) * det = Det(Adjoint(m) * m)
//              So for real nonsquare matrices, the sign is arbitrary,
//              and for complex nonsquare matrices, it is multiplied
//              by an arbitrary phase.
//
//    real_type m.logDet(value_type* sign=0) const   or LogDet(m,sign)
//        Returns the logarithm of the absolute value of the determinant.
//        For many large matrices, the determinant yields to overflow.
//        Hence, this function is provided, which stably calculates the
//        natural logarithm of the absolute value of the determinant.
//        The optional sign argument returns the sign of the determinant
//        if value_type is real, or the factor exp(it) by which exp(logdet) 
//        would have to be multiplied to get the actual determinant.
//
//    value_type m.trace() const    or Trace(m)
//        Returns the trace of a matrix.
//        = sum_i ( a_ii )
//
//    real_type m.norm() const    or Norm(m)
//    real_type m.normF() const    or NormF(m)
//        Return the Frobenius norm of a matrix.
//        = sqrt( sum_ij |a_ij|^2 )
//
//    real_type m.normSq() const    or NormSq()
//    real_type m.normSq(real_type scale) const
//        Returns the square of Norm().
//        In the method version, you can provide an optional scale, in
//        which case the output is equal to NormSq(scale*m).
//
//    real_type m.norm1() const    or Norm1(m)
//        Returns the 1-norm of a matrix.
//        = max_j (sum_i |a_ij|)
//
//    real_type m.norm2() const    or Norm2(m)
//        Returns the 2-norm of a matrix.
//        = sqrt( Max Eigenvalue of (A.adjoint * A) )
//        = Maximum singular value
//        Note: This norm is costly to calculate if one is not 
//              otherwise doing a singular value decomposition
//              of the matrix.
//
//    real_type m.normInf() const    or NormInf(m)
//        Returns the infinity-norm of a matrix.
//        = max_i (sum_j |a_ij|)
//
//    value_type sumElements() const    or SumElements(m) 
//        Returns the sum of all elements in the matrix.
//
//    real_type sumAbsElements() const    or SumAbsElements(m) 
//        Returns the sum of absolute values of elements in the matrix.
//
//    real_type sumAbs2Elements() const    or SumAbs2Elements(m) 
//        Returns the sum of absolute values of the real and imaginary
//        components of the elements in the matrix.
//        i.e. realPart().sumAbsElements() + imagPart().sumAbsElements()
//
//    value_type maxElement() const    or MaxElement(m) 
//        Returns the maximum value of any element in the matrix.
//        As "max" doesn't make sense for complex values, for these
//        we use just the real components.
//
//    real_type maxAbsElement() const    or MaxAbsElement(m) 
//        The same as MaxElement, except absolute values are used.
//
//    real_type maxAbs2Element() const    or MaxAbsElement(m) 
//        The same as MaxElement, except abs(m.real()) + abs(m.imag()) is used.
//
//    void m.makeInverse(minv) const
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
//    void m.makeInverseATA(invata) const
//        Sets invata to the Inverse of (A.adjoint * A) for matrix m = A
//        If Ax=b is solved for x, then (AT A)^-1 is the 
//        covariance matrix of the least-squares solution x
//
//    inverse_type m.inverse() const    or Inverse(m)
//        Returns an auxiliary object that delays the calculation of the 
//        inverse until there is appropriate storage for it.
//        m.Inverse(minv) is equivalent to minv = m.Inverse().
//
//
// Modifying Functions
//
//    type& setZero()
//        Sets all elements to 0
//
//    type& setAllTo(value_type x)
//        Sets all elements to x
//
//    type& addToAll(value_type x)
//        Add x to each element of m
//
//    type& clip(real_type thresh)
//        Set to 0 all elements whose abolute value is < thresh
//
//    type& conjugateSelf()
//        Sets all elements to its conjugate
//
//    type& transposeSelf() 
//        Transposes the elements of a square matrix
//
//    type& setToIdentity(value_type x = 1)
//        Set to the identity matrix, or 
//        with a parameter, set to x times identity matrix
//
//    type& swapRows(int i1, int i2)
//        Swap two rows
//
//    type& swapCols(int j1, int j2)
//        Swap two columns
//
//    type& permuteRows(const int* p)
//    type& permuteRows(const int* p, int i1, int i2)
//        Perform a series of row swaps (0,p[0]), (1,p[1])...
//        In the second case, only do (i1,p[i1)...(i2-1,p[i2-1])
//    type& reversePermuteRows(const int* p)
//    type& reversePermuteRows(const int* p, int i1, int i2)
//        The same, but perform the swaps in reverse order
//
//    type& permuteCols(const int* p)
//    type& permuteCols(const int* p, int j1, int j2)
//        Perform a series of column swaps (0,p[0]), (1,p[1])...
//        In the second case, only do (j1,p[j1)...(j2-1,p[j2-1])
//    type& reversePermuteCols(const int* p)
//    type& reversePermuteCols(const int* p, int j1, int j2)
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
//    MatrixView<T,Si,Sj,C,I>(T* p, size_t m, size_t n, int si=Si, int sj=Sj)
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
//    submatrix_type subMatrix(int i1, int i2, int j1, int j2,
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
//           A.subMatrix(0,2,2,4).setToIdentity();
//           A.subMatrix(2,4,0,2).setAllTo(2);
//        The substep values allow you to space the elements of 
//        the submatrix at steps larger than 1.
//        eg. To make an 8x8 checkerboard of 1's and 0's, you could write:
//           Matrix<int> A(8,8,0);
//           A.subMatrix(0,8,1,9,2,2) = 1;
//           A.subMatrix(1,9,0,8,2,2) = 1;
//        Note that in this case the i2,j2 values need to be the 
//        "one past the end" value given the step size, so 9 here when
//        starting at 1.
//
//    subvector_type subVector(int i, int j, int istep, int jstep, int size)
//        Returns a subvector which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.
//        For example, the cross-diagonal from the lower left to the upper
//        right of a 6x6 matrix could be accessed using:
//        m.subVector(5,0,-1,1,6)
//
//    colpair_type colPair(int j1, int j2)
//        This returns an Mx2 submatrix which consists of the 
//        columns j1 and j2.  This is useful for multiplying two 
//        (not necessarily adjacent) columns of a matrix by a 2x2 matrix.
//
//    rowpair_type rowPair(int i1, int i2)
//        Same as above, but two rows.
//
//    colrange_type colRangeint j1, int j2)
//        This is shorthand for subMatrix(0,m.colsize(),j1,j2)
//
//    rowrange_type rowRange(int i1, int i2)
//        This is shorthand for subMatrix(i1,i2,0,m.rowsize())
//
//    realpart_type realPart()
//    imagpart_type imagPart()
//        Returns a view to the real/imag elements of a complex matrix.
//
//    view_type view()
//        Returns a view of a matrix.
//
//    conjugate_type conjugate()
//        Returns a view to the conjugate of a matrix.
//
//    transpose_type transpose()
//        Returns a view to the transpose of a matrix.
//
//    adjoint_type adjoint()
//        Returns a view to the adjoint (conjugate transpose) of a matrix.
//        Note: Some people define the adjoint as the cofactor matrix.
//              This is not the same as our definition of the Adjoint.
//
//    uppertri_type upperTri()
//        Returns a view to the upper triangle portion of the matrix.
//
//    unit_uppertri_type unitUpperTri()
//        Returns a view to the upper triangle portion of the matrix
//        where the diagonal elements are all taken to be equal to 1.
//
//    unknown_uppertri_type upperTri(dt)
//        Returns a view to the upper triangle portion of the matrix
//        with the specification of UnitDiag or NonUnitDiag given by dt.
//
//    lowertri_type lowerTri()
//        Returns a view to the lower triangle portion of the matrix.
//
//    unit_lowertri_type unitLowerTri()
//        Returns a view to the lower triangle portion of the matrix
//        where the diagonal elements are all taken to be equal to 1.
//
//    unknown_lowertri_type lowerTri(dt)
//        Returns a view to the lower triangle portion of the matrix
//        with the specification of UnitDiag or NonUnitDiag given by dt.
//
//    linearview_type linearView()
//        Returns a VectorView with all the elements of the matrix.
//        This is mostly used internally for things like MaxElement
//        and ConjugateSelf, where the matrix structure is irrelevant,
//        and we just want to do something to all the elements.
//        The correlary function CanLinearize() returns whether this is 
//        allowed.
//
//    nonconj_type nonConj()
//        Returns a view of the underlying memory elements, removing
//        any isconj that might be set in the current view.
//        This is sometimes useful when an operation has the same
//        effect regardless of isconj, so it is better to just ignore
//        any isconj value.
//
//    nonconst_type nonConst()
//        Returns a mutable view of a const Matrix.
//
//    const_view_type constView()
//        Returns a const view of a mutable Matrix.
//
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
//       and m += 3 means m += 3*I (or equivalently, m.diag().addToAll(3)).
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
//    void write(std::ostream& os) const    or os << m
//        Writes m to ostream os in the following format:
//          colsize rowsize
//          ( m(0,0) m(0,1) m(0,2) ... m(0,rowsize) )
//          ( m(1,0) m(1,1) m(1,2) ... m(1,rowsize) )
//          ...
//          ( m(colsize,0)  ...  m(colsize,rowsize) )
//
//
//    m.write(ostream& os, real_type thresh)
//        Write m to os as above, but if |m(i,j)| < thresh, write 0 instead
//
//    void read(std::istream& is)    or is >> m
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
//    m.divideUsing(dt)
//    where dt is LU, QR, QRP, or SV
//
//    Each of these also has an in-place version whcih overwrites the
//    current matrix memory with the decomposition needed for 
//    doing the division.  Obviously, if you try to use the matrix
//    after doing setDiv (implicitly or explicitly), the results will
//    be invalid.
//
//    The default method is LU (LU Decomposition) for square
//    matrices and QR (QR Decomposition) for non-square.
//   
//    See the header comment to TMV_Divider.h for more info about
//    the different algorithms.
//
//    There are also shorthands for accessing the decomposition.
//    If dt = LU, then lud() returns the LUD<T> class, performing
//    the decomposition if it hasn't be done yet.
//
//    Likewise:
//    qrd(), qrpd(), svd() return the corresponding Divider classes for
//    those algorithms.
//


#ifndef TMV_Matrix_H
#define TMV_Matrix_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Vector.h"
#include "TMV_Divider.h"
#include <vector>

namespace tmv {

    template <class M> class LUD;
    template <class M> class QRD;
    template <class M> class QRPD;
    template <class M> class SVD;

    // This contains all the stuff for doing division according to different
    // DivType values: divideUsing, setDiv, saveDiv, etc.
    // SmallMatrix varieties don't have it, so we put it all here, and let
    // Matrix, ConstMatrixView and MatrixView inherit from this.
    template <class M, bool hasdivider>
    class MatrixDivHelper;

    template <class M>
    class MatrixDivHelper<M,true>
    {
    public:

        typedef const Divider<typename Traits<M>::value_type> div_type;
        typedef const div_type* getdiv_type;
        typedef typename Traits<M>::lud_type lud_type;
        typedef typename Traits<M>::qrd_type qrd_type;
        typedef typename Traits<M>::qrpd_type qrpd_type;
        typedef typename Traits<M>::svd_type svd_type;

        // Constructor starts with a default of LU or QR depending on 
        // whether matrix is square.
        MatrixDivHelper(bool isSquare) : 
            divtype(isSquare ? tmv::LU : tmv::QR) {}

        ~MatrixDivHelper() {}

        void divideInPlace() const
        { divtype |= tmv::DivInPlaceFlag; saveDiv(); }

        bool divIsInPlace() const 
        { return divtype & tmv::DivInPlaceFlag; }

        void saveDiv() const
        { divtype |= tmv::SaveDivFlag; }

        void unsaveDiv() const
        { divtype &= ~tmv::SaveDivFlag; }

        bool divIsSaved() const 
        { return divtype & tmv::SaveDivFlag; }

        void divideUsing(DivType dt) const
        {
            TMVAssert(dt == tmv::LU || dt == tmv::QR || 
                      dt == tmv::QRP || dt == tmv::SV);
            if (!(divtype & dt)) {
                unsetDiv();
                divtype &= ~tmv::DivTypeFlags;
                divtype |= dt;
            }
        }

        DivType getDivType() const 
        { return divtype & tmv::DivTypeFlags; }

        void setDiv() const
        {
            TMVStaticAssert(!Traits<typename M::real_type>::isinteger);
            if (!divIsSet()) {
                DivType dt = getDivType();
                TMVAssert(dt == tmv::LU /*|| dt == tmv::QR || 
                          dt == tmv::QRP || dt == tmv::SV*/);
                switch (dt) {
                  case tmv::LU : 
                       divider.reset(makeLUD(mat2(),divIsInPlace()));
                       break;
                  default :
                       // The above assert should have already failed.
                       // So go ahead and fall through.
                       break;
                }
            }
        }

        void unsetDiv() const
        { divider.reset(); }

        bool divIsSet() const
        { return getDiv(); }

        void resetDiv() const
        { unsetDiv(); setDiv(); }

        void doneDiv() const
        { if (!divIsSaved()) unsetDiv(); }

        getdiv_type getDiv() const
        { return divider.get(); }

        lud_type lud() const
        {
            divideUsing(LU);
            setDiv();
            TMVAssert(dynamic_cast<const tmv::LUD<M>*>(getDiv()));
            return static_cast<lud_type>(*getDiv());
        }

        qrd_type qrd() const
        {
            divideUsing(QR);
            setDiv();
            TMVAssert(dynamic_cast<const tmv::QRD<M>*>(getDiv()));
            return static_cast<qrd_type>(*getDiv());
        }

        qrpd_type qrpd() const
        {
            divideUsing(QRP);
            setDiv();
            TMVAssert(dynamic_cast<const tmv::QRPD<M>*>(getDiv()));
            return static_cast<qrpd_type>(*getDiv());
        }

        svd_type svd() const
        {
            divideUsing(SV);
            setDiv();
            TMVAssert(dynamic_cast<const tmv::SVD<M>*>(getDiv()));
            return static_cast<svd_type>(*getDiv());
        }

        // use name mat2 rather than mat to avoid ambiguating mat() from
        // BaseMatrix when inheriting from both BaseMatrix and MatrixDivHelper.
        const M& mat2() const
        { return *static_cast<const M*>(this); }

        bool checkDecomp(std::ostream* fout=0) const
        { return getDiv()->checkDecomp(mat2(),fout); }

        template <class M2>
        bool checkDecomp(
            const BaseMatrix_Rec<M2>& m, std::ostream* fout=0) const
        { return getDiv()->checkDecomp(m,fout); }

    private:

        mutable std::auto_ptr<div_type> divider;
        mutable DivType divtype;

    }; // MatrixDivHelper

    template <class M>
    class MatrixDivHelper<M,false>
    {
    public:

        MatrixDivHelper(bool ) {}
        ~MatrixDivHelper() {}
        void divideInPlace() const {}
        bool divIsInPlace() const { return false; }
        void saveDiv() const {}
        void unsaveDiv() const {}
        bool divIsSaved() const { return false; }
        void divideUsing(DivType dt) const {}
        DivType getDivType() const { return XXX; }
        void setDiv() const {}
        void unsetDiv() const {}
        bool divIsSet() const { return false; }
        void resetDiv() const {}
        void doneDiv() const {}
        void* getDiv() const{ return 0; }
        LUD<M> lud() const { return LUD<M>(); }
        void qrd() const {}
        void qrpd() const {}
        void svd() const {}
        bool checkDecomp(std::ostream* =0) const { return false; }

        template <class M2>
        bool checkDecomp(const BaseMatrix_Rec<M2>& , std::ostream* =0) const 
        { return false; }
    };

    template <class T, StorageType S, IndexStyle I>
    struct Traits<Matrix<T,S,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef Matrix<T,S,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef type copy_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _shape = Rec };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (S == RowMajor) };
        enum { _colmajor = (S == ColMajor) };
        enum { _stor = S };
        enum { _stepi = (S==ColMajor ? 1 : UNKNOWN) };
        enum { _stepj = (S==RowMajor ? 1 : UNKNOWN) };
        enum { _diagstep = UNKNOWN };
        enum { _conj = false };
        enum { _canlin = true };
        enum { twoSi = isreal ? int(_stepi) : int(IntTraits<_stepi>::twoS) };
        enum { twoSj = isreal ? int(_stepj) : int(IntTraits<_stepj>::twoS) };
        enum { notC = iscomplex };

        enum { _hasdivider = Traits<real_type>::isinst && !Traits<real_type>::isinteger };
        typedef typename TypeSelect<_hasdivider ,
                const LUD<type>& , LUD<type> >::type lud_type;
#if 1
        typedef InvalidType qrd_type;
        typedef InvalidType qrpd_type;
        typedef InvalidType svd_type;
#else
        typedef const QRD<type>& qrd_type;
        typedef const QRPD<type>& qrpd_type;
        typedef const SVD<type>& svd_type;
#endif

        typedef ConstVectorView<T,_stepi,false,I> const_col_type;
        typedef ConstVectorView<T,_stepi,false,I> const_col_sub_type;
        typedef ConstVectorView<T,_stepj,false,I> const_row_type;
        typedef ConstVectorView<T,_stepj,false,I> const_row_sub_type;
        typedef ConstVectorView<T,_diagstep,false,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,false,I> const_diag_sub_type;

        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,false,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,false,I> const_subvector_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,false,I> 
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,false,I> 
            const_rowpair_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_colrange_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_rowrange_type;

        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_view_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,CStyle> const_cview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,FortranStyle> 
            const_fview_type;
        typedef ConstMatrixView<T> const_xview_type;
        typedef ConstMatrixView<T,1,_stepj,false,I> const_cmview_type;
        typedef ConstMatrixView<T,_stepi,1,false,I> const_rmview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,notC,I> const_conjugate_type;
        typedef ConstMatrixView<T,_stepj,_stepi,false,I> const_transpose_type;
        typedef ConstMatrixView<T,_stepj,_stepi,notC,I> const_adjoint_type;
        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            const_uppertri_type;
        typedef ConstUpperTriMatrixView<T,UnitDiag,_stepi,_stepj,false,I> 
            const_unit_uppertri_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,false,I> 
            const_unknown_uppertri_type;
        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            const_lowertri_type;
        typedef ConstLowerTriMatrixView<T,UnitDiag,_stepi,_stepj,false,I> 
            const_unit_lowertri_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,false,I> 
            const_unknown_lowertri_type;
        typedef ConstMatrixView<real_type,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,1,false,I> const_linearview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_nonconj_type;
        typedef MatrixView<T,_stepi,_stepj,false,I> nonconst_type;

        typedef T& reference;

        typedef VectorView<T,_stepi,false,I> col_type;
        typedef VectorView<T,_stepi,false,I> col_sub_type;
        typedef VectorView<T,_stepj,false,I> row_type;
        typedef VectorView<T,_stepj,false,I> row_sub_type;
        typedef VectorView<T,_diagstep,false,I> diag_type;
        typedef VectorView<T,_diagstep,false,I> diag_sub_type;

        typedef MatrixView<T,_stepi,_stepj,false,I> submatrix_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,false,I> submatrix_step_type;
        typedef VectorView<T,UNKNOWN,false,I> subvector_type;
        typedef SmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,false,I> 
            colpair_type;
        typedef SmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,false,I> 
            rowpair_type;
        typedef MatrixView<T,_stepi,_stepj,false,I> colrange_type;
        typedef MatrixView<T,_stepi,_stepj,false,I> rowrange_type;

        typedef MatrixView<T,_stepi,_stepj,false,I> view_type;
        typedef MatrixView<T,_stepi,_stepj,false,CStyle> cview_type;
        typedef MatrixView<T,_stepi,_stepj,false,FortranStyle> fview_type;
        typedef MatrixView<T> xview_type;
        typedef MatrixView<T,1,_stepj,false,I> cmview_type;
        typedef MatrixView<T,_stepi,1,false,I> rmview_type;
        typedef MatrixView<T,_stepi,_stepj,notC,I> conjugate_type;
        typedef MatrixView<T,_stepj,_stepi,false,I> transpose_type;
        typedef MatrixView<T,_stepj,_stepi,notC,I> adjoint_type;
        typedef UpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            uppertri_type;
        typedef UpperTriMatrixView<T,UnitDiag,_stepi,_stepj,false,I> 
            unit_uppertri_type;
        typedef UpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,false,I> 
            unknown_uppertri_type;
        typedef LowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,false,I> 
            lowertri_type;
        typedef LowerTriMatrixView<T,UnitDiag,_stepi,_stepj,false,I> 
            unit_lowertri_type;
        typedef LowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,false,I> 
            unknown_lowertri_type;
        typedef MatrixView<real_type,twoSi,twoSj,false,I> realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<T,1,false,I> linearview_type;
        typedef MatrixView<T,_stepi,_stepj,false,I> nonconj_type;
    };

    template <class T, StorageType S, IndexStyle I>
    class Matrix : 
        public BaseMatrix_Rec_Mutable<Matrix<T,S,I> >,
        public MatrixDivHelper<Matrix<T,S,I> , Traits<Matrix<T,S,I> >::_hasdivider >
    {
    public:

        typedef Matrix<T,S,I> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        typedef BaseMatrix_Rec_Mutable<type> base_div;
        typedef MatrixDivHelper<type,Traits<type>::_hasdivider> divhelper;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _canlin = Traits<type>::_canlin };


        //
        // Constructors
        //

        Matrix() : 
            divhelper(true),
            itscs(0), itsrs(0), linsize(0), itsm(0)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
        }

        Matrix(size_t cs, size_t rs) :
            divhelper(cs==rs),
            itscs(cs), itsrs(rs), linsize(cs*rs), itsm(linsize)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

        Matrix(size_t cs, size_t rs, T x) :
            divhelper(cs==rs),
            itscs(cs), itsrs(rs), linsize(cs*rs), itsm(linsize)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            this->setAllTo(x);
        }

        Matrix(size_t cs, size_t rs, const T* vv) :
            divhelper(cs==rs),
            itscs(cs), itsrs(rs), linsize(cs*rs), itsm(linsize)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            typename type::linearview_type lv = this->linearView();
            ConstVectorView<T,1>(vv,linsize).newAssignTo(lv);
        }

        Matrix(size_t cs, size_t rs, const std::vector<T>& vv) : 
            divhelper(cs==rs),
            itscs(cs), itsrs(rs), linsize(cs*rs), itsm(linsize)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVAssert(vv.size() == linsize);
            typename type::linearview_type lv = this->linearView();
            ConstVectorView<T,1>(&vv[0],linsize).newAssignTo(lv);
        }

        Matrix(const std::vector<std::vector<T> >& vv) :
            divhelper(false),
            itscs(vv.size()), itsrs(0), linsize(0), itsm(0)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            if (itscs > 0) {
                resize(itscs,vv[0].size());
                for(size_t i=0;i<itscs;++i) {
                    TMVAssert(vv[i].size() == itsrs);
                    typename type::row_type mi = this->row(i);
                    ConstVectorView<T,1>(&vv[i][0],itsrs).newAssignTo(mi);
                }
            }
        }

        Matrix(const type& m2) :
            divhelper(m2.isSquare()),
            itscs(m2.itscs), itsrs(m2.itsrs),
            linsize(m2.linsize), itsm(linsize)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            m2.newAssignTo(*this);
        }

        template <class M2>
        Matrix(const BaseMatrix<M2>& m2) :
            divhelper(m2.isSquare()),
            itscs(m2.colsize()), itsrs(m2.rowsize()),
            linsize(itscs * itsrs), itsm(linsize)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            m2.newAssignTo(*this);
        }

        template <class M2>
        Matrix(const BaseMatrix_Rec<M2>& m2) :
            divhelper(m2.isSquare()),
            itscs(m2.colsize()), itsrs(m2.rowsize()),
            linsize(m2.ls()), itsm(linsize)
        {
            TMVStaticAssert(S==RowMajor || S==ColMajor);
            m2.newAssignTo(*this);
        }

        ~Matrix() 
        {
#ifdef TMV_DEBUG
            this->setAllTo(T(999));
#endif
        }

        //
        // Op=
        //

        type& operator=(const type& m2)
        {
            if (&m2 != this) base_mut::operator=(m2);
            return *this; 
        }

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2);
            return *this; 
        }

        type& operator=(T x)
        {
            base_mut::operator=(x);
            return *this; 
        }


        //
        // Auxilliary Functions
        //

        const T* cptr() const { return itsm; }
        T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return itsm[S==RowMajor ? i*stepi()+j : i+j*stepj()]; }

        T& ref(int i, int j)
        { return itsm[S==RowMajor ? i*stepi()+j : i+j*stepj()]; }

        void swapWith(type& m2)
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            if (itsm.get() == m2.itsm.get()) return;
            itsm.swapWith(m2.itsm);
        }
        
        void resize(const size_t cs, const size_t rs)
        {
            itscs = cs;
            itsrs = rs;
            if (cs == rs) divhelper::divideUsing(tmv::LU);
            else divhelper::divideUsing(tmv::QR);
            linsize = cs*rs;
            itsm.resize(linsize);
#ifdef TMV_DEBUG
            this->setAllTo(T(888));
#endif
        }

        size_t ls() const { return linsize; }
        size_t colsize() const { return itscs; }
        size_t rowsize() const { return itsrs; }
        int nElements() const { return itscs*itsrs; }
        int stepi() const { return S == RowMajor ? itsrs : 1; }
        int stepj() const { return S == RowMajor ? 1 : itscs; }
        bool isconj() const { return false; }
        bool isrm() const { return S == RowMajor; }
        bool iscm() const { return S == ColMajor; }
        StorageType stor() const { return S; }

    private:

        size_t itscs;
        size_t itsrs;
        size_t linsize;
        AlignedArray<T> itsm;

    }; // Matrix

    template <class T, StorageType S>
    class MatrixF : public Matrix<T,S,FortranStyle>
    {
    public:
        typedef MatrixF<T,S> type;
        typedef Matrix<T,S,FortranStyle> mtype;

        MatrixF(size_t cs, size_t rs) : mtype(cs,rs) {}
        MatrixF(size_t cs, size_t rs, T x) : mtype(cs,rs,x) {}
        MatrixF(size_t cs, size_t rs, const T* vv) : mtype(cs,rs,vv) {}
        MatrixF(size_t cs, size_t rs, const std::vector<T>& vv) : 
            mtype(cs,rs,vv) {}
        MatrixF(const std::vector<std::vector<T> >& vv) : mtype(vv) {}
        MatrixF(const type& m2) : mtype(m2) {}
        template <class M2>
        MatrixF(const BaseMatrix<M2>& m2) : mtype(m2) {}
        ~MatrixF() {}

        type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        type& operator=(T x)
        { mtype::operator=(x); return *this; }

    }; // MatrixF

    template <class T, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<ConstMatrixView<T,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstMatrixView<T,Si,Sj,C,I> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef Matrix<T,Sj==1?RowMajor:ColMajor,I> copy_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _shape = Rec };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (Sj == 1) };
        enum { _colmajor = (Si == 1) };
        enum { _stor = (_rowmajor ? RowMajor : ColMajor) };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = C };
        enum { _canlin = false };
        enum { twoSi = isreal ? Si : IntTraits<Si>::twoS };
        enum { twoSj = isreal ? Sj : IntTraits<Sj>::twoS };
        enum { notC = !C && iscomplex };

        enum { _hasdivider = Traits<real_type>::isinst && !Traits<real_type>::isinteger };
        typedef typename TypeSelect<_hasdivider ,
                const LUD<type>& , LUD<type> >::type lud_type;
#if 1
        typedef InvalidType qrd_type;
        typedef InvalidType qrpd_type;
        typedef InvalidType svd_type;
#else
        typedef const QRD<type>& qrd_type;
        typedef const QRPD<type>& qrpd_type;
        typedef const SVD<type>& svd_type;
#endif

        typedef ConstVectorView<T,_stepi,C,I> const_col_type;
        typedef ConstVectorView<T,_stepi,C,I> const_col_sub_type;
        typedef ConstVectorView<T,_stepj,C,I> const_row_type;
        typedef ConstVectorView<T,_stepj,C,I> const_row_sub_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_sub_type;

        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,C,I> 
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,C,I> 
            const_rowpair_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_colrange_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_rowrange_type;

        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_view_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,CStyle> const_cview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,FortranStyle> 
            const_fview_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C> const_xview_type;
        typedef ConstMatrixView<T,1,_stepj,C,I> const_cmview_type;
        typedef ConstMatrixView<T,_stepi,1,C,I> const_rmview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,notC,I> const_conjugate_type;
        typedef ConstMatrixView<T,_stepj,_stepi,C,I> const_transpose_type;
        typedef ConstMatrixView<T,_stepj,_stepi,notC,I> const_adjoint_type;
        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            const_uppertri_type;
        typedef ConstUpperTriMatrixView<T,UnitDiag,_stepi,_stepj,C,I> 
            const_unit_uppertri_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,C,I> 
            const_unknown_uppertri_type;
        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            const_lowertri_type;
        typedef ConstLowerTriMatrixView<T,UnitDiag,_stepi,_stepj,C,I> 
            const_unit_lowertri_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,C,I> 
            const_unknown_lowertri_type;
        typedef ConstMatrixView<real_type,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,1,C,I> const_linearview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_nonconj_type;
        typedef MatrixView<T,_stepi,_stepj,C,I> nonconst_type;
    };

    template <class T, int Si, int Sj, bool C, IndexStyle I>
    class ConstMatrixView : 
        public BaseMatrix_Rec<ConstMatrixView<T,Si,Sj,C,I> >,
        public MatrixDivHelper<ConstMatrixView<T,Si,Sj,C,I>,Traits<ConstMatrixView<T,Si,Sj,C,I> >::_hasdivider >
    {
    public:
        typedef ConstMatrixView<T,Si,Sj,C,I> type;
        typedef MatrixDivHelper<type,Traits<type>::_hasdivider> divhelper;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _canlin = Traits<type>::_canlin };

        //
        // Constructors
        //

        ConstMatrixView(
            const T* m, size_t cs, size_t rs, int si, int sj) :
            divhelper(cs==rs),
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) 
        {}

        ConstMatrixView(const T* m, size_t cs, size_t rs, int si) :
            divhelper(cs==rs),
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj) 
        { TMVStaticAssert(Sj != UNKNOWN); }

        ConstMatrixView(const T* m, size_t cs, size_t rs) :
            divhelper(cs==rs),
            itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        ConstMatrixView(const type& m2) :
            divhelper(m2.isSquare()),
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {}

        ConstMatrixView(const MatrixView<T,Si,Sj,C,I>& m2) :
            divhelper(m2.isSquare()),
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {}

        template <int Si2, int Sj2, IndexStyle I2>
        ConstMatrixView(const ConstMatrixView<T,Si2,Sj2,C,I2>& m2) :
            divhelper(m2.isSquare()),
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {}

        template <int Si2, int Sj2, IndexStyle I2>
        ConstMatrixView(const MatrixView<T,Si2,Sj2,C,I2>& m2) :
            divhelper(m2.isSquare()),
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {}

        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        ConstMatrixView(
            const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) :
            divhelper(m2.isSquare()),
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {}

        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        ConstMatrixView(
            const SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) :
            divhelper(m2.isSquare()),
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        {}

        ~ConstMatrixView() { 
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }

    private :
        void operator=(const type& m2);
    public :

        //
        // Auxilliary Functions
        //

        const T* cptr() const { return itsm; }

        T cref(int i, int j) const
        { return DoConj<C>(itsm[i*stepi()+j*stepj()]); }

        size_t ls() const { return itscs*itsrs; }
        size_t colsize() const { return itscs; }
        size_t rowsize() const { return itsrs; }
        int nElements() const { return itscs*itsrs; }
        int stepi() const { return itssi; }
        int stepj() const { return itssj; }
        bool isconj() const { return C; }
        bool isrm() const 
        { return _rowmajor || (!_colmajor &&  stepj() == 1); }
        bool iscm() const 
        { return _colmajor || (!_rowmajor &&  stepi() == 1); }

    private :

        const T* itsm;
        const size_t itscs;
        const size_t itsrs;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // ConstMatrixView

    template <class T, int Si, int Sj, bool C>
    class ConstMatrixViewF : public ConstMatrixView<T,Si,Sj,C,FortranStyle>
    {
    public:

        typedef ConstMatrixViewF<T,Si,Sj,C> type;
        typedef ConstMatrixView<T,Si,Sj,C,FortranStyle> mtype;

        ConstMatrixViewF(
            const T* m, size_t cs, size_t rs, int si, int sj) :
            mtype(m,cs,rs,si,sj) 
        {}
        ConstMatrixViewF(const T* m, size_t cs, size_t rs, int si) :
            mtype(m,cs,rs,si) 
        {}
        ConstMatrixViewF(const T* m, size_t cs, size_t rs) : 
            mtype(m,cs,rs) 
        {}
        ConstMatrixViewF(const type& m2) : mtype(m2) {}
        template <int Si2, int Sj2, IndexStyle I2>
        ConstMatrixViewF(const ConstMatrixView<T,Si2,Sj2,C,I2>& m2) :
            mtype(m2) 
        {}
        template <int Si2, int Sj2, IndexStyle I2>
        ConstMatrixViewF(const MatrixView<T,Si2,Sj2,C,I2>& m2) :
            mtype(m2) 
        {}
        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        ConstMatrixViewF(
            const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) : mtype(m2) 
        {}
        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        ConstMatrixViewF(
            const SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2>& m2) : mtype(m2) 
        {}
        ~ConstMatrixViewF() 
        {}

    private :
        void operator=(const type& m2);

    }; // ConstMatrixViewF

    template <class T, int Si, int Sj, bool C, IndexStyle I>
    struct Traits<MatrixView<T,Si,Sj,C,I> >
    {
        typedef T value_type;

        typedef typename Traits<T>::real_type real_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef MatrixView<T,Si,Sj,C,I> type;
        typedef const ConstMatrixView<T,Si,Sj,C,I> calc_type;
        typedef calc_type eval_type;
        typedef Matrix<T,Sj==1?RowMajor:ColMajor,I> copy_type;

        typedef QuotXM<1,real_type,type> inverse_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _shape = Rec };
        enum { _fort = (I == FortranStyle) };
        enum { _calc = true };
        enum { _rowmajor = (Sj == 1) };
        enum { _colmajor = (Si == 1) };
        enum { _stor = (_rowmajor ? RowMajor : ColMajor) };
        enum { _stepi = Si };
        enum { _stepj = Sj };
        enum { _diagstep = IntTraits2<Si,Sj>::sum };
        enum { _conj = C };
        enum { _canlin = false };
        enum { twoSi = isreal ? Si : IntTraits<Si>::twoS };
        enum { twoSj = isreal ? Sj : IntTraits<Sj>::twoS };
        enum { notC = !C && iscomplex };

        enum { _hasdivider = Traits<real_type>::isinst && !Traits<real_type>::isinteger };
        typedef typename TypeSelect<_hasdivider ,
                const LUD<type>& , LUD<type> >::type lud_type;
#if 1
        typedef InvalidType qrd_type;
        typedef InvalidType qrpd_type;
        typedef InvalidType svd_type;
#else
        typedef const QRD<type>& qrd_type;
        typedef const QRPD<type>& qrpd_type;
        typedef const SVD<type>& svd_type;
#endif

        typedef ConstVectorView<T,_stepi,C,I> const_col_type;
        typedef ConstVectorView<T,_stepi,C,I> const_col_sub_type;
        typedef ConstVectorView<T,_stepj,C,I> const_row_type;
        typedef ConstVectorView<T,_stepj,C,I> const_row_sub_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_type;
        typedef ConstVectorView<T,_diagstep,C,I> const_diag_sub_type;

        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_submatrix_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C,I> 
            const_submatrix_step_type;
        typedef ConstVectorView<T,UNKNOWN,C,I> const_subvector_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,C,I> 
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,C,I> 
            const_rowpair_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_colrange_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_rowrange_type;

        typedef ConstMatrixView<T,_stepi,_stepj,C,I> const_view_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,CStyle> const_cview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,C,FortranStyle> 
            const_fview_type;
        typedef ConstMatrixView<T,UNKNOWN,UNKNOWN,C> const_xview_type;
        typedef ConstMatrixView<T,1,_stepj,C,I> const_cmview_type;
        typedef ConstMatrixView<T,_stepi,1,C,I> const_rmview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,notC,I> const_conjugate_type;
        typedef ConstMatrixView<T,_stepj,_stepi,C,I> const_transpose_type;
        typedef ConstMatrixView<T,_stepj,_stepi,notC,I> const_adjoint_type;
        typedef ConstUpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            const_uppertri_type;
        typedef ConstUpperTriMatrixView<T,UnitDiag,_stepi,_stepj,C,I> 
            const_unit_uppertri_type;
        typedef ConstUpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,C,I> 
            const_unknown_uppertri_type;
        typedef ConstLowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            const_lowertri_type;
        typedef ConstLowerTriMatrixView<T,UnitDiag,_stepi,_stepj,C,I> 
            const_unit_lowertri_type;
        typedef ConstLowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,C,I> 
            const_unknown_lowertri_type;
        typedef ConstMatrixView<real_type,twoSi,twoSj,false,I> 
            const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,1,C,I> const_linearview_type;
        typedef ConstMatrixView<T,_stepi,_stepj,false,I> const_nonconj_type;
        typedef MatrixView<T,_stepi,_stepj,C,I> nonconst_type;

        typedef typename AuxRef<T,C>::reference reference;

        typedef VectorView<T,_stepi,C,I> col_type;
        typedef VectorView<T,_stepi,C,I> col_sub_type;
        typedef VectorView<T,_stepj,C,I> row_type;
        typedef VectorView<T,_stepj,C,I> row_sub_type;
        typedef VectorView<T,_diagstep,C,I> diag_type;
        typedef VectorView<T,_diagstep,C,I> diag_sub_type;

        typedef MatrixView<T,_stepi,_stepj,C,I> submatrix_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,C,I> submatrix_step_type;
        typedef VectorView<T,UNKNOWN,C,I> subvector_type;
        typedef SmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,C,I> colpair_type;
        typedef SmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,C,I> rowpair_type;
        typedef MatrixView<T,_stepi,_stepj,C,I> colrange_type;
        typedef MatrixView<T,_stepi,_stepj,C,I> rowrange_type;

        typedef MatrixView<T,_stepi,_stepj,C,I> view_type;
        typedef MatrixView<T,_stepi,_stepj,C,CStyle> cview_type;
        typedef MatrixView<T,_stepi,_stepj,C,FortranStyle> fview_type;
        typedef MatrixView<T,UNKNOWN,UNKNOWN,C> xview_type;
        typedef MatrixView<T,1,_stepj,C,I> cmview_type;
        typedef MatrixView<T,_stepi,1,C,I> rmview_type;
        typedef MatrixView<T,_stepi,_stepj,notC,I> conjugate_type;
        typedef MatrixView<T,_stepj,_stepi,C,I> transpose_type;
        typedef MatrixView<T,_stepj,_stepi,notC,I> adjoint_type;
        typedef UpperTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            uppertri_type;
        typedef UpperTriMatrixView<T,UnitDiag,_stepi,_stepj,C,I> 
            unit_uppertri_type;
        typedef UpperTriMatrixView<T,UnknownDiag,_stepi,_stepj,C,I> 
            unknown_uppertri_type;
        typedef LowerTriMatrixView<T,NonUnitDiag,_stepi,_stepj,C,I> 
            lowertri_type;
        typedef LowerTriMatrixView<T,UnitDiag,_stepi,_stepj,C,I> 
            unit_lowertri_type;
        typedef LowerTriMatrixView<T,UnknownDiag,_stepi,_stepj,C,I> 
            unknown_lowertri_type;
        typedef MatrixView<real_type,twoSi,twoSj,false,I> realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<T,1,C,I> linearview_type;
        typedef MatrixView<T,_stepi,_stepj,false,I> nonconj_type;
    };

    template <class T, int Si, int Sj, bool C, IndexStyle I>
    class MatrixView : 
        public BaseMatrix_Rec_Mutable<MatrixView<T,Si,Sj,C,I> >,
        public MatrixDivHelper<MatrixView<T,Si,Sj,C,I>,Traits<MatrixView<T,Si,Sj,C,I> >::_hasdivider >
    {
    public:

        typedef MatrixView<T,Si,Sj,C,I> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;
        typedef MatrixDivHelper<type,Traits<type>::_hasdivider> divhelper;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stor = Traits<type>::_stor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _canlin = Traits<type>::_canlin };

        //
        // Constructors
        //

        MatrixView(T* m, size_t cs, size_t rs, int si, int sj) :
            divhelper(cs==rs),
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) {}

        MatrixView(T* m, size_t cs, size_t rs, int si) :
            divhelper(cs==rs),
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(Sj) 
        { TMVStaticAssert(Sj != UNKNOWN); }

        MatrixView(T* m, size_t cs, size_t rs) :
            divhelper(cs==rs),
            itsm(m), itscs(cs), itsrs(rs), itssi(Si), itssj(Sj)
        { TMVStaticAssert(Si != UNKNOWN); TMVStaticAssert(Sj != UNKNOWN); }

        MatrixView(const type& m2) :
            divhelper(m2.isSquare()),
            itsm(m2.itsm), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int Si2, int Sj2, IndexStyle I2>
        MatrixView(MatrixView<T,Si2,Sj2,C,I2> m2) :
            divhelper(m2.isSquare()),
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        MatrixView(SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2> m2) :
            divhelper(m2.isSquare()),
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) {}

        ~MatrixView() { 
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }


        //
        // Op = 
        //

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        {
            base_mut::operator=(m2); 
            return *this;
        }

        type& operator=(const type& m2)
        {
            base_mut::operator=(m2); 
            return *this;
        }

        type& operator=(const T x)
        {
            base_mut::operator=(x); 
            return *this;
        }


        //
        // Auxilliary Functions
        //

        const T* cptr() const { return itsm; }
        T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return DoConj<C>(itsm[i*stepi()+j*stepj()]); }

        reference ref(int i, int j) 
        { return reference(itsm[i*stepi()+j*stepj()]); }

        size_t ls() const { return itscs*itsrs; }
        size_t colsize() const { return itscs; }
        size_t rowsize() const { return itsrs; }
        int nElements() const { return itscs*itsrs; }
        int stepi() const { return itssi; }
        int stepj() const { return itssj; }
        bool isconj() const { return C; }
        bool isrm() const 
        { return _rowmajor || (!_colmajor &&  stepj() == 1); }
        bool iscm() const 
        { return _colmajor || (!_rowmajor &&  stepi() == 1); }

    private :

        T* itsm;
        const size_t itscs;
        const size_t itsrs;
        const CheckedInt<Si> itssi;
        const CheckedInt<Sj> itssj;

    }; // MatrixView

    template <class T, int Si, int Sj, bool C>
    class MatrixViewF : public MatrixView<T,Si,Sj,C,FortranStyle>
    {
    public:

        typedef MatrixViewF<T,Si,Sj,C> type;
        typedef MatrixView<T,Si,Sj,C,FortranStyle> mtype;

        MatrixViewF(T* m, size_t cs, size_t rs, int si, int sj) :
            mtype(m,cs,rs,si,sj) 
        {}
        MatrixViewF(T* m, size_t cs, size_t rs, int si) :
            mtype(m,cs,rs,si) 
        {}
        MatrixViewF(T* m, size_t cs, size_t rs) : mtype(m,cs,rs) 
        {}
        MatrixViewF(const type& m2) : mtype(m2) 
        {}
        template <int Si2, int Sj2, IndexStyle I2>
        MatrixViewF(MatrixView<T,Si2,Sj2,C,I2> m2) : mtype(m2) 
        {}
        template <int M2, int N2, int Si2, int Sj2, IndexStyle I2>
        MatrixViewF(SmallMatrixView<T,M2,N2,Si2,Sj2,C,I2> m2) :
            mtype(m2) 
        {}
        ~MatrixViewF() 
        {}

        template <class M2>
        type& operator=(const BaseMatrix<M2>& m2)
        { mtype::operator=(m2); return *this; }
        type& operator=(const type& m2)
        { mtype::operator=(m2); return *this; }
        type& operator=(const T x)
        { mtype::operator=(x); return *this; }

    }; // MatrixViewF



    //-------------------------------------------------------------------------

    //
    // Special Creators: 
    //   MatrixViewOf(m,colsize,rowsize,stor) = MatrixView of m 
    //   MatrixViewOf(m,colsize,rowsize,stepi,stepj) = MatrixView of m 
    //

    // MatrixView of raw memory:
    template <class T> 
    static MatrixView<T> MatrixViewOf(
        T* m, size_t colsize, size_t rowsize, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor) 
            return MatrixView<T>(m,colsize,rowsize,rowsize,1);
        else 
            return MatrixView<T>(m,colsize,rowsize,1,colsize);
    }

    template <class T> 
    static ConstMatrixView<T> MatrixViewOf(
        const T* m, size_t colsize, size_t rowsize, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstMatrixView<T>(m,colsize,rowsize,rowsize,1);
        else 
            return ConstMatrixView<T>(m,colsize,rowsize,1,colsize);
    }

    template <class T> 
    static MatrixView<T> MatrixViewOf(
        T* m, size_t colsize, size_t rowsize, int stepi, int stepj)
    { return MatrixView<T>(m,colsize,rowsize,stepi,stepj); }

    template <class T> 
    static ConstMatrixView<T> MatrixViewOf(
        const T* m, size_t colsize, size_t rowsize, int stepi, int stepj)
    { return ConstMatrixView<T>(m,colsize,rowsize,stepi,stepj); }


    //
    // Swap
    //

    template <class T, StorageType S, IndexStyle I>
    static void Swap(Matrix<T,S,I>& m1, Matrix<T,S,I>& m2)
    { m1.swapWith(m2); }
    template <class M, class T, int Si, int Sj, bool C, IndexStyle I>
    static void Swap(BaseMatrix_Rec_Mutable<M>& m1, MatrixView<T,Si,Sj,C,I> m2)
    { DoSwap(m1,m2); }
    template <class M, class T, int Si, int Sj, bool C, IndexStyle I>
    static void Swap(MatrixView<T,Si,Sj,C,I> m1, BaseMatrix_Rec_Mutable<M>& m2)
    { DoSwap(m1,m2); }
    template <class T, int Si1, int Sj1, bool C1, IndexStyle I1, int Si2, int Sj2, bool C2, IndexStyle I2>
    static void Swap(
        MatrixView<T,Si1,Sj1,C1,I1> m1, MatrixView<T,Si2,Sj2,C2,I2> m2)
    { DoSwap(m1,m2); }


    //
    // Conjugate, Transpose, Adjoint
    //

    template <class T, StorageType S, IndexStyle I>
    static typename Matrix<T,S,I>::conjugate_type Conjugate(Matrix<T,S,I>& m)
    { return m.conjugate(); }
    template <class T, int Si, int Sj, bool C, IndexStyle I>
    static typename MatrixView<T,Si,Sj,C,I>::conjugate_type Conjugate(
        MatrixView<T,Si,Sj,C,I>& m)
    { return m.conjugate(); }

    template <class T, StorageType S, IndexStyle I>
    static typename Matrix<T,S,I>::transpose_type Transpose(Matrix<T,S,I>& m)
    { return m.transpose(); }
    template <class T, int Si, int Sj, bool C, IndexStyle I>
    static typename MatrixView<T,Si,Sj,C,I>::transpose_type Transpose(
        MatrixView<T,Si,Sj,C,I>& m)
    { return m.transpose(); }

    template <class T, StorageType S, IndexStyle I>
    static typename Matrix<T,S,I>::adjoint_type Adjoint(Matrix<T,S,I>& m)
    { return m.adjoint(); }
    template <class T, int Si, int Sj, bool C, IndexStyle I>
    static typename MatrixView<T,Si,Sj,C,I>::adjoint_type Adjoint(
        MatrixView<T,Si,Sj,C,I>& m)
    { return m.adjoint(); }


    //
    // TMV_Text 
    //

    template <class T, StorageType S, IndexStyle I>
    static std::string TMV_Text(const Matrix<T,S,I>& )
    {
        std::ostringstream s;
        s << "Matrix<"<<TMV_Text(T())<<","<<TMV_Text(S)<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int Si, int Sj, bool C, IndexStyle I>
    static std::string TMV_Text(const ConstMatrixView<T,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "ConstMatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

    template <class T, int Si, int Sj, bool C, IndexStyle I>
    static std::string TMV_Text(const MatrixView<T,Si,Sj,C,I>& m)
    {
        std::ostringstream s;
        s << "MatrixView<"<<TMV_Text(T());
        s << ","<<IntTraits<Si>::text();
        if (Si == UNKNOWN) s << "("<<m.stepi()<<")";
        s << ","<<IntTraits<Sj>::text();
        if (Sj == UNKNOWN) s << "("<<m.stepj()<<")";
        s << ","<<C<<","<<TMV_Text(I)<<">";
        return s.str();
    }

} // namespace tmv

#endif
