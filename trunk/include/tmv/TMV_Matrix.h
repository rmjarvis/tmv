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
// An optional second template parameter specifies the known attributes
// of the matrix.  The valid options are:
// ColMajor or RowMajor
// CStyle or FortranStyle
// WithDivider or NoDivider
// The default is (ColMajor,CStyle,WithDivider).
//
// The options are treated as a bit field, so you | them together to 
// get the complete value.  e.g.
// Matrix<double,ColMajor | NoDivider | FortranStyle>
//
// Also, you only need to specify things that differ from the default, so
// Matrix<T,RowMajor> means Matrix<T,RowMajor|CStyle|WithDivider>
//
// The *Style options indicate whether to use C-style or Fortran-style indexing.
// With C-style (the default), the upper left corner of an MxN matrix is 
// m(0,0), and the lower right is m(M-1,N-1).  
// With Fortran-style, these are m(1,1) and m(M,N) respectively.  
// Also, when a function takes an index range, i1,i2, 
// then with C-style, this means elements from i1...i2-1 inclusive. 
// With Fortran-style, this means i1..i2 inclusive.
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
//    Matrix<T,A>(size_t colsize, size_t rowsize)
//        Makes a matrix with column size = colsize and row size = rowsize
//        with _uninitialized_ values
//
//    Matrix<T,A>(size_t colsize, size_t rowsize, T x)
//        Makes a matrix of size n with all values = x
//
//    Matrix<T,A>(const vector<vector<T> >& m)
//        Makes a matrix with a_ij = m[i][j]
//
//    Matrix<T,A>(size_t colsize, size_t rowsize, const T* vv)
//    Matrix<T,A>(size_t colsize, size_t rowsize, 
//            const std::vector<T>& vv)
//        Make a matrix which copies the elements of vv.
//        If A is tmv::RowMajor then the elements are taken in row major
//        order (m00,m01,..m0n,m10,m11...).  If A is tmv::ColMajor
//        then the elements are taken in column major order.
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
//    Both of these take two template arguments:
//    T = the underlying data type.
//    A = the known attributes of the matrix/
//        Options are:
//        NonConj || Conj
//        CStyle || FortranStyle
//        ColMajor || RowMajor || NonMajor
//        WithDivider || NoDivider
// 
//    The default attributes if you don't specify them are:
//    NonConj, CStyle, NonMajor, NoDivider
//
//    ConstMatrixView<T,A>(const T* p, size_t m, size_t n, int si, int sj)
//    MatrixView<T,A>(T* p, size_t m, size_t n, int si, int sj)
//        Make a matrix view starting at memory location p, 
//        with m rows and n cols, stepping over the data with step sizes 
//        si and sj along the columns and rows respectively.
//
//    You can also convert a MatrixView into a ConstMatrixView or
//    convert to a view with different attributes, so long as the 
//    conversion is logically valid.  e.g. you cannot just cast away
//    a Conj attribute, but you can change the indexing style or 
//    specify a storage attribute if you know that the matrix really
//    is row-major or column-major.
//
//
// Views:
//
//    All of these methods return some kind of view into the matrix
//    data.  The return types are either ConstMatrixView or MatrixView
//    with the appropriate attributes.
//    (Except for SubVector and LinearView, which return either 
//    ConstVectorView or VectorView.)
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
//    Note: These functions are only available for matrices that 
//    have the WithDivider attribute specified.
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

#include <vector>
#include "TMV_BaseMatrix_Rec.h"
#include "TMV_BaseMatrix_Tri.h"
#include "TMV_Divider.h"
#include "TMV_Array.h"
#include "TMV_VIt.h"

namespace tmv {


    //
    // Define the DivHelper class for rectangular matrices:
    //

    template <class M> class LUD;
    template <class T> class InstLUD;
    template <class M> class QRD;
    template <class T> class InstQRD;
    template <class M> class QRPD;
    template <class T> class InstQRPD;
    template <class M> class SVD;
    template <class T> class InstSVD;

    // The MatrixDivHelper class breaks up into two other classes:
    // MatrixDivHelper1<M,true> and MatrixDivHelper1<M,false>
    // depending on whether the Matrix class M is set to have a divider
    // or not.  
    //
    // If it does not, then MatrixDivHelper1<M,false> just sets up
    // the lud(), etc. functions to return a new object on the spot.
    //
    // If it does, then MatrixDivHelper1<M,true> derives from yet another
    // class that is a generic DivHelper<T> class and thus has all the
    // functions like divideUsing, setDiv, resetDiv, etc.
    
    template <class T>
    class MatrixDivHelper2 :
        public DivHelper<T>
    {
    public:
        // There is no inline version of this.  These functions are defined
        // in TMV_Matrix.cpp for the instantiated types.

        typedef const InstLUD<T>& lud_type;
        typedef const InstQRD<T>& qrd_type;
        typedef void qrpd_type;
        typedef void svd_type;
        //typedef InstQRPD<T> qrpd_type;
        //typedef InstSVD<T> svd_type;

        MatrixDivHelper2();
        ~MatrixDivHelper2();

        void setDiv() const;
        Matrix<T> getM() const;
        bool mIsSquare() const;

        lud_type lud() const;
        qrd_type qrd() const;
        qrpd_type qrpd() const;
        svd_type svd() const;

        virtual ConstMatrixView<T> getConstView() const=0;

    }; // MatrixDivHelper2

    template <class M, bool hasdivider>
    class MatrixDivHelper1;

    template <class M>
    class MatrixDivHelper1<M,true> :
        public MatrixDivHelper2<typename Traits<M>::value_type>
    {
    public:

        typedef typename Traits<M>::value_type T;

        TMV_INLINE ConstMatrixView<T> getConstView() const
        { return mat2().constView(); }

        // Use mat2() rather than mat() to avoid having to disambiguate
        // it with the normal BaseMatrix mat() function.
        TMV_INLINE const M& mat2() const
        { return static_cast<const M&>(*this); }
    };

    template <class M>
    class MatrixDivHelper1<M,false>
    {
    public:
        typedef typename Traits<M>::lud_type lud_type;
        typedef typename Traits<M>::qrd_type qrd_type;

        TMV_INLINE void resetDivType() const {} 
        TMV_INLINE lud_type lud() const { return lud_type(mat2(),false); }
        TMV_INLINE qrd_type qrd() const { return qrd_type(mat2(),false); }
        //TMV_INLINE qrd_type qrpd() const { return qrpd_type(mat2(),false); }
        //TMV_INLINE qrd_type svd() const { return svd_type(mat2(),false); }
        TMV_INLINE void qrpd() const {}
        TMV_INLINE void svd() const {}

        TMV_INLINE const M& mat2() const
        { return static_cast<const M&>(*this); }
    };

    template <class M>
    class MatrixDivHelper :
        public MatrixDivHelper1<M,Traits<M>::_hasdivider>
    {};

    template <class T, int A0, int A1>
    struct Traits<Matrix<T,A0,A1> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A01 = A0 | A1 };
        enum { A = A01 | (
                (Attrib<A01>::rowmajor ? 0 : ColMajor) |
                ( !Traits<real_type>::isinst ? NoDivider :
                  Traits<real_type>::isinteger ? NoDivider :
                  Attrib<A01>::nodivider ? 0 : WithDivider ) )};
        enum { okA = (
                !Attrib<A>::conj &&
                (Attrib<A>::rowmajor || Attrib<A>::colmajor) &&
                (Attrib<A>::rowmajor != int(Attrib<A>::colmajor) ) &&
                !Attrib<A>::diagmajor &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                !Attrib<A>::upper &&
                (Attrib<A>::nodivider || Attrib<A>::withdivider) &&
                (Attrib<A>::nodivider != int(Attrib<A>::withdivider) ) &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst && 
                    !Traits<real_type>::isinteger) ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef Matrix<T,A0,A1> type;
        typedef const type& calc_type;
        typedef const type& eval_type;
        typedef Matrix<T,A01> copy_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _shape = Rec };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = (_colmajor ? 1 : UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : UNKNOWN) };
        enum { _diagstep = UNKNOWN };
        enum { _conj = false };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _canlin = true };
        enum { twoSi = isreal ? int(_stepi) : IntTraits<_stepi>::twoS };
        enum { twoSj = isreal ? int(_stepj) : IntTraits<_stepj>::twoS };

        enum { _hasdivider = Attrib<A>::withdivider };
        typedef QuotXM<1,real_type,type> inverse_type;

        typedef typename TypeSelect<_hasdivider ,
                const InstLUD<value_type>& , LUD<type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRD<value_type>& , QRD<type> >::type qrd_type;
        typedef InvalidType qrpd_type;
        typedef InvalidType svd_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_checkalias ? 0 : NoAlias) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA };
        enum { twosA = isreal ? int(ndA) : (ndA & ~Conj & ~AllStorageType) };
        enum { Ac = _checkalias ? (ndA | CheckAlias) : (ndA & ~NoAlias) };
        enum { colpairAc = Ac & ~RowMajor };
        enum { rowpairAc = Ac & ~ColMajor };
        enum { twosAc = isreal ? int(Ac) : (Ac & ~Conj & ~AllStorageType) };

        typedef ConstVectorView<T,colA> const_col_type;
        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,vecA> const_diag_type;
        typedef ConstVectorView<T,vecA> const_diag_sub_type;

        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,nmA> const_submatrix_step_type;
        typedef ConstVectorView<T,vecA> const_subvector_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,colpairAc>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,rowpairAc>
            const_rowpair_type;
        typedef ConstMatrixView<T,ndA> const_colrange_type;
        typedef ConstMatrixView<T,ndA> const_rowrange_type;

        typedef ConstMatrixView<T,ndA> const_view_type;
        typedef ConstMatrixView<T,cstyleA> const_cview_type;
        typedef ConstMatrixView<T,fstyleA> const_fview_type;
        typedef ConstMatrixView<T> const_xview_type;
        typedef ConstMatrixView<T,cmA> const_cmview_type;
        typedef ConstMatrixView<T,rmA> const_rmview_type;
        typedef ConstMatrixView<T,conjA> const_conjugate_type;
        typedef ConstMatrixView<T,trA> const_transpose_type;
        typedef ConstMatrixView<T,adjA> const_adjoint_type;
        typedef ConstUpperTriMatrixView<T,(ndA|NonUnitDiag)>
            const_uppertri_type;
        typedef ConstUpperTriMatrixView<T,(ndA|UnitDiag)>
            const_unit_uppertri_type;
        typedef ConstUpperTriMatrixView<T,(ndA|UnknownDiag)>
            const_unknown_uppertri_type;
        typedef ConstLowerTriMatrixView<T,(ndA|NonUnitDiag)>
            const_lowertri_type;
        typedef ConstLowerTriMatrixView<T,(ndA|UnitDiag)>
            const_unit_lowertri_type;
        typedef ConstLowerTriMatrixView<T,(ndA|UnknownDiag)>
            const_unknown_lowertri_type;
        typedef typename TypeSelect< iscomplex ,
                ConstSmallMatrixView<real_type,UNKNOWN,UNKNOWN,twoSi,twoSj,twosAc> ,
                ConstMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,(vecA|Unit)> const_linearview_type;
        typedef ConstMatrixView<T,nonconjA> const_nonconj_type;
        typedef MatrixView<T,ndA> nonconst_type;

        typedef T& reference;
        typedef VIt<T,1,false> linear_iterator;

        typedef VectorView<T,colA> col_type;
        typedef VectorView<T,colA> col_sub_type;
        typedef VectorView<T,rowA> row_type;
        typedef VectorView<T,rowA> row_sub_type;
        typedef VectorView<T,vecA> diag_type;
        typedef VectorView<T,vecA> diag_sub_type;

        typedef MatrixView<T,ndA> submatrix_type;
        typedef MatrixView<T,nmA> submatrix_step_type;
        typedef VectorView<T,vecA> subvector_type;
        typedef SmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,colpairAc>
            colpair_type;
        typedef SmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,rowpairAc>
            rowpair_type;
        typedef MatrixView<T,ndA> colrange_type;
        typedef MatrixView<T,ndA> rowrange_type;

        typedef MatrixView<T,ndA> view_type;
        typedef MatrixView<T,cstyleA> cview_type;
        typedef MatrixView<T,fstyleA> fview_type;
        typedef MatrixView<T> xview_type;
        typedef MatrixView<T,cmA> cmview_type;
        typedef MatrixView<T,rmA> rmview_type;
        typedef MatrixView<T,conjA> conjugate_type;
        typedef MatrixView<T,trA> transpose_type;
        typedef MatrixView<T,adjA> adjoint_type;
        typedef UpperTriMatrixView<T,(ndA|NonUnitDiag)> uppertri_type;
        typedef UpperTriMatrixView<T,(ndA|UnitDiag)> unit_uppertri_type;
        typedef UpperTriMatrixView<T,(ndA|UnknownDiag)>
            unknown_uppertri_type;
        typedef LowerTriMatrixView<T,(ndA|NonUnitDiag)> lowertri_type;
        typedef LowerTriMatrixView<T,(ndA|UnitDiag)> unit_lowertri_type;
        typedef LowerTriMatrixView<T,(ndA|UnknownDiag)>
            unknown_lowertri_type;
        typedef typename TypeSelect< iscomplex ,
                SmallMatrixView<real_type,UNKNOWN,UNKNOWN,twoSi,twoSj,twosAc> ,
                MatrixView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<T,(vecA|Unit)> linearview_type;
        typedef MatrixView<T,nonconjA> nonconj_type;
    };

    template <class T, int A0, int A1>
    class Matrix : 
        public BaseMatrix_Rec_Mutable<Matrix<T,A0,A1> >,
        public MatrixDivHelper<Matrix<T,A0,A1> >
    {
    public:

        typedef Matrix<T,A0,A1> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        typedef BaseMatrix_Rec_Mutable<type> base_div;
        typedef MatrixDivHelper<type> divhelper;

        typedef typename Traits<T>::real_type real_type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _attrib = Traits<type>::_attrib };


        //
        // Constructors
        //

        Matrix() : 
            itscs(0), itsrs(0), linsize(0), itsm(0)
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        Matrix(size_t cs, size_t rs) :
            itscs(cs), itsrs(rs), linsize(cs*rs), itsm(linsize)
        {
            TMVStaticAssert(Traits<type>::okA);
#ifdef TMV_DEBUG
            this->linearView().flatten().setAllTo(Traits<real_type>::constr_value());
#endif
        }

        Matrix(size_t cs, size_t rs, T x) :
            itscs(cs), itsrs(rs), linsize(cs*rs), itsm(linsize)
        {
            TMVStaticAssert(Traits<type>::okA);
            this->setAllTo(x);
        }

        Matrix(size_t cs, size_t rs, const T* vv) :
            itscs(cs), itsrs(rs), linsize(cs*rs), itsm(linsize)
        {
            TMVStaticAssert(Traits<type>::okA);
            typename type::linearview_type lv = this->linearView();
            ConstVectorView<T,Unit>(vv,linsize).newAssignTo(lv);
        }

        Matrix(size_t cs, size_t rs, const std::vector<T>& vv) : 
            itscs(cs), itsrs(rs), linsize(cs*rs), itsm(linsize)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVAssert(vv.size() == linsize);
            typename type::linearview_type lv = this->linearView();
            ConstVectorView<T,Unit>(&vv[0],linsize).newAssignTo(lv);
        }

        Matrix(const std::vector<std::vector<T> >& vv) :
            itscs(vv.size()), itsrs(0), linsize(0), itsm(0)
        {
            TMVStaticAssert(Traits<type>::okA);
            if (itscs > 0) {
                resize(itscs,vv[0].size());
                for(size_t i=0;i<itscs;++i) {
                    TMVAssert(vv[i].size() == itsrs);
                    typename type::row_type mi = this->row(i);
                    ConstVectorView<T,Unit>(&vv[i][0],itsrs).newAssignTo(mi);
                }
            }
        }

        Matrix(const type& m2) :
            itscs(m2.itscs), itsrs(m2.itsrs),
            linsize(m2.linsize), itsm(linsize)
        {
            TMVStaticAssert(Traits<type>::okA);
            m2.newAssignTo(*this);
        }

        template <class M2>
        Matrix(const BaseMatrix<M2>& m2) :
            itscs(m2.colsize()), itsrs(m2.rowsize()),
            linsize(itscs * itsrs), itsm(linsize)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert((ShapeTraits2<M2::_shape,_shape>::assignable));
            m2.newAssignTo(*this);
        }

        template <class M2>
        Matrix(const BaseMatrix_Rec<M2>& m2) :
            itscs(m2.colsize()), itsrs(m2.rowsize()),
            linsize(m2.ls()), itsm(linsize)
        {
            TMVStaticAssert(Traits<type>::okA);
            m2.newAssignTo(*this);
        }

        ~Matrix() 
        {
#ifdef TMV_DEBUG
            this->linearView().flatten().setAllTo(Traits<real_type>::destr_value());
#endif
        }

        //
        // Op=
        //

        TMV_INLINE type& operator=(const type& m2)
        { if (this != &m2) base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(T x)
        { base_mut::operator=(x); return *this; }


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return itsm[_rowmajor ? i*stepi()+j : i+j*stepj()]; }

        T& ref(int i, int j)
        { return itsm[_rowmajor ? i*stepi()+j : i+j*stepj()]; }

        void swapWith(type& m2)
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            if (itsm.get() == m2.itsm.get()) return;
            itsm.swapWith(m2.itsm);
        }

        void resize(const size_t cs, const size_t rs)
        {
#ifdef TMV_DEBUG
            this->linearView().flatten().setAllTo(Traits<real_type>::destr_value());
#endif
            itscs = cs;
            itsrs = rs;
            divhelper::resetDivType();
            linsize = cs*rs;
            itsm.resize(linsize);
#ifdef TMV_DEBUG
            this->linearView().flatten().setAllTo(Traits<real_type>::constr_value());
#endif
        }

        TMV_INLINE size_t ls() const { return linsize; }
        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        int nElements() const { return itscs*itsrs; }
        TMV_INLINE int stepi() const { return _rowmajor ? itsrs : 1; }
        TMV_INLINE int stepj() const { return _rowmajor ? 1 : itscs; }
        TMV_INLINE bool isconj() const { return false; }
        TMV_INLINE bool isrm() const { return _rowmajor; }
        TMV_INLINE bool iscm() const { return _colmajor; }

    private:

        size_t itscs;
        size_t itsrs;
        size_t linsize;
        AlignedArray<T> itsm;

    }; // Matrix

    template <class T, int A0>
    struct Traits<ConstMatrixView<T,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider) };
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst && 
                    !Traits<real_type>::isinteger ) ) &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef ConstMatrixView<T,A0> type;
        typedef const type& calc_type;
        typedef const type& eval_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _shape = Rec };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = (_colmajor ? 1 : UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : UNKNOWN) };
        enum { _diagstep = UNKNOWN };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _canlin = false };
        enum { twoSi = isreal ? int(_stepi) : IntTraits<_stepi>::twoS };
        enum { twoSj = isreal ? int(_stepj) : IntTraits<_stepj>::twoS };

        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) |
                NoDivider ) };
        typedef Matrix<T,copyA> copy_type;
        enum { _hasdivider = Attrib<A>::withdivider };
        typedef QuotXM<1,real_type,type> inverse_type;

        typedef typename TypeSelect<_hasdivider ,
                const InstLUD<value_type>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRD<value_type>& , QRD<copy_type> >::type qrd_type;
        typedef InvalidType qrpd_type;
        typedef InvalidType svd_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) |
                (_checkalias ? 0 : NoAlias) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA & ~Conj };
        enum { twosA = isreal ? int(ndA) : (ndA & ~Conj & ~AllStorageType) };
        enum { Ac = _checkalias ? (ndA | CheckAlias) : (ndA & ~NoAlias) };
        enum { colpairAc = Ac & ~RowMajor };
        enum { rowpairAc = Ac & ~ColMajor };
        enum { twosAc = isreal ? int(Ac) : (Ac & ~Conj & ~AllStorageType) };

        typedef ConstVectorView<T,colA> const_col_type;
        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,vecA> const_diag_type;
        typedef ConstVectorView<T,vecA> const_diag_sub_type;

        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,nmA> const_submatrix_step_type;
        typedef ConstVectorView<T,vecA> const_subvector_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,colpairAc>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,rowpairAc>
            const_rowpair_type;
        typedef ConstMatrixView<T,ndA> const_colrange_type;
        typedef ConstMatrixView<T,ndA> const_rowrange_type;

        typedef ConstMatrixView<T,ndA> const_view_type;
        typedef ConstMatrixView<T,cstyleA> const_cview_type;
        typedef ConstMatrixView<T,fstyleA> const_fview_type;
        typedef ConstMatrixView<T,(_conj ? Conj : NonConj)> const_xview_type;
        typedef ConstMatrixView<T,cmA> const_cmview_type;
        typedef ConstMatrixView<T,rmA> const_rmview_type;
        typedef ConstMatrixView<T,conjA> const_conjugate_type;
        typedef ConstMatrixView<T,trA> const_transpose_type;
        typedef ConstMatrixView<T,adjA> const_adjoint_type;
        typedef ConstUpperTriMatrixView<T,(ndA|NonUnitDiag)>
            const_uppertri_type;
        typedef ConstUpperTriMatrixView<T,(ndA|UnitDiag)>
            const_unit_uppertri_type;
        typedef ConstUpperTriMatrixView<T,(ndA|UnknownDiag)>
            const_unknown_uppertri_type;
        typedef ConstLowerTriMatrixView<T,(ndA|NonUnitDiag)>
            const_lowertri_type;
        typedef ConstLowerTriMatrixView<T,(ndA|UnitDiag)>
            const_unit_lowertri_type;
        typedef ConstLowerTriMatrixView<T,(ndA|UnknownDiag)>
            const_unknown_lowertri_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                ConstSmallMatrixView<real_type,UNKNOWN,UNKNOWN,twoSi,twoSj,twosAc> ,
                ConstMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,(vecA|Unit)> const_linearview_type;
        typedef ConstMatrixView<T,nonconjA> const_nonconj_type;
        typedef MatrixView<T,ndA> nonconst_type;
    };

    template <class T, int A>
    class ConstMatrixView : 
        public BaseMatrix_Rec<ConstMatrixView<T,A> >,
        public MatrixDivHelper<ConstMatrixView<T,A> >
    {
    public:
        typedef ConstMatrixView<T,A> type;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        ConstMatrixView(
            const T* m, size_t cs, size_t rs, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) 
        { 
            TMVStaticAssert(Traits<type>::okA);
        }

        ConstMatrixView(const T* m, size_t cs, size_t rs, int si) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(_stepj) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepj != UNKNOWN); 
        }

        ConstMatrixView(const T* m, size_t cs, size_t rs) :
            itsm(m), itscs(cs), itsrs(rs), itssi(_stepi), itssj(_stepj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepi != UNKNOWN);
            TMVStaticAssert(_stepj != UNKNOWN); 
        }

        ConstMatrixView(const type& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        ConstMatrixView(const ConstMatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int A2>
        ConstMatrixView(const MatrixView<T,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int Si2, int Sj2, int A2>
        ConstMatrixView(
            const ConstSmallMatrixView<T,M2,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int Si2, int Sj2, int A2>
        ConstMatrixView(
            const SmallMatrixView<T,M2,N2,Si2,Sj2,A2>& m2) :
            itsm(m2.cptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj()) 
        { 
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

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

        TMV_INLINE const T* cptr() const { return itsm; }

        T cref(int i, int j) const
        { return DoConj<_conj>(itsm[i*stepi()+j*stepj()]); }

        TMV_INLINE size_t ls() const { return itscs*itsrs; }
        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        int nElements() const { return itscs*itsrs; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const 
        { return _rowmajor || (!_colmajor &&  stepj() == 1); }
        TMV_INLINE bool iscm() const 
        { return _colmajor || (!_rowmajor &&  stepi() == 1); }

    private :

        const T* itsm;
        const size_t itscs;
        const size_t itsrs;
        const CheckedInt<_stepi> itssi;
        const CheckedInt<_stepj> itssj;

    }; // ConstMatrixView

    template <class T, int A0>
    struct Traits<MatrixView<T,A0> >
    {
        typedef typename Traits<T>::real_type real_type;

        enum { A = (A0 & ~NoDivider) };
        enum { okA = (
                !Attrib<A>::diagmajor &&
                !Attrib<A>::nonunitdiag &&
                !Attrib<A>::unitdiag &&
                !Attrib<A>::zerodiag &&
                !Attrib<A>::packed &&
                !Attrib<A>::lower &&
                ( !Attrib<A>::withdivider ||
                  ( Traits<real_type>::isinst && 
                    !Traits<real_type>::isinteger ) ) &&
                ( Traits<T>::iscomplex || !Attrib<A>::conj ) )};
        enum { _attrib = A };

        typedef T value_type;
        typedef typename Traits<T>::complex_type complex_type;
        enum { isreal = Traits<T>::isreal };
        enum { iscomplex = Traits<T>::iscomplex };

        typedef MatrixView<T,A0> type;
        typedef ConstMatrixView<T,(A&~WithDivider)> calc_type;
        typedef const type& eval_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _shape = Rec };
        enum { _fort = Attrib<A>::fort };
        enum { _calc = true };
        enum { _rowmajor = Attrib<A>::rowmajor };
        enum { _colmajor = Attrib<A>::colmajor };
        enum { _stepi = (_colmajor ? 1 : UNKNOWN) };
        enum { _stepj = (_rowmajor ? 1 : UNKNOWN) };
        enum { _diagstep = UNKNOWN };
        enum { _conj = Attrib<A>::conj };
        enum { _checkalias = !Attrib<A>::noalias };
        enum { _canlin = false };
        enum { twoSi = isreal ? int(_stepi) : IntTraits<_stepi>::twoS };
        enum { twoSj = isreal ? int(_stepj) : IntTraits<_stepj>::twoS };

        enum { copyA = (
                (_rowmajor ? RowMajor : ColMajor) |
                (_fort ? FortranStyle : CStyle) |
                NoDivider ) };
        typedef Matrix<T,copyA> copy_type;
        enum { _hasdivider = Attrib<A>::withdivider };
        typedef QuotXM<1,real_type,type> inverse_type;

        typedef typename TypeSelect<_hasdivider ,
                const InstLUD<value_type>& , LUD<copy_type> >::type lud_type;
        typedef typename TypeSelect<_hasdivider ,
                const InstQRD<value_type>& , QRD<copy_type> >::type qrd_type;
        typedef InvalidType qrpd_type;
        typedef InvalidType svd_type;

        enum { vecA = (
                (_fort ? FortranStyle : 0) |
                (_conj ? Conj : 0 ) |
                (_checkalias ? 0 : NoAlias) )};
        enum { colA = vecA | (_colmajor ? Unit : 0) };
        enum { rowA = vecA | (_rowmajor ? Unit : 0) };
        enum { ndA = (A & ~AllDivStatus) };
        enum { cstyleA = ndA & ~FortranStyle };
        enum { fstyleA = ndA | FortranStyle };
        enum { nmA = (ndA & ~AllStorageType) };
        enum { cmA = nmA | ColMajor };
        enum { rmA = nmA | RowMajor };
        enum { conjA = iscomplex ? (ndA ^ Conj) : int(ndA) };
        enum { trA = _colmajor ? rmA : _rowmajor ? cmA : int(ndA) };
        enum { adjA = iscomplex ? (trA ^ Conj) : int(trA) };
        enum { nonconjA = ndA & ~Conj };
        enum { twosA = isreal ? int(ndA) : (ndA & ~Conj & ~AllStorageType) };
        enum { Ac = _checkalias ? (ndA | CheckAlias) : (ndA & ~NoAlias) };
        enum { colpairAc = Ac & ~RowMajor };
        enum { rowpairAc = Ac & ~ColMajor };
        enum { twosAc = isreal ? int(Ac) : (Ac & ~Conj & ~AllStorageType) };

        typedef ConstVectorView<T,colA> const_col_type;
        typedef ConstVectorView<T,colA> const_col_sub_type;
        typedef ConstVectorView<T,rowA> const_row_type;
        typedef ConstVectorView<T,rowA> const_row_sub_type;
        typedef ConstVectorView<T,vecA> const_diag_type;
        typedef ConstVectorView<T,vecA> const_diag_sub_type;

        typedef ConstMatrixView<T,ndA> const_submatrix_type;
        typedef ConstMatrixView<T,nmA> const_submatrix_step_type;
        typedef ConstVectorView<T,vecA> const_subvector_type;
        typedef ConstSmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,colpairAc>
            const_colpair_type;
        typedef ConstSmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,rowpairAc>
            const_rowpair_type;
        typedef ConstMatrixView<T,ndA> const_colrange_type;
        typedef ConstMatrixView<T,ndA> const_rowrange_type;

        typedef ConstMatrixView<T,ndA> const_view_type;
        typedef ConstMatrixView<T,cstyleA> const_cview_type;
        typedef ConstMatrixView<T,fstyleA> const_fview_type;
        typedef ConstMatrixView<T,(_conj ? Conj : NonConj)> const_xview_type;
        typedef ConstMatrixView<T,cmA> const_cmview_type;
        typedef ConstMatrixView<T,rmA> const_rmview_type;
        typedef ConstMatrixView<T,conjA> const_conjugate_type;
        typedef ConstMatrixView<T,trA> const_transpose_type;
        typedef ConstMatrixView<T,adjA> const_adjoint_type;
        typedef ConstUpperTriMatrixView<T,(ndA|NonUnitDiag)>
            const_uppertri_type;
        typedef ConstUpperTriMatrixView<T,(ndA|UnitDiag)>
            const_unit_uppertri_type;
        typedef ConstUpperTriMatrixView<T,(ndA|UnknownDiag)>
            const_unknown_uppertri_type;
        typedef ConstLowerTriMatrixView<T,(ndA|NonUnitDiag)>
            const_lowertri_type;
        typedef ConstLowerTriMatrixView<T,(ndA|UnitDiag)>
            const_unit_lowertri_type;
        typedef ConstLowerTriMatrixView<T,(ndA|UnknownDiag)>
            const_unknown_lowertri_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                ConstSmallMatrixView<real_type,UNKNOWN,UNKNOWN,twoSi,twoSj,twosAc> ,
                ConstMatrixView<real_type,twosA> >::type const_realpart_type;
        typedef const_realpart_type const_imagpart_type;
        typedef ConstVectorView<T,(vecA|Unit)> const_linearview_type;
        typedef ConstMatrixView<T,nonconjA> const_nonconj_type;
        typedef MatrixView<T,ndA> nonconst_type;

        typedef typename AuxRef<T,_conj>::reference reference;
        typedef VIt<T,1,_conj> linear_iterator;

        typedef VectorView<T,colA> col_type;
        typedef VectorView<T,colA> col_sub_type;
        typedef VectorView<T,rowA> row_type;
        typedef VectorView<T,rowA> row_sub_type;
        typedef VectorView<T,vecA> diag_type;
        typedef VectorView<T,vecA> diag_sub_type;

        typedef MatrixView<T,ndA> submatrix_type;
        typedef MatrixView<T,nmA> submatrix_step_type;
        typedef VectorView<T,vecA> subvector_type;
        typedef SmallMatrixView<T,UNKNOWN,2,_stepi,UNKNOWN,colpairAc>
            colpair_type;
        typedef SmallMatrixView<T,2,UNKNOWN,UNKNOWN,_stepj,rowpairAc>
            rowpair_type;
        typedef MatrixView<T,ndA> colrange_type;
        typedef MatrixView<T,ndA> rowrange_type;

        typedef MatrixView<T,ndA> view_type;
        typedef MatrixView<T,cstyleA> cview_type;
        typedef MatrixView<T,fstyleA> fview_type;
        typedef MatrixView<T,(_conj ? Conj : NonConj)> xview_type;
        typedef MatrixView<T,cmA> cmview_type;
        typedef MatrixView<T,rmA> rmview_type;
        typedef MatrixView<T,conjA> conjugate_type;
        typedef MatrixView<T,trA> transpose_type;
        typedef MatrixView<T,adjA> adjoint_type;
        typedef UpperTriMatrixView<T,(ndA|NonUnitDiag)> uppertri_type;
        typedef UpperTriMatrixView<T,(ndA|UnitDiag)> unit_uppertri_type;
        typedef UpperTriMatrixView<T,(ndA|UnknownDiag)>
            unknown_uppertri_type;
        typedef LowerTriMatrixView<T,(ndA|NonUnitDiag)> lowertri_type;
        typedef LowerTriMatrixView<T,(ndA|UnitDiag)> unit_lowertri_type;
        typedef LowerTriMatrixView<T,(ndA|UnknownDiag)>
            unknown_lowertri_type;
        typedef typename TypeSelect< (iscomplex && (_colmajor||_rowmajor)) ,
                SmallMatrixView<real_type,UNKNOWN,UNKNOWN,twoSi,twoSj,twosAc> ,
                MatrixView<real_type,twosA> >::type realpart_type;
        typedef realpart_type imagpart_type;
        typedef VectorView<T,(vecA|Unit)> linearview_type;
        typedef MatrixView<T,nonconjA> nonconj_type;
    };

    template <class T, int A>
    class MatrixView : 
        public BaseMatrix_Rec_Mutable<MatrixView<T,A> >,
        public MatrixDivHelper<MatrixView<T,A> >
    {
    public:

        typedef MatrixView<T,A> type;
        typedef BaseMatrix_Rec_Mutable<type> base_mut;
        typedef typename base_mut::reference reference;

        enum { _colsize = Traits<type>::_colsize };
        enum { _rowsize = Traits<type>::_rowsize };
        enum { _shape = Traits<type>::_shape };
        enum { _fort = Traits<type>::_fort };
        enum { _calc = Traits<type>::_calc };
        enum { _rowmajor = Traits<type>::_rowmajor };
        enum { _colmajor = Traits<type>::_colmajor };
        enum { _stepi = Traits<type>::_stepi };
        enum { _stepj = Traits<type>::_stepj };
        enum { _diagstep = Traits<type>::_diagstep };
        enum { _conj = Traits<type>::_conj };
        enum { _checkalias = Traits<type>::_checkalias };
        enum { _canlin = Traits<type>::_canlin };
        enum { _attrib = Traits<type>::_attrib };

        //
        // Constructors
        //

        MatrixView(T* m, size_t cs, size_t rs, int si, int sj) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(sj) 
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        MatrixView(T* m, size_t cs, size_t rs, int si) :
            itsm(m), itscs(cs), itsrs(rs), itssi(si), itssj(_stepj) 
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepj != UNKNOWN); 
        }

        MatrixView(T* m, size_t cs, size_t rs) :
            itsm(m), itscs(cs), itsrs(rs), itssi(_stepi), itssj(_stepj)
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(_stepi != UNKNOWN);
            TMVStaticAssert(_stepj != UNKNOWN); 
        }

        MatrixView(const type& m2) :
            itsm(m2.itsm), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj())
        {
            TMVStaticAssert(Traits<type>::okA);
        }

        template <int A2>
        MatrixView(MatrixView<T,A2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        template <int M2, int N2, int Si2, int Sj2, int A2>
        MatrixView(SmallMatrixView<T,M2,N2,Si2,Sj2,A2> m2) :
            itsm(m2.ptr()), itscs(m2.colsize()), itsrs(m2.rowsize()),
            itssi(m2.stepi()), itssj(m2.stepj())
        {
            TMVStaticAssert(Traits<type>::okA);
            TMVStaticAssert(Attrib<A>::conj == int(Attrib<A2>::conj)); 
        }

        ~MatrixView() {
#ifdef TMV_DEBUG
            itsm = 0; 
#endif
        }


        //
        // Op = 
        //

        TMV_INLINE type& operator=(const type& m2)
        { if (this != &m2) base_mut::operator=(m2); return *this; }

        template <class M2>
        TMV_INLINE type& operator=(const BaseMatrix<M2>& m2)
        { base_mut::operator=(m2); return *this; }

        TMV_INLINE type& operator=(const T x)
        { base_mut::operator=(x); return *this; }


        //
        // Auxilliary Functions
        //

        TMV_INLINE const T* cptr() const { return itsm; }
        TMV_INLINE T* ptr() { return itsm; }

        T cref(int i, int j) const
        { return DoConj<_conj>(itsm[i*stepi()+j*stepj()]); }

        reference ref(int i, int j) 
        { return reference(itsm[i*stepi()+j*stepj()]); }

        TMV_INLINE size_t ls() const { return itscs*itsrs; }
        TMV_INLINE size_t colsize() const { return itscs; }
        TMV_INLINE size_t rowsize() const { return itsrs; }
        int nElements() const { return itscs*itsrs; }
        TMV_INLINE int stepi() const { return itssi; }
        TMV_INLINE int stepj() const { return itssj; }
        TMV_INLINE bool isconj() const { return _conj; }
        TMV_INLINE bool isrm() const 
        { return _rowmajor || (!_colmajor &&  stepj() == 1); }
        TMV_INLINE bool iscm() const 
        { return _colmajor || (!_rowmajor &&  stepi() == 1); }

    private :

        T* itsm;
        const size_t itscs;
        const size_t itsrs;
        const CheckedInt<_stepi> itssi;
        const CheckedInt<_stepj> itssj;

    }; // MatrixView


    //-------------------------------------------------------------------------

    //
    // Special Creators: 
    //   MatrixViewOf(m,colsize,rowsize,stor) = MatrixView of m 
    //   MatrixViewOf(m,colsize,rowsize,stepi,stepj) = MatrixView of m 
    //

    // MatrixView of raw memory:
    template <class T>
    static inline MatrixView<T> MatrixViewOf(
        T* m, size_t colsize, size_t rowsize, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor) 
            return MatrixView<T>(m,colsize,rowsize,rowsize,1);
        else 
            return MatrixView<T>(m,colsize,rowsize,1,colsize);
    }

    template <class T>
    static inline ConstMatrixView<T> MatrixViewOf(
        const T* m, size_t colsize, size_t rowsize, StorageType stor)
    {
        TMVAssert(stor == RowMajor || stor == ColMajor);
        if (stor == RowMajor)
            return ConstMatrixView<T>(m,colsize,rowsize,rowsize,1);
        else 
            return ConstMatrixView<T>(m,colsize,rowsize,1,colsize);
    }

    template <class T>
    static TMV_INLINE MatrixView<T> MatrixViewOf(
        T* m, size_t colsize, size_t rowsize, int stepi, int stepj)
    { return MatrixView<T>(m,colsize,rowsize,stepi,stepj); }

    template <class T>
    static TMV_INLINE ConstMatrixView<T> MatrixViewOf(
        const T* m, size_t colsize, size_t rowsize, int stepi, int stepj)
    { return ConstMatrixView<T>(m,colsize,rowsize,stepi,stepj); }


    //
    // Swap
    //

    template <class T, int A0, int A1>
    static TMV_INLINE void Swap(Matrix<T,A0,A1>& m1, Matrix<T,A0,A1>& m2)
    { m1.swapWith(m2); }
    template <class M, class T, int A>
    static TMV_INLINE void Swap(
        BaseMatrix_Rec_Mutable<M>& m1, MatrixView<T,A> m2)
    { DoSwap(m1,m2); }
    template <class M, class T, int A>
    static TMV_INLINE void Swap(
        MatrixView<T,A> m1, BaseMatrix_Rec_Mutable<M>& m2)
    { DoSwap(m1,m2); }
    template <class T, int A1, int A2>
    static TMV_INLINE void Swap(MatrixView<T,A1> m1, MatrixView<T,A2> m2)
    { DoSwap(m1,m2); }


    //
    // Conjugate, Transpose, Adjoint
    //

    template <class T, int A0, int A1>
    static TMV_INLINE typename Matrix<T,A0,A1>::conjugate_type Conjugate(
        Matrix<T,A0,A1>& m)
    { return m.conjugate(); }
    template <class T, int A>
    static TMV_INLINE typename MatrixView<T,A>::conjugate_type Conjugate(
        MatrixView<T,A>& m)
    { return m.conjugate(); }

    template <class T, int A0, int A1>
    static TMV_INLINE typename Matrix<T,A0,A1>::transpose_type Transpose(
        Matrix<T,A0,A1>& m)
    { return m.transpose(); }
    template <class T, int A>
    static TMV_INLINE typename MatrixView<T,A>::transpose_type Transpose(
        MatrixView<T,A>& m)
    { return m.transpose(); }

    template <class T, int A0, int A1>
    static TMV_INLINE typename Matrix<T,A0,A1>::adjoint_type Adjoint(
        Matrix<T,A0,A1>& m)
    { return m.adjoint(); }
    template <class T, int A>
    static TMV_INLINE typename MatrixView<T,A>::adjoint_type Adjoint(
        MatrixView<T,A>& m)
    { return m.adjoint(); }


    //
    // TMV_Text 
    //

#ifdef TMV_TEXT
    template <class T, int A0, int A1>
    static inline std::string TMV_Text(const Matrix<T,A0,A1>& m)
    {
        const int A = A0 | A1;
        std::ostringstream s;
        s << "Matrix<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int A>
    static inline std::string TMV_Text(const ConstMatrixView<T,A>& m)
    {
        std::ostringstream s;
        s << "ConstMatrixView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }

    template <class T, int A>
    static inline std::string TMV_Text(const MatrixView<T,A>& m)
    {
        std::ostringstream s;
        s << "MatrixView<"<<TMV_Text(T());
        s << ","<<Attrib<A>::text()<<">";
        s << "("<<m.colsize()<<","<<m.rowsize()<<",";
        s << m.stepi()<<","<<m.stepj()<<")";
        return s.str();
    }
#endif

} // namespace tmv

#endif
