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
// The Matrix class and all associated functions are contained
// in the namespace tmv.  Alse, the Matrix class is a template, 
// so for a Matrix of long doubles, one would write 
// tmv::Matrix<long double>.  
//
// There is also an optional second template parameter which is 
// either RowMajor or ColMajor, referring to the storage pattern of 
// the matrix elements.
//
// Finally, a third optional template parameter is either CStyle or
// FortranStyle to indicate whether to use C-style or Fortran-style indexing.  
// With C-style (the default), the upper left corner of an MxN matrix is m(0,0), 
// and the lower right is m(M-1,N-1).  With Fortran-style, these are m(1,1) and
// m(M,N) respectively.  Also, when a function takes an index range, i1,i2, 
// then with C-style, this means elements from i1...i2-1 inclusive. 
// With Fortran-style, this means i1..i2 inclusive.
//
// Note that if you want to specify Fortran-style, then you must also 
// specify the storage pattern, since that parameter precedes it.  For example:
// Matrix<T,RowMajor,FortranStyle> m(10,10);
//
// Also, GenMatrix and MatrixComposite always use C-style indexing.
// However, they are both castable as a Fortran-style matrix.  Just remember
// to do so if you want to use Fortran-style indexing.
//
// Constructors:
//
//    Matrix<T,stor,I>(int colsize, int rowsize)
//        Makes a Matrix with column size = colsize and row size = rowsize
//        with _uninitialized_ values
//
//    Matrix<T,stor,I>(int colsize, int rowsize, T x)
//        Makes a Matrix of size n with all values = x
//
//
// Special Constructors
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
//    MatrixView MatrixViewOf(T* m, int colsize, int rowsize,
//            StorageType stor)
//    MatrixView MatrixViewOf(T* m, int colsize, int rowsize,
//            int stepi, int stepj)
//    ConstMatrixView MatrixViewOf(const T* m, int colsize, int rowsize,
//            StorageType stor)
//    ConstMatrixView MatrixViewOf(const T* m, int colsize, int rowsize,
//            int stepi, int stepj)
//        Returns a MatrixView of the elements in m, using the actual
//        elements m for the storage.  This is essentially the same as the 
//        constructor with (const T*m), except that the data isn't duplicated.
//
// Access Functions
//
//    int colsize() const
//    int rowsize() const
//        Return the dimensions of the Matrix
//
//    T& operator()(int i, int j)
//    T operator()(int i, int j) const
//        Return the (i,j) element of the Matrix
//
//    VectorView& operator[](int i)
//    ConstVectorView& operator[](int i) const
//        Return the ith row of the Matrix as a Vector.
//        This allows for m[i][j] style access.
//
//    VectorView& row(int i, int j1, int j2)
//    ConstVectorView& row(int i, int j1, int j2) const
//        Return the ith row of the Matrix as a Vector
//        If j1,j2 are given, it returns the subVector from j1 to j2 
//        (not including j2) within the row.
//
//    VectorView& col(int j, int i1, int i2)
//    ConstVectorView& col(int j) const
//        Return the jth column of the Matrix as a Vector
//        If i1,i2 are given, it returns the subVector from i1 to i2 
//        (not including i2) within the column.
//
//    VectorView& diag()
//    ConstVectorView& diag() const
//        Return the diagonal of the Matrix as a Vector
//
//    VectorView& diag(int i, int j1, int j2)
//    ConstVectorView& diag(i, j1, j2) const
//        Return the super- or sub-diagonal i
//        If i > 0 return the super diagonal starting at m_0i
//        If i < 0 return the sub diagonal starting at m_|i|0
//        If j1,j2 are given, it returns the diagonal subVector 
//        either from m_j1,i+j1 to m_j2,i+j2 (for i>0) 
//        or from m_|i|+j1,j1 to m_|i|+j2,j2 (for i<0)
//
// Modifying Functions
//
//    Matrix& setZero()
//        Sets all elements to 0
//
//    Matrix& setAllTo(T x)
//        Sets all elements to x
//
//    Matrix& addToAll(T x)
//        Adds x to all elements
//
//    Matrix& clip(RT thresh)
//        Set to 0 all elements whose abolute value is < thresh
//
//    Matrix<T>& transposeSelf() 
//        Transposes the elements of a square Matrix or subMatrix
//
//    Matrix& conjugateSelf()
//        Sets all elements to its conjugate
//
//    Matrix& setToIdentity(x = 1)
//        Set to Identity Matrix, or 
//        with a parameter, set to x times Identity Matrix
//
//    Matrix& swapRows(i1, i2)
//        Swap two rows
//
//    Matrix& swapCols(j1, j2)
//        Swap two columns
//
//    Matrix& permuteRows(const int* p)
//    Matrix& permuteRows(const int* p, int i1, int i2)
//        Perform a series of row swaps (0,p[0]), (1,p[1])...
//        In the second case, only do (i1,p[i1)...(i2-1,p[i2-1])
//    Matrix& reversePermuteRows(const int* p)
//    Matrix& reversePermuteRows(const int* p, int i1, int i2)
//        The same, but perform the swaps in reverse order
//
//    Matrix& permuteCols(const int* p)
//    Matrix& permuteCols(const int* p, int j1, int j2)
//        Perform a series of column swaps (0,p[0]), (1,p[1])...
//        In the second case, only do (j1,p[j1)...(j2-1,p[j2-1])
//    Matrix& reversePermuteCols(const int* p)
//    Matrix& reversePermuteCols(const int* p, int j1, int j2)
//        The same, but perform the swaps in reverse order
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
//    MatrixView subMatrix(int i1, int i2, int j1, int j2,
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
//    VectorView subVector(int i, int j, int istep, int jstep, int size)
//        Returns a subVector which starts at position (i,j) in the 
//        matrix, moves in the directions (istep,jstep) and has a length
//        of size.
//        For example, the cross-diagonal from the lower left to the upper
//        right of a 6x6 matrix could be accessed using:
//        m.subVector(5,0,-1,1,6)
//
//    UpperTriMatrixView upperTri()
//    LowerTriMatrixView lowerTri()
//        Returns a view of the upper or lower triangle portion of the Matrix.
//
//    UpperTriMatrixView unitUpperTri()
//    LowerTriMatrixView unitLowerTri()
//        Returns a view of the upper or lower triangle portion of the Matrix
//        with unit-diagonal elements, rather than what is in the matrix.
//
//    MatrixView colPair(int j1, int j2)
//        This returns an Mx2 submatrix which consists of the 
//        columns j1 and j2.  This is useful for multiplying two 
//        (not necessarily adjacent) columns of a matrix by a 2x2 matrix.
//
//    MatrixView rowPair(int i1, int i2)
//        Same as above, but two rows.
//
//    MatrixView colRange(int j1, int j2)
//        This is shorthand for subMatrix(0,m.colsize(),j1,j2)
//
//    MatrixView rowRange(int i1, int i2)
//        This is shorthand for subMatrix(i1,i2,0,m.rowsize())
//
//    MatrixView realPart()
//    MatrixView imagPart()
//        Returns a view to the real/imag elements of a complex Matrix.
//
//    MatrixView view()
//        Returns a view of a matrix.
//
//    MatrixView conjugate()
//        Returns a view to the conjugate of a Matrix.
//
//    MatrixView transpose()
//        Returns a view to the transpose of a Matrix.
//
//    MatrixView adjoint()
//        Returns a view to the adjoint (conjugate transpose) of a Matrix.
//        Note: Some people define the adjoint as the cofactor matrix.
//              This is not the same as our definition of the adjoint.
//
//    ConstVectorView constLinearView()
//    VectorView linearView()
//        Returns a VectorView with all the elements of the Matrix.
//        This is mostly used internally for things like MaxElement
//        and conjugateSelf, where the matrix structure is irrelevant,
//        and we just want to do something to all the elements.
//        The corrolary function canLinearize() returns whether this is 
//        allowed.
//
// Functions of Matrixs:
//
//    Det(m)
//    m.det()
//        Returns the determinant of a Matrix.
//        Note: If the matrix is not square, the determinant is not
//              well defined.  The returned value is such that
//              conj(det) * det = det(adjoint(m) * m)
//              So for real nonsquare matrices, the sign is arbitrary,
//              and for complex nonsquare matrices, it is multiplied
//              by an arbitrary phase.
//
//    LogDet(m)
//    m.logDet(T* sign=0)
//        Returns the logarithm of the absolute value of the determinant.
//        For many large matrices, the determinant yields to overflow.
//        Hence, this function is provided, which stably calculates the
//        natural logarithm of the absolute value of the determinant.
//        The optional sign argument returns the sign of the determinant
//        if T is real, or the factor exp(it) factor by which exp(logdet) 
//        would have to be multiplied to get the actual determinant.
//
//    Trace(m)
//    m.trace()
//        Returns the trace of a Matrix.
//        = sum_i ( a_ii )
//
//    SumElements(m)
//    m.sumElements()
//        Returns the sum of the elements of a Matrix.
//        = sum_ij ( a_ij )
//
//    SumAbsElements(m)
//    m.sumAbsElements()
//        Returns the sum of the absolute values of the elements of a Matrix.
//        = sum_ij |a_ij|
//
//    SumAbs2Elements(m)
//    m.sumAbs2Elements()
//        Returns the sum of the absolute values of the elements of a Matrix.
//        = sum_ij |real(a_ij)| + |imag(a_ij)|
//
//    Norm(m) or NormF(m)
//    m.norm() or m.normF()
//        Return the Frobenius norm of a Matrix.
//        = sqrt( sum_ij |a_ij|^2 )
//
//    NormSq(m)
//    m.normSq(RT scale = 1.)
//        Returns the square of norm().
//        = sum_ij |a_ij|^2
//
//    Norm1(m) 
//    m.norm1() 
//        Returns the 1-norm of a Matrix.
//        = max_j (sum_i |a_ij|)
//
//    Norm2(m) 
//    m.norm2() 
//        Returns the 2-norm of a Matrix.
//        = sqrt( Max Eigenvalue of (A.adjoint * A) )
//        = Maximum singular value
//        Note: This norm is costly to calculate if one is not 
//              otherwise doing a singular value decomposition
//              of the Matrix.
//
//    NormInf(m) 
//    m.normInf() 
//        Returns the infinity-norm of a Matrix.
//        = max_i (sum_j |a_ij|)
//
//    m.makeInverse(minv)
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
//    m.makeInverseATA(invata)
//        Sets invata to the Inverse of (A.adjoint * A) for matrix m = A
//        If Ax=b is solved for x, then (AT A)^-1 is the 
//        covariance matrix of the least-squares solution x
//
//    m.inverse()
//    Inverse(m)
//        Returns an auxiliary object that delays the calculation of the 
//        inverse until there is appropriate storage for it.
//        m.makeInverse(minv) is equivalent to minv = m.inverse().
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
//       Most of these behave in the logical way for dealing with matrices.
//       Two comments about behavior that might not be obvious:
//
//       1) Vectors are either row or column Vectors as appropriate.
//          eg. For m*v, v is a column Vector, but for v*m, v is a row Vector
//
//       2) Sometimes x should be thought of a x*I, where I is an appropriately
//          sized identity matrix.  For example m-1 means m-I,
//          and m += 3 means m += 3*I (or equivalently, m.diag().addToAll(3)).
//
//       3) Division by a matrix can be from either the left or the right.
//          ie. v/m can mean either m^-1 v or v m^-1
//          The first case is the solution of mx=v, the second is of xm = v.
//          Since the first case is the way one generally poses a problem
//          for solving a set of equations, we take v/m to be left-division.
//          If you want right-division (v m^-1), then we supply the % operator
//          to do so.  
//          ie. v%m means v m^-1
//          If you want to be more explicit, you can write:
//          v/m as m.inverse() * v and v%m as v * m.inverse().
//          In all cases, the actual calculation is delayed until there is
//          storage to put it.  (Unless you string too many calculations 
//          together, in which case it will use a temporary.)
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
//
//
// Division Control Functions:
//
//    There are a number of algorithms available for dividing
//    matrices.  We provide functions to allow you to 
//    change the algorithm used by the code on the fly.
//    In particular, you can write:
//    m.divideUsing(dt)
//    where dt is LU, QR, QRP, or SV
//    (ie. anything but CH)
//    Each of these also has an in-place version whcih overwrites the
//    current Matrix memory with the decomposition needed for 
//    doing the division.  Obviously, if you try to use the Matrix
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
//    If dt = LU, and the decomposition has been set, then 
//    lud() returns the LUDiv<T> class
//    If it is not set yet, or divideUsing was called with some other dt,
//    then lud() sets the divider to LU for you and returns the 
//    LU decomposition class.
//
//    Likewise:
//    qrd(), qrpd(), svd() return the corresponding Divider classes.
//


#ifndef TMV_Matrix_H
#define TMV_Matrix_H

#include "tmv/TMV_BaseMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_Permutation.h"
#include "tmv/TMV_Array.h"
#include "tmv/TMV_MIt.h"
#include <vector>

namespace tmv {

    template <class T1, class T2> 
    inline void Copy(const GenMatrix<T1>& m1, const MatrixView<T2>& m2);

    template <class T, class Tm> 
    class QuotXM;

    template <class T> 
    class LUDiv;

    template <class T> 
    class QRDiv;

    template <class T> 
    class QRPDiv;

    template <class T> 
    class SVDiv;

    template <class T> 
    class GenMatrix : 
        public BaseMatrix<T>,
        public DivHelper<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef T value_type;
        typedef RT real_type;
        typedef CT complex_type;
        typedef GenMatrix<T> type;
        typedef Matrix<T> copy_type;
        typedef ConstVectorView<T> const_vec_type;
        typedef ConstMatrixView<T> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstMatrixView<RT> const_realpart_type;
        typedef ConstUpperTriMatrixView<T> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T> const_lowertri_type;
        typedef MatrixView<T> nonconst_type;
        typedef CRMIt<type> const_rowmajor_iterator;
        typedef CCMIt<type> const_colmajor_iterator;

        //
        // Constructors
        //

        inline GenMatrix() {}
        inline GenMatrix(const type&) {}
        virtual inline ~GenMatrix() {}

        //
        // Access Functions
        //

        using AssignableToMatrix<T>::colsize;
        using AssignableToMatrix<T>::rowsize;

        inline const_vec_type row(int i) const 
        { 
            TMVAssert(i>=0 && i<colsize());
            return const_vec_type(
                cptr()+i*stepi(),rowsize(),stepj(),ct()); 
        }

        inline const_vec_type col(int j) const
        {
            TMVAssert(j>=0 && j<rowsize());
            return const_vec_type(
                cptr()+j*stepj(),colsize(),stepi(),ct()); 
        }

        inline const_vec_type diag() const
        {
            return const_vec_type(
                cptr(),TMV_MIN(colsize(),rowsize()),stepi()+stepj(),ct());
        }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const int diagstep = stepi() + stepj();
            if (i >= 0) {
                const int diagsize = TMV_MIN(colsize(),rowsize()-i);
                return const_vec_type(
                    cptr()+i*stepj(),diagsize,diagstep,ct());
            } else {
                const int diagsize = TMV_MIN(colsize()+i,rowsize());
                return const_vec_type(
                    cptr()-i*stepi(),diagsize,diagstep,ct());
            }
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            return const_vec_type(
                cptr()+i*stepi()+j1*stepj(),j2-j1,stepj(),ct()); 
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<rowsize());
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            return const_vec_type(
                cptr()+i1*stepi()+j*stepj(),i2-i1,stepi(),ct()); 
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const int diagstep = stepi() + stepj();
            if (i >= 0) {
                TMVAssert(j1>=0 && j1-j2<=0 && 
                          j2<=TMV_MIN(colsize(),rowsize()-i));
                return const_vec_type(
                    cptr()+i*stepj()+j1*diagstep,j2-j1, diagstep,ct());
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && 
                          j2<=TMV_MIN(colsize()+i,rowsize()));
                return const_vec_type(
                    cptr()-i*stepi()+j1*diagstep,j2-j1, diagstep,ct());
            }
        }

        inline T operator()(int i, int j) const
        { 
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j>=0 && j<rowsize());
            return cref(i,j);
        }

        inline const_vec_type operator[](int i) const
        { 
            TMVAssert(i>=0 && i<colsize());
            return row(i); 
        }

        template <class T2> 
        inline bool isSameAs(const BaseMatrix<T2>&) const
        { return false; }

        inline bool isSameAs(const type& m2) const
        {
            return 
                this==&m2 || 
                ( cptr()==m2.cptr() && 
                  rowsize()==m2.rowsize() && colsize()==m2.colsize() &&
                  stepi()==m2.stepi() && stepj()==m2.stepj() && ct()==m2.ct()
                );
        }

        inline void assignToM(const MatrixView<RT>& m2) const
        { 
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            if (!isSameAs(m2)) Copy(*this,m2);
        }

        inline void assignToM(const MatrixView<CT>& m2) const
        { 
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            if (!isSameAs(m2)) Copy(*this,m2);
        }

        //
        // subMatrix
        //

        bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const;

        inline const_view_type cSubMatrix(int i1, int i2, int j1, int j2) const
        {
            return const_view_type(
                cptr()+i1*stepi()+j1*stepj(),
                i2-i1,j2-j1,stepi(),stepj(),stor(),ct());
        }

        inline const_view_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline const_view_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            StorageType newstor = 
                iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            return const_view_type(
                cptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),
                newstor,ct());
        }

        inline const_view_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        bool hasSubVector(int i, int j, int istep, int jstep, int s) const;

        inline const_vec_type cSubVector(
            int i, int j, int istep, int jstep, int s) const
        {
            return const_vec_type(
                cptr()+i*stepi()+j*stepj(),s,
                istep*stepi()+jstep*stepj(),ct());
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int s) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,s));
            return cSubVector(i,j,istep,jstep,s);
        }

        inline const_uppertri_type unitUpperTri() const
        {
            TMVAssert(rowsize() <= colsize());
            return const_uppertri_type(
                cptr(),rowsize(),stepi(),stepj(),UnitDiag,stor(),ct() );
        }

        inline const_uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        {
            TMVAssert(rowsize() <= colsize());
            return const_uppertri_type(
                cptr(),rowsize(),stepi(),stepj(),dt,stor(),ct() );
        }

        inline const_lowertri_type unitLowerTri() const
        {
            TMVAssert(colsize() <= rowsize());
            return const_lowertri_type(
                cptr(),colsize(),stepi(),stepj(),UnitDiag,stor(),ct() );
        }

        inline const_lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        {
            TMVAssert(colsize() <= rowsize());
            return const_lowertri_type(
                cptr(),colsize(),stepi(),stepj(),dt,stor(),ct() );
        }

        inline const_view_type cColPair(int j1, int j2) const
        {
            StorageType newstor = 
                iscm() ? ColMajor : 
                isrm() ? (j2==j1+1 ? RowMajor : NoMajor) : NoMajor;
            return const_view_type(
                cptr()+j1*stepj(),colsize(),2,
                stepi(),(j2-j1)*stepj(),newstor,ct());
        }

        inline const_view_type colPair(int j1, int j2) const
        {
            TMVAssert(j1>=0 && j1<rowsize() && j2>=0 && j2<rowsize());
            return cColPair(j1,j2);
        }

        inline const_view_type cRowPair(int i1, int i2) const
        {
            StorageType newstor = 
                isrm() ? RowMajor : 
                iscm() ? (i2==i1+1 ? ColMajor : NoMajor) : NoMajor;
            return const_view_type(
                cptr()+i1*stepi(),2,rowsize(),
                (i2-i1)*stepi(),stepj(),newstor,ct());
        }

        inline const_view_type rowPair(int i1, int i2) const
        {
            TMVAssert(i1>=0 && i1<colsize() && i2>=0 && i2<colsize());
            return cRowPair(i1,i2);
        }

        inline const_view_type cColRange(int j1, int j2) const
        {
            return const_view_type(
                cptr()+j1*stepj(),colsize(),j2-j1,
                stepi(),stepj(),stor(),ct(),(iscm()&&ls())?1:0);
        }

        inline const_view_type colRange(int j1, int j2) const
        {
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            return cColRange(j1,j2);
        }

        inline const_view_type cRowRange(int i1, int i2) const
        {
            return const_view_type(
                cptr()+i1*stepi(),i2-i1,rowsize(),
                stepi(),stepj(),stor(),ct(),(isrm()&&ls())?1:0);
        }

        inline const_view_type rowRange(int i1, int i2) const
        {
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            return cRowRange(i1,i2);
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr()),colsize(),rowsize(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? stor() : NoMajor,NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr())+1,
                colsize(),rowsize(),2*stepi(),2*stepj(),NoMajor,NonConj);
        }

        //
        // Views
        //

        inline const_view_type view() const
        { 
            return const_view_type(
                cptr(),colsize(),rowsize(),stepi(),stepj(),stor(),ct(),ls());
        }

        inline const_view_type transpose() const
        { 
            return const_view_type(
                cptr(),rowsize(),colsize(),
                stepj(),stepi(),TMV_TransOf(stor()),ct(),ls());
        }

        inline const_view_type conjugate() const
        { 
            return const_view_type(
                cptr(),colsize(),rowsize(),
                stepi(),stepj(),stor(),TMV_ConjOf(T,ct()),ls());
        }

        inline const_view_type adjoint() const
        { 
            return const_view_type(
                cptr(),rowsize(),colsize(),
                stepj(),stepi(),TMV_TransOf(stor()),TMV_ConjOf(T,ct()),ls());
        }

        inline const_vec_type constLinearView() const
        {
            TMVAssert(ls() != 1 || (rowsize() == 1 && colsize() == 1));
            // To assure that next Assert has no effect

            TMVAssert(canLinearize());
            TMVAssert(ls() == colsize()*rowsize());
            return const_vec_type(cptr(),ls(),1,ct());
        }

        inline nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<T*>(cptr()),colsize(),rowsize(),
                stepi(),stepj(),stor(),ct(),ls()
                TMV_FIRSTLAST1(cptr(),row(colsize()-1).end().getP()));
        }

        //
        // Functions of Matrix
        //

        T det() const;

        RT logDet(T* sign=0) const;

        bool isSingular() const;

        inline T trace() const
        { return diag().sumElements(); }

        T sumElements() const;

        RT sumAbsElements() const;

        RT sumAbs2Elements() const;

        inline RT norm() const 
        { return normF(); }

        RT normF() const;

        // normF()^2
        RT normSq(const RT scale = RT(1)) const;

        // 1-norm = max_j (sum_i |a_ij|)
        RT norm1() const;

        // inf-norm = max_i (sum_j |a_ij|)
        inline RT normInf() const
        { return transpose().norm1(); }

        // = max_i,j (|a_ij|)
        RT maxAbsElement() const;

        // = max_i,j (|real(a_ij)|+|imag(a_ij)|)
        RT maxAbs2Element() const;

        RT doNorm2() const;
        inline RT norm2() const
        {
            if (this->divIsSet() && this->getDivType() == SV)
                return DivHelper<T>::norm2();
            else return doNorm2();
        }

        RT doCondition() const;
        inline RT condition() const
        {
            if (this->divIsSet() && this->getDivType() == SV)
                return DivHelper<T>::condition();
            else return doCondition();
        }

        // icpc gives a strange compiler error if I don't do this 
        // throwugh QInverse.  I think I should be able to just write:
        // QuotXM<T,T> inverse() const;
        // and then define this in TMV_Matrix.cpp, but that doesn't work.
        QuotXM<T,T> QInverse() const;
        inline QuotXM<T,T> inverse() const
        { return QInverse(); }

        // 
        // Division Control
        //

        void setDiv() const;

        inline void divideUsing(DivType dt) const
        {
            TMVAssert(dt == LU || dt == QR || dt == QRP || dt == SV);
            DivHelper<T>::divideUsing(dt); 
        }

        inline const LUDiv<T>& lud() const
        {
            divideUsing(LU);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(divIsLUDiv());
            return static_cast<const LUDiv<T>&>(*this->getDiv());
        }

        inline const QRDiv<T>& qrd() const
        {
            divideUsing(QR);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(divIsQRDiv());
            return static_cast<const QRDiv<T>&>(*this->getDiv());
        }

        inline const QRPDiv<T>& qrpd() const
        {
            divideUsing(QRP);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(divIsQRPDiv());
            return static_cast<const QRPDiv<T>&>(*this->getDiv());
        }

        inline const SVDiv<T>& svd() const
        {
            divideUsing(SV);
            setDiv();
            TMVAssert(this->getDiv());
            TMVAssert(divIsSVDiv());
            return static_cast<const SVDiv<T>&>(*this->getDiv());
        }


        //
        // I/O
        //

        void write(const TMV_Writer& writer) const;

        virtual const T* cptr() const = 0;
        virtual int stepi() const = 0;
        virtual int stepj() const = 0;
        virtual int ls() const = 0;
        virtual inline bool isrm() const { return stor() == RowMajor; }
        virtual inline bool iscm() const { return stor() == ColMajor; }
        virtual inline bool isconj() const 
        { return isComplex(T()) && ct()==Conj; }
        virtual StorageType stor() const =0;
        virtual ConjType ct() const =0;

        virtual bool canLinearize() const = 0;
        virtual T cref(int i, int j) const;

        inline int rowstart(int ) const { return 0; }
        inline int rowend(int ) const { return rowsize(); }
        inline int colstart(int ) const { return 0; }
        inline int colend(int ) const { return colsize(); }

        inline const_rowmajor_iterator rowmajor_begin() const 
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const 
        { return const_rowmajor_iterator(this,colsize(),0); }

        inline const_colmajor_iterator colmajor_begin() const 
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const 
        { return const_colmajor_iterator(this,0,rowsize()); }

    protected :

        inline const BaseMatrix<T>& getMatrix() const { return *this; }

    private :

        type& operator=(const type&);

        bool divIsLUDiv() const;
        bool divIsQRDiv() const;
        bool divIsQRPDiv() const;
        bool divIsSVDiv() const;

    }; // GenMatrix

    template <class T, IndexStyle I> 
    class ConstMatrixView : 
        public GenMatrix<T>
    {
    public :

        typedef GenMatrix<T> base;
        typedef ConstMatrixView<T,I> type;

        inline ConstMatrixView(const type& rhs) :
            itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs), 
            itssi(rhs.itssi), itssj(rhs.itssj),
            itsstor(rhs.itsstor), itsct(rhs.itsct), linsize(rhs.linsize) {}

        inline ConstMatrixView(const base& rhs) :
            itsm(rhs.cptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()), 
            itssi(rhs.stepi()), itssj(rhs.stepj()),
            itsstor(rhs.stor()), itsct(rhs.ct()), linsize(rhs.ls()) {}

        inline ConstMatrixView(
            const T* _m, int _cs, int _rs, int _si, int _sj, 
            StorageType _stor, ConjType _ct, int _ls=0) : 
            itsm(_m), itscs(_cs), itsrs(_rs), itssi(_si), itssj(_sj),
            itsstor(_stor), itsct(_ct), linsize(_ls)
        { 
            TMVAssert(_stor==RowMajor ? _sj==1 :
                      _stor==ColMajor ? _si==1 : true);
        }

        virtual inline ~ConstMatrixView()
        {
#ifdef TMV_EXTRA_DEBUG
            const_cast<const T*&>(itsm) = 0;
#endif
        }

        virtual inline int colsize() const { return itscs; }
        virtual inline int rowsize() const { return itsrs; }
        virtual inline const T* cptr() const { return itsm; }
        virtual inline int stepi() const { return itssi; }
        virtual inline int stepj() const { return itssj; }
        virtual inline StorageType stor() const { return itsstor; }
        virtual inline ConjType ct() const { return itsct; }
        virtual inline int ls() const { return linsize; }
        using GenMatrix<T>::isrm;
        using GenMatrix<T>::iscm;

        virtual inline bool canLinearize() const 
        { 
            if (linsize == 1 && !(rowsize() == 1 && colsize() == 1))
                linsize = rowsize() * colsize();
            TMVAssert(linsize == 0 || isrm() || iscm());
            TMVAssert(linsize == 0 || !isrm() || stepi() == rowsize());
            TMVAssert(linsize == 0 || !iscm() || stepj() == colsize());
            return linsize > 0; 
        }

    protected :

        const T*const itsm;
        const int itscs;
        const int itsrs;
        const int itssi;
        const int itssj;
        const StorageType itsstor;
        const ConjType itsct;

        mutable int linsize;

    private :

        type& operator=(const type&);

    }; // ConstMatrixView

    template <class T> 
    class ConstMatrixView<T,FortranStyle> : 
        public ConstMatrixView<T,CStyle>
    {
    public :

        typedef TMV_RealType(T) RT;
        typedef GenMatrix<T> base;
        typedef ConstMatrixView<T,FortranStyle> type;
        typedef ConstMatrixView<T,CStyle> c_type;
        typedef ConstMatrixView<T,FortranStyle> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstVectorView<T,FortranStyle> const_vec_type;
        typedef ConstMatrixView<RT,FortranStyle> const_realpart_type;
        typedef ConstUpperTriMatrixView<T,FortranStyle> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,FortranStyle> const_lowertri_type;
        typedef MatrixView<T,FortranStyle> nonconst_type;

        inline ConstMatrixView(const ConstMatrixView<T>& rhs) : c_type(rhs) {}

        inline ConstMatrixView(const GenMatrix<T>& rhs) : c_type(rhs) {}

        inline ConstMatrixView(
            const T* _m, int _cs, int _rs, int _si, int _sj, 
            StorageType instor, ConjType inct, int linsize=0) : 
            c_type(_m,_cs,_rs,_si,_sj,instor,inct,linsize) {}

        virtual inline ~ConstMatrixView() {}


        // 
        // Access
        //

        inline T operator()(int i, int j)
        { 
            TMVAssert(i>0 && i<=this->colsize());
            TMVAssert(j>0 && j<=this->rowsize());
            return base::cref(i-1,j-1);
        }

        inline const_vec_type row(int i) const 
        { 
            TMVAssert(i>0 && i<=this->colsize());
            return base::row(i-1);
        }

        inline const_vec_type col(int j) const
        {
            TMVAssert(j>0 && j<=this->rowsize());
            return base::col(j-1);
        }

        inline const_vec_type diag() const
        { return base::diag(); }

        inline const_vec_type diag(int i) const
        { return base::diag(i); }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>0 && i<=this->colsize());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->rowsize());
            return base::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=this->rowsize());
            TMVAssert(i1 > 0 && i1 <= i2 && i2 <= this->colsize());
            return base::col(j-1,i1-1,i2);
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(j1 > 0);
            return base::diag(i,j1-1,j2);
        }

        inline const_vec_type operator[](int i) const
        { return row(i); }

        //
        // subMatrix
        //

        bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const;

        inline const_view_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return base::cSubMatrix(i1-1,i2,j1-1,j2);
        }

        inline const_view_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return base::cSubMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        bool hasSubVector(int i, int j, int istep, int jstep, int s) const;

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int s) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,s));
            return base::cSubVector(i-1,j-1,istep,jstep,s);
        }

        inline const_uppertri_type unitUpperTri() const
        { return base::upperTri(UnitDiag); }

        inline const_uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        { return base::upperTri(dt); }

        inline const_lowertri_type unitLowerTri() const
        { return base::lowerTri(UnitDiag); }

        inline const_lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        { return base::lowerTri(dt); }

        inline const_view_type colPair(int j1, int j2) const
        {
            TMVAssert(j1 > 0 && j1 <= this->rowsize());
            TMVAssert(j2 > 0 && j2 <= this->rowsize());
            return base::cColPair(j1-1,j2-1);
        }

        inline const_view_type rowPair(int i1, int i2) const
        {
            TMVAssert(i1 > 0 && i1 <= this->colsize());
            TMVAssert(i2 > 0 && i2 <= this->colsize());
            return base::cRowPair(i1-1,i2-1);
        }

        inline const_view_type colRange(int j1, int j2) const
        {
            TMVAssert(j1 > 0 && j1 <= j2 && j2 <= this->rowsize());
            return base::cColRange(j1-1,j2);
        }

        inline const_view_type rowRange(int i1, int i2) const
        {
            TMVAssert(i1 > 0 && i1 <= i2 && i2 <= this->colsize());
            return base::cRowRange(i1-1,i2);
        }

        inline const_realpart_type realPart() const
        { return base::realPart(); }

        inline const_realpart_type imagPart() const
        { return base::imagPart(); }

        //
        // Views
        //

        inline const_view_type view() const
        { return base::view(); }

        inline const_view_type transpose() const
        { return base::transpose(); }

        inline const_view_type conjugate() const
        { return base::conjugate(); }

        inline const_view_type adjoint() const
        { return base::adjoint(); }

        inline const_vec_type constLinearView() const
        { return base::constLinearView(); }

        inline nonconst_type nonConst() const
        { return base::nonConst(); }

    private :

        type& operator=(const type&);

    }; // FortranStyle ConstMatrixView

    template <class T, IndexStyle I> 
    class MatrixView : 
        public GenMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenMatrix<T> base;
        typedef MatrixView<T,I> type;
        typedef MatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef MatrixView<RT,I> realpart_type;
        typedef VectorView<T,I> vec_type;
        typedef UpperTriMatrixView<T,I> uppertri_type;
        typedef LowerTriMatrixView<T,I> lowertri_type;
        typedef TMV_RefType(T) reference;
        typedef RMIt<const type> rowmajor_iterator;
        typedef CMIt<const type> colmajor_iterator;

        //
        // Constructors
        //

        inline MatrixView(const type& rhs) : 
            itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs),
            itssi(rhs.itssi), itssj(rhs.itssj),
            itsstor(rhs.itsstor), itsct(rhs.itsct),
            linsize(rhs.linsize) 
            TMV_DEFFIRSTLAST(rhs._first,rhs._last) {}

        inline MatrixView(
            T* _m, int _cs, int _rs, int _si, int _sj,
            StorageType _stor, ConjType _ct, int _ls 
            TMV_PARAMFIRSTLAST(T) ) :
            itsm(_m), itscs(_cs), itsrs(_rs), itssi(_si), itssj(_sj),
            itsstor(_stor), itsct(_ct), linsize(_ls) 
            TMV_DEFFIRSTLAST(_first,_last)
        {
            TMVAssert(_stor==RowMajor ? _sj==1 :
                      _stor==ColMajor ? _si==1 : true);
        }

        inline MatrixView(
            T* _m, int _cs, int _rs, int _si, int _sj,
            StorageType _stor, ConjType _ct TMV_PARAMFIRSTLAST(T) ) :
            itsm(_m), itscs(_cs), itsrs(_rs), itssi(_si), itssj(_sj),
            itsstor(_stor), itsct(_ct), linsize(0) 
            TMV_DEFFIRSTLAST(_first,_last)
        {
            TMVAssert(_stor==RowMajor ? _sj==1 :
                      _stor==ColMajor ? _si==1 : true);
        }

        virtual inline ~MatrixView() 
        {
            TMV_SETFIRSTLAST(0,0);
#ifdef TMV_EXTRA_DEBUG
            const_cast<T*&>(itsm) = 0;
#endif
        }

        //
        // Op=
        //

        inline const type& operator=(const MatrixView<T,I>& m2) const
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            m2.assignToM(*this); 
            return *this; 
        }

        // This next one is to make sure the compiler doesn't try to make
        // an implicit non-const assignment operator.  I think the above 
        // const version should be sufficient, but some compilers need the 
        // non-const one as well.  (e.g. icc version 10)
        inline const type& operator=(const MatrixView<T,I>& m2)
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            m2.assignToM(*this); 
            return *this; 
        }

        inline const type& operator=(const GenMatrix<RT>& m2) const
        { 
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            m2.assignToM(*this); 
            return *this; 
        }

        inline const type& operator=(const GenMatrix<CT>& m2) const
        { 
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVAssert(isComplex(T()));
            m2.assignToM(*this); 
            return *this; 
        }

        template <class T2> 
        inline const type& operator=(const GenMatrix<T2>& m2) const
        { 
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(m2,*this);
            return *this; 
        }

        inline const type& operator=(const T& x) const 
        { return setToIdentity(x); }

        inline const type& operator=(
            const AssignableToMatrix<RT>& m2) const
        { 
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignToM(*this);
            return *this;
        }

        inline const type& operator=(
            const AssignableToMatrix<CT>& m2) const
        { 
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(isComplex(T()));
            m2.assignToM(*this);
            return *this;
        }

        inline const type& operator=(const Permutation& m2) const
        {
            m2.assignToM(*this);
            return *this; 
        }

        template <class T2, int M, int N, StorageType S2, IndexStyle I2> 
        inline const type& operator=(const SmallMatrix<T2,M,N,S2,I2>& m2) const
        { 
            TMVAssert(colsize() == M && rowsize() == N);
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(m2.view(),*this);
            return *this; 
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x) const
        { return MyListAssigner(rowmajor_begin(),colsize()*rowsize(),x); }

        //
        // Access
        //

        inline reference operator()(int i,int j) const 
        { 
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j>=0 && j<rowsize());
            return ref(i,j); 
        }

        inline vec_type operator[](int i) const 
        { 
            TMVAssert(i>=0 && i<colsize());
            return row(i); 
        }

        inline vec_type row(int i) const
        {
            TMVAssert(i>=0 && i<colsize());
            return vec_type(
                ptr()+i*stepi(), rowsize(),stepj(),ct() TMV_FIRSTLAST ); 
        }

        inline vec_type col(int j) const
        {
            TMVAssert(j>=0 && j<rowsize());
            return vec_type(
                ptr()+j*stepj(), colsize(),stepi(),ct() TMV_FIRSTLAST ); 
        }

        inline vec_type diag() const
        {
            return vec_type(
                ptr(), TMV_MIN(colsize(),rowsize()),stepi()+stepj(),ct() 
                TMV_FIRSTLAST);
        }

        inline vec_type diag(int i) const
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const int diagstep = stepi() + stepj();
            if (i >= 0) {
                const int diagsize = TMV_MIN(colsize(),rowsize()-i);
                return vec_type(
                    ptr()+i*stepj(), diagsize,diagstep,ct() TMV_FIRSTLAST );
            } else {
                const int diagsize = TMV_MIN(colsize()+i,rowsize());
                return vec_type(
                    ptr()-i*stepi(), diagsize,diagstep,ct() TMV_FIRSTLAST );
            }
        }

        inline vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            return vec_type(
                ptr()+i*stepi()+j1*stepj(), j2-j1,stepj(),ct() TMV_FIRSTLAST ); 
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>=0 && j<rowsize());
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            return vec_type(
                ptr()+i1*stepi()+j*stepj(), i2-i1,stepi(),ct() TMV_FIRSTLAST ); 
        }

        inline vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const int diagstep = stepi() + stepj();
            if (i >= 0) {
                TMVAssert(
                    j1>=0 && j1-j2<=0 && 
                    j2<=TMV_MIN(colsize(),rowsize()-i));
                return vec_type(
                    ptr()+i*stepj()+j1*diagstep, j2-j1,diagstep,ct() 
                    TMV_FIRSTLAST );
            } else {
                TMVAssert(
                    j1>=0 && j1-j2<=0 && 
                    j2<=TMV_MIN(colsize()+i,rowsize()));
                return vec_type(
                    ptr()-i*stepi()+j1*diagstep, j2-j1,diagstep,ct() 
                    TMV_FIRSTLAST );
            }
        }

        //
        // Modifying Functions
        //

        const type& setZero() const;

        const type& setAllTo(const T& x) const;

        const type& addToAll(const T& x) const;

        const type& clip(RT thresh) const;

        const type& transposeSelf() const;

        const type& conjugateSelf() const;

        const type& setToIdentity(const T& x=T(1)) const;

        inline const type& swapRows(int i1, int i2) const
        {
            TMVAssert(i1>=0 && i1 < colsize() && 
                      i2>=0 && i2 < colsize());
            if (i1!=i2) Swap(row(i1),row(i2));
            return *this;
        }

        inline const type& swapCols(int j1, int j2) const
        {
            TMVAssert(j1>=0 && j1 < rowsize() && 
                      j2>=0 && j2 < rowsize());
            if (j1!=j2) Swap(col(j1),col(j2));
            return *this;
        }

        const type& permuteRows(const int* p, int i1, int i2) const;

        inline const type& permuteRows(const int* p) const
        { return permuteRows(p,0,colsize()); }

        inline const type& permuteCols(
            const int* p, int j1, int j2) const
        { transpose().permuteRows(p,j1,j2); return *this; }

        inline const type& permuteCols(const int* p) const
        { return permuteCols(p,0,rowsize()); }

        const type& reversePermuteRows(
            const int* p, int i1, int i2) const;

        inline const type& reversePermuteRows(const int* p) const
        { return reversePermuteRows(p,0,colsize()); }

        inline const type& reversePermuteCols(
            const int* p, int j1, int j2) const
        { transpose().reversePermuteRows(p,j1,j2); return *this; }

        inline const type& reversePermuteCols(const int* p) const
        { return reversePermuteCols(p,0,rowsize()); }

        //
        // subMatrix
        //

        inline view_type cSubMatrix(int i1, int i2, int j1, int j2) const
        {
            return type(
                ptr()+i1*stepi()+j1*stepj(),
                i2-i1,j2-j1,stepi(),stepj(),stor(),ct() TMV_FIRSTLAST );
        }

        inline view_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,1,1));
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline view_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            StorageType newstor = 
                this->iscm() ? (istep == 1 ? ColMajor : NoMajor) :
                this->isrm() ? (jstep == 1 ? RowMajor : NoMajor) : NoMajor;
            return type(
                ptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),
                newstor,ct() 
                TMV_FIRSTLAST );
        }

        inline view_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline vec_type cSubVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            return vec_type(
                ptr()+i*stepi()+j*stepj(),size,
                istep*stepi()+jstep*stepj(),ct() 
                TMV_FIRSTLAST );
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int size) const
        {
            TMVAssert(size >= 0);
            TMVAssert(base::hasSubVector(i,j,istep,jstep,size));
            return cSubVector(i,j,istep,jstep,size);
        }

        inline uppertri_type unitUpperTri() const
        {
            TMVAssert(rowsize() <= colsize());
            return uppertri_type(
                ptr(),rowsize(),stepi(),stepj(),UnitDiag,stor(),ct() 
                TMV_FIRSTLAST);
        }

        inline uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        {
            TMVAssert(rowsize() <= colsize());
            return uppertri_type(
                ptr(),rowsize(),stepi(),stepj(), dt,stor(),ct() TMV_FIRSTLAST);
        }

        inline lowertri_type unitLowerTri() const
        {
            TMVAssert(colsize() <= rowsize());
            return lowertri_type(
                ptr(),colsize(),stepi(),stepj(),UnitDiag,stor(),ct() 
                TMV_FIRSTLAST);
        }

        inline lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        {
            TMVAssert(colsize() <= rowsize());
            return lowertri_type(
                ptr(),colsize(),stepi(),stepj(), dt,stor(),ct() TMV_FIRSTLAST);
        }

        inline view_type cColPair(int j1, int j2) const
        {
            StorageType newstor = 
                this->iscm() ? ColMajor : 
                this->isrm() ? (j2==j1+1 ? RowMajor : NoMajor) : NoMajor;
            return type(
                ptr()+j1*stepj(),colsize(),2,
                stepi(),(j2-j1)*stepj(),newstor,ct() 
                TMV_FIRSTLAST );
        }

        inline view_type colPair(int j1, int j2) const
        {
            TMVAssert(j1>=0 && j1<rowsize() && j2>=0 && j2<rowsize());
            return cColPair(j1,j2);
        }

        inline view_type cRowPair(int i1, int i2) const
        {
            StorageType newstor = 
                this->isrm() ? RowMajor : 
                this->iscm() ? (i2==i1+1 ? ColMajor : NoMajor) : NoMajor;
            return type(
                ptr()+i1*stepi(),2,rowsize(),
                (i2-i1)*stepi(),stepj(),newstor,ct() 
                TMV_FIRSTLAST );
        }

        inline view_type rowPair(int i1, int i2) const
        {
            TMVAssert(i1>=0 && i1<colsize() && i2>=0 && i2<colsize());
            return cRowPair(i1,i2);
        }

        inline view_type cColRange(int j1, int j2) const
        {
            return type(
                ptr()+j1*stepj(),colsize(),j2-j1,
                stepi(),stepj(),stor(),ct(),
                (this->iscm()&&ls())?1:0 
                TMV_FIRSTLAST);
        }

        inline view_type colRange(int j1, int j2) const
        {
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            return cColRange(j1,j2);
        }

        inline view_type cRowRange(int i1, int i2) const
        {
            return type(
                ptr()+i1*stepi(),i2-i1,rowsize(),
                stepi(),stepj(),stor(),ct(),
                (this->isrm()&&ls())?1:0 
                TMV_FIRSTLAST);
        }

        inline view_type rowRange(int i1, int i2) const
        {
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            return cRowRange(i1,i2);
        }

        inline realpart_type realPart() const
        {
            return realpart_type(
                reinterpret_cast<RT*>(ptr()),colsize(),rowsize(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? stor() : NoMajor,NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return realpart_type(
                reinterpret_cast<RT*>(ptr())+1,
                colsize(),rowsize(),2*stepi(),2*stepj(), NoMajor,NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        //
        // Views
        //

        inline view_type view() const
        { return *this; }

        inline view_type transpose() const
        { 
            return type(
                ptr(),rowsize(),colsize(),
                stepj(),stepi(),TMV_TransOf(stor()),ct(),ls() TMV_FIRSTLAST);
        }

        inline view_type conjugate() const
        { 
            return type(
                ptr(),colsize(),rowsize(),
                stepi(),stepj(),stor(),TMV_ConjOf(T,ct()),ls() TMV_FIRSTLAST);
        }

        inline view_type adjoint() const
        { 
            return type(
                ptr(),rowsize(),colsize(),
                stepj(),stepi(),TMV_TransOf(stor()),TMV_ConjOf(T,ct()),ls() 
                TMV_FIRSTLAST);
        }

        inline vec_type linearView() const
        {
            TMVAssert(ls() != 1 || (rowsize() == 1 && colsize() == 1));
            // To assure that next Assert has no effect

            TMVAssert(canLinearize());
            TMVAssert(ls() == colsize()*rowsize());
            return vec_type(ptr(),ls(),1,ct() TMV_FIRSTLAST );
        }

        //
        // I/O
        //

        void read(const TMV_Reader& reader) const;

        virtual inline int colsize() const { return itscs; }
        virtual inline int rowsize() const { return itsrs; }
        virtual inline const T* cptr() const { return itsm; }
        inline T* ptr() const { return itsm; }
        virtual inline int stepi() const { return itssi; }
        virtual inline int stepj() const { return itssj; }
        virtual inline StorageType stor() const { return itsstor; }
        virtual inline ConjType ct() const { return itsct; }
        virtual inline int ls() const { return linsize; }

        virtual inline bool canLinearize() const 
        { 
            if (linsize == 1 && !(rowsize() == 1 && colsize() == 1))
                linsize = rowsize() * colsize();
            TMVAssert(linsize == 0 || this->isrm() || this->iscm());
            TMVAssert(linsize == 0 || !this->isrm() || 
                      stepi() == rowsize());
            TMVAssert(linsize == 0 || !this->iscm() || 
                      stepj() == colsize());
            return linsize > 0; 
        }

        reference ref(int i, int j) const;

        inline rowmajor_iterator rowmajor_begin() const
        { return rowmajor_iterator(this,0,0); }
        inline rowmajor_iterator rowmajor_end() const
        { return rowmajor_iterator(this,colsize(),0); }

        inline colmajor_iterator colmajor_begin() const
        { return colmajor_iterator(this,0,0); }
        inline colmajor_iterator colmajor_end() const
        { return colmajor_iterator(this,0,rowsize()); }

    protected:

        T*const itsm;
        const int itscs;
        const int itsrs;
        const int itssi;
        const int itssj;
        const StorageType itsstor;
        const ConjType itsct;

        mutable int linsize;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
#endif

    }; // MatrixView

    template <class T> 
    class MatrixView<T,FortranStyle> : 
        public MatrixView<T,CStyle>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef MatrixView<T,FortranStyle> type;
        typedef ConstMatrixView<T,FortranStyle> const_type;
        typedef MatrixView<T,CStyle> c_type;
        typedef MatrixView<T,FortranStyle> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef MatrixView<RT,FortranStyle> realpart_type;
        typedef VectorView<T,FortranStyle> vec_type;
        typedef UpperTriMatrixView<T,FortranStyle> uppertri_type;
        typedef LowerTriMatrixView<T,FortranStyle> lowertri_type;

        //
        // Constructors
        //

        inline MatrixView(const type& rhs) : c_type(rhs) {}

        inline MatrixView(const c_type& rhs) : c_type(rhs) {}

        inline MatrixView(
            T* _m, int _cs, int _rs, int _si, int _sj,
            StorageType instor, ConjType inct, int linsize
            TMV_PARAMFIRSTLAST(T) ) :
            c_type(_m,_cs,_rs,_si,_sj,instor,inct,linsize 
                   TMV_FIRSTLAST1(_first,_last) ) 
        {}

        inline MatrixView(
            T* _m, int _cs, int _rs, int _si, int _sj,
            StorageType instor, ConjType inct TMV_PARAMFIRSTLAST(T) ) :
            c_type(_m,_cs,_rs,_si,_sj,instor,inct TMV_FIRSTLAST1(_first,_last))
        {}

        virtual inline ~MatrixView() {} 

        //
        // Op=
        //

        inline const type& operator=(const type& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const type& m2)
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const c_type& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const c_type& m2)
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        template <class T2> 
        inline const type& operator=(const GenMatrix<T2>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const T& x) const 
        { c_type::operator=(x); return *this; }

        inline const type& operator=(const AssignableToMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const AssignableToMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        template <class T2, int M, int N, StorageType S2, IndexStyle I2> 
        inline const type& operator=(const SmallMatrix<T2,M,N,S2,I2>& m2) const
        { c_type::operator=(m2); return *this; }

        typedef typename c_type::MyListAssigner MyListAssigner;
        inline MyListAssigner operator<<(const T& x) const
        { return c_type::operator<<(x); }

        //
        // Access
        //

        inline TMV_RefType(T) operator()(int i,int j) const 
        { 
            TMVAssert(i > 0 && i <= this->colsize());
            TMVAssert(j > 0 && j <= this->rowsize());
            return ref(i-1,j-1); 
        }

        inline vec_type operator[](int i) const 
        { 
            TMVAssert(i>0 && i<=this->colsize());
            return row(i); 
        }

        inline vec_type row(int i) const
        {
            TMVAssert(i>0 && i<=this->colsize());
            return c_type::row(i-1);
        }

        inline vec_type col(int j) const
        {
            TMVAssert(j>0 && j<=this->rowsize());
            return c_type::col(j-1);
        }

        inline vec_type diag() const
        { return c_type::diag(); }

        inline vec_type diag(int i) const
        { return c_type::diag(i); }

        inline vec_type row(int i, int j1, int j2) const 
        { 
            TMVAssert(i>0 && i<=this->colsize());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->rowsize());
            return c_type::row(i-1,j1-1,j2);
        }

        inline vec_type col(int j, int i1, int i2) const
        {
            TMVAssert(j>0 && j<=this->rowsize());
            TMVAssert(i1>0 && i1-i2<=0 && i2<=this->colsize());
            return c_type::col(j-1,i1-1,i2);
        }

        inline vec_type diag(int i, int j1, int j2) const
        {
            TMVAssert(j1 > 0);
            return c_type::diag(i,j1-1,j2);
        }

        //
        // Modifying Functions
        //

        inline const type& setZero() const 
        { c_type::setZero(); return *this; }

        inline const type& setAllTo(const T& x) const
        { c_type::setAllTo(x); return *this; }

        inline const type& addToAll(const T& x) const
        { c_type::addToAll(x); return *this; }

        inline const type& clip(RT thresh) const
        { c_type::clip(thresh); return *this; }

        inline const type& transposeSelf() const
        { c_type::transposeSelf(); return *this; }

        inline const type& conjugateSelf() const
        { c_type::conjugateSelf(); return *this; }

        inline const type& setToIdentity(const T& x=T(1)) const
        { c_type::setToIdentity(x); return *this; }

        inline const type& swapRows(int i1, int i2) const
        { 
            TMVAssert(i1 > 0 && i1 <= this->colsize());
            TMVAssert(i2 > 0 && i2 <= this->colsize());
            if (i1 != i2)
                c_type::swapRows(i1-1,i2-1); 
            return *this; 
        }

        inline const type& swapCols(int j1, int j2) const
        { 
            TMVAssert(j1 > 0 && j1 <= this->rowsize());
            TMVAssert(j2 > 0 && j2 <= this->rowsize());
            if (j1 != j2)
                c_type::swapCols(j1-1,j2-1); 
            return *this; 
        }

        inline const type& permuteRows(const int* p, int i1, int i2) const
        {
            TMVAssert(i1>0);
            c_type::permuteRows(p,i1-1,i2);
            return *this; 
        }

        inline const type& permuteRows(const int* p) const
        { c_type::permuteRows(p); return *this; }

        inline const type& permuteCols(const int* p, int j1, int j2) const
        { transpose().permuteRows(p,j1,j2); return *this; }

        inline const type& permuteCols(const int* p) const
        { transpose().permuteRows(p); return *this; }

        inline const type& reversePermuteRows(
            const int* p, int i1, int i2) const
        { 
            TMVAssert(i1>0);
            c_type::reversePermuteRows(p,i1-1,i2);  
            return *this;
        }

        inline const type& reversePermuteRows(const int* p) const
        { c_type::reversePermuteRows(p); return *this; }

        inline const type& reversePermuteCols(
            const int* p, int j1, int j2) const
        { transpose().reversePermuteRows(p,j1,j2); return *this; }

        inline const type& reversePermuteCols(const int* p) const
        { transpose().reversePermuteRows(p); return *this; }

        //
        // subMatrix
        //

        inline bool hasSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            return const_type(*this).hasSubMatrix(
                i1,i2,j1,j2,istep,jstep); 
        }

        inline bool hasSubVector(
            int i, int j, int istep, int jstep, int s) const
        { return const_type(*this).hasSubVector(i,j,istep,jstep,s); }

        inline view_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return c_type::cSubMatrix(i1-1,i2,j1-1,j2);
        }

        inline view_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return c_type::cSubMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline vec_type subVector(
            int i, int j, int istep, int jstep, int s) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,s));
            return c_type::cSubVector(i-1,j-1,istep,jstep,s);
        }

        inline uppertri_type unitUpperTri() const
        { return c_type::upperTri(UnitDiag); }

        inline uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        { return c_type::upperTri(dt); }

        inline lowertri_type unitLowerTri() const
        { return c_type::lowerTri(UnitDiag); }

        inline lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        { return c_type::lowerTri(dt); }

        inline view_type colPair(int j1, int j2) const
        {
            TMVAssert(j1 > 0 && j1 <= this->rowsize());
            TMVAssert(j2 > 0 && j2 <= this->rowsize());
            return c_type::cColPair(j1-1,j2-1);
        }

        inline view_type rowPair(int i1, int i2) const
        {
            TMVAssert(i1 > 0 && i1 <= this->rowsize());
            TMVAssert(i2 > 0 && i2 <= this->rowsize());
            return c_type::cRowPair(i1-1,i2-1);
        }

        inline view_type colRange(int j1, int j2) const
        {
            TMVAssert(j1 > 0 && j1 <= j2 && j2 <= this->rowsize());
            return c_type::cColRange(j1-1,j2);
        }

        inline view_type rowRange(int i1, int i2) const
        {
            TMVAssert(i1 > 0 && i1 <= i2 && i2 <= this->colsize());
            return c_type::cRowRange(i1-1,i2);
        }

        inline realpart_type realPart() const
        { return c_type::realPart(); }

        inline realpart_type imagPart() const
        { return c_type::imagPart(); }

        //
        // Views
        //

        inline view_type view() const
        { return c_type::view(); }

        inline view_type transpose() const
        { return c_type::transpose(); }

        inline view_type conjugate() const
        { return c_type::conjugate(); }

        inline view_type adjoint() const
        { return c_type::adjoint(); }

        inline vec_type linearView() const
        { return c_type::linearView(); }

    protected:

        using c_type::ref;

    }; // FortranStyle MatrixView

    template <StorageType S, class M>
    struct MatrixIterHelper;

    template <class M>
    struct MatrixIterHelper<RowMajor,M>
    {
        typedef typename M::value_type T;
        typedef VIt<T,Unit,NonConj> rowmajor_iterator;
        typedef CVIt<T,Unit,NonConj> const_rowmajor_iterator;

        static rowmajor_iterator rowmajor_begin(M* m)
        { return rowmajor_iterator(m->ptr(),1); }
        static rowmajor_iterator rowmajor_end(M* m)
        { return rowmajor_iterator(m->ptr()+m->ls(),1); }

        static const_rowmajor_iterator rowmajor_begin(const M* m)
        { return const_rowmajor_iterator(m->cptr(),1); }
        static const_rowmajor_iterator rowmajor_end(const M* m)
        { return const_rowmajor_iterator(m->cptr()+m->ls(),1); }

        typedef CMIt<M> colmajor_iterator;
        typedef CCMIt<M> const_colmajor_iterator;

        static colmajor_iterator colmajor_begin(M* m)
        { return colmajor_iterator(m,0,0); }
        static colmajor_iterator colmajor_end(M* m)
        { return colmajor_iterator(m,0,m->rowsize()); }

        static const_colmajor_iterator colmajor_begin(const M* m)
        { return const_colmajor_iterator(m,0,0); }
        static const_colmajor_iterator colmajor_end(const M* m)
        { return const_colmajor_iterator(m,0,m->rowsize()); }

    };

    template <class M>
    struct MatrixIterHelper<ColMajor,M>
    {
        typedef RMIt<M> rowmajor_iterator;
        typedef CRMIt<M> const_rowmajor_iterator;

        static rowmajor_iterator rowmajor_begin(M* m)
        { return rowmajor_iterator(m,0,0); }
        static rowmajor_iterator rowmajor_end(M* m)
        { return rowmajor_iterator(m,m->colsize(),0); }

        static const_rowmajor_iterator rowmajor_begin(const M* m)
        { return const_rowmajor_iterator(m,0,0); }
        static const_rowmajor_iterator rowmajor_end(const M* m)
        { return const_rowmajor_iterator(m,m->colsize(),0); }

        typedef typename M::value_type T;
        typedef VIt<T,Unit,NonConj> colmajor_iterator;
        typedef CVIt<T,Unit,NonConj> const_colmajor_iterator;

        static colmajor_iterator colmajor_begin(M* m)
        { return colmajor_iterator(m->ptr(),1); }
        static colmajor_iterator colmajor_end(M* m)
        { return colmajor_iterator(m->ptr()+m->ls(),1); }

        static const_colmajor_iterator colmajor_begin(const M* m)
        { return const_colmajor_iterator(m->cptr(),1); }
        static const_colmajor_iterator colmajor_end(const M* m)
        { return const_colmajor_iterator(m->cptr()+m->ls(),1); }

    };

    template <class T, StorageType S, IndexStyle I> 
    class Matrix : 
        public GenMatrix<T> 
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef Matrix<T,S,I> type;
        typedef ConstMatrixView<T,I> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstMatrixView<RT,I> const_realpart_type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef ConstUpperTriMatrixView<T,I> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,I> const_lowertri_type;
        typedef MatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef MatrixView<RT,I> realpart_type;
        typedef VectorView<T,I> vec_type;
        typedef UpperTriMatrixView<T,I> uppertri_type;
        typedef LowerTriMatrixView<T,I> lowertri_type;
        typedef T& reference;
        typedef typename MatrixIterHelper<S,type>::rowmajor_iterator 
            rowmajor_iterator;
        typedef typename MatrixIterHelper<S,type>::const_rowmajor_iterator 
            const_rowmajor_iterator;
        typedef typename MatrixIterHelper<S,type>::colmajor_iterator 
            colmajor_iterator;
        typedef typename MatrixIterHelper<S,type>::const_colmajor_iterator 
            const_colmajor_iterator;

        //
        // Constructors
        //

#define NEW_SIZE(cs,rs) \
        linsize((cs)*(rs)), \
        itsm(linsize), itscs(cs), itsrs(rs) \
        TMV_DEFFIRSTLAST(itsm.get(),itsm.get()+linsize)

        inline Matrix() : NEW_SIZE(0,0)
        {
            TMVAssert(S==RowMajor || S==ColMajor);
        }

        inline Matrix(int _colsize, int _rowsize) :
            NEW_SIZE(_colsize,_rowsize)
        {
            TMVAssert(_colsize >= 0 && _rowsize >= 0);
            TMVAssert(S==RowMajor || S==ColMajor);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline Matrix(int _colsize, int _rowsize, const T& x) :
            NEW_SIZE(_colsize,_rowsize)
        {
            TMVAssert(_colsize >= 0 && _rowsize >= 0);
            TMVAssert(S==RowMajor || S==ColMajor);
            setAllTo(x);
        }

        inline Matrix(const type& rhs) : NEW_SIZE(rhs.colsize(),rhs.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(rhs.cptr(),rhs.cptr()+linsize,itsm.get());
        }

        template <IndexStyle I2> 
        inline Matrix(const Matrix<T,S,I2>& rhs) :
            NEW_SIZE(rhs.colsize(),rhs.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            std::copy(rhs.cptr(),rhs.cptr()+linsize,itsm.get());
        }

        inline Matrix(const GenMatrix<RT>& rhs) :
            NEW_SIZE(rhs.colsize(),rhs.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            rhs.assignToM(view());
        }

        inline Matrix(const GenMatrix<CT>& rhs) :
            NEW_SIZE(rhs.colsize(),rhs.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(isComplex(T()));
            rhs.assignToM(view());
        }

        template <class T2> 
        inline Matrix(const GenMatrix<T2>& rhs) :
            NEW_SIZE(rhs.colsize(),rhs.rowsize()) 
        { 
            TMVAssert(isComplex(T()) || isReal(T2()));
            TMVAssert(S==RowMajor || S==ColMajor);
            Copy(rhs,view()); 
        }

        inline Matrix(const AssignableToMatrix<RT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            m2.assignToM(view());
        }

        inline Matrix(const AssignableToMatrix<CT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            TMVAssert(isComplex(T()));
            m2.assignToM(view());
        }

        inline Matrix(const Permutation& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize())
        {
            TMVAssert(S==RowMajor || S==ColMajor);
            m2.assignToM(view());
        }

        template <class T2, int M, int N, StorageType S2, IndexStyle I2> 
        inline Matrix(const SmallMatrix<T2,M,N,S2,I2>& rhs) :
            NEW_SIZE(rhs.colsize(),rhs.rowsize()) 
        { 
            TMVAssert(isComplex(T()) || isReal(T2()));
            TMVAssert(S==RowMajor || S==ColMajor);
            Copy(rhs.view(),view()); 
        }


#undef NEW_SIZE

        virtual inline ~Matrix() 
        {
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(999));
#endif
        }

        //
        // Op=
        //

        inline type& operator=(const Matrix<T,S,I>& m2)
        { 
            TMVAssert(m2.colsize() == colsize() && m2.rowsize() == rowsize());
            if (&m2 != this)
                std::copy(m2.cptr(),m2.cptr()+linsize,itsm.get());
            return *this; 
        }

        template <IndexStyle I2> 
        inline type& operator=(const Matrix<T,S,I2>& m2)
        { 
            TMVAssert(m2.colsize() == colsize() && m2.rowsize() == rowsize());
            std::copy(m2.cptr(),m2.cptr()+linsize,itsm.get());
            return *this; 
        }

        inline type& operator=(const GenMatrix<RT>& m2)
        { 
            TMVAssert(m2.colsize() == colsize() && m2.rowsize() == rowsize());
            m2.assignToM(view());
            return *this; 
        }

        inline type& operator=(const GenMatrix<CT>& m2)
        { 
            TMVAssert(m2.colsize() == colsize() && m2.rowsize() == rowsize());
            TMVAssert(isComplex(T()));
            m2.assignToM(view());
            return *this; 
        }

        template <class T2> 
        inline type& operator=(const GenMatrix<T2>& m2)
        {
            TMVAssert(m2.colsize() == colsize() && m2.rowsize() == rowsize());
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(m2,view()); 
            return *this; 
        }

        inline type& operator=(const T& x) 
        { return setToIdentity(x); }

        inline type& operator=(const AssignableToMatrix<RT>& m2)
        { 
            TMVAssert(m2.colsize() == colsize() && m2.rowsize() == rowsize());
            m2.assignToM(view()); 
            return *this; 
        }

        inline type& operator=(const AssignableToMatrix<CT>& m2)
        { 
            TMVAssert(m2.colsize() == colsize() && m2.rowsize() == rowsize());
            TMVAssert(isComplex(T()));
            m2.assignToM(view()); 
            return *this; 
        }

        inline type& operator=(const Permutation& m2)
        {
            TMVAssert(m2.colsize() == colsize() && m2.rowsize() == rowsize());
            m2.assignToM(view());
            return *this; 
        }

        template <class T2, int M, int N, StorageType S2, IndexStyle I2> 
        inline type& operator=(const SmallMatrix<T2,M,N,S2,I2>& m2)
        { 
            TMVAssert(m2.colsize() == colsize() && m2.rowsize() == rowsize());
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(m2.view(),view()); 
            return *this; 
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { 
            TMVAssert(linsize == colsize() * rowsize());
            return MyListAssigner(rowmajor_begin(),linsize,x);
        }

        //
        // Access
        //

        inline T operator()(int i,int j) const
        {
            if (I == CStyle) {
                TMVAssert(i>=0 && i < colsize());
                TMVAssert(j>=0 && j < rowsize());
                return cref(i,j); 
            } else {
                TMVAssert(i > 0 && i <= colsize());
                TMVAssert(j > 0 && j <= rowsize());
                return cref(i-1,j-1); 
            }
        }

        inline T& operator()(int i,int j) 
        {
            if (I == CStyle) {
                TMVAssert(i>=0 && i < colsize());
                TMVAssert(j>=0 && j < rowsize());
                return ref(i,j); 
            } else {
                TMVAssert(i > 0 && i <= colsize());
                TMVAssert(j > 0 && j <= rowsize());
                return ref(i-1,j-1); 
            }
        }

        inline const_vec_type row(int i) const 
        { 
            if (I == FortranStyle) { TMVAssert(i>0 && i<=colsize()); --i; }
            else TMVAssert(i>=0 && i<colsize());
            return const_vec_type(
                itsm.get()+i*stepi(),rowsize(),stepj(),NonConj); 
        }

        inline const_vec_type row(int i, int j1, int j2) const 
        { 
            if (I == FortranStyle) {
                TMVAssert(i>0 && i<=colsize()); --i; 
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize()); --j1; 
            } else {
                TMVAssert(i>=0 && i<colsize()); 
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }
            return const_vec_type(
                itsm.get()+i*stepi()+j1*stepj(),j2-j1,stepj(),NonConj); 
        }

        inline const_vec_type operator[](int i) const
        { return row(i); }

        inline const_vec_type col(int j) const
        {
            if (I == FortranStyle) { TMVAssert(j>0 && j<=rowsize()); --j; }
            else TMVAssert(j>=0 && j<rowsize());
            return const_vec_type(
                itsm.get()+j*stepj(),colsize(),stepi(),NonConj); 
        }

        inline const_vec_type col(int j, int i1, int i2) const
        {
            if (I == FortranStyle) {
                TMVAssert(j>0 && j<=rowsize()); --j; 
                TMVAssert(i1>0 && i1-i2<=0 && i2<=colsize()); --i1; 
            } else {
                TMVAssert(j>=0 && j<rowsize()); 
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            }
            return const_vec_type(
                itsm.get()+i1*stepi()+j*stepj(),i2-i1,stepi(),NonConj); 
        }

        inline const_vec_type diag() const
        {
            return const_vec_type(
                itsm.get(),TMV_MIN(colsize(),rowsize()),
                stepi()+stepj(),NonConj);
        }

        inline const_vec_type diag(int i) const
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const int diagstep = stepi() + stepj();
            if (i >= 0) {
                const int diagsize = TMV_MIN(colsize(),rowsize()-i);
                return const_vec_type(
                    itsm.get()+i*stepj(),diagsize,diagstep,NonConj);
            } else {
                const int diagsize = TMV_MIN(colsize()+i,rowsize());
                return const_vec_type(
                    itsm.get()-i*stepi(),diagsize,diagstep,NonConj);
            }
        }

        inline const_vec_type diag(int i, int j1, int j2) const
        {
            if (I == FortranStyle) {
                TMVAssert(j1 > 0 && j1-j2<=0);
                --j1;
            } else { 
                TMVAssert( j1>=0 && j1-j2<=0);
            }
            TMVAssert(i>=-colsize() && i<=rowsize());
            const int diagstep = stepi() + 1;
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(colsize(),rowsize()-i));
                return const_vec_type(
                    itsm.get()+i*stepj()+j1*diagstep,j2-j1,diagstep,NonConj);
            } else {
                TMVAssert(j2<=TMV_MIN(colsize()+i,rowsize()));
                return const_vec_type(
                    itsm.get()-i*stepi()+j1*diagstep,j2-j1,diagstep,NonConj);
            }
        }

        inline vec_type row(int i)
        { 
            if (I == FortranStyle) { 
                TMVAssert(i>0 && i<=colsize());
                --i;
            } else {
                TMVAssert(i>=0 && i<colsize());
            }
            return vec_type(
                ptr()+i*stepi(),rowsize(),stepj(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type row(int i, int j1, int j2)
        { 
            if (I == FortranStyle) {
                TMVAssert(i>0 && i<=colsize());
                --i; 
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize());
                --j1; 
            } else {
                TMVAssert(i>=0 && i<colsize()); 
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }
            return vec_type(
                ptr()+i*stepi()+j1*stepj(),
                j2-j1,stepj(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type operator[](int i)
        { return row(i); }

        inline vec_type col(int j)
        {
            if (I == FortranStyle) { 
                TMVAssert(j>0 && j<=rowsize());
                --j; 
            } else {
                TMVAssert(j>=0 && j<rowsize());
            }
            return vec_type(
                ptr()+j*stepj(),colsize(),stepi(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type col(int j, int i1, int i2)
        {
            if (I == FortranStyle) {
                TMVAssert(j>0 && j<=rowsize()); --j; 
                TMVAssert(i1>0 && i1-i2<=0 && i2<=colsize()); --i1; 
            } else {
                TMVAssert(j>=0 && j<rowsize()); 
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            }
            return vec_type(
                ptr()+i1*stepi()+j*stepj(),
                i2-i1,stepi(),NonConj TMV_FIRSTLAST); 
        }

        inline vec_type diag()
        {
            return vec_type(
                ptr(),TMV_MIN(colsize(),rowsize()),stepi()+stepj(),NonConj 
                TMV_FIRSTLAST);
        }

        inline vec_type diag(int i)
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const int diagstep = stepi() + stepj();
            if (i >= 0) {
                const int diagsize = TMV_MIN(colsize(),rowsize()-i);
                return vec_type(
                    ptr()+i*stepj(), diagsize,diagstep,NonConj 
                    TMV_FIRSTLAST);
            } else {
                const int diagsize = TMV_MIN(colsize()+i,rowsize());
                return vec_type(
                    ptr()-i*stepi(), diagsize,diagstep,NonConj 
                    TMV_FIRSTLAST);
            }
        }

        inline vec_type diag(int i, int j1, int j2) 
        {
            if (I == FortranStyle) { TMVAssert(j1 > 0 && j1-j2<=0); --j1; }
            else { TMVAssert(j1>=0 && j1-j2<=0); }
            TMVAssert(i>=-colsize() && i<=rowsize());
            const int diagstep = stepi() + stepj();
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(colsize(),rowsize()-i));
                return vec_type(
                    ptr()+i*stepj() + j1*diagstep, j2-j1, diagstep, NonConj 
                    TMV_FIRSTLAST);
            } else {
                TMVAssert(j2<=TMV_MIN(colsize(),rowsize()-i));
                return vec_type(
                    ptr()-i*stepi() + j1*diagstep, j2-j1, diagstep, NonConj 
                    TMV_FIRSTLAST);
            }
        }

        //
        // Modifying Functions
        //

        inline type& setZero() 
        { linearView().setZero(); return *this; }

        inline type& setAllTo(const T& x) 
        { linearView().setAllTo(x); return *this; }

        inline type& addToAll(const T& x) 
        { linearView().addToAll(x); return *this; }

        inline type& clip(RT thresh)
        { linearView().clip(thresh); return *this; }

        inline type& transposeSelf() 
        { view().transposeSelf(); return *this; }

        inline type& conjugateSelf() 
        { linearView().conjugateSelf(); return *this; }

        inline type& setToIdentity(const T& x=T(1)) 
        { 
            TMVAssert(colsize() == rowsize());
            setZero(); diag().setAllTo(x);
            return *this; 
        }

        inline type& swapRows(int i1, int i2)
        { 
            if (I == CStyle) { 
                TMVAssert(i1>=0 && i1<colsize());
                TMVAssert(i2>=0 && i2<colsize()); 
            } else { 
                TMVAssert(i1>0 && i1<=colsize());
                TMVAssert(i2>0 && i2<=colsize()); 
            }
            if (i1!=i2) Swap(row(i1),row(i2));
            return *this; 
        }

        inline type& swapCols(int j1, int j2)
        { 
            if (I == CStyle) { 
                TMVAssert(j1>=0 && j1<rowsize());
                TMVAssert(j2>=0 && j2<rowsize()); 
            } else { 
                TMVAssert(j1>0 && j1<=rowsize());
                TMVAssert(j2>0 && j2<=rowsize()); 
            }
            if (j1!=j2) Swap(col(j1),col(j2));
            return *this; 
        }

        inline type& permuteRows(const int* p, int i1, int i2)
        { view().permuteRows(p,i1,i2); return *this; }

        inline type& permuteRows(const int* p)
        { view().permuteRows(p); return *this; }

        inline type& permuteCols(const int* p, int j1, int j2)
        { view().permuteCols(p,j1,j2); return *this; }

        inline type& permuteCols(const int* p)
        { view().permuteCols(p); return *this; }

        inline type& reversePermuteRows(const int* p, int i1, int i2)
        { view().reversePermuteRows(p,i1,i2); return *this; }

        inline type& reversePermuteRows(const int* p)
        { view().reversePermuteRows(p); return *this; }

        inline type& reversePermuteCols(const int* p, int j1, int j2)
        { view().reversePermuteCols(p,j1,j2); return *this; }

        inline type& reversePermuteCols(const int* p)
        { view().reversePermuteCols(p); return *this; }

        //
        // subMatrix
        //

        inline const_view_type cSubMatrix(int i1, int i2, int j1, int j2) const
        {
            return const_view_type(
                itsm.get()+i1*stepi()+j1*stepj(),
                i2-i1,j2-j1,stepi(),stepj(),S,NonConj);
        }

        inline const_view_type subMatrix(int i1, int i2, int j1, int j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I == FortranStyle) { --i1; --j1; }
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline const_view_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            StorageType newstor =
                S == RowMajor ?
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            return const_view_type(
                itsm.get()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep,(j2-j1)/jstep,
                istep*stepi(),jstep*stepj(), 
                newstor,NonConj);
        }

        inline const_view_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I == FortranStyle) {
                --i1; --j1; 
                i2 += istep-1; j2 += jstep-1; 
            }
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline const_vec_type cSubVector(
            int i, int j, int istep, int jstep, int s) const
        {
            return const_vec_type(
                itsm.get()+i*stepi()+j*stepj(),s,
                istep*stepi()+jstep*stepj(),NonConj);
        }

        inline const_vec_type subVector(
            int i, int j, int istep, int jstep, int s) const
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,s));
            if (I==FortranStyle) { --i; --j; }
            return cSubVector(i,j,istep,jstep,s);
        }

        inline const_uppertri_type unitUpperTri() const
        {
            TMVAssert(rowsize() <= colsize());
            return const_uppertri_type(
                cptr(),rowsize(), stepi(),stepj(),UnitDiag,stor(),ct());
        }

        inline const_uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        {
            TMVAssert(rowsize() <= colsize());
            return const_uppertri_type(
                cptr(),rowsize(), stepi(),stepj(),dt,stor(),ct());
        }

        inline const_lowertri_type unitLowerTri() const
        {
            TMVAssert(colsize() <= rowsize());
            return const_lowertri_type(
                cptr(),colsize(), stepi(),stepj(),UnitDiag,stor(),ct());
        }

        inline const_lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        {
            TMVAssert(colsize() <= rowsize());
            return const_lowertri_type(
                cptr(),colsize(), stepi(),stepj(),dt,stor(),ct());
        }

        inline const_view_type cColPair(int j1, int j2) const
        {
            StorageType newstor =
                S == RowMajor ?
                j2==j1+1 ? RowMajor : NoMajor : ColMajor;
            return const_view_type(
                itsm.get()+j1*stepj(),colsize(),2,
                stepi(),(j2-j1)*stepj(),newstor,NonConj);
        }

        inline const_view_type colPair(int j1, int j2) const
        {
            if (I == CStyle) {
                TMVAssert(j1>=0 && j1<rowsize() && 
                          j2>=0 && j2<rowsize()); 
            } else  {
                TMVAssert(j1>0 && j1<=rowsize() && 
                          j2>0 && j2<=rowsize()); 
                --j1; --j2;
            }
            return cColPair(j1,j2);
        }

        inline const_view_type cRowPair(int i1, int i2) const
        {
            StorageType newstor = 
                S == RowMajor ?  RowMajor : i2==i1+1 ? ColMajor : NoMajor;
            return const_view_type(
                itsm.get()+i1*stepi(),2,rowsize(),
                (i2-i1)*stepi(),stepj(),newstor,NonConj);
        }

        inline const_view_type rowPair(int i1, int i2) const
        {
            if (I == CStyle)  {
                TMVAssert(i1>=0 && i1<colsize() && 
                          i2>=0 && i2<colsize());
            } else  {
                TMVAssert(i1>0 && i1<=colsize() && 
                          i2>0 && i2<=colsize()); 
                --i1; --i2;
            }
            return cRowPair(i1,i2);
        }

        inline const_view_type cColRange(int j1, int j2) const
        {
            return const_view_type(
                itsm.get()+j1*stepj(),colsize(),j2-j1,
                stepi(),stepj(),S,NonConj,iscm()?1:0);
        }

        inline const_view_type colRange(int j1, int j2) const
        {
            if (I==FortranStyle) { 
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize()); 
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize()); 
            }
            return cColRange(j1,j2);
        }

        inline const_view_type cRowRange(int i1, int i2) const
        {
            return const_view_type(
                itsm.get()+i1*stepi(),i2-i1,rowsize(),
                stepi(),stepj(),S,NonConj,isrm()?1:0);
        }

        inline const_view_type rowRange(int i1, int i2) const
        {
            if (I==FortranStyle) {
                TMVAssert(i1>0 && i1-i2<=0 && i2<=colsize()); 
                --i1;
            } else {
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize()); 
            }
            return cRowRange(i1,i2);
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm.get()),
                colsize(),rowsize(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? S : NoMajor,NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm.get())+1,
                colsize(),rowsize(),2*stepi(),2*stepj(),NoMajor,NonConj);
        }

        inline view_type cSubMatrix(int i1, int i2, int j1, int j2)
        {
            return view_type(
                ptr()+i1*stepi()+j1*stepj(),
                i2-i1,j2-j1,stepi(),stepj(),S,NonConj 
                TMV_FIRSTLAST);
        }

        inline view_type subMatrix(int i1, int i2, int j1, int j2)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==FortranStyle) { --i1; --j1; }
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline view_type cSubMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) 
        {
            StorageType newstor =
                S == RowMajor ? 
                jstep == 1 ? RowMajor : NoMajor :
                istep == 1 ? ColMajor : NoMajor;
            return view_type(
                ptr()+i1*stepi()+j1*stepj(),
                (i2-i1)/istep,(j2-j1)/jstep,
                istep*stepi(),jstep*stepj(), 
                newstor,NonConj TMV_FIRSTLAST);
        }

        inline view_type subMatrix(
            int i1, int i2, int j1, int j2, int istep, int jstep) 
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I == FortranStyle) {
                --i1; --j1;
                i2 += istep-1; j2 += jstep-1; 
            }
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline vec_type cSubVector(int i, int j, int istep, int jstep, int s) 
        {
            return vec_type(
                ptr()+i*stepi()+j*stepj(),s,
                istep*stepi()+jstep*stepj(),NonConj 
                TMV_FIRSTLAST);
        }

        inline vec_type subVector(int i, int j, int istep, int jstep, int s) 
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,s));
            if (I == FortranStyle) { --i; --j; }
            return cSubVector(i,j,istep,jstep,s);
        }

        inline uppertri_type unitUpperTri()
        {
            TMVAssert(rowsize() <= colsize());
            return uppertri_type(
                ptr(),rowsize(), stepi(),stepj(),UnitDiag,stor(),ct() 
                TMV_FIRSTLAST);
        }

        inline uppertri_type upperTri(DiagType dt=NonUnitDiag)
        {
            TMVAssert(rowsize() <= colsize());
            return uppertri_type(
                ptr(),rowsize(), stepi(),stepj(),dt,stor(),ct() 
                TMV_FIRSTLAST);
        }

        inline lowertri_type unitLowerTri()
        { 
            TMVAssert(colsize() <= rowsize());
            return lowertri_type(
                ptr(),colsize(), stepi(),stepj(),UnitDiag,stor(),ct() 
                TMV_FIRSTLAST);
        }

        inline lowertri_type lowerTri(DiagType dt=NonUnitDiag)
        { 
            TMVAssert(colsize() <= rowsize());
            return lowertri_type(
                ptr(),colsize(), stepi(),stepj(),dt,stor(),ct() 
                TMV_FIRSTLAST);
        }

        inline view_type cColPair(int j1, int j2) 
        {
            StorageType newstor = 
                S == RowMajor ?  j2==j1+1 ? RowMajor : NoMajor : ColMajor;
            return view_type(
                ptr()+j1*stepj(),colsize(),2,
                stepi(),(j2-j1)*stepj(),newstor,NonConj 
                TMV_FIRSTLAST);
        }

        inline view_type colPair(int j1, int j2) 
        {
            if (I == CStyle) 
                TMVAssert(j1>=0 && j1<rowsize() && 
                          j2>=0 && j2<rowsize());
            else {
                TMVAssert(j1>0 && j1<=rowsize() &&
                          j2>0 && j2<=rowsize());
                --j1; --j2;
            }
            return cColPair(j1,j2);
        }

        inline view_type cRowPair(int i1, int i2) 
        {
            StorageType newstor = 
                S == RowMajor ?  RowMajor : i2==i1+1 ? ColMajor : NoMajor;
            return view_type(
                ptr()+i1*stepi(),2,rowsize(),
                (i2-i1)*stepi(),stepj(),newstor,NonConj 
                TMV_FIRSTLAST);
        }

        inline view_type rowPair(int i1, int i2) 
        {
            if (I == CStyle)
                TMVAssert(i1>=0 && i1<colsize() && 
                          i2>=0 && i2<colsize());
            else {
                TMVAssert(i1>0 && i1<=colsize() && 
                          i2>0 && i2<=colsize());
                --i1; --i2;
            }
            return cRowPair(i1,i2);
        }

        inline view_type cColRange(int j1, int j2) 
        {
            return view_type(
                ptr()+j1*stepj(),colsize(),j2-j1,
                stepi(),stepj(),S,NonConj,iscm()?1:0 
                TMV_FIRSTLAST);
        }

        inline view_type colRange(int j1, int j2) 
        {
            if (I==FortranStyle) { 
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize());
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }
            return cColRange(j1,j2);
        }

        inline view_type cRowRange(int i1, int i2) 
        {
            return view_type(
                ptr()+i1*stepi(),i2-i1,rowsize(),
                stepi(),stepj(),S,NonConj,isrm()?1:0 
                TMV_FIRSTLAST);
        }

        inline view_type rowRange(int i1, int i2) 
        {
            if (I==FortranStyle) {
                TMVAssert(i1>0 && i1-i2<=0 && i2<=colsize());
                --i1;
            } else {
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            }
            return cRowRange(i1,i2);
        }

        inline realpart_type realPart() 
        {
            return realpart_type(
                reinterpret_cast<RT*>(ptr()),
                colsize(),rowsize(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),
                isReal(T()) ? S : NoMajor,NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        inline realpart_type imagPart() 
        {
            TMVAssert(isComplex(T()));
            return realpart_type(
                reinterpret_cast<RT*>(ptr())+1, 
                colsize(),rowsize(),
                2*stepi(),2*stepj(),NoMajor,NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        //
        // Views
        //

        inline const_view_type view() const
        { 
            return const_view_type(
                itsm.get(),colsize(),rowsize(),
                stepi(),stepj(),S,NonConj,linsize);
        }

        inline const_view_type transpose() const
        { 
            return const_view_type(
                itsm.get(),rowsize(),colsize(),
                stepj(),stepi(),TMV_TransOf(S),NonConj,linsize);
        }

        inline const_view_type conjugate() const
        { 
            return const_view_type(
                itsm.get(),colsize(),rowsize(),
                stepi(),stepj(),S,TMV_ConjOf(T,NonConj),linsize);
        }

        inline const_view_type adjoint() const
        { 
            return const_view_type(
                itsm.get(),rowsize(),colsize(),
                stepj(),stepi(),TMV_TransOf(S),TMV_ConjOf(T,NonConj),linsize);
        }

        inline const_vec_type constLinearView() const
        {
            TMVAssert(linsize == colsize()*rowsize());
            return const_vec_type(itsm.get(),linsize,1,NonConj); 
        }

        inline view_type view()
        { 
            return view_type(
                ptr(),colsize(),rowsize(),
                stepi(),stepj(),S,NonConj,linsize 
                TMV_FIRSTLAST);
        }

        inline view_type transpose()
        { 
            return view_type(
                ptr(),rowsize(),colsize(),
                stepj(),stepi(),TMV_TransOf(S),NonConj,linsize 
                TMV_FIRSTLAST);
        }

        inline view_type conjugate()
        { 
            return view_type(
                ptr(),colsize(),rowsize(),
                stepi(),stepj(),S,TMV_ConjOf(T,NonConj),linsize 
                TMV_FIRSTLAST);
        }

        inline view_type adjoint()
        { 
            return view_type(
                ptr(),rowsize(),colsize(),
                stepj(),stepi(),TMV_TransOf(S),TMV_ConjOf(T,NonConj),linsize 
                TMV_FIRSTLAST);
        }

        inline vec_type linearView()
        {
            TMVAssert(linsize == colsize()*rowsize());
            return vec_type(ptr(),linsize,1,NonConj TMV_FIRSTLAST); 
        }

        //
        // I/O
        //
        
        void read(const TMV_Reader& reader);

        virtual inline int ls() const { return linsize; }
        virtual inline int colsize() const { return itscs; }
        virtual inline int rowsize() const { return itsrs; }
        virtual inline const T* cptr() const { return itsm.get(); }
        inline T* ptr() { return itsm.get(); }
        virtual inline int stepi() const { return S == RowMajor ? itsrs : 1; }
        virtual inline int stepj() const { return S == RowMajor ? 1 : itscs; }
        inline bool isrm() const { return S==RowMajor; }
        inline bool iscm() const { return S==ColMajor; }
        inline bool isconj() const { return false; }
        virtual inline StorageType stor() const { return S; }
        virtual inline ConjType ct() const { return NonConj; }

        virtual inline bool canLinearize() const { return true; }

        virtual inline T cref(int i, int j) const
        { return itsm.get()[S==RowMajor ? i*stepi()+j : i+j*stepj()]; }

        inline T& ref(int i, int j)
        { return itsm.get()[S==RowMajor ? i*stepi()+j : i+j*stepj()]; }

        inline void resize(int cs, int rs)
        {
            TMVAssert(cs >= 0 && rs >= 0);
            linsize = cs*rs;
            itsm.resize(linsize);
            itscs = cs;
            itsrs = rs;
            DivHelper<T>::resetDivType();
#ifdef TMVFLDEBUG
            _first = itsm.get();
            _last = _first + linsize;
#endif
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline rowmajor_iterator rowmajor_begin() 
        { return MatrixIterHelper<S,type>::rowmajor_begin(this); }
        inline rowmajor_iterator rowmajor_end() 
        { return MatrixIterHelper<S,type>::rowmajor_end(this); }

        inline const_rowmajor_iterator rowmajor_begin() const 
        { return MatrixIterHelper<S,type>::rowmajor_begin(this); }
        inline const_rowmajor_iterator rowmajor_end() const 
        { return MatrixIterHelper<S,type>::rowmajor_end(this); }

        inline colmajor_iterator colmajor_begin()
        { return MatrixIterHelper<S,type>::colmajor_begin(this); }
        inline colmajor_iterator colmajor_end()
        { return MatrixIterHelper<S,type>::colmajor_end(this); }

        inline const_colmajor_iterator colmajor_begin() const 
        { return MatrixIterHelper<S,type>::colmajor_begin(this); }
        inline const_colmajor_iterator colmajor_end() const 
        { return MatrixIterHelper<S,type>::colmajor_end(this); }


    protected :

        int linsize;
        AlignedArray<T> itsm;
        int itscs;
        int itsrs;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
    protected:
#endif

        // If two matrices are the same size and storage, then 
        // swap can be much faster by copying the pointers to the data.
        template <IndexStyle I2>
        friend void Swap(Matrix<T,S,I>& m1, Matrix<T,S,I2>& m2)
        {
            TMVAssert(m1.colsize() == m2.colsize());
            TMVAssert(m1.rowsize() == m2.rowsize());
            m1.itsm.swapWith(m2.itsm);
#ifdef TMVFLDEBUG
            TMV_SWAP(m1._first,m2._first);
            TMV_SWAP(m1._last,m2._last);
#endif
        }

    }; // Matrix

    //-------------------------------------------------------------------------

    //
    // Special Creators: 
    //   RowVectorViewOf(v) = 1xn Matrix with v in only row - Same Storage
    //   ColVectorViewOf(v) = nx1 Matrix with v in only col - Same Storage
    //   MatrixViewOf(m,colsize,rowsize,storage) = MatrixView of m 
    //   MatrixViewOf(m,colsize,rowsize,stepi,stepj) = MatrixView of m 
    //

    template <class T> 
    inline ConstMatrixView<T> RowVectorViewOf(const GenVector<T>& v)
    {
        return ConstMatrixView<T>(
            v.cptr(),1,v.size(),v.size(),v.step(),
            v.step()==1?RowMajor:NoMajor,v.ct(),v.step()==1?v.size():0);
    }

    template <class T, IndexStyle I> 
    inline ConstMatrixView<T,I> RowVectorViewOf(
        const ConstVectorView<T,I>& v)
    {
        return ConstMatrixView<T,I>(
            v.cptr(),1,v.size(),v.size(),v.step(),
            v.step()==1?RowMajor:NoMajor,v.ct(),v.step()==1?v.size():0);
    }

    template <class T, IndexStyle I> 
    inline ConstMatrixView<T,I> RowVectorViewOf(const Vector<T,I>& v)
    {
        return ConstMatrixView<T,I>(
            v.cptr(),1,v.size(),v.size(),1,RowMajor,v.ct(),v.size());
    }

    template <class T, IndexStyle I> 
    inline MatrixView<T,I> RowVectorViewOf(const VectorView<T,I>& v)
    {
        return MatrixView<T,I>(
            v.ptr(),1,v.size(),v.size(),v.step(),
            v.step()==1?RowMajor:NoMajor,v.ct(),v.step()==1?v.size():0
            TMV_FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

    template <class T, IndexStyle I> 
    inline MatrixView<T,I> RowVectorViewOf(Vector<T,I>& v)
    {
        return MatrixView<T,I>(
            v.ptr(),1,v.size(),v.size(),1,RowMajor,v.ct(),v.size()
            TMV_FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

    template <class T> 
    inline ConstMatrixView<T> ColVectorViewOf(const GenVector<T>& v)
    { 
        return ConstMatrixView<T>(
            v.cptr(),v.size(),1,v.step(),v.size(),
            v.step()==1?ColMajor:NoMajor,v.ct(),v.step()==1?v.size():0);
    }

    template <class T, IndexStyle I> 
    inline ConstMatrixView<T,I> ColVectorViewOf(
        const ConstVectorView<T,I>& v)
    { 
        return ConstMatrixView<T,I>(
            v.cptr(),v.size(),1,v.step(),v.size(),
            v.step()==1?ColMajor:NoMajor,v.ct(),v.step()==1?v.size():0);
    }

    template <class T, IndexStyle I> 
    inline ConstMatrixView<T,I> ColVectorViewOf(const Vector<T,I>& v)
    { 
        return ConstMatrixView<T,I>(
            v.cptr(),v.size(),1,1,v.size(),ColMajor,v.ct(),v.size());
    }

    template <class T, IndexStyle I> 
    inline MatrixView<T,I> ColVectorViewOf(const VectorView<T,I>& v)
    { 
        return MatrixView<T,I>(
            v.ptr(),v.size(),1,v.step(),v.size(),
            v.step()==1?ColMajor:NoMajor,v.ct(),v.step()==1?v.size():0
            TMV_FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

    template <class T, IndexStyle I> 
    inline MatrixView<T,I> ColVectorViewOf(Vector<T,I>& v)
    { 
        return MatrixView<T,I>(
            v.ptr(),v.size(),1,1,v.size(),ColMajor,v.ct(),v.size() 
            TMV_FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

    template <class T> 
    inline MatrixView<T> MatrixViewOf(
        T* m, int colsize, int rowsize, StorageType stor)
    {
        TMVAssert(colsize >= 0 && rowsize >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        const int linsize = colsize * rowsize;
        const int stepi = stor == RowMajor ? rowsize : 1;
        const int stepj = stor == RowMajor ? 1 : colsize;
        return MatrixView<T>(
            m,colsize,rowsize,stepi,stepj,stor,NonConj,linsize 
            TMV_FIRSTLAST1(m,m+linsize));
    }

    template <class T> 
    inline ConstMatrixView<T> MatrixViewOf(
        const T* m, int colsize, int rowsize, StorageType stor)
    {
        TMVAssert(colsize >= 0 && rowsize >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        const int linsize = colsize*rowsize;
        const int stepi = stor == RowMajor ? rowsize : 1;
        const int stepj = stor == RowMajor ? 1 : colsize;
        return ConstMatrixView<T>(
            m,colsize,rowsize,stepi,stepj,stor,NonConj,linsize);
    }

    template <class T> 
    inline MatrixView<T> MatrixViewOf(
        T* m, int colsize, int rowsize, int stepi, int stepj)
    {
        TMVAssert(colsize >= 0 && rowsize >= 0);
        const StorageType stor = (
            stepi==1 ? ColMajor :
            stepj==1 ? RowMajor : NoMajor );
        const int linsize = (
            (stepi==1 && stepj==colsize) ? colsize * rowsize :
            (stepj==1 && stepi==rowsize) ? colsize * rowsize : 
            0 );
        return MatrixView<T>(
            m,colsize,rowsize,stepi,stepj,stor,NonConj,linsize 
            TMV_FIRSTLAST1(m,m+stepi*(colsize-1)+stepj*(rowsize-1)+1));
    }

    template <class T> 
    inline ConstMatrixView<T> MatrixViewOf(
        const T* m, int colsize, int rowsize, int stepi, int stepj)
    {
        TMVAssert(colsize >= 0 && rowsize >= 0);
        const StorageType stor = (
            stepi==1 ? ColMajor :
            stepj==1 ? RowMajor : NoMajor );
        const int linsize = (
            (stepi==1 && stepj==colsize) ? colsize * rowsize :
            (stepj==1 && stepi==rowsize) ? colsize * rowsize : 
            0 );
        return ConstMatrixView<T>(
            m,colsize,rowsize,stepi,stepj,stor,NonConj,linsize);
    }



    //
    // Copy Matrices
    //

    template <class T> 
    void DoCopySameType(const GenMatrix<T>& m1, const MatrixView<T>& m2);

    template <class T> 
    inline void DoCopy(const GenMatrix<T>& m1, const MatrixView<T>& m2)
    { DoCopySameType(m1,m2); }

    template <class T, class T1> 
    inline void DoCopyDiffType(const GenMatrix<T1>& m1, const MatrixView<T>& m2)
    {
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m2.rowsize() > 0);
        TMVAssert(m2.colsize() > 0);
        TMVAssert(m1.ct()==NonConj);
        TMVAssert(m2.ct()==NonConj);
        TMVAssert(!m2.isrm());
        TMVAssert(!m2.isSameAs(m1));
        TMVAssert(isComplex(T()) || isReal(T1()));

        if (m1.iscm() && m2.iscm() && m1.colsize() > 1) 
            for(int j=0;j<m2.rowsize();++j)
                DoCopyDiffType(m1.col(j),m2.col(j));
        else if (m2.colsize() < m2.rowsize())
            if (shouldReverse(m1.stepj(),m2.stepj()))
                for(int i=0;i<m2.colsize();++i) 
                    DoCopyDiffType(m1.row(i).reverse(),m2.row(i).reverse());
            else
                for(int i=0;i<m2.colsize();++i) 
                    DoCopyDiffType(m1.row(i),m2.row(i));
        else
            if (shouldReverse(m1.stepi(),m2.stepi()))
                for(int j=0;j<m2.rowsize();++j) 
                    DoCopyDiffType(m1.col(j).reverse(),m2.col(j).reverse());
            else
                for(int j=0;j<m2.rowsize();++j) 
                    DoCopyDiffType(m1.col(j),m2.col(j));
    }

    template <class T, class T1> 
    inline void DoCopy(const GenMatrix<T1>& m1, const MatrixView<T>& m2)
    { DoCopyDiffType(m1,m2); }

    template <class T> 
    inline void DoCopy(const GenMatrix<std::complex<T> >&, const MatrixView<T>&)
    { TMVAssert(TMV_FALSE); }

    template <class T1, class T2> 
    inline void nonconjCopy(const GenMatrix<T1>& m1, const MatrixView<T2>& m2)
    {
        TMVAssert(isComplex(T2()) || isReal(T1()));
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m1.ct() == NonConj);
        TMVAssert(m2.ct() == NonConj);

        if (m2.isrm() || (m1.isrm() && !m2.iscm()))
            DoCopy(m1.transpose(),m2.transpose());
        else
            DoCopy(m1,m2);
    }

    template <class T1, class T2> 
    inline void Copy(const GenMatrix<T1>& m1, const MatrixView<T2>& m2)
    {
        TMVAssert(isComplex(T2()) || isReal(T1()));
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        if (m2.rowsize() > 0 && m2.colsize() > 0) {
            if (SameStorage(m1,m2)) {
                if (m2.isSameAs(m1)) {
                    // Do Nothing
                } else if (m2.transpose().isSameAs(m1)) {
                    m2.transposeSelf(); 
                } else if (m1.isrm()) {
                    Matrix<T1,RowMajor> m1x = m1;
                    m2 = m1x;
                } else {
                    Matrix<T1,ColMajor> m1x = m1;
                    m2 = m1x;
                }
            } else if (m1.stor()==m2.stor() && m1.canLinearize() && 
                       m2.canLinearize()) {
                TMVAssert(m1.constLinearView().size()==m2.linearView().size());
                TMVAssert((m1.stepi()==m2.stepi() && m1.stepj()==m2.stepj()));
                m2.linearView() = m1.constLinearView();
            } else {
                if (m1.isconj()) {
                    if (m2.isconj()) {
                        nonconjCopy(m1.conjugate(),m2.conjugate());
                    } else {
                        nonconjCopy(m1.conjugate(),m2);
                        m2.conjugateSelf();
                    }
                } else {
                    if (m2.isconj()) {
                        nonconjCopy(m1,m2.conjugate());
                        m2.conjugateSelf();
                    } else  {
                        nonconjCopy(m1,m2);
                    }
                }
            }
        }
    }

    //
    // Swap Matrices
    //

    template <class T> 
    void Swap(const MatrixView<T>& m1, const MatrixView<T>& m2);

    template <class T, StorageType S, IndexStyle I> 
    inline void Swap(const MatrixView<T>& m1, Matrix<T,S,I>& m2)
    { Swap(m1,m2.view()); }

    template <class T, StorageType S, IndexStyle I> 
    inline void Swap(Matrix<T,S,I>& m1, const MatrixView<T>& m2)
    { Swap(m1.view(),m2); }

    template <class T, StorageType S1, IndexStyle I1, StorageType S2, IndexStyle I2> 
    inline void Swap(Matrix<T,S1,I1>& m1, Matrix<T,S2,I2>& m2)
    { Swap(m1.view(),m2.view()); }


    //
    // Views of a Matrix:
    //

    template <class T> 
    inline ConstMatrixView<T> Transpose(const GenMatrix<T>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline ConstMatrixView<T,I> Transpose(const ConstMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstMatrixView<T,I> Transpose(const Matrix<T,S,I>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline MatrixView<T,I> Transpose(const MatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, StorageType S, IndexStyle I> 
    inline MatrixView<T,I> Transpose(Matrix<T,S,I>& m)
    { return m.transpose(); }

    template <class T> 
    inline ConstMatrixView<T> Conjugate(const GenMatrix<T>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline ConstMatrixView<T,I> Conjugate(const ConstMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstMatrixView<T,I> Conjugate(const Matrix<T,S,I>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline MatrixView<T,I> Conjugate(const MatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, StorageType S, IndexStyle I> 
    inline MatrixView<T,I> Conjugate(Matrix<T,S,I>& m)
    { return m.conjugate(); }

    template <class T> 
    inline ConstMatrixView<T> Adjoint(const GenMatrix<T>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline ConstMatrixView<T,I> Adjoint(const ConstMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, StorageType S, IndexStyle I> 
    inline ConstMatrixView<T,I> Adjoint(const Matrix<T,S,I>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline MatrixView<T,I> Adjoint(const MatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, StorageType S, IndexStyle I> 
    inline MatrixView<T,I> Adjoint(Matrix<T,S,I>& m)
    { return m.adjoint(); }

    template <class T> 
    inline QuotXM<T,T> Inverse(const GenMatrix<T>& m)
    { return m.inverse(); }


    //
    // Matrix ==, != Matrix
    //

    template <class T1, class T2> 
    bool operator==(
        const GenMatrix<T1>& m1, const GenMatrix<T2>& m2);
    template <class T1, class T2> 
    inline bool operator!=(
        const GenMatrix<T1>& m1, const GenMatrix<T2>& m2)
    { return !(m1 == m2); }


    //
    // I/O
    //

    template <class T>
    inline std::istream& operator>>(
        const TMV_Reader& reader, const MatrixView<T>& m)
    { m.read(reader); return reader.getis(); }

    template <class T, StorageType S, IndexStyle I> 
    inline std::istream& operator>>(
        const TMV_Reader& reader, Matrix<T,S,I>& m)
    { m.read(reader); return reader.getis(); }

    template <class T>
    std::istream& operator>>(std::istream& is, const MatrixView<T>& m)
    { return is >> IOStyle() >> m; }

    template <class T, StorageType S, IndexStyle I> 
    std::istream& operator>>(std::istream& is, Matrix<T,S,I>& m)
    { return is >> IOStyle() >> m; }

} // namespace tmv

#endif
