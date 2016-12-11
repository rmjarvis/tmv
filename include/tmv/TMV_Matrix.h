///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// The Template Matrix/Vector Library for C++ was created by Mike Jarvis     //
// Copyright (C) 1998 - 2016                                                 //
// All rights reserved                                                       //
//                                                                           //
// The project is hosted at https://code.google.com/p/tmv-cpp/               //
// where you can find the current version and current documention.           //
//                                                                           //
// For concerns or problems with the software, Mike may be contacted at      //
// mike_jarvis17 [at] gmail.                                                 //
//                                                                           //
// This software is licensed under a FreeBSD license.  The file              //
// TMV_LICENSE should have bee included with this distribution.              //
// It not, you can get a copy from https://code.google.com/p/tmv-cpp/.       //
//                                                                           //
// Essentially, you can use this software however you want provided that     //
// you include the TMV_LICENSE file in any distribution that uses it.        //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


//---------------------------------------------------------------------------
//
// This file defines the TMV Matrix class.
//
// The Matrix class and all associated functions are contained
// in the namespace tmv.  Alse, the Matrix class is a template,
// so for a Matrix of doubles, one would write tmv::Matrix<double>.
//
// An optional second template parameter specifies the known attributes
// of the matrix.  The valid attributes for a Matrix are:
// - ColMajor or RowMajor
// - CStyle or FortranStyle
// The defaults are ColMajor and CStyle if you do not specify otherwise.
//
// The attributes are treated as a bit field, so you | them together to
// get the complete value.  e.g.
// Matrix<double, ColMajor | FortranStyle>
//
// Also, you only need to specify things that differ from the default, so
// Matrix<T,RowMajor> means Matrix<T,RowMajor|CStyle>, and
// Matrix<T> means Matrix<T,ColMajor|CStyle>.
//
// The *Style attributes indicate whether to use C-style or Fotran-style
// indexing.
// With C-style (the default), the upper left corner of an MxN matrix is
// m(0,0), and the lower right is m(M-1,N-1).
// With Fortran-style, these are m(1,1) and m(M,N) respectively.
// Also, when a function takes an index range, i1,i2,
// then with C-style, this means elements from i1...i2-1 inclusive.
// With Fortran-style, this means i1..i2 inclusive.
//
// Constructors:
//
//    Matrix<T,A>(int colsize, int rowsize)
//        Makes a Matrix with column size = colsize and row size = rowsize
//        with _uninitialized_ values
//
//    Matrix<T,A>(int colsize, int rowsize, T x)
//        Makes a Matrix of size n with all values = x
//
//
// Special Constructors
//
//    MatrixView RowVectorViewOf(Vector& v)
//    MatrixView RowVectorViewOf(VectorView v)
//    ConstMatrixView RowVectorViewOf(const Vector& v)
//    ConstMatrixView RowVectorViewOf(const ConstVectorView& v)
//        Returns a 1xn MatrixView with v in the (only) row.
//        Unlike creating a Matrix with RowVector(v), this uses the same
//        data storage as the actual Vector v.
//        For a const Vector or a ConstVectorView, this returns a
//        ConstMatrixView.
//
//    MatrixView ColVectorViewOf(Vector& v)
//    MatrixView ColVectorViewOf(VectorView v)
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

    template <typename T1, typename T2>
    inline void Copy(const GenMatrix<T1>& m1, MatrixView<T2> m2);

    template <typename T, typename Tm>
    class QuotXM;

    template <typename T>
    class LUDiv;

    template <typename T>
    class QRDiv;

    template <typename T>
    class QRPDiv;

    template <typename T>
    class SVDiv;

    template <typename T>
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

        inline const_vec_type row(ptrdiff_t i) const
        {
            TMVAssert(i>=0 && i<colsize());
            return const_vec_type(
                cptr()+i*ptrdiff_t(stepi()),rowsize(),stepj(),ct());
        }

        inline const_vec_type col(ptrdiff_t j) const
        {
            TMVAssert(j>=0 && j<rowsize());
            return const_vec_type(
                cptr()+j*ptrdiff_t(stepj()),colsize(),stepi(),ct());
        }

        inline const_vec_type diag() const
        {
            return const_vec_type(
                cptr(),TMV_MIN(colsize(),rowsize()),stepi()+stepj(),ct());
        }

        inline const_vec_type diag(ptrdiff_t i) const
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const ptrdiff_t diagstep = stepi() + stepj();
            if (i >= 0) {
                const ptrdiff_t diagsize = TMV_MIN(colsize(),rowsize()-i);
                return const_vec_type(
                    cptr()+i*ptrdiff_t(stepj()),diagsize,diagstep,ct());
            } else {
                const ptrdiff_t diagsize = TMV_MIN(colsize()+i,rowsize());
                return const_vec_type(
                    cptr()-i*ptrdiff_t(stepi()),diagsize,diagstep,ct());
            }
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            return const_vec_type(
                cptr()+i*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),j2-j1,stepj(),ct());
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(j>=0 && j<rowsize());
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            return const_vec_type(
                cptr()+i1*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),i2-i1,stepi(),ct());
        }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const ptrdiff_t diagstep = stepi() + stepj();
            if (i >= 0) {
                TMVAssert(j1>=0 && j1-j2<=0 &&
                          j2<=TMV_MIN(colsize(),rowsize()-i));
                return const_vec_type(
                    cptr()+i*ptrdiff_t(stepj())+j1*diagstep,j2-j1, diagstep,ct());
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 &&
                          j2<=TMV_MIN(colsize()+i,rowsize()));
                return const_vec_type(
                    cptr()-i*ptrdiff_t(stepi())+j1*diagstep,j2-j1, diagstep,ct());
            }
        }

        inline T operator()(ptrdiff_t i, ptrdiff_t j) const
        {
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j>=0 && j<rowsize());
            return cref(i,j);
        }

        inline const_vec_type operator[](ptrdiff_t i) const
        {
            TMVAssert(i>=0 && i<colsize());
            return row(i);
        }

        template <typename T2>
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

        inline void assignToM(MatrixView<RT> m2) const
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVAssert(isReal(T()));
            if (!isSameAs(m2)) Copy(*this,m2);
        }

        inline void assignToM(MatrixView<CT> m2) const
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            if (!isSameAs(m2)) Copy(*this,m2);
        }

        //
        // subMatrix
        //

        bool hasSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const;

        inline const_view_type cSubMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            return const_view_type(
                cptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1,j2-j1,stepi(),stepj(),ct());
        }

        inline const_view_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline const_view_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            return const_view_type(
                cptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep,(j2-j1)/jstep,istep*stepi(),jstep*stepj(),ct());
        }

        inline const_view_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        bool hasSubVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const;

        inline const_vec_type cSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        {
            return const_vec_type(
                cptr()+i*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),s,
                istep*stepi()+jstep*stepj(),ct());
        }

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,s));
            return cSubVector(i,j,istep,jstep,s);
        }

        inline const_uppertri_type unitUpperTri() const
        {
            TMVAssert(rowsize() <= colsize());
            return const_uppertri_type(
                cptr(),rowsize(),stepi(),stepj(),UnitDiag,ct() );
        }

        inline const_uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        {
            TMVAssert(rowsize() <= colsize());
            return const_uppertri_type(
                cptr(),rowsize(),stepi(),stepj(),dt,ct() );
        }

        inline const_lowertri_type unitLowerTri() const
        {
            TMVAssert(colsize() <= rowsize());
            return const_lowertri_type(
                cptr(),colsize(),stepi(),stepj(),UnitDiag,ct() );
        }

        inline const_lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        {
            TMVAssert(colsize() <= rowsize());
            return const_lowertri_type(
                cptr(),colsize(),stepi(),stepj(),dt,ct() );
        }

        inline const_view_type cColPair(ptrdiff_t j1, ptrdiff_t j2) const
        {
            return const_view_type(
                cptr()+j1*ptrdiff_t(stepj()),colsize(),2,stepi(),(j2-j1)*stepj(),ct());
        }

        inline const_view_type colPair(ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1>=0 && j1<rowsize() && j2>=0 && j2<rowsize());
            return cColPair(j1,j2);
        }

        inline const_view_type cRowPair(ptrdiff_t i1, ptrdiff_t i2) const
        {
            return const_view_type(
                cptr()+i1*ptrdiff_t(stepi()),2,rowsize(),(i2-i1)*stepi(),stepj(),ct());
        }

        inline const_view_type rowPair(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(i1>=0 && i1<colsize() && i2>=0 && i2<colsize());
            return cRowPair(i1,i2);
        }

        inline const_view_type cColRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            return const_view_type(
                cptr()+j1*ptrdiff_t(stepj()),colsize(),j2-j1,
                stepi(),stepj(),ct(),(iscm()&&ls())?-1:0);
        }

        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            return cColRange(j1,j2);
        }

        inline const_view_type cRowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            return const_view_type(
                cptr()+i1*ptrdiff_t(stepi()),i2-i1,rowsize(),
                stepi(),stepj(),ct(),(isrm()&&ls())?-1:0);
        }

        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            return cRowRange(i1,i2);
        }

        inline const_realpart_type realPart() const
        {
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr()),colsize(),rowsize(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(),NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(cptr())+1,
                colsize(),rowsize(),2*stepi(),2*stepj(),NonConj);
        }

        //
        // Views
        //

        inline const_view_type view() const
        {
            return const_view_type(
                cptr(),colsize(),rowsize(),stepi(),stepj(),ct(),ls());
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                cptr(),rowsize(),colsize(),stepj(),stepi(),ct(),ls());
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                cptr(),colsize(),rowsize(),
                stepi(),stepj(),TMV_ConjOf(T,ct()),ls());
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                cptr(),rowsize(),colsize(),
                stepj(),stepi(),TMV_ConjOf(T,ct()),ls());
        }

        inline const_vec_type constLinearView() const
        {
            TMVAssert(ls() != -1);
            // To assure that next Assert has no effect

            TMVAssert(canLinearize());
            TMVAssert(ls() == colsize()*rowsize());
            return const_vec_type(cptr(),ls(),1,ct());
        }

        inline nonconst_type nonConst() const
        {
            return nonconst_type(
                const_cast<T*>(cptr()),colsize(),rowsize(),
                stepi(),stepj(),ct(),ls()
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
        virtual ptrdiff_t stepi() const = 0;
        virtual ptrdiff_t stepj() const = 0;
        virtual ptrdiff_t ls() const = 0;
        virtual inline bool isrm() const { return stepj() == 1; }
        virtual inline bool iscm() const { return stepi() == 1; }
        virtual inline bool isconj() const
        { return isComplex(T()) && ct()==Conj; }
        virtual ConjType ct() const =0;

        virtual bool canLinearize() const = 0;
        virtual T cref(ptrdiff_t i, ptrdiff_t j) const;

        inline ptrdiff_t rowstart(ptrdiff_t ) const { return 0; }
        inline ptrdiff_t rowend(ptrdiff_t ) const { return rowsize(); }
        inline ptrdiff_t colstart(ptrdiff_t ) const { return 0; }
        inline ptrdiff_t colend(ptrdiff_t ) const { return colsize(); }

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

    template <typename T, int A>
    class ConstMatrixView : public GenMatrix<T>
    {
    public :

        typedef GenMatrix<T> base;
        typedef ConstMatrixView<T,A> type;

        inline ConstMatrixView(const type& rhs) :
            itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs),
            itssi(rhs.itssi), itssj(rhs.itssj),
            itsct(rhs.itsct), linsize(rhs.linsize)
        { TMVAssert(Attrib<A>::viewok); }

        inline ConstMatrixView(const base& rhs) :
            itsm(rhs.cptr()), itscs(rhs.colsize()), itsrs(rhs.rowsize()),
            itssi(rhs.stepi()), itssj(rhs.stepj()),
            itsct(rhs.ct()), linsize(rhs.ls())
        { TMVAssert(Attrib<A>::viewok); }

        inline ConstMatrixView(
            const T* _m, ptrdiff_t _cs, ptrdiff_t _rs, ptrdiff_t _si, ptrdiff_t _sj,
            ConjType _ct, ptrdiff_t _ls=0) :
            itsm(_m), itscs(_cs), itsrs(_rs), itssi(_si), itssj(_sj),
            itsct(_ct), linsize(_ls)
        {
            TMVAssert(Attrib<A>::viewok);
            TMVAssert(linsize==0 || linsize==-1 ||
                      ((itssi==1 || itssj==1) && linsize == itscs*itsrs));
        }

        virtual inline ~ConstMatrixView()
        {
#ifdef TMV_EXTRA_DEBUG
            const_cast<const T*&>(itsm) = 0;
#endif
        }

        virtual inline ptrdiff_t colsize() const { return itscs; }
        virtual inline ptrdiff_t rowsize() const { return itsrs; }
        virtual inline const T* cptr() const { return itsm; }
        virtual inline ptrdiff_t stepi() const { return itssi; }
        virtual inline ptrdiff_t stepj() const { return itssj; }
        virtual inline ConjType ct() const { return itsct; }
        virtual inline ptrdiff_t ls() const { return linsize; }
        using GenMatrix<T>::isrm;
        using GenMatrix<T>::iscm;

        virtual inline bool canLinearize() const
        {
            if (linsize == -1) {
                if ( (stepi() == 1 && stepj() == colsize()) ||
                     (stepj() == 1 && stepi() == rowsize()) )
                    linsize = rowsize() * colsize();
                else
                    linsize = 0;
            }
            TMVAssert(linsize == 0 ||
                      (this->isrm() && stepi() == rowsize()) ||
                      (this->iscm() && stepj() == colsize()));
            return linsize > 0;
        }

    protected :

        const T*const itsm;
        const ptrdiff_t itscs;
        const ptrdiff_t itsrs;
        const ptrdiff_t itssi;
        const ptrdiff_t itssj;
        const ConjType itsct;

        mutable ptrdiff_t linsize;

    private :

        type& operator=(const type&);

    }; // ConstMatrixView

    template <typename T>
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
            const T* _m, ptrdiff_t _cs, ptrdiff_t _rs, ptrdiff_t _si, ptrdiff_t _sj,
            ConjType _ct, ptrdiff_t linsize=0) :
            c_type(_m,_cs,_rs,_si,_sj,_ct,linsize) {}

        virtual inline ~ConstMatrixView() {}


        //
        // Access
        //

        inline T operator()(ptrdiff_t i, ptrdiff_t j)
        {
            TMVAssert(i>0 && i<=this->colsize());
            TMVAssert(j>0 && j<=this->rowsize());
            return base::cref(i-1,j-1);
        }

        inline const_vec_type row(ptrdiff_t i) const
        {
            TMVAssert(i>0 && i<=this->colsize());
            return base::row(i-1);
        }

        inline const_vec_type col(ptrdiff_t j) const
        {
            TMVAssert(j>0 && j<=this->rowsize());
            return base::col(j-1);
        }

        inline const_vec_type diag() const
        { return base::diag(); }

        inline const_vec_type diag(ptrdiff_t i) const
        { return base::diag(i); }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>0 && i<=this->colsize());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->rowsize());
            return base::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(j>0 && j<=this->rowsize());
            TMVAssert(i1 > 0 && i1 <= i2 && i2 <= this->colsize());
            return base::col(j-1,i1-1,i2);
        }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1 > 0);
            return base::diag(i,j1-1,j2);
        }

        inline const_vec_type operator[](ptrdiff_t i) const
        { return row(i); }

        //
        // subMatrix
        //

        bool hasSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const;

        inline const_view_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return base::cSubMatrix(i1-1,i2,j1-1,j2);
        }

        inline const_view_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return base::cSubMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        bool hasSubVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const;

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
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

        inline const_view_type colPair(ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1 > 0 && j1 <= this->rowsize());
            TMVAssert(j2 > 0 && j2 <= this->rowsize());
            return base::cColPair(j1-1,j2-1);
        }

        inline const_view_type rowPair(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(i1 > 0 && i1 <= this->colsize());
            TMVAssert(i2 > 0 && i2 <= this->colsize());
            return base::cRowPair(i1-1,i2-1);
        }

        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1 > 0 && j1 <= j2 && j2 <= this->rowsize());
            return base::cColRange(j1-1,j2);
        }

        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
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

    template <typename T, int A>
    class MatrixView : public GenMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenMatrix<T> base;
        typedef MatrixView<T,A> type;
        typedef MatrixView<T,A> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef MatrixView<RT,A> realpart_type;
        typedef VectorView<T,A> vec_type;
        typedef UpperTriMatrixView<T,A> uppertri_type;
        typedef LowerTriMatrixView<T,A> lowertri_type;
        typedef typename RefHelper<T>::reference reference;
        typedef ConstMatrixView<T,A> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstMatrixView<RT,A> const_realpart_type;
        typedef ConstVectorView<T,A> const_vec_type;
        typedef ConstUpperTriMatrixView<T,A> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,A> const_lowertri_type;
        typedef RMIt<type> rowmajor_iterator;
        typedef CMIt<type> colmajor_iterator;
        typedef RMIt<const type> const_rowmajor_iterator;
        typedef CMIt<const type> const_colmajor_iterator;

        //
        // Constructors
        //

        inline MatrixView(const type& rhs) :
            itsm(rhs.itsm), itscs(rhs.itscs), itsrs(rhs.itsrs),
            itssi(rhs.itssi), itssj(rhs.itssj),
            itsct(rhs.itsct), linsize(rhs.linsize)
            TMV_DEFFIRSTLAST(rhs._first,rhs._last)
        { TMVAssert(Attrib<A>::viewok); }

        inline MatrixView(
            T* _m, ptrdiff_t _cs, ptrdiff_t _rs, ptrdiff_t _si, ptrdiff_t _sj,
            ConjType _ct, ptrdiff_t _ls=0 TMV_PARAMFIRSTLAST(T) ) :
            itsm(_m), itscs(_cs), itsrs(_rs), itssi(_si), itssj(_sj),
            itsct(_ct), linsize(_ls)
            TMV_DEFFIRSTLAST(_first,_last)
        {
            TMVAssert(Attrib<A>::viewok);
            TMVAssert(linsize==0 || linsize==-1 ||
                      ((itssi==1 || itssj==1) && linsize == itscs*itsrs));
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

        inline type& operator=(const MatrixView<T,A>& m2)
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            m2.assignToM(*this);
            return *this;
        }

        inline type& operator=(const GenMatrix<RT>& m2)
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            m2.assignToM(*this);
            return *this;
        }

        inline type& operator=(const GenMatrix<CT>& m2)
        {
            TMVAssert(m2.colsize() == colsize());
            TMVAssert(m2.rowsize() == rowsize());
            TMVAssert(isComplex(T()));
            m2.assignToM(*this);
            return *this;
        }

        template <typename T2>
        inline type& operator=(const GenMatrix<T2>& m2)
        {
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(m2,*this);
            return *this;
        }

        inline type& operator=(const T& x)
        { return setToIdentity(x); }

        inline type& operator=(const AssignableToMatrix<RT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            m2.assignToM(*this);
            return *this;
        }

        inline type& operator=(const AssignableToMatrix<CT>& m2)
        {
            TMVAssert(colsize() == m2.colsize());
            TMVAssert(rowsize() == m2.rowsize());
            TMVAssert(isComplex(T()));
            m2.assignToM(*this);
            return *this;
        }

        inline type& operator=(const Permutation& m2)
        {
            m2.assignToM(*this);
            return *this;
        }

        template <typename T2, ptrdiff_t M, ptrdiff_t N, int A2>
        inline type& operator=(const SmallMatrix<T2,M,N,A2>& m2)
        {
            TMVAssert(colsize() == M && rowsize() == N);
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(m2.view(),*this);
            return *this;
        }

        typedef ListAssigner<T,rowmajor_iterator> MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { return MyListAssigner(rowmajor_begin(),colsize()*rowsize(),x); }

        //
        // Access
        //

        inline reference operator()(ptrdiff_t i,ptrdiff_t j)
        {
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j>=0 && j<rowsize());
            return ref(i,j);
        }

        inline vec_type operator[](ptrdiff_t i)
        {
            TMVAssert(i>=0 && i<colsize());
            return row(i);
        }

        inline vec_type row(ptrdiff_t i)
        {
            TMVAssert(i>=0 && i<colsize());
            return vec_type(
                ptr()+i*ptrdiff_t(stepi()), rowsize(),stepj(),ct() TMV_FIRSTLAST );
        }

        inline vec_type col(ptrdiff_t j)
        {
            TMVAssert(j>=0 && j<rowsize());
            return vec_type(
                ptr()+j*ptrdiff_t(stepj()), colsize(),stepi(),ct() TMV_FIRSTLAST );
        }

        inline vec_type diag()
        {
            return vec_type(
                ptr(), TMV_MIN(colsize(),rowsize()),stepi()+stepj(),ct()
                TMV_FIRSTLAST);
        }

        inline vec_type diag(ptrdiff_t i)
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const ptrdiff_t diagstep = stepi() + stepj();
            if (i >= 0) {
                const ptrdiff_t diagsize = TMV_MIN(colsize(),rowsize()-i);
                return vec_type(
                    ptr()+i*ptrdiff_t(stepj()), diagsize,diagstep,ct() TMV_FIRSTLAST );
            } else {
                const ptrdiff_t diagsize = TMV_MIN(colsize()+i,rowsize());
                return vec_type(
                    ptr()-i*ptrdiff_t(stepi()), diagsize,diagstep,ct() TMV_FIRSTLAST );
            }
        }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>=0 && i<colsize());
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            return vec_type(
                ptr()+i*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()), j2-j1,stepj(),ct() TMV_FIRSTLAST );
        }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(j>=0 && j<rowsize());
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            return vec_type(
                ptr()+i1*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()), i2-i1,stepi(),ct() TMV_FIRSTLAST );
        }

        inline vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const ptrdiff_t diagstep = stepi() + stepj();
            if (i >= 0) {
                TMVAssert(
                    j1>=0 && j1-j2<=0 &&
                    j2<=TMV_MIN(colsize(),rowsize()-i));
                return vec_type(
                    ptr()+i*ptrdiff_t(stepj())+j1*diagstep, j2-j1,diagstep,ct()
                    TMV_FIRSTLAST );
            } else {
                TMVAssert(
                    j1>=0 && j1-j2<=0 &&
                    j2<=TMV_MIN(colsize()+i,rowsize()));
                return vec_type(
                    ptr()-i*ptrdiff_t(stepi())+j1*diagstep, j2-j1,diagstep,ct()
                    TMV_FIRSTLAST );
            }
        }

        // Repeat the const versions:
        inline T operator()(ptrdiff_t i,ptrdiff_t j) const
        { return base::operator()(i,j); }
        inline const_vec_type operator[](ptrdiff_t i) const
        { return base::operator[](i); }
        inline const_vec_type row(ptrdiff_t i) const
        { return base::row(i); }
        inline const_vec_type col(ptrdiff_t j) const
        { return base::col(j); }
        inline const_vec_type diag() const
        { return base::diag(); }
        inline const_vec_type diag(ptrdiff_t i) const
        { return base::diag(i); }
        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::row(i,j1,j2); }
        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        { return base::col(j,i1,i2); }
        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::diag(i,j1,j2); }

        //
        // Modifying Functions
        //

        type& setZero();

        type& setAllTo(const T& x);

        type& addToAll(const T& x);

        type& clip(RT thresh);

        type& transposeSelf();

        type& conjugateSelf();

        type& setToIdentity(const T& x=T(1));

        inline type& swapRows(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>=0 && i1 < colsize() &&
                      i2>=0 && i2 < colsize());
            if (i1!=i2) Swap(row(i1),row(i2));
            return *this;
        }

        inline type& swapCols(ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1>=0 && j1 < rowsize() &&
                      j2>=0 && j2 < rowsize());
            if (j1!=j2) Swap(col(j1),col(j2));
            return *this;
        }

        type& permuteRows(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2);

        inline type& permuteRows(const ptrdiff_t* p)
        { return permuteRows(p,0,colsize()); }

        inline type& permuteCols(const ptrdiff_t* p, ptrdiff_t j1, ptrdiff_t j2)
        { transpose().permuteRows(p,j1,j2); return *this; }

        inline type& permuteCols(const ptrdiff_t* p)
        { return permuteCols(p,0,rowsize()); }

        type& reversePermuteRows(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2);

        inline type& reversePermuteRows(const ptrdiff_t* p)
        { return reversePermuteRows(p,0,colsize()); }

        inline type& reversePermuteCols(const ptrdiff_t* p, ptrdiff_t j1, ptrdiff_t j2)
        { transpose().reversePermuteRows(p,j1,j2); return *this; }

        inline type& reversePermuteCols(const ptrdiff_t* p)
        { return reversePermuteCols(p,0,rowsize()); }

        //
        // subMatrix
        //

        inline view_type cSubMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            return type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1,j2-j1,stepi(),stepj(),ct() TMV_FIRSTLAST );
        }

        inline view_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,1,1));
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline view_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            return type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep, (j2-j1)/jstep, istep*stepi(), jstep*stepj(),ct()
                TMV_FIRSTLAST );
        }

        inline view_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(base::hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline vec_type cSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size)
        {
            TMVAssert(size >= 0);
            return vec_type(
                ptr()+i*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),size,
                istep*stepi()+jstep*stepj(),ct()
                TMV_FIRSTLAST );
        }

        inline vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size)
        {
            TMVAssert(size >= 0);
            TMVAssert(base::hasSubVector(i,j,istep,jstep,size));
            return cSubVector(i,j,istep,jstep,size);
        }

        inline uppertri_type unitUpperTri()
        {
            TMVAssert(rowsize() <= colsize());
            return uppertri_type(
                ptr(),rowsize(),stepi(),stepj(),UnitDiag,ct()
                TMV_FIRSTLAST);
        }

        inline uppertri_type upperTri(DiagType dt=NonUnitDiag)
        {
            TMVAssert(rowsize() <= colsize());
            return uppertri_type(
                ptr(),rowsize(),stepi(),stepj(), dt,ct() TMV_FIRSTLAST);
        }

        inline lowertri_type unitLowerTri()
        {
            TMVAssert(colsize() <= rowsize());
            return lowertri_type(
                ptr(),colsize(),stepi(),stepj(),UnitDiag,ct()
                TMV_FIRSTLAST);
        }

        inline lowertri_type lowerTri(DiagType dt=NonUnitDiag)
        {
            TMVAssert(colsize() <= rowsize());
            return lowertri_type(
                ptr(),colsize(),stepi(),stepj(),dt,ct() TMV_FIRSTLAST);
        }

        inline view_type cColPair(ptrdiff_t j1, ptrdiff_t j2)
        {
            return type(
                ptr()+j1*ptrdiff_t(stepj()),colsize(),2,stepi(),(j2-j1)*stepj(),ct()
                TMV_FIRSTLAST );
        }

        inline view_type colPair(ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1>=0 && j1<rowsize() && j2>=0 && j2<rowsize());
            return cColPair(j1,j2);
        }

        inline view_type cRowPair(ptrdiff_t i1, ptrdiff_t i2)
        {
            return type(
                ptr()+i1*ptrdiff_t(stepi()),2,rowsize(),(i2-i1)*stepi(),stepj(),ct()
                TMV_FIRSTLAST );
        }

        inline view_type rowPair(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>=0 && i1<colsize() && i2>=0 && i2<colsize());
            return cRowPair(i1,i2);
        }

        inline view_type cColRange(ptrdiff_t j1, ptrdiff_t j2)
        {
            return type(
                ptr()+j1*ptrdiff_t(stepj()),colsize(),j2-j1,
                stepi(),stepj(),ct(),(this->iscm()&&ls())?-1:0
                TMV_FIRSTLAST);
        }

        inline view_type colRange(ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            return cColRange(j1,j2);
        }

        inline view_type cRowRange(ptrdiff_t i1, ptrdiff_t i2)
        {
            return type(
                ptr()+i1*ptrdiff_t(stepi()),i2-i1,rowsize(),
                stepi(),stepj(),ct(),(this->isrm()&&ls())?-1:0
                TMV_FIRSTLAST);
        }

        inline view_type rowRange(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            return cRowRange(i1,i2);
        }

        inline realpart_type realPart()
        {
            return realpart_type(
                reinterpret_cast<RT*>(ptr()),colsize(),rowsize(),
                isReal(T()) ? stepi() : 2*stepi(),
                isReal(T()) ? stepj() : 2*stepj(), NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)
                ,reinterpret_cast<const RT*>(_last)
#endif
            );
        }

        inline realpart_type imagPart()
        {
            TMVAssert(isComplex(T()));
            return realpart_type(
                reinterpret_cast<RT*>(ptr())+1,
                colsize(),rowsize(),2*stepi(),2*stepj(),NonConj
#ifdef TMVFLDEBUG
                ,reinterpret_cast<const RT*>(_first)+1
                ,reinterpret_cast<const RT*>(_last)+1
#endif
            );
        }

        //
        // Views
        //

        inline view_type view()
        { return *this; }

        inline view_type transpose()
        {
            return type(
                ptr(),rowsize(),colsize(),
                stepj(),stepi(),ct(),ls() TMV_FIRSTLAST);
        }

        inline view_type conjugate()
        {
            return type(
                ptr(),colsize(),rowsize(),
                stepi(),stepj(),TMV_ConjOf(T,ct()),ls() TMV_FIRSTLAST);
        }

        inline view_type adjoint()
        {
            return type(
                ptr(),rowsize(),colsize(),
                stepj(),stepi(),TMV_ConjOf(T,ct()),ls()
                TMV_FIRSTLAST);
        }

        inline vec_type linearView()
        {
            TMVAssert(ls() != -1);
            // To assure that next Assert has no effect

            TMVAssert(canLinearize());
            TMVAssert(ls() == colsize()*rowsize());
            return vec_type(ptr(),ls(),1,ct() TMV_FIRSTLAST );
        }

        // Repeat the const versions
        inline const_view_type cSubMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::cSubMatrix(i1,i2,j1,j2); }
        inline const_view_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        { return base::subMatrix(i1,i2,j1,j2); }
        inline const_view_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        { return base::cSubMatrix(i1,i2,j1,j2,istep,jstep); }
        inline const_view_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        { return base::subMatrix(i1,i2,j1,j2,istep,jstep); }
        inline const_vec_type cSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const
        { return base::cSubVector(i,j,istep,jstep,size); }
        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t size) const
        { return base::subVector(i,j,istep,jstep,size); }
        inline const_uppertri_type unitUpperTri() const
        { return base::unitUpperTri(); }
        inline const_uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        { return base::upperTri(dt); }
        inline const_lowertri_type unitLowerTri() const
        { return base::unitLowerTri(); }
        inline const_lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        { return base::lowerTri(dt); }
        inline const_view_type cColPair(ptrdiff_t j1, ptrdiff_t j2) const
        { return base::cColPair(j1,j2); }
        inline const_view_type colPair(ptrdiff_t j1, ptrdiff_t j2) const
        { return base::colPair(j1,j2); }
        inline const_view_type cRowPair(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::cRowPair(i1,i2); }
        inline const_view_type rowPair(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::rowPair(i1,i2); }
        inline const_view_type cColRange(ptrdiff_t j1, ptrdiff_t j2) const
        { return base::cColRange(j1,j2); }
        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        { return base::colRange(j1,j2); }
        inline const_view_type cRowRange(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::cRowRange(i1,i2); }
        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        { return base::rowRange(i1,i2); }
        inline const_realpart_type realPart() const
        { return base::realPart(); }
        inline const_realpart_type imagPart() const
        { return base::imagPart(); }
        inline const_view_type view() const
        { return base::view(); }
        inline const_view_type transpose() const
        { return base::transpose(); }
        inline const_view_type conjugate() const
        { return base::conjugate(); }
        inline const_view_type adjoint() const
        { return base::adjoint(); }
        inline const_vec_type linearView() const
        { return base::linearView(); }

        //
        // I/O
        //

        void read(const TMV_Reader& reader);

        virtual inline ptrdiff_t colsize() const { return itscs; }
        virtual inline ptrdiff_t rowsize() const { return itsrs; }
        virtual inline const T* cptr() const { return itsm; }
        inline T* ptr() const { return itsm; }
        virtual inline ptrdiff_t stepi() const { return itssi; }
        virtual inline ptrdiff_t stepj() const { return itssj; }
        virtual inline ConjType ct() const { return itsct; }
        virtual inline ptrdiff_t ls() const { return linsize; }

        virtual inline bool canLinearize() const
        {
            if (linsize == -1) {
                if ( (stepi() == 1 && stepj() == colsize()) ||
                     (stepj() == 1 && stepi() == rowsize()) )
                    linsize = rowsize() * colsize();
                else
                    linsize = 0;
            }
            TMVAssert(linsize == 0 ||
                      (this->isrm() && stepi() == rowsize()) ||
                      (this->iscm() && stepj() == colsize()));
            return linsize > 0;
        }

        reference ref(ptrdiff_t i, ptrdiff_t j);

        inline rowmajor_iterator rowmajor_begin()
        { return rowmajor_iterator(this,0,0); }
        inline rowmajor_iterator rowmajor_end()
        { return rowmajor_iterator(this,colsize(),0); }

        inline colmajor_iterator colmajor_begin()
        { return colmajor_iterator(this,0,0); }
        inline colmajor_iterator colmajor_end()
        { return colmajor_iterator(this,0,rowsize()); }

        inline const_rowmajor_iterator rowmajor_begin() const
        { return const_rowmajor_iterator(this,0,0); }
        inline const_rowmajor_iterator rowmajor_end() const
        { return const_rowmajor_iterator(this,colsize(),0); }

        inline const_colmajor_iterator colmajor_begin() const
        { return const_colmajor_iterator(this,0,0); }
        inline const_colmajor_iterator colmajor_end() const
        { return const_colmajor_iterator(this,0,rowsize()); }

    protected:

        T*const itsm;
        const ptrdiff_t itscs;
        const ptrdiff_t itsrs;
        const ptrdiff_t itssi;
        const ptrdiff_t itssj;
        const ConjType itsct;

        mutable ptrdiff_t linsize;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
#endif

    }; // MatrixView

    template <typename T>
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
        typedef ConstMatrixView<T,FortranStyle> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstMatrixView<RT,FortranStyle> const_realpart_type;
        typedef ConstVectorView<T,FortranStyle> const_vec_type;
        typedef ConstUpperTriMatrixView<T,FortranStyle> const_uppertri_type;
        typedef ConstLowerTriMatrixView<T,FortranStyle> const_lowertri_type;
        typedef typename RefHelper<T>::reference reference;

        //
        // Constructors
        //

        inline MatrixView(const type& rhs) : c_type(rhs) {}
        inline MatrixView(const c_type& rhs) : c_type(rhs) {}

        inline MatrixView(
            T* _m, ptrdiff_t _cs, ptrdiff_t _rs, ptrdiff_t _si, ptrdiff_t _sj,
            ConjType _ct, ptrdiff_t _ls=0 TMV_PARAMFIRSTLAST(T) ) :
            c_type(_m,_cs,_rs,_si,_sj,_ct,_ls TMV_FIRSTLAST1(_first,_last) )
        {}

        virtual inline ~MatrixView() {}

        //
        // Op=
        //

        inline type& operator=(const type& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const c_type& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenMatrix<RT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const GenMatrix<CT>& m2)
        { c_type::operator=(m2); return *this; }

        template <typename T2>
        inline type& operator=(const GenMatrix<T2>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const T& x)
        { c_type::operator=(x); return *this; }

        inline type& operator=(const AssignableToMatrix<RT>& m2)
        { c_type::operator=(m2); return *this; }

        inline type& operator=(const AssignableToMatrix<CT>& m2)
        { c_type::operator=(m2); return *this; }

        template <typename T2, ptrdiff_t M, ptrdiff_t N, int A2>
        inline type& operator=(const SmallMatrix<T2,M,N,A2>& m2)
        { c_type::operator=(m2); return *this; }

        typedef typename c_type::MyListAssigner MyListAssigner;
        inline MyListAssigner operator<<(const T& x)
        { return c_type::operator<<(x); }

        //
        // Access
        //

        inline reference operator()(ptrdiff_t i,ptrdiff_t j)
        {
            TMVAssert(i > 0 && i <= this->colsize());
            TMVAssert(j > 0 && j <= this->rowsize());
            return c_type::ref(i-1,j-1);
        }

        inline vec_type operator[](ptrdiff_t i)
        {
            TMVAssert(i>0 && i<=this->colsize());
            return row(i);
        }

        inline vec_type row(ptrdiff_t i)
        {
            TMVAssert(i>0 && i<=this->colsize());
            return c_type::row(i-1);
        }

        inline vec_type col(ptrdiff_t j)
        {
            TMVAssert(j>0 && j<=this->rowsize());
            return c_type::col(j-1);
        }

        inline vec_type diag()
        { return c_type::diag(); }

        inline vec_type diag(ptrdiff_t i)
        { return c_type::diag(i); }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(i>0 && i<=this->colsize());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->rowsize());
            return c_type::row(i-1,j1-1,j2);
        }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(j>0 && j<=this->rowsize());
            TMVAssert(i1>0 && i1-i2<=0 && i2<=this->colsize());
            return c_type::col(j-1,i1-1,i2);
        }

        inline vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1 > 0);
            return c_type::diag(i,j1-1,j2);
        }

        inline T operator()(ptrdiff_t i,ptrdiff_t j) const
        {
            TMVAssert(i > 0 && i <= this->colsize());
            TMVAssert(j > 0 && j <= this->rowsize());
            return c_type::cref(i-1,j-1);
        }

        inline const_vec_type operator[](ptrdiff_t i) const
        {
            TMVAssert(i>0 && i<=this->colsize());
            return row(i);
        }

        inline const_vec_type row(ptrdiff_t i) const
        {
            TMVAssert(i>0 && i<=this->colsize());
            return c_type::row(i-1);
        }

        inline const_vec_type col(ptrdiff_t j) const
        {
            TMVAssert(j>0 && j<=this->rowsize());
            return c_type::col(j-1);
        }

        inline const_vec_type diag() const
        { return c_type::diag(); }

        inline const_vec_type diag(ptrdiff_t i) const
        { return c_type::diag(i); }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(i>0 && i<=this->colsize());
            TMVAssert(j1>0 && j1-j2<=0 && j2<=this->rowsize());
            return c_type::row(i-1,j1-1,j2);
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(j>0 && j<=this->rowsize());
            TMVAssert(i1>0 && i1-i2<=0 && i2<=this->colsize());
            return c_type::col(j-1,i1-1,i2);
        }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1 > 0);
            return c_type::diag(i,j1-1,j2);
        }

        //
        // Modifying Functions
        //

        inline type& setZero()
        { c_type::setZero(); return *this; }

        inline type& setAllTo(const T& x)
        { c_type::setAllTo(x); return *this; }

        inline type& addToAll(const T& x)
        { c_type::addToAll(x); return *this; }

        inline type& clip(RT thresh)
        { c_type::clip(thresh); return *this; }

        inline type& transposeSelf()
        { c_type::transposeSelf(); return *this; }

        inline type& conjugateSelf()
        { c_type::conjugateSelf(); return *this; }

        inline type& setToIdentity(const T& x=T(1))
        { c_type::setToIdentity(x); return *this; }

        inline type& swapRows(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1 > 0 && i1 <= this->colsize());
            TMVAssert(i2 > 0 && i2 <= this->colsize());
            if (i1 != i2)
                c_type::swapRows(i1-1,i2-1);
            return *this;
        }

        inline type& swapCols(ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1 > 0 && j1 <= this->rowsize());
            TMVAssert(j2 > 0 && j2 <= this->rowsize());
            if (j1 != j2)
                c_type::swapCols(j1-1,j2-1);
            return *this;
        }

        inline type& permuteRows(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>0);
            c_type::permuteRows(p,i1-1,i2);
            return *this;
        }

        inline type& permuteRows(const ptrdiff_t* p)
        { c_type::permuteRows(p); return *this; }

        inline type& permuteCols(const ptrdiff_t* p, ptrdiff_t j1, ptrdiff_t j2)
        { transpose().permuteRows(p,j1,j2); return *this; }

        inline type& permuteCols(const ptrdiff_t* p)
        { transpose().permuteRows(p); return *this; }

        inline type& reversePermuteRows(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1>0);
            c_type::reversePermuteRows(p,i1-1,i2);
            return *this;
        }

        inline type& reversePermuteRows(const ptrdiff_t* p)
        { c_type::reversePermuteRows(p); return *this; }

        inline type& reversePermuteCols(const ptrdiff_t* p, ptrdiff_t j1, ptrdiff_t j2)
        { transpose().reversePermuteRows(p,j1,j2); return *this; }

        inline type& reversePermuteCols(const ptrdiff_t* p)
        { transpose().reversePermuteRows(p); return *this; }

        //
        // subMatrix
        //

        inline bool hasSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            return const_type(*this).hasSubMatrix(
                i1,i2,j1,j2,istep,jstep);
        }

        inline bool hasSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        { return const_type(*this).hasSubVector(i,j,istep,jstep,s); }

        inline view_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return c_type::cSubMatrix(i1-1,i2,j1-1,j2);
        }

        inline view_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return c_type::cSubMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s)
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,s));
            return c_type::cSubVector(i-1,j-1,istep,jstep,s);
        }

        inline uppertri_type unitUpperTri()
        { return c_type::upperTri(UnitDiag); }

        inline uppertri_type upperTri(DiagType dt=NonUnitDiag)
        { return c_type::upperTri(dt); }

        inline lowertri_type unitLowerTri()
        { return c_type::lowerTri(UnitDiag); }

        inline lowertri_type lowerTri(DiagType dt=NonUnitDiag)
        { return c_type::lowerTri(dt); }

        inline view_type colPair(ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1 > 0 && j1 <= this->rowsize());
            TMVAssert(j2 > 0 && j2 <= this->rowsize());
            return c_type::cColPair(j1-1,j2-1);
        }

        inline view_type rowPair(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1 > 0 && i1 <= this->rowsize());
            TMVAssert(i2 > 0 && i2 <= this->rowsize());
            return c_type::cRowPair(i1-1,i2-1);
        }

        inline view_type colRange(ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(j1 > 0 && j1 <= j2 && j2 <= this->rowsize());
            return c_type::cColRange(j1-1,j2);
        }

        inline view_type rowRange(ptrdiff_t i1, ptrdiff_t i2)
        {
            TMVAssert(i1 > 0 && i1 <= i2 && i2 <= this->colsize());
            return c_type::cRowRange(i1-1,i2);
        }

        inline realpart_type realPart()
        { return c_type::realPart(); }

        inline realpart_type imagPart()
        { return c_type::imagPart(); }

        inline const_view_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,1,1));
            return c_type::cSubMatrix(i1-1,i2,j1-1,j2);
        }

        inline const_view_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            return c_type::cSubMatrix(
                i1-1,i2-1+istep,j1-1,j2-1+jstep,istep,jstep);
        }

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        {
            TMVAssert(hasSubVector(i,j,istep,jstep,s));
            return c_type::cSubVector(i-1,j-1,istep,jstep,s);
        }

        inline const_uppertri_type unitUpperTri() const
        { return c_type::upperTri(UnitDiag); }

        inline const_uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        { return c_type::upperTri(dt); }

        inline const_lowertri_type unitLowerTri() const
        { return c_type::lowerTri(UnitDiag); }

        inline const_lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        { return c_type::lowerTri(dt); }

        inline const_view_type colPair(ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1 > 0 && j1 <= this->rowsize());
            TMVAssert(j2 > 0 && j2 <= this->rowsize());
            return c_type::cColPair(j1-1,j2-1);
        }

        inline const_view_type rowPair(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(i1 > 0 && i1 <= this->rowsize());
            TMVAssert(i2 > 0 && i2 <= this->rowsize());
            return c_type::cRowPair(i1-1,i2-1);
        }

        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(j1 > 0 && j1 <= j2 && j2 <= this->rowsize());
            return c_type::cColRange(j1-1,j2);
        }

        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            TMVAssert(i1 > 0 && i1 <= i2 && i2 <= this->colsize());
            return c_type::cRowRange(i1-1,i2);
        }

        inline const_realpart_type realPart() const
        { return c_type::realPart(); }

        inline const_realpart_type imagPart() const
        { return c_type::imagPart(); }


        //
        // Views
        //

        inline view_type view()
        { return c_type::view(); }

        inline view_type transpose()
        { return c_type::transpose(); }

        inline view_type conjugate()
        { return c_type::conjugate(); }

        inline view_type adjoint()
        { return c_type::adjoint(); }

        inline vec_type linearView()
        { return c_type::linearView(); }

        inline const_view_type view() const
        { return c_type::view(); }

        inline const_view_type transpose() const
        { return c_type::transpose(); }

        inline const_view_type conjugate() const
        { return c_type::conjugate(); }

        inline const_view_type adjoint() const
        { return c_type::adjoint(); }

        inline const_vec_type linearView() const
        { return c_type::linearView(); }

    }; // FortranStyle MatrixView

    template <ptrdiff_t S, class M>
    struct MatrixIterHelper;

    template <class M>
    struct MatrixIterHelper<RowMajor,M>
    {
        typedef typename M::value_type T;
        typedef VIt<T,1,NonConj> rowmajor_iterator;
        typedef CVIt<T,1,NonConj> const_rowmajor_iterator;

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
        typedef VIt<T,1,NonConj> colmajor_iterator;
        typedef CVIt<T,1,NonConj> const_colmajor_iterator;

        static colmajor_iterator colmajor_begin(M* m)
        { return colmajor_iterator(m->ptr(),1); }
        static colmajor_iterator colmajor_end(M* m)
        { return colmajor_iterator(m->ptr()+m->ls(),1); }

        static const_colmajor_iterator colmajor_begin(const M* m)
        { return const_colmajor_iterator(m->cptr(),1); }
        static const_colmajor_iterator colmajor_end(const M* m)
        { return const_colmajor_iterator(m->cptr()+m->ls(),1); }

    };

    template <typename T, int A>
    class Matrix :
        public GenMatrix<T>
    {
    public:

        enum { S = A & AllStorageType };
        enum { I = A & FortranStyle };

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef Matrix<T,A> type;
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
            TMVAssert(Attrib<A>::matrixok);
        }

        inline Matrix(ptrdiff_t _colsize, ptrdiff_t _rowsize) :
            NEW_SIZE(_colsize,_rowsize)
        {
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(_colsize >= 0 && _rowsize >= 0);
#ifdef TMV_EXTRA_DEBUG
            setAllTo(T(888));
#endif
        }

        inline Matrix(ptrdiff_t _colsize, ptrdiff_t _rowsize, const T& x) :
            NEW_SIZE(_colsize,_rowsize)
        {
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(_colsize >= 0 && _rowsize >= 0);
            setAllTo(x);
        }

        inline Matrix(const type& rhs) : NEW_SIZE(rhs.colsize(),rhs.rowsize())
        {
            TMVAssert(Attrib<A>::matrixok);
            std::copy(rhs.cptr(),rhs.cptr()+linsize,itsm.get());
        }

        inline Matrix(const GenMatrix<RT>& rhs) :
            NEW_SIZE(rhs.colsize(),rhs.rowsize())
        {
            TMVAssert(Attrib<A>::matrixok);
            rhs.assignToM(view());
        }

        inline Matrix(const GenMatrix<CT>& rhs) :
            NEW_SIZE(rhs.colsize(),rhs.rowsize())
        {
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(isComplex(T()));
            rhs.assignToM(view());
        }

        template <typename T2>
        inline Matrix(const GenMatrix<T2>& rhs) :
            NEW_SIZE(rhs.colsize(),rhs.rowsize())
        {
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(isComplex(T()) || isReal(T2()));
            Copy(rhs,view());
        }

        inline Matrix(const AssignableToMatrix<RT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize())
        {
            TMVAssert(Attrib<A>::matrixok);
            m2.assignToM(view());
        }

        inline Matrix(const AssignableToMatrix<CT>& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize())
        {
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(isComplex(T()));
            m2.assignToM(view());
        }

        inline Matrix(const Permutation& m2) :
            NEW_SIZE(m2.colsize(),m2.rowsize())
        {
            TMVAssert(Attrib<A>::matrixok);
            m2.assignToM(view());
        }

        template <typename T2, ptrdiff_t M, ptrdiff_t N, int A2>
        inline Matrix(const SmallMatrix<T2,M,N,A2>& rhs) :
            NEW_SIZE(rhs.colsize(),rhs.rowsize())
        {
            TMVAssert(Attrib<A>::matrixok);
            TMVAssert(isComplex(T()) || isReal(T2()));
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

        inline type& operator=(const Matrix<T,A>& m2)
        {
            TMVAssert(m2.colsize() == colsize() && m2.rowsize() == rowsize());
            if (&m2 != this)
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

        template <typename T2>
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

        template <typename T2, ptrdiff_t M, ptrdiff_t N, int A2>
        inline type& operator=(const SmallMatrix<T2,M,N,A2>& m2)
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

        inline T operator()(ptrdiff_t i,ptrdiff_t j) const
        {
            if (I == int(CStyle)) {
                TMVAssert(i>=0 && i < colsize());
                TMVAssert(j>=0 && j < rowsize());
                return cref(i,j);
            } else {
                TMVAssert(i > 0 && i <= colsize());
                TMVAssert(j > 0 && j <= rowsize());
                return cref(i-1,j-1);
            }
        }

        inline T& operator()(ptrdiff_t i,ptrdiff_t j)
        {
            if (I == int(CStyle)) {
                TMVAssert(i>=0 && i < colsize());
                TMVAssert(j>=0 && j < rowsize());
                return ref(i,j);
            } else {
                TMVAssert(i > 0 && i <= colsize());
                TMVAssert(j > 0 && j <= rowsize());
                return ref(i-1,j-1);
            }
        }

        inline const_vec_type row(ptrdiff_t i) const
        {
            if (I == int(FortranStyle)) { TMVAssert(i>0 && i<=colsize()); --i; }
            else TMVAssert(i>=0 && i<colsize());
            return const_vec_type(
                itsm.get()+i*ptrdiff_t(stepi()),rowsize(),stepj(),NonConj);
        }

        inline const_vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            if (I == int(FortranStyle)) {
                TMVAssert(i>0 && i<=colsize()); --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize()); --j1;
            } else {
                TMVAssert(i>=0 && i<colsize());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }
            return const_vec_type(
                itsm.get()+i*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),j2-j1,stepj(),NonConj);
        }

        inline const_vec_type operator[](ptrdiff_t i) const
        { return row(i); }

        inline const_vec_type col(ptrdiff_t j) const
        {
            if (I == int(FortranStyle)) { TMVAssert(j>0 && j<=rowsize()); --j; }
            else TMVAssert(j>=0 && j<rowsize());
            return const_vec_type(
                itsm.get()+j*ptrdiff_t(stepj()),colsize(),stepi(),NonConj);
        }

        inline const_vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2) const
        {
            if (I == int(FortranStyle)) {
                TMVAssert(j>0 && j<=rowsize()); --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=colsize()); --i1;
            } else {
                TMVAssert(j>=0 && j<rowsize());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            }
            return const_vec_type(
                itsm.get()+i1*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),i2-i1,stepi(),NonConj);
        }

        inline const_vec_type diag() const
        {
            return const_vec_type(
                itsm.get(),TMV_MIN(colsize(),rowsize()),
                stepi()+stepj(),NonConj);
        }

        inline const_vec_type diag(ptrdiff_t i) const
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const ptrdiff_t diagstep = stepi() + stepj();
            if (i >= 0) {
                const ptrdiff_t diagsize = TMV_MIN(colsize(),rowsize()-i);
                return const_vec_type(
                    itsm.get()+i*ptrdiff_t(stepj()),diagsize,diagstep,NonConj);
            } else {
                const ptrdiff_t diagsize = TMV_MIN(colsize()+i,rowsize());
                return const_vec_type(
                    itsm.get()-i*ptrdiff_t(stepi()),diagsize,diagstep,NonConj);
            }
        }

        inline const_vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2) const
        {
            if (I == int(FortranStyle)) {
                TMVAssert(j1 > 0 && j1-j2<=0);
                --j1;
            } else {
                TMVAssert( j1>=0 && j1-j2<=0);
            }
            TMVAssert(i>=-colsize() && i<=rowsize());
            const ptrdiff_t diagstep = stepi() + 1;
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(colsize(),rowsize()-i));
                return const_vec_type(
                    itsm.get()+i*ptrdiff_t(stepj())+j1*diagstep,j2-j1,diagstep,NonConj);
            } else {
                TMVAssert(j2<=TMV_MIN(colsize()+i,rowsize()));
                return const_vec_type(
                    itsm.get()-i*ptrdiff_t(stepi())+j1*diagstep,j2-j1,diagstep,NonConj);
            }
        }

        inline vec_type row(ptrdiff_t i)
        {
            if (I == int(FortranStyle)) {
                TMVAssert(i>0 && i<=colsize());
                --i;
            } else {
                TMVAssert(i>=0 && i<colsize());
            }
            return vec_type(
                ptr()+i*ptrdiff_t(stepi()),rowsize(),stepj(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type row(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I == int(FortranStyle)) {
                TMVAssert(i>0 && i<=colsize());
                --i;
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize());
                --j1;
            } else {
                TMVAssert(i>=0 && i<colsize());
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }
            return vec_type(
                ptr()+i*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                j2-j1,stepj(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type operator[](ptrdiff_t i)
        { return row(i); }

        inline vec_type col(ptrdiff_t j)
        {
            if (I == int(FortranStyle)) {
                TMVAssert(j>0 && j<=rowsize());
                --j;
            } else {
                TMVAssert(j>=0 && j<rowsize());
            }
            return vec_type(
                ptr()+j*ptrdiff_t(stepj()),colsize(),stepi(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type col(ptrdiff_t j, ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I == int(FortranStyle)) {
                TMVAssert(j>0 && j<=rowsize()); --j;
                TMVAssert(i1>0 && i1-i2<=0 && i2<=colsize()); --i1;
            } else {
                TMVAssert(j>=0 && j<rowsize());
                TMVAssert(i1>=0 && i1-i2<=0 && i2<=colsize());
            }
            return vec_type(
                ptr()+i1*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),
                i2-i1,stepi(),NonConj TMV_FIRSTLAST);
        }

        inline vec_type diag()
        {
            return vec_type(
                ptr(),TMV_MIN(colsize(),rowsize()),stepi()+stepj(),NonConj
                TMV_FIRSTLAST);
        }

        inline vec_type diag(ptrdiff_t i)
        {
            TMVAssert(i>=-colsize() && i<=rowsize());
            const ptrdiff_t diagstep = stepi() + stepj();
            if (i >= 0) {
                const ptrdiff_t diagsize = TMV_MIN(colsize(),rowsize()-i);
                return vec_type(
                    ptr()+i*ptrdiff_t(stepj()), diagsize,diagstep,NonConj
                    TMV_FIRSTLAST);
            } else {
                const ptrdiff_t diagsize = TMV_MIN(colsize()+i,rowsize());
                return vec_type(
                    ptr()-i*ptrdiff_t(stepi()), diagsize,diagstep,NonConj
                    TMV_FIRSTLAST);
            }
        }

        inline vec_type diag(ptrdiff_t i, ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I == int(FortranStyle)) { TMVAssert(j1 > 0 && j1-j2<=0); --j1; }
            else { TMVAssert(j1>=0 && j1-j2<=0); }
            TMVAssert(i>=-colsize() && i<=rowsize());
            const ptrdiff_t diagstep = stepi() + stepj();
            if (i >= 0) {
                TMVAssert(j2<=TMV_MIN(colsize(),rowsize()-i));
                return vec_type(
                    ptr()+i*ptrdiff_t(stepj()) + j1*diagstep, j2-j1, diagstep, NonConj
                    TMV_FIRSTLAST);
            } else {
                TMVAssert(j2<=TMV_MIN(colsize(),rowsize()-i));
                return vec_type(
                    ptr()-i*ptrdiff_t(stepi()) + j1*diagstep, j2-j1, diagstep, NonConj
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

        inline type& swapRows(ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I == int(CStyle)) {
                TMVAssert(i1>=0 && i1<colsize());
                TMVAssert(i2>=0 && i2<colsize());
            } else {
                TMVAssert(i1>0 && i1<=colsize());
                TMVAssert(i2>0 && i2<=colsize());
            }
            if (i1!=i2) Swap(row(i1),row(i2));
            return *this;
        }

        inline type& swapCols(ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I == int(CStyle)) {
                TMVAssert(j1>=0 && j1<rowsize());
                TMVAssert(j2>=0 && j2<rowsize());
            } else {
                TMVAssert(j1>0 && j1<=rowsize());
                TMVAssert(j2>0 && j2<=rowsize());
            }
            if (j1!=j2) Swap(col(j1),col(j2));
            return *this;
        }

        inline type& permuteRows(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        { view().permuteRows(p,i1,i2); return *this; }

        inline type& permuteRows(const ptrdiff_t* p)
        { view().permuteRows(p); return *this; }

        inline type& permuteCols(const ptrdiff_t* p, ptrdiff_t j1, ptrdiff_t j2)
        { view().permuteCols(p,j1,j2); return *this; }

        inline type& permuteCols(const ptrdiff_t* p)
        { view().permuteCols(p); return *this; }

        inline type& reversePermuteRows(const ptrdiff_t* p, ptrdiff_t i1, ptrdiff_t i2)
        { view().reversePermuteRows(p,i1,i2); return *this; }

        inline type& reversePermuteRows(const ptrdiff_t* p)
        { view().reversePermuteRows(p); return *this; }

        inline type& reversePermuteCols(const ptrdiff_t* p, ptrdiff_t j1, ptrdiff_t j2)
        { view().reversePermuteCols(p,j1,j2); return *this; }

        inline type& reversePermuteCols(const ptrdiff_t* p)
        { view().reversePermuteCols(p); return *this; }

        //
        // subMatrix
        //

        inline const_view_type cSubMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            return const_view_type(
                itsm.get()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1,j2-j1,stepi(),stepj(),NonConj);
        }

        inline const_view_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I == int(FortranStyle)) { --i1; --j1; }
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline const_view_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            return const_view_type(
                itsm.get()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep,(j2-j1)/jstep,
                istep*stepi(),jstep*stepj(),NonConj);
        }

        inline const_view_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep) const
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I == int(FortranStyle)) {
                --i1; --j1;
                i2 += istep-1; j2 += jstep-1;
            }
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline const_vec_type cSubVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        {
            return const_vec_type(
                itsm.get()+i*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),s,
                istep*stepi()+jstep*stepj(),NonConj);
        }

        inline const_vec_type subVector(
            ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s) const
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,s));
            if (I==int(FortranStyle)) { --i; --j; }
            return cSubVector(i,j,istep,jstep,s);
        }

        inline const_uppertri_type unitUpperTri() const
        {
            TMVAssert(rowsize() <= colsize());
            return const_uppertri_type(
                cptr(),rowsize(), stepi(),stepj(),UnitDiag,ct());
        }

        inline const_uppertri_type upperTri(DiagType dt=NonUnitDiag) const
        {
            TMVAssert(rowsize() <= colsize());
            return const_uppertri_type(
                cptr(),rowsize(), stepi(),stepj(),dt,ct());
        }

        inline const_lowertri_type unitLowerTri() const
        {
            TMVAssert(colsize() <= rowsize());
            return const_lowertri_type(
                cptr(),colsize(), stepi(),stepj(),UnitDiag,ct());
        }

        inline const_lowertri_type lowerTri(DiagType dt=NonUnitDiag) const
        {
            TMVAssert(colsize() <= rowsize());
            return const_lowertri_type(
                cptr(),colsize(), stepi(),stepj(),dt,ct());
        }

        inline const_view_type cColPair(ptrdiff_t j1, ptrdiff_t j2) const
        {
            return const_view_type(
                itsm.get()+j1*ptrdiff_t(stepj()),colsize(),2,
                stepi(),(j2-j1)*stepj(),NonConj);
        }

        inline const_view_type colPair(ptrdiff_t j1, ptrdiff_t j2) const
        {
            if (I == int(CStyle)) {
                TMVAssert(j1>=0 && j1<rowsize() &&
                          j2>=0 && j2<rowsize());
            } else  {
                TMVAssert(j1>0 && j1<=rowsize() &&
                          j2>0 && j2<=rowsize());
                --j1; --j2;
            }
            return cColPair(j1,j2);
        }

        inline const_view_type cRowPair(ptrdiff_t i1, ptrdiff_t i2) const
        {
            return const_view_type(
                itsm.get()+i1*ptrdiff_t(stepi()),2,rowsize(),
                (i2-i1)*stepi(),stepj(),NonConj);
        }

        inline const_view_type rowPair(ptrdiff_t i1, ptrdiff_t i2) const
        {
            if (I == int(CStyle))  {
                TMVAssert(i1>=0 && i1<colsize() &&
                          i2>=0 && i2<colsize());
            } else  {
                TMVAssert(i1>0 && i1<=colsize() &&
                          i2>0 && i2<=colsize());
                --i1; --i2;
            }
            return cRowPair(i1,i2);
        }

        inline const_view_type cColRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            return const_view_type(
                itsm.get()+j1*ptrdiff_t(stepj()),colsize(),j2-j1,
                stepi(),stepj(),NonConj,iscm()?-1:0);
        }

        inline const_view_type colRange(ptrdiff_t j1, ptrdiff_t j2) const
        {
            if (I==int(FortranStyle)) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize());
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }
            return cColRange(j1,j2);
        }

        inline const_view_type cRowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            return const_view_type(
                itsm.get()+i1*ptrdiff_t(stepi()),i2-i1,rowsize(),
                stepi(),stepj(),NonConj,isrm()?-1:0);
        }

        inline const_view_type rowRange(ptrdiff_t i1, ptrdiff_t i2) const
        {
            if (I==int(FortranStyle)) {
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
                isReal(T()) ? stepj() : 2*stepj(),NonConj);
        }

        inline const_realpart_type imagPart() const
        {
            TMVAssert(isComplex(T()));
            return const_realpart_type(
                reinterpret_cast<const RT*>(itsm.get())+1,
                colsize(),rowsize(),2*stepi(),2*stepj(),NonConj);
        }

        inline view_type cSubMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            return view_type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                i2-i1,j2-j1,stepi(),stepj(),NonConj
                TMV_FIRSTLAST);
        }

        inline view_type subMatrix(ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,1,1));
            if (I==int(FortranStyle)) { --i1; --j1; }
            return cSubMatrix(i1,i2,j1,j2);
        }

        inline view_type cSubMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            return view_type(
                ptr()+i1*ptrdiff_t(stepi())+j1*ptrdiff_t(stepj()),
                (i2-i1)/istep,(j2-j1)/jstep,
                istep*stepi(),jstep*stepj(), NonConj TMV_FIRSTLAST);
        }

        inline view_type subMatrix(
            ptrdiff_t i1, ptrdiff_t i2, ptrdiff_t j1, ptrdiff_t j2, ptrdiff_t istep, ptrdiff_t jstep)
        {
            TMVAssert(view().hasSubMatrix(i1,i2,j1,j2,istep,jstep));
            if (I == int(FortranStyle)) {
                --i1; --j1;
                i2 += istep-1; j2 += jstep-1;
            }
            return cSubMatrix(i1,i2,j1,j2,istep,jstep);
        }

        inline vec_type cSubVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s)
        {
            return vec_type(
                ptr()+i*ptrdiff_t(stepi())+j*ptrdiff_t(stepj()),s,
                istep*stepi()+jstep*stepj(),NonConj
                TMV_FIRSTLAST);
        }

        inline vec_type subVector(ptrdiff_t i, ptrdiff_t j, ptrdiff_t istep, ptrdiff_t jstep, ptrdiff_t s)
        {
            TMVAssert(view().hasSubVector(i,j,istep,jstep,s));
            if (I == int(FortranStyle)) { --i; --j; }
            return cSubVector(i,j,istep,jstep,s);
        }

        inline uppertri_type unitUpperTri()
        {
            TMVAssert(rowsize() <= colsize());
            return uppertri_type(
                ptr(),rowsize(), stepi(),stepj(),UnitDiag,ct()
                TMV_FIRSTLAST);
        }

        inline uppertri_type upperTri(DiagType dt=NonUnitDiag)
        {
            TMVAssert(rowsize() <= colsize());
            return uppertri_type(
                ptr(),rowsize(), stepi(),stepj(),dt,ct()
                TMV_FIRSTLAST);
        }

        inline lowertri_type unitLowerTri()
        {
            TMVAssert(colsize() <= rowsize());
            return lowertri_type(
                ptr(),colsize(), stepi(),stepj(),UnitDiag,ct()
                TMV_FIRSTLAST);
        }

        inline lowertri_type lowerTri(DiagType dt=NonUnitDiag)
        {
            TMVAssert(colsize() <= rowsize());
            return lowertri_type(
                ptr(),colsize(), stepi(),stepj(),dt,ct()
                TMV_FIRSTLAST);
        }

        inline view_type cColPair(ptrdiff_t j1, ptrdiff_t j2)
        {
            return view_type(
                ptr()+j1*ptrdiff_t(stepj()),colsize(),2,stepi(),(j2-j1)*stepj(),NonConj
                TMV_FIRSTLAST);
        }

        inline view_type colPair(ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I == int(CStyle))
                TMVAssert(j1>=0 && j1<rowsize() &&
                          j2>=0 && j2<rowsize());
            else {
                TMVAssert(j1>0 && j1<=rowsize() &&
                          j2>0 && j2<=rowsize());
                --j1; --j2;
            }
            return cColPair(j1,j2);
        }

        inline view_type cRowPair(ptrdiff_t i1, ptrdiff_t i2)
        {
            return view_type(
                ptr()+i1*ptrdiff_t(stepi()),2,rowsize(),(i2-i1)*stepi(),stepj(),NonConj
                TMV_FIRSTLAST);
        }

        inline view_type rowPair(ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I == int(CStyle))
                TMVAssert(i1>=0 && i1<colsize() &&
                          i2>=0 && i2<colsize());
            else {
                TMVAssert(i1>0 && i1<=colsize() &&
                          i2>0 && i2<=colsize());
                --i1; --i2;
            }
            return cRowPair(i1,i2);
        }

        inline view_type cColRange(ptrdiff_t j1, ptrdiff_t j2)
        {
            return view_type(
                ptr()+j1*ptrdiff_t(stepj()),colsize(),j2-j1,
                stepi(),stepj(),NonConj,iscm()?-1:0 TMV_FIRSTLAST);
        }

        inline view_type colRange(ptrdiff_t j1, ptrdiff_t j2)
        {
            if (I==int(FortranStyle)) {
                TMVAssert(j1>0 && j1-j2<=0 && j2<=rowsize());
                --j1;
            } else {
                TMVAssert(j1>=0 && j1-j2<=0 && j2<=rowsize());
            }
            return cColRange(j1,j2);
        }

        inline view_type cRowRange(ptrdiff_t i1, ptrdiff_t i2)
        {
            return view_type(
                ptr()+i1*ptrdiff_t(stepi()),i2-i1,rowsize(),
                stepi(),stepj(),NonConj,isrm()?-1:0 TMV_FIRSTLAST);
        }

        inline view_type rowRange(ptrdiff_t i1, ptrdiff_t i2)
        {
            if (I==int(FortranStyle)) {
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
                isReal(T()) ? stepj() : 2*stepj(), NonConj
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
                colsize(),rowsize(),2*stepi(),2*stepj(),NonConj
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
                stepi(),stepj(),NonConj,linsize);
        }

        inline const_view_type transpose() const
        {
            return const_view_type(
                itsm.get(),rowsize(),colsize(),
                stepj(),stepi(),NonConj,linsize);
        }

        inline const_view_type conjugate() const
        {
            return const_view_type(
                itsm.get(),colsize(),rowsize(),
                stepi(),stepj(),TMV_ConjOf(T,NonConj),linsize);
        }

        inline const_view_type adjoint() const
        {
            return const_view_type(
                itsm.get(),rowsize(),colsize(),
                stepj(),stepi(),TMV_ConjOf(T,NonConj),linsize);
        }

        inline const_vec_type constLinearView() const
        {
            TMVAssert(linsize == colsize()*rowsize());
            return const_vec_type(itsm.get(),linsize,1,NonConj);
        }

        inline view_type view()
        {
            return view_type(
                ptr(),colsize(),rowsize(), stepi(),stepj(),NonConj,linsize
                TMV_FIRSTLAST);
        }

        inline view_type transpose()
        {
            return view_type(
                ptr(),rowsize(),colsize(), stepj(),stepi(),NonConj,linsize
                TMV_FIRSTLAST);
        }

        inline view_type conjugate()
        {
            return view_type(
                ptr(),colsize(),rowsize(),
                stepi(),stepj(),TMV_ConjOf(T,NonConj),linsize TMV_FIRSTLAST);
        }

        inline view_type adjoint()
        {
            return view_type(
                ptr(),rowsize(),colsize(),
                stepj(),stepi(),TMV_ConjOf(T,NonConj),linsize TMV_FIRSTLAST);
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

        virtual inline ptrdiff_t ls() const { return linsize; }
        virtual inline ptrdiff_t colsize() const { return itscs; }
        virtual inline ptrdiff_t rowsize() const { return itsrs; }
        virtual inline const T* cptr() const { return itsm.get(); }
        inline T* ptr() { return itsm.get(); }
        virtual inline ptrdiff_t stepi() const
        { return S == int(RowMajor) ? itsrs : 1; }
        virtual inline ptrdiff_t stepj() const
        { return S == int(RowMajor) ? 1 : itscs; }
        inline bool isrm() const { return S==int(RowMajor); }
        inline bool iscm() const { return S==int(ColMajor); }
        inline bool isconj() const { return false; }
        virtual inline ConjType ct() const { return NonConj; }

        virtual inline bool canLinearize() const { return true; }

        virtual inline T cref(ptrdiff_t i, ptrdiff_t j) const
        { return itsm.get()[S==int(RowMajor) ? i*ptrdiff_t(stepi())+j : i+j*ptrdiff_t(stepj())]; }

        inline T& ref(ptrdiff_t i, ptrdiff_t j)
        { return itsm.get()[S==int(RowMajor) ? i*ptrdiff_t(stepi())+j : i+j*ptrdiff_t(stepj())]; }

        inline void resize(ptrdiff_t cs, ptrdiff_t rs)
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

        ptrdiff_t linsize;
        AlignedArray<T> itsm;
        ptrdiff_t itscs;
        ptrdiff_t itsrs;

#ifdef TMVFLDEBUG
    public:
        const T* _first;
        const T* _last;
    protected:
#endif

        // If two matrices are the same size and storage, then
        // swap can be much faster by copying the pointers to the data.
        friend void Swap(Matrix<T,A>& m1, Matrix<T,A>& m2)
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

    template <typename T>
    inline ConstMatrixView<T> RowVectorViewOf(const GenVector<T>& v)
    {
        return ConstMatrixView<T>(
            v.cptr(),1,v.size(),v.size(),v.step(),
            v.ct(),v.step()==1?v.size():0);
    }

    template <typename T, int A>
    inline ConstMatrixView<T,A> RowVectorViewOf(const ConstVectorView<T,A>& v)
    {
        return ConstMatrixView<T,A>(
            v.cptr(),1,v.size(),v.size(),v.step(),
            v.ct(),v.step()==1?v.size():0);
    }

    template <typename T, int A>
    inline ConstMatrixView<T,A> RowVectorViewOf(const Vector<T,A>& v)
    {
        return ConstMatrixView<T,A>(
            v.cptr(),1,v.size(),v.size(),1,v.ct(),v.size());
    }

    template <typename T, int A>
    inline MatrixView<T,A> RowVectorViewOf(VectorView<T,A> v)
    {
        return MatrixView<T,A>(
            v.ptr(),1,v.size(),v.size(),v.step(),
            v.ct(),v.step()==1?v.size():0
            TMV_FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

    template <typename T, int A>
    inline MatrixView<T,A> RowVectorViewOf(Vector<T,A>& v)
    {
        return MatrixView<T,A>(
            v.ptr(),1,v.size(),v.size(),1,v.ct(),v.size()
            TMV_FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

    template <typename T>
    inline ConstMatrixView<T> ColVectorViewOf(const GenVector<T>& v)
    {
        return ConstMatrixView<T>(
            v.cptr(),v.size(),1,v.step(),v.size(),
            v.ct(),v.step()==1?v.size():0);
    }

    template <typename T, int A>
    inline ConstMatrixView<T,A> ColVectorViewOf(const ConstVectorView<T,A>& v)
    {
        return ConstMatrixView<T,A>(
            v.cptr(),v.size(),1,v.step(),v.size(),
            v.ct(),v.step()==1?v.size():0);
    }

    template <typename T, int A>
    inline ConstMatrixView<T,A> ColVectorViewOf(const Vector<T,A>& v)
    {
        return ConstMatrixView<T,A>(
            v.cptr(),v.size(),1,1,v.size(),v.ct(),v.size());
    }

    template <typename T, int A>
    inline MatrixView<T,A> ColVectorViewOf(VectorView<T,A> v)
    {
        return MatrixView<T,A>(
            v.ptr(),v.size(),1,v.step(),v.size(),
            v.ct(),v.step()==1?v.size():0
            TMV_FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

    template <typename T, int A>
    inline MatrixView<T,A> ColVectorViewOf(Vector<T,A>& v)
    {
        return MatrixView<T,A>(
            v.ptr(),v.size(),1,1,v.size(),v.ct(),v.size()
            TMV_FIRSTLAST1(v.ptr(),v.ptr()+v.size()) );
    }

    template <typename T>
    inline MatrixView<T> MatrixViewOf(
        T* m, ptrdiff_t colsize, ptrdiff_t rowsize, StorageType stor)
    {
        TMVAssert(colsize >= 0 && rowsize >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        const ptrdiff_t linsize = colsize * rowsize;
        const ptrdiff_t stepi = stor == RowMajor ? rowsize : 1;
        const ptrdiff_t stepj = stor == RowMajor ? 1 : colsize;
        return MatrixView<T>(
            m,colsize,rowsize,stepi,stepj,NonConj,linsize
            TMV_FIRSTLAST1(m,m+linsize));
    }

    template <typename T>
    inline ConstMatrixView<T> MatrixViewOf(
        const T* m, ptrdiff_t colsize, ptrdiff_t rowsize, StorageType stor)
    {
        TMVAssert(colsize >= 0 && rowsize >= 0);
        TMVAssert(stor == RowMajor || stor == ColMajor);
        const ptrdiff_t linsize = colsize*rowsize;
        const ptrdiff_t stepi = stor == RowMajor ? rowsize : 1;
        const ptrdiff_t stepj = stor == RowMajor ? 1 : colsize;
        return ConstMatrixView<T>(
            m,colsize,rowsize,stepi,stepj,NonConj,linsize);
    }

    template <typename T>
    inline MatrixView<T> MatrixViewOf(
        T* m, ptrdiff_t colsize, ptrdiff_t rowsize, ptrdiff_t stepi, ptrdiff_t stepj)
    {
        TMVAssert(colsize >= 0 && rowsize >= 0);
        const ptrdiff_t linsize = (
            (stepi==1 && stepj==colsize) ? colsize * rowsize :
            (stepj==1 && stepi==rowsize) ? colsize * rowsize :
            0 );
        return MatrixView<T>(
            m,colsize,rowsize,stepi,stepj,NonConj,linsize
            TMV_FIRSTLAST1(m,m+stepi*(colsize-1)+stepj*(rowsize-1)+1));
    }

    template <typename T>
    inline ConstMatrixView<T> MatrixViewOf(
        const T* m, ptrdiff_t colsize, ptrdiff_t rowsize, ptrdiff_t stepi, ptrdiff_t stepj)
    {
        TMVAssert(colsize >= 0 && rowsize >= 0);
        const ptrdiff_t linsize = (
            (stepi==1 && stepj==colsize) ? colsize * rowsize :
            (stepj==1 && stepi==rowsize) ? colsize * rowsize :
            0 );
        return ConstMatrixView<T>(
            m,colsize,rowsize,stepi,stepj,NonConj,linsize);
    }



    //
    // Copy Matrices
    //

    template <typename T>
    void DoCopySameType(const GenMatrix<T>& m1, MatrixView<T> m2);

    template <typename T>
    inline void DoCopy(const GenMatrix<T>& m1, MatrixView<T> m2)
    { DoCopySameType(m1,m2); }

    template <typename T, typename T1>
    inline void DoCopyDiffType(const GenMatrix<T1>& m1, MatrixView<T> m2)
    {
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m2.rowsize() > 0);
        TMVAssert(m2.colsize() > 0);
        TMVAssert(m1.ct()==NonConj);
        TMVAssert(m2.ct()==NonConj);
        TMVAssert(!m2.isSameAs(m1));
        TMVAssert(isComplex(T()) || isReal(T1()));

        if (m1.iscm() && m2.iscm() && m1.colsize() > 1)  {
            for(ptrdiff_t j=0;j<m2.rowsize();++j)
                DoCopyDiffType(m1.col(j),m2.col(j));
        } else if (m2.colsize() < m2.rowsize()) {
            if (shouldReverse(m1.stepj(),m2.stepj())) {
                for(ptrdiff_t i=0;i<m2.colsize();++i)
                    DoCopyDiffType(m1.row(i).reverse(),m2.row(i).reverse());
            } else {
                for(ptrdiff_t i=0;i<m2.colsize();++i)
                    DoCopyDiffType(m1.row(i),m2.row(i));
            }
        } else {
            if (shouldReverse(m1.stepi(),m2.stepi())) {
                for(ptrdiff_t j=0;j<m2.rowsize();++j)
                    DoCopyDiffType(m1.col(j).reverse(),m2.col(j).reverse());
            } else {
                for(ptrdiff_t j=0;j<m2.rowsize();++j)
                    DoCopyDiffType(m1.col(j),m2.col(j));
            }
        }
    }

    template <typename T, typename T1>
    inline void DoCopy(const GenMatrix<T1>& m1, MatrixView<T> m2)
    { DoCopyDiffType(m1,m2); }

    template <typename T>
    inline void DoCopy(const GenMatrix<std::complex<T> >&, MatrixView<T>)
    { TMVAssert(TMV_FALSE); }

    template <typename T1, typename T2>
    inline void nonconjCopy(const GenMatrix<T1>& m1, MatrixView<T2> m2)
    {
        TMVAssert(isComplex(T2()) || isReal(T1()));
        TMVAssert(m2.rowsize() == m1.rowsize());
        TMVAssert(m2.colsize() == m1.colsize());
        TMVAssert(m1.ct() == NonConj);
        TMVAssert(m2.ct() == NonConj);

        if (!m2.iscm() && (m2.isrm() || m1.isrm()))
            DoCopy(m1.transpose(),m2.transpose());
        else
            DoCopy(m1,m2);
    }

    template <typename T1, typename T2>
    inline void Copy(const GenMatrix<T1>& m1, MatrixView<T2> m2)
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
            } else if (m1.canLinearize() && m2.canLinearize() &&
                       m1.stepi() == m2.stepi() && m1.stepj() == m2.stepj()) {
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

    template <typename T>
    void Swap(MatrixView<T> m1, MatrixView<T> m2);

    template <typename T, int A>
    inline void Swap(MatrixView<T> m1, Matrix<T,A>& m2)
    { Swap(m1,m2.view()); }

    template <typename T, int A>
    inline void Swap(Matrix<T,A>& m1, MatrixView<T> m2)
    { Swap(m1.view(),m2); }

    template <typename T, int A1, int A2>
    inline void Swap(Matrix<T,A1>& m1, Matrix<T,A2>& m2)
    { Swap(m1.view(),m2.view()); }


    //
    // Views of a Matrix:
    //

    template <typename T>
    inline ConstMatrixView<T> Transpose(const GenMatrix<T>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline ConstMatrixView<T,A> Transpose(const ConstMatrixView<T,A>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline ConstMatrixView<T,A&FortranStyle> Transpose(const Matrix<T,A>& m)
    { return m.transpose(); }

    template <typename T, int A>
    inline MatrixView<T,A> Transpose(MatrixView<T,A> m)
    { return m.transpose(); }

    template <typename T, int A>
    inline MatrixView<T,A&FortranStyle> Transpose(Matrix<T,A>& m)
    { return m.transpose(); }

    template <typename T>
    inline ConstMatrixView<T> Conjugate(const GenMatrix<T>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline ConstMatrixView<T,A> Conjugate(const ConstMatrixView<T,A>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline ConstMatrixView<T,A&FortranStyle> Conjugate(const Matrix<T,A>& m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline MatrixView<T,A> Conjugate(MatrixView<T,A> m)
    { return m.conjugate(); }

    template <typename T, int A>
    inline MatrixView<T,A&FortranStyle> Conjugate(Matrix<T,A>& m)
    { return m.conjugate(); }

    template <typename T>
    inline ConstMatrixView<T> Adjoint(const GenMatrix<T>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline ConstMatrixView<T,A> Adjoint(const ConstMatrixView<T,A>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline ConstMatrixView<T,A&FortranStyle> Adjoint(const Matrix<T,A>& m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline MatrixView<T,A> Adjoint(MatrixView<T,A> m)
    { return m.adjoint(); }

    template <typename T, int A>
    inline MatrixView<T,A&FortranStyle> Adjoint(Matrix<T,A>& m)
    { return m.adjoint(); }

    template <typename T>
    inline QuotXM<T,T> Inverse(const GenMatrix<T>& m)
    { return m.inverse(); }


    //
    // Matrix ==, != Matrix
    //

    template <typename T1, typename T2>
    bool operator==(
        const GenMatrix<T1>& m1, const GenMatrix<T2>& m2);
    template <typename T1, typename T2>
    inline bool operator!=(
        const GenMatrix<T1>& m1, const GenMatrix<T2>& m2)
    { return !(m1 == m2); }


    //
    // I/O
    //

    template <typename T>
    inline std::istream& operator>>(
        const TMV_Reader& reader, MatrixView<T> m)
    { m.read(reader); return reader.getis(); }

    template <typename T, int A>
    inline std::istream& operator>>(const TMV_Reader& reader, Matrix<T,A>& m)
    { m.read(reader); return reader.getis(); }

    template <typename T>
    inline std::istream& operator>>(std::istream& is, MatrixView<T> m)
    { return is >> IOStyle() >> m; }

    template <typename T, int A>
    inline std::istream& operator>>(std::istream& is, Matrix<T,A>& m)
    { return is >> IOStyle() >> m; }

} // namespace tmv

#endif
