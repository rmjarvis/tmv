//---------------------------------------------------------------------------
//
// This file defines the TMV Permutation class.
//
// A Permutation is an operator that rearranges the elements of 
// a Vector or the rows/columns of a matrix.  
//
// When viewed as a matrix each row and column has exactly one 1 and the 
// rest of the elements are 0.
//
// However, internally, this is not how a Permutation is stored.  It is
// stored as an ordered set of pairwise interchanges.  This is because
// it is very fast to apply a permutation in this form to a matrix.
//
// As a result, Permutation is derived from BaseMatrix, rather than
// BaseMatrix_Calc, since it takes a bit of calculation to figure out
// the values in matrix form.
//
// Also, we do not provide any public non-const methods for a permutation.
// The only way to create one is by passing it the vector of int's with
// the indices of the interchanges.  The intent is that the end user will
// never create a permutation.  This class is primarily designed to 
// make it easier to use the permutations that various TMV algorithms
// calculate.  (e.g. LU decomposition, QRP decomposition, and vector sort.)
//
// Note that Transpose(P) x P = I, so P is an orthogonal Matrix.
// Thus, inverse() and transpose() are the same thing.  Both return a 
// new permutation that knows to apply the permutations in the reverse
// order.  The indices are not copied.  It just passes a pointer to the 
// same array of int's that the first one uses.  So writing something like
// m = p.inverse() * m is efficient.
//
// Constructors:
//
//    Permutation(size_t n, const int* pp, bool isinv, int det)
//        Make a Permutation P from array pp.
//        If isinv = false, then the swaps are applied in order from the 
//        left to the rows of m.  If isinv = true, then in reverse order.
//        If det is given, it is taken to be the determinant of P. 
//        If it is omitted, the determinant will be calculated.
//
// Access Functions
//
//    size_t size() const
//    size_t colsize() const
//    size_t rowsize() const
//        Return the size of the Permutation.  
//        Since P is a square matrix, all three of these return the same value.
//
//    const int* getP() const
//        Return the integer array.
//
//    bool isInverse() const
//        Return whether the storage order is inverted.
//
//    size_t operator()(size_t i, size_t j) const
//        Returns the (i,j) element of P (always either 1 or 0).
//
//
// Functions of Permutations:
//        (These are all both member functions and functions of a Permutation,
//         so Det(p) and p.Det() for example are equivalent.)
//
//    p.sumElements() const
//    ... and other similarly trivial functions required by BaseMatrix.
//
//    p.trace() const
//       This one is required, but not trivial.  Requires a bit of calculation
//       to figure out.  O(N).
//
//    p.det() const
//       This is also O(N) to calculate, but it is calculated during
//       the creation of the permutation, so actually O(1) here.
//    
//
// Inverse:
//
//    p.inverse() const
//    p.transpose() const
//        Returns the inverse Permutation of p 
//
//
// I/O: 
//
//    p.writeCompact(os) const;
//        Writes p to ostream os in vector format:
//          size ( p(0) p(1) p(2) ... p(n-1) )
//
//    p.Write(os)
//        Write p as a Matrix<int>
//
// Operators:
//       Here we use p for a Permutation, m for a Matrix and v for a Vector
//
//
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


#ifndef TMV_Permutation_H
#define TMV_Permutation_H

#include "TMV_Base.h"
#include "TMV_Matrix.h"
#include <algorithm>

namespace tmv {

    class Permutation;

    template <>
    struct Traits<Permutation>
    {
        typedef int value_type;
        typedef int real_type;
        typedef std::complex<int> complex_type;
        enum { isreal = true };
        enum { iscomplex = false };

        typedef Permutation type;
        typedef Matrix<int> copy_type;
        typedef copy_type calc_type;
        typedef copy_type eval_type;
        typedef Permutation inverse_type;

        enum { _colsize = UNKNOWN };
        enum { _rowsize = UNKNOWN };
        enum { _shape = SquareRec };
        enum { _fort = false };
        enum { _calc = false };
    };

    class Permutation : public BaseMatrix<Permutation>
    {

    public:

        //
        // Constructors
        //

        Permutation(size_t n, const int* p, bool _isinv, int d) :
            itsn(n), itsp(p), isinv(_isinv), itsdet(d) {}

        // if det is unknown, calculate it now.
        Permutation(size_t n, const int* p, bool _isinv=false) :
            itsn(n), itsp(p), isinv(_isinv), itsdet(1) 
        { for(int i=0;i<itsn;++i) if (p[i] != i) itsdet = -itsdet; }

        Permutation(const Permutation& rhs) :
            itsn(rhs.itsn), itsp(rhs.itsp),
            isinv(rhs.isinv), itsdet(rhs.itsdet) {}

        ~Permutation() {}

        //
        // Access 
        //

        inline size_t size() const { return itsn; }

        inline size_t colsize() const { return itsn; }

        inline size_t rowsize() const { return itsn; }

        // TODO: This is fairly slow.  I don't know if there is a faster
        // algorithm, but I can certainly do this without the temporary
        // vector, since the only information about temp is the position
        // of the 1.  So this can be kept track of with a single int, rather
        // than a full vector.  And with that, I could return early whenever
        // it becomes impossible for the 1 to swap into the element that
        // we will return.
        inline int cref(size_t i,size_t j) const
        { 
            // Two options:
            // 1) P = P * I = I.permuteRows(p)
            //    P.col(j) = I.col(j).permute(p)
            //             = e_j.permute(p)
            //    P(i,j) = e_j.permute(p)(i)
            //           = e_j.permute(p,0,i+1)(i)
            //
            // 2) P = I * P = I.reversePermuteCols(p)
            //    P.row(i) = I.row(i).reversePermute(p)
            //             = e_i.reversePermute(p)
            //    P(i,j) = e_i.reversePermute(p)(j)
            //           = e_i.reversePermute(p,0,i+1)(j)
            // 
            // If isinv = true, these two options become:
            // 1) P(i,j) = e_j.reversePermute(p,0,j+1)(i)
            // 2) P(i,j) = e_i.permute(p,0,j+1)(j)
            //
            // I don't see any reason to pick one over the other, since they
            // both can stop at the same number of permutations, so the 
            // ops count is the same in either case.  
            // So I arbitrarily pick the permute version rather than the 
            // reversePermute version.
            // 
            tmv::Vector<int> temp(itsn);
            if (isinv) {
                temp.cMakeBasis(i);
                temp.cPermute(itsp,0,j+1);
                return temp.cref(j);
            } else {
                temp.cMakeBasis(j);
                temp.cPermute(itsp,0,i+1);
                return temp.cref(i);
            }
        }

        inline const int* getValues() const { return itsp; }

        inline bool isInverse() const { return isinv; }

        //
        // Functions of Permutations
        //

        inline Permutation inverse() const
        { return Permutation(itsn,itsp,!isinv,itsdet); }

        inline Permutation transpose() const
        { return inverse(); }

        inline int det() const
        { return itsdet; }

        inline int trace() const
        {
            // The trace of a permutation is the number of 1's on the diagonal.
            // This corresponds to the number of elements whose position
            // is unchanged after the permutation is applied.
            // So the algorithm used here makes a vector with each value
            // equal to its index.
            // Then apply the permutation and count how many are still in
            // the same position.
            Vector<int> v(itsn);
            for(int i=0;i<itsn;++i) v[i] = i;
            v.permute(itsp);  // istrans is irrelevant here.
            int t = 0;
            for(int i=0;i<itsn;++i) if (v[i] == i) ++t;
            return t;
        }

        inline int sumElements() const
        { return itsn; }

        inline int sumAbsElements() const
        { return itsn; }

        inline int maxAbsElement() const
        { return 1; }

        inline int normSq() const
        { return itsn; }

        inline int normSq(const int scale) const
        { return itsn * scale * scale; }

        // Normally these two would return int, but since there is 
        // a sqrt involved, better to upgrade to double.
        inline double normF() const
        { return TMV_SQRT(double(itsn)); }

        inline double norm() const
        { return normF(); }

        inline int norm1() const
        { return 1; }

        inline int norm2() const
        { return 1; }

        inline int normInf() const
        { return 1; }


        //
        // Create matrix version
        //
        
        template <class M2>
        inline void assignTo(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            m2.setToIdentity();
            applyOnLeft(m2);
        }

        template <class M2>
        inline void newAssignTo(BaseMatrix_Rec_Mutable<M2>& m2) const
        { assignTo(m2); }

        //
        // Apply permutation to a vector
        //

        template <class V2>
        inline void apply(BaseVector_Mutable<V2>& v2) const
        {
            if (isinv) v2.reversePermute(itsp);
            else v2.permute(itsp);
        }

        //
        // Apply permutation to a matrix
        //

        template <class M2>
        inline void applyOnLeft(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            if (isinv) m2.reversePermuteRows(itsp);
            else m2.permuteRows(itsp);
        }

        template <class M2>
        inline void applyOnRight(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            if (isinv) m2.permuteCols(itsp);
            else m2.reversePermuteCols(itsp);
        }

    protected:

        int itsn;
        const int*const itsp;
        bool isinv;
        int itsdet; // det = 1 or -1

    }; // Permutation

    inline Permutation Transpose(const Permutation& m)
    { return m.transpose(); }
    inline Permutation Adjoint(const Permutation& m)
    { return m.transpose(); }
    inline const Permutation& Conjugate(const Permutation& m)
    { return m; }
    inline Permutation Inverse(const Permutation& m)
    { return m.transpose(); }

    // 
    // Vector::sort 
    // Wait unil here to define the version that returns a Permutation.
    //
    
    template <class V>
    inline Permutation BaseVector_Mutable<V>::sort(
        int*const P, ADType ad, CompType comp)
    { tmv::Sort(*this,P,ad,comp); return Permutation(size(),P); }


    //
    // TMV_Text
    //
    
    inline std::string TMV_Text(const Permutation& )
    { return "Permutation"; }

 
} // namespace tmv

#endif
