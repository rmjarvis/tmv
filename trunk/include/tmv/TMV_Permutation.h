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
// There are only two ways to create a Permutation object.  The first
// constructor creates an identity permutation of a given size. 
// This can be passed to some functions which are able to modify it.
// Or, you can create one from a vector of int's with the indices
// of the interchanges.  Probably most users won't ever use this 
// constructor, but it is public to allow that possibility.
//
// This class is primarily designed to make it easier to use the 
// permutations that various TMV algorithms calculate.  
// (e.g. LU decomposition, QRP decomposition, and vector sort.)
//
// The copy semantics are similar to matrix Views.  The copies do not
// own their data.  However, if you perform a mutable action on the
// Permutation, then it will allocate memory and take ownership of 
// its values at that point.
//
// Note that Transpose(P) x P = I, so P is an orthogonal Matrix.
// Thus, inverse() and transpose() are the same thing.  Both return a 
// new permutation that knows to apply the permutations in the reverse
// order.  As with copies, the indices are not copied.  It just passes 
// a pointer to the same data that the first one uses.  So writing 
// something like m = p.inverse() * m is efficient.
//
// Constructors:
//
//    Permutation(size_t n)
//        Make an identity Permutation of size n.
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
//    const int* getValues() const
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
// Mutable functions
//
//   p.setToIdentity() 
//       Make p the identity permutation.
//   
//   p.transposeSelf()
//   p.invertSelf()
//       Make p into its inverse.
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
//    v / p  = p^-1 * v
//    v % p  = v * p^-1
//
//    p * m
//    m * p
//    m / p  = p^-1 * m
//    m % p  = m * p^-1
//
//    p == p
//    p != p
//
//       These all act in the way you would expect for the Matrix
//       formulation of p.


#ifndef TMV_Permutation_H
#define TMV_Permutation_H

#include "TMV_BaseMatrix_Rec.h"
#include "TMV_Array.h"

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

        enum { _colsize = TMV_UNKNOWN };
        enum { _rowsize = TMV_UNKNOWN };
        enum { _nlo = TMV_UNKNOWN };
        enum { _nhi = TMV_UNKNOWN };
        enum { _shape = Rec };
        enum { _fort = false };
        enum { _calc = false };
    };

    class Permutation : 
        public BaseMatrix<Permutation>
    {

    public:

        //
        // Constructors
        //

        explicit Permutation(size_t n=0) :
            itsn(n), itsmem(n), itsp(itsmem), isinv(false), itsdet(1) 
        { for(int i=0;i<itsn;++i) itsmem[i] = i; }

        Permutation(size_t n, const int* p, bool _isinv, int d) :
            itsn(n), itsp(p), isinv(_isinv), itsdet(d) {}

        // if det is unknown, calculate it now.
        Permutation(size_t n, const int* p, bool _isinv=false) :
            itsn(n), itsp(p), isinv(_isinv), itsdet(1) 
        { calcDet(); }

        Permutation(const Permutation& rhs) :
            itsn(rhs.itsn), itsp(rhs.itsp),
            isinv(rhs.isinv), itsdet(rhs.itsdet) {}

        ~Permutation() {}

        Permutation& operator=(const Permutation& rhs)
        {
            TMVAssert(size() == rhs.size());
            itsmem.resize(0); // also deallocates the memory
            itsp = rhs.itsp;
            isinv = rhs.isinv;
            itsdet = rhs.itsdet;
            return *this;
        }

        //
        // Access 
        //

        int cref(int i, int j) const
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
            // which is equivalent to just swapping the i,j values
            // (and using the other algorithm).
            //
            // The two options use the same number of permutations.
            // However, with the forward permutation, it is easy 
            // to see where we can stop the loop early.  Specifically,
            // once the location of the 1 is smaller than k,
            // it cannot be modified further, since itsp[k] >= k.
            // So we choose to use the forward loop option.

            if (isinv) TMV_SWAP(i,j);
            int temp = j;
            for(int k=0;k<=i && k<=temp;++k) if (itsp[k]!=k) {
                if (temp == k) temp = itsp[k];
                else if (temp == itsp[k]) temp = k;
            }
            return (temp == i) ? 1 : 0;
        }


        //
        // Functions of Permutations
        //

        TMV_INLINE Permutation inverse() const
        { return Permutation(itsn,itsp,!isinv,itsdet); }

        TMV_INLINE Permutation transpose() const
        { return inverse(); }

        TMV_INLINE int det() const
        { return itsdet; }

        int logDet(int* sign=0) const
        { if (sign) *sign = itsdet; return 0; }

        int trace() const
        {
            // The trace of a permutation is the number of 1's on the diagonal.
            // This corresponds to the number of elements whose position
            // is unchanged after the permutation is applied.
            // So the algorithm used here makes a vector with each value
            // equal to its index.
            // Then apply the permutation and count how many are still in
            // the same position.
            AlignedArray<int> temp(itsn);
            makeIndex(temp.get());
            int t = 0;
            for(int k=0;k<itsn;++k) if (temp[k] == k) ++t;
            return t;
        }

        TMV_INLINE int sumElements() const
        { return itsn; }

        TMV_INLINE int sumAbsElements() const
        { return itsn; }

        TMV_INLINE int sumAbs2Elements() const
        { return itsn; }

        TMV_INLINE int maxAbsElement() const
        { return 1; }

        TMV_INLINE int maxAbs2Element() const
        { return 1; }

        TMV_INLINE int normSq() const
        { return itsn; }

        int normSq(const int scale) const
        { return itsn * scale * scale; }

        // Normally these two would return int, but since there is 
        // a sqrt involved, better to upgrade to double.
        double normF() const
        { return TMV_SQRT(double(itsn)); }

        TMV_INLINE double norm() const
        { return normF(); }

        TMV_INLINE int norm1() const
        { return 1; }

        TMV_INLINE int norm2() const
        { return 1; }

        TMV_INLINE int normInf() const
        { return 1; }

        TMV_INLINE int condition() const
        { return 1; }

        TMV_INLINE bool isSinular() const
        { return false; }


        //
        // op==
        //
        friend bool operator==(
            const Permutation& p1, const Permutation& p2)
        {
            TMVAssert(p1.size() == p2.size());
            const int n = p1.itsn;
            if (p1.isinv == p2.isinv) {
                for(int i=0;i<n;++i) {
                    if (p1.itsp[i] != p2.itsp[i]) return false;
                }
                return true;
            } else {
                // If not the same storage, then this requires a bit of work
                // to see if they effect the same permutation.
                AlignedArray<int> temp1(n);
                AlignedArray<int> temp2(n);
                p1.makeIndex(temp1.get());
                p2.makeIndex(temp2.get());
                for(int i=0;i<n;++i) {
                    if (temp1[i] != temp2[i]) return false;
                }
                return true;
            }
        }

        friend bool operator!=(
            const Permutation& p1, const Permutation& p2)
        { return !(p1==p2); }

        //
        // Mutable Functions
        // (Not too many here, since most normal functions are invalid
        //  for permutations.)
        //

        Permutation& setToIdentity()
        {
            allocateMem();
            for(int i=0;i<itsn;++i) itsmem[i] = i;
            isinv = false;
            itsdet = 1;
            return *this;
        }

        Permutation& transposeSelf()
        { isinv = !isinv; return *this; }

        Permutation& invertSelf()
        { return transposeSelf(); }


        //
        // Create matrix version
        //

        template <class M2>
        void assignTo(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            m2.setToIdentity();
            applyOnLeft(m2);
        }

        template <class M2>
        void newAssignTo(BaseMatrix_Rec_Mutable<M2>& m2) const
        { assignTo(m2); }

        //
        // MakeInverse
        //

        template <class T2>
        void makeInverse(BaseMatrix_Rec_Mutable<T2>& minv) const
        { inverse().assignTo(minv); }

        // (PtP)^-1 = P^-1 Pt^-1 = Pt P = I
        template <class T2>
        void makeInverseATA(BaseMatrix_Rec_Mutable<T2>& ata) const
        { ata.setToIdentity(); }


        //
        // I/O
        //

        void write(std::ostream& os) const
        {
            AlignedArray<int> temp(itsn);
            makeIndex(temp.get());
            os<<itsn<<"  "<<itsn<<std::endl;
            for(int i=0;i<itsn;++i) {
                os << "( ";
                for(int j=0;j<itsn;++j)
                    os<<' '<<(temp[i]==j ? 1 : 0)<<' ';
                os << " )\n";
            }
        }

        void writeCompact(std::ostream& os) const
        {
            os<<"P "<<itsn<<" "<<isinv<<" (";
            for(int i=0;i<itsn;++i) os<<' '<<itsp[i]<<' ';
            os << ")";
        }

        void read(std::istream& is);

        //
        // Apply permutation to a vector
        //
        
        template <class V2>
        void applyOnLeft(BaseVector_Mutable<V2>& v2) const
        {
            if (isinv) v2.reversePermute(itsp);
            else v2.permute(itsp);
        }
 
        template <class V2>
        void applyOnRight(BaseVector_Mutable<V2>& v2) const
        {
            if (isinv) v2.permute(itsp);
            else v2.reversePermute(itsp);
        }

        //
        // Apply permutation to a matrix
        //

        template <class M2>
        void applyOnLeft(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            if (isinv) m2.reversePermuteRows(itsp);
            else m2.permuteRows(itsp);
        }

        template <class M2>
        void applyOnRight(BaseMatrix_Rec_Mutable<M2>& m2) const
        {
            if (isinv) m2.permuteCols(itsp);
            else m2.reversePermuteCols(itsp);
        }

        //
        // Friend functions that can act on a mutable Permutation.
        //

        friend void Swap(Permutation& p1, Permutation& p2)
        {
            TMVAssert(p1.size() == p2.size());
            p1.itsmem.swapWith(p2.itsmem);
            TMV_SWAP(p1.itsp,p2.itsp);
            TMV_SWAP(p1.isinv,p2.isinv);
            TMV_SWAP(p1.itsdet,p2.itsdet);
        }

        // Defined below.
        template <class V>
        friend V& BaseVector_Mutable<V>::sort(
            Permutation& P, ADType ad, CompType comp);

        // In TMV_LUDecompose.h
        template <class M>
        friend void LU_Decompose(
            BaseMatrix_Rec_Mutable<M>& m, Permutation& P);

        // In TMV_QRPDecompose.h
        template <class M, class V>
        friend void QRP_Decompose(
            BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta,
            Permutation& P, bool strict=false);

#if 0
        template <class T>
        friend void LU_Decompose(
            const GenBandMatrix<T>& A, const LowerTriMatrixView<T>& L,
            const BandMatrixView<T>& U, Permutation& P)
        {
            TMVAssert(P.size() == A.colsize());
            P.allocateMem();
            LU_Decompose(A,L,U,P.getMem());
            P.isinv = true;
            P.calcDet();
        }

        template <class T>
        friend void LU_Decompose(
            const BandMatrixView<T>& A, Permutation& P, int nhi)
        {
            TMVAssert(P.size() == A.colsize());
            P.allocateMem();
            LU_Decompose(A,P.getMem(),P.itsdet=1,nhi);
            P.isinv = true;
        }

        template <class T>
        friend void LDL_Decompose(
            const SymMatrixView<T>& A, const SymBandMatrixView<T>& D,
            Permutation& P)
        {
            TMVAssert(P.size() == A.colsize());
            P.allocateMem();
            LDL_Decompose(A,D,P.getMem());
            P.isinv = true;
            P.calcDet();
        }

        template <class T>
        friend void LDL_Decompose(
            const SymMatrixView<T>& A, const VectorView<T>& xD,
            Permutation& P, TMV_RealType(T)& logdet, T& signdet)
        {
            TMVAssert(P.size() == A.size());
            P.allocateMem();
            LDL_Decompose(A,xD,P.getMem(),logdet,signdet);
            P.isinv = true;
            P.calcDet();
        }
#endif

        template <class V>
        friend void DoVectorSort(
            BaseVector_Mutable<V>& v, Permutation& P, ADType ad, CompType comp)
        {
            TMVAssert(P.size() == v.size());
            P.allocateMem();
            tmv::Sort(v,P.getMem(),ad,comp);
            P.calcDet();
            P.isinv = false;
        }

        //
        // Auxilliary functions
        //

        void resize(size_t n)
        {
            if (int(n) > itsn) {
                itsmem.resize(n);
                itsp = itsmem.get();
            }
            itsn = n;
            isinv = false;
            itsdet = 1;
        }

        TMV_INLINE size_t size() const { return itsn; }
        TMV_INLINE size_t colsize() const { return itsn; }
        TMV_INLINE int nlo() const { return TMV_MAX(itsn-1,0); }
        TMV_INLINE int nhi() const { return TMV_MAX(itsn-1,0); }
        TMV_INLINE size_t rowsize() const { return itsn; }
        TMV_INLINE const int* getValues() const { return itsp; }
        TMV_INLINE bool isInverse() const { return isinv; }


    protected:

        int itsn;
        AlignedArray<int> itsmem;
        const int* itsp;
        bool isinv;
        int itsdet; // det = 1 or -1

        void makeIndex(int* index) const
        {
            for(int k=0;k<itsn;++k) index[k] = k;
            if (isinv) {
                for(int k=itsn-1;k>=0;--k)
                    if (itsp[k]!=k) TMV_SWAP(index[k],index[itsp[k]]);
            } else {
                for(int k=0;k<itsn;++k)
                    if (itsp[k]!=k) TMV_SWAP(index[k],index[itsp[k]]);
            }
        }

        // 
        // Helper functions for mutable actions.
        //

        void allocateMem()
        {
            if (!itsmem.get()) {
                itsmem.resize(itsn);
                itsp = itsmem.get();
            }
        }

        void saveToMem()
        {
            if (!itsmem.get()) {
                itsmem.resize(itsn);
                for(int i=0;i<itsn;++i) itsmem[i] = itsp[i];
                itsp = itsmem.get();
            }
        }

        void makeCopyOf(const Permutation& orig)
        {
            itsn = orig.itsn;
            itsmem.resize(itsn);
            for(int i=0;i<itsn;++i) itsmem[i] = orig.itsp[i];
            itsp = itsmem.get();
            isinv = orig.isinv;
            itsdet = orig.itsdet;
        }

        void calcDet() 
        {
            itsdet = 1; 
            for(int i=0;i<itsn;++i) 
                if (itsp[i] != i) itsdet = -itsdet; 
        }

        TMV_INLINE_ND int* getMem() 
        {
            // Make sure P owns its memory:
            TMVAssert(itsn==0 || itsmem.get());
            // This next one shoudl be true if the previous one passes.
            TMVAssert(itsmem.get() == itsp);
            return itsmem;
        }

    }; // Permutation

    static TMV_INLINE Permutation Transpose(const Permutation& m)
    { return m.transpose(); }
    static TMV_INLINE Permutation Adjoint(const Permutation& m)
    { return m.transpose(); }
    static TMV_INLINE const Permutation& Conjugate(const Permutation& m)
    { return m; }
    static TMV_INLINE Permutation Inverse(const Permutation& m)
    { return m.transpose(); }
    static TMV_INLINE double Norm(const Permutation& m)
    { return m.norm(); }
    static TMV_INLINE double NormF(const Permutation& m)
    { return m.normF(); }
    static TMV_INLINE int NormSq(const Permutation& m)
    { return m.normSq(); }
    static TMV_INLINE int Norm1(const Permutation& m)
    { return m.norm1(); }
    static TMV_INLINE int Norm2(const Permutation& m)
    { return m.norm2(); }
    static TMV_INLINE int NormInf(const Permutation& m)
    { return m.normInf(); }
    static TMV_INLINE int MaxAbsElement(const Permutation& m)
    { return m.maxAbsElement(); }
    static TMV_INLINE int MaxAbs2Element(const Permutation& m)
    { return m.maxAbs2Element(); }
    static TMV_INLINE int Trace(const Permutation& m)
    { return m.trace(); }
    static TMV_INLINE int Det(const Permutation& m)
    { return m.det(); }
    static TMV_INLINE int LogDet(const Permutation& m)
    { return m.logDet(); }

    static inline std::ostream& operator<<(
        std::ostream& os, const Permutation& p)
    { p.write(os); return os; }


#ifndef TMV_NO_THROWW
    class PermutationReadError : 
        public ReadError
    {
    public :
        Permutation m;
        int i;
        char exp,got;
        size_t n;
        bool is, iseof, isbad;

        PermutationReadError(
            int _i, const Permutation& _m, std::istream& _is
        ) throw() :
            ReadError("Permutation."),
            i(_i), m(_m), exp(0), got(0),
            n(_m.size()), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        PermutationReadError(std::istream& _is) throw() :
            ReadError("Permutation."),
            i(0), exp(0), got(0), n(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        PermutationReadError(
            int _i, const Permutation& _m,
            std::istream& _is, char _e, char _g
        ) throw() :
            ReadError("Permutation."),
            i(_i), m(_m), exp(_e), got(_g),
            n(_m.size()), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        PermutationReadError(std::istream& _is, char _e, char _g) throw() :
            ReadError("Permutation."),
            i(0), exp(_e), got(_g),
            n(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        PermutationReadError(
            const Permutation& _m, std::istream& _is, size_t _n
        ) throw() :
            ReadError("Permutation."),
            i(0), m(_m), exp(0), got(0),
            n(_n), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        PermutationReadError(const PermutationReadError& rhs) :
            m(rhs.m), i(rhs.i), exp(rhs.exp), got(rhs.got),
            n(rhs.n), is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        virtual ~PermutationReadError() throw() {}

        virtual void write(std::ostream& os) const throw()
        {
            os<<"TMV Read Error: Reading istream input for Permutation\n";
            if (exp != got) {
                os<<"Wrong format: expected '"<<exp<<"', got '"<<got<<"'.\n";
            }
            if (n != m.size()) {
                os<<"Wrong size: expected "<<m.size()<<", got "<<n<<".\n";
            }
            if (!is) {
                if (iseof) {
                    os<<"Input stream reached end-of-file prematurely.\n";
                } else if (isbad) {
                    os<<"Input stream is corrupted.\n";
                } else {
                    os<<"Input stream cannot read next character.\n";
                }
            }
            if (m.size() > 0) {
                os<<"The portion of the Permutation which was successfully "
                    "read is: \n";
                os<<"( ";
                for(int k=0;k<i;++k)
                    os<<' '<<m.getValues()[k]<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    inline void Permutation::read(std::istream& is)
    {
        allocateMem();
        int _inv;
        is >> _inv;
        if (!is) {
#ifdef TMV_NO_THROWW
            std::cerr<<"Permutation Read Error !is \n";
            exit(1);
#else
            throw PermutationReadError(is);
#endif
        }
        isinv = _inv;
        char paren;
        is >> paren;
        if (!is || paren != '(') {
#ifdef TMV_NO_THROWW
            std::cerr<<"Permutation ReadError: "<<paren<<" != (\n";
            exit(1);
#else
            throw PermutationReadError(0,*this,is,'(',is?paren:'(');
#endif
        }
        for(int i=0;i<itsn;++i) {
            is >> itsmem[i];
            if (!is) {
#ifdef TMV_NO_THROWW
                std::cerr<<"Permutation ReadError: !is \n";
                exit(1);
#else
                throw PermutationReadError(i,*this,is);
#endif
            }
        }
        is >> paren;
        if (!is || paren != ')') {
#ifdef TMV_NO_THROWW
            std::cerr<<"Permutation ReadError: "<<paren<<" != )\n";
            exit(1);
#else
            throw PermutationReadError(size(),*this,is,')',is?paren:')');
#endif
        }
        calcDet();
    }

    static inline std::istream& operator>>(std::istream& is, Permutation& m)
    {
        char p;
        is >> p;
        if (!is) {
#ifdef TMV_NO_THROWW
            std::cerr<<"Permutation Read Error !is \n";
            exit(1);
#else
            throw PermutationReadError(is);
#endif
        }
        if (p != 'P') {
#ifdef TMV_NO_THROWW
            std::cerr<<"Permutation Read Error "<<p<<" != P\n";
            exit(1);
#else
            throw PermutationReadError(is,'P',p);
#endif
        }
        size_t n;
        is >> n;
        if (!is) {
#ifdef TMV_NO_THROWW
            std::cerr<<"Permutation Read Error !is \n";
            exit(1);
#else
            throw PermutationReadError(is);
#endif
        }
        m.resize(n);
        m.read(is);
        return is;
    }


    // 
    // Vector::sort 
    // Wait unil here to define the version with Permutation.
    //

    template <class V>
    TMV_INLINE V& BaseVector_Mutable<V>::sort(
        Permutation& P, ADType ad, CompType comp)
    {
        DoVectorSort(vec(),P,ad,comp);
        return vec();
    }


    //
    // TMV_Text
    //
    
#ifdef TMV_TEXT
    static inline std::string TMV_Text(const Permutation& )
    { return "Permutation"; }
#endif

 
} // namespace tmv

#endif
