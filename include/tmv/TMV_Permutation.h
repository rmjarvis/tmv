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
//    Permutation(int n)
//        Make an identity Permutation of size n.
//
//    Permutation(int n, const int* pp, bool isinv)
//        Make a Permutation P from array pp.
//        If isinv = false, then the swaps are applied in order from the 
//        left to the rows of m.  If isinv = true, then in reverse order.
//
// Access Functions
//
//    int size() const
//    int colsize() const
//    int rowsize() const
//        Return the size of the Permutation.  
//        Since P is a square matrix, all three of these return the same value.
//
//    const int* getValues() const
//        Return the integer array.
//
//    bool isInverse() const
//        Return whether the storage order is inverted.
//
//    int operator()(int i, int j) const
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
//       This is also O(N) to calculate
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
//    os << p;
//        Write p in the same format as a Matrix<int>
//
//    os << tmv::CompactIO() << p;
//        Writes p to ostream os in vector format:
//          P size inv ( p(0) p(1) p(2) ... p(n-1) )
//
//    is >> tmv::CompactIO() >> p;
//        Read p in the compact format.
//        It is not possible to read in a permutation that was written
//        using the normal format.
//
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
#include "TMV_IOStyle.h"

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

        enum { _colsize = Unknown };
        enum { _rowsize = Unknown };
        enum { _nlo = Unknown };
        enum { _nhi = Unknown };
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

        explicit Permutation(ptrdiff_t n=0) :
            itsn(n), itsmem(n), itsp(itsmem), isinv(false)
        {
            TMVAssert(n >= 0);
            for(ptrdiff_t i=0;i<itsn;++i) itsmem[i] = i; 
        }

        Permutation(ptrdiff_t n, const ptrdiff_t* p, bool _isinv=false) :
            itsn(n), itsp(p), isinv(_isinv)
        { TMVAssert(n >= 0); }

        Permutation(const Permutation& rhs) :
            itsn(rhs.itsn), itsp(rhs.itsp), isinv(rhs.isinv) {}

        ~Permutation() {}

        Permutation& operator=(const Permutation& rhs)
        {
            TMVAssert(size() == rhs.size());
            itsmem.resize(0); // also deallocates the memory
            itsp = rhs.itsp;
            isinv = rhs.isinv;
            return *this;
        }

        //
        // Access 
        //

        int cref(ptrdiff_t i, ptrdiff_t j) const
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
            ptrdiff_t temp = j;
            for(ptrdiff_t k=0;k<=i && k<=temp;++k) if (itsp[k]!=k) {
                if (temp == k) temp = itsp[k];
                else if (temp == itsp[k]) temp = k;
            }
            return (temp == i) ? 1 : 0;
        }


        //
        // Functions of Permutations
        //

        TMV_INLINE Permutation inverse() const
        { return Permutation(itsn,itsp,!isinv); }

        TMV_INLINE Permutation transpose() const
        { return inverse(); }

        TMV_INLINE int det() const
        {
            int d = 1; 
            for(ptrdiff_t i=0;i<itsn;++i) if (itsp[i] != i) d = -d; 
            return d;
        }

        int logDet(int* sign=0) const
        { if (sign) *sign = det(); return 0; }

        int trace() const
        {
            // The trace of a permutation is the number of 1's on the diagonal.
            // This corresponds to the number of elements whose position
            // is unchanged after the permutation is applied.
            // So the algorithm used here makes a vector with each value
            // equal to its index.
            // Then apply the permutation and count how many are still in
            // the same position.
            AlignedArray<ptrdiff_t> temp(itsn);
            makeIndex(temp.get());
            int t = 0;
            for(ptrdiff_t k=0;k<itsn;++k) if (temp[k] == k) ++t;
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
            const ptrdiff_t n = p1.itsn;
            if (p1.isinv == p2.isinv) {
                for(ptrdiff_t i=0;i<n;++i) {
                    if (p1.itsp[i] != p2.itsp[i]) return false;
                }
                return true;
            } else {
                // If not the same storage, then this requires a bit of work
                // to see if they effect the same permutation.
                AlignedArray<ptrdiff_t> temp1(n);
                AlignedArray<ptrdiff_t> temp2(n);
                p1.makeIndex(temp1.get());
                p2.makeIndex(temp2.get());
                for(ptrdiff_t i=0;i<n;++i) {
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
            for(ptrdiff_t i=0;i<itsn;++i) itsmem[i] = i;
            isinv = false;
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
    
        inline void write(const TMV_Writer& writer) const
        {
            const ptrdiff_t N = size();
            writer.begin();
            writer.writeCode("P");
            writer.writeSize(N);
            writer.writeSimpleSize(N);
            writer.writeFullSize(isinv);
            writer.writeStart();

            if (writer.isCompact()) {
                writer.writeLParen();
                for(ptrdiff_t i=0;i<N;++i) {
                    if (i > 0) writer.writeSpace();
                    writer.writeValue(itsp[i]);
                }
                writer.writeRParen();
            } else {
                AlignedArray<ptrdiff_t> temp(N);
                makeIndex(temp.get());
                for(ptrdiff_t i=0;i<N;++i) {
                    writer.writeLParen();
                    for(ptrdiff_t j=0;j<N;++j) {
                        if (j>0) writer.writeSpace();
                        writer.writeValue(temp[i]==j ? 1 : 0);
                    }
                    writer.writeRParen();
                    if (i < N-1) writer.writeRowEnd();
                }
            }

            writer.writeFinal();
            writer.end();
        }

        inline void read(const TMV_Reader& reader);

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
        }

        // Defined below.
        template <class V>
        friend V& BaseVector_Mutable<V>::sort(Permutation& P, ADType ad, CompType comp);

        // In TMV_LUDecompose.h
        template <class M>
        friend void LU_Decompose(BaseMatrix_Rec_Mutable<M>& m, Permutation& P);

        // In TMV_QRPDecompose.h
        template <class M, class V>
        friend void QRP_Decompose(
            BaseMatrix_Rec_Mutable<M>& m, BaseVector_Mutable<V>& beta,
            Permutation& P, bool strict=false);

        // In TMV_BandLUDecompose.h
        template <class M>
        friend void BandLU_Decompose(BaseMatrix_Band_Mutable<M>& m, Permutation& P);

#if 0
        template <class T>
        friend void LDL_Decompose(
            const SymMatrixView<T>& A, const SymBandMatrixView<T>& D,
            Permutation& P)
        {
            TMVAssert(P.size() == A.colsize());
            P.allocateMem();
            LDL_Decompose(A,D,P.getMem());
            P.isinv = true;
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
        }
#endif

        template <class V>
        friend void DoVectorSort(
            BaseVector_Mutable<V>& v, Permutation& P, ADType ad, CompType comp)
        {
            TMVAssert(P.size() == v.size());
            P.allocateMem();
            tmv::Sort(v,P.getMem(),ad,comp);
            P.isinv = false;
        }

        //
        // Auxilliary functions
        //

        void resize(ptrdiff_t n)
        {
            TMVAssert(n >= 0);
            if (n > itsn) {
                itsmem.resize(n);
                itsp = itsmem.get();
            }
            itsn = n;
            isinv = false;
        }

        TMV_INLINE ptrdiff_t size() const { return itsn; }
        TMV_INLINE ptrdiff_t colsize() const { return itsn; }
        TMV_INLINE ptrdiff_t nlo() const { return TMV_MAX(itsn-1,ptrdiff_t(0)); }
        TMV_INLINE ptrdiff_t nhi() const { return TMV_MAX(itsn-1,ptrdiff_t(0)); }
        TMV_INLINE ptrdiff_t rowsize() const { return itsn; }
        TMV_INLINE const ptrdiff_t* getValues() const { return itsp; }
        TMV_INLINE bool isInverse() const { return isinv; }


    protected:

        ptrdiff_t itsn;
        AlignedArray<ptrdiff_t> itsmem;
        const ptrdiff_t* itsp;
        bool isinv;

        inline void makeIndex(ptrdiff_t* index) const
        {
            for(ptrdiff_t k=0;k<itsn;++k) index[k] = k;
            if (isinv) {
                for(ptrdiff_t k=itsn-1;k>=0;--k)
                    if (itsp[k]!=k) TMV_SWAP(index[k],index[itsp[k]]);
            } else {
                for(ptrdiff_t k=0;k<itsn;++k)
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
                for(ptrdiff_t i=0;i<itsn;++i) itsmem[i] = itsp[i];
                itsp = itsmem.get();
            }
        }

        void makeCopyOf(const Permutation& orig)
        {
            itsn = orig.itsn;
            itsmem.resize(itsn);
            for(ptrdiff_t i=0;i<itsn;++i) itsmem[i] = orig.itsp[i];
            itsp = itsmem.get();
            isinv = orig.isinv;
        }

        TMV_INLINE_ND ptrdiff_t* getMem() 
        {
            // Make sure P owns its memory:
            TMVAssert(itsn==0 || itsmem.get());
            // This next one shoudl be true if the previous one passes.
            TMVAssert(itsmem.get() == itsp);
            return itsmem;
        }

    }; // Permutation

    TMV_INLINE Permutation Transpose(const Permutation& m)
    { return m.transpose(); }
    TMV_INLINE Permutation Adjoint(const Permutation& m)
    { return m.transpose(); }
    TMV_INLINE const Permutation& Conjugate(const Permutation& m)
    { return m; }
    TMV_INLINE Permutation Inverse(const Permutation& m)
    { return m.transpose(); }
    TMV_INLINE double Norm(const Permutation& m)
    { return m.norm(); }
    TMV_INLINE double NormF(const Permutation& m)
    { return m.normF(); }
    TMV_INLINE int NormSq(const Permutation& m)
    { return m.normSq(); }
    TMV_INLINE int Norm1(const Permutation& m)
    { return m.norm1(); }
    TMV_INLINE int Norm2(const Permutation& m)
    { return m.norm2(); }
    TMV_INLINE int NormInf(const Permutation& m)
    { return m.normInf(); }
    TMV_INLINE int MaxAbsElement(const Permutation& m)
    { return m.maxAbsElement(); }
    TMV_INLINE int MaxAbs2Element(const Permutation& m)
    { return m.maxAbs2Element(); }
    TMV_INLINE int Trace(const Permutation& m)
    { return m.trace(); }
    TMV_INLINE int Det(const Permutation& m)
    { return m.det(); }
    TMV_INLINE int LogDet(const Permutation& m)
    { return m.logDet(); }

    inline std::ostream& operator<<(const TMV_Writer& writer, const Permutation& p)
    { p.write(writer); return writer.getos(); }

    inline std::ostream& operator<<(std::ostream& os, const Permutation& p)
    { return os << IOStyle() << p; }


#ifndef TMV_NO_THROWW
    class PermutationReadError : 
        public ReadError
    {
    public :
        Permutation m;
        ptrdiff_t i;
        std::string exp,got;
        ptrdiff_t n;
        bool is, iseof, isbad;

        PermutationReadError(std::istream& _is) throw() :
            ReadError("Permutation."),
            i(0), n(0), is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        PermutationReadError(
            std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("Permutation."),
            i(0), exp(_e), got(_g), n(0),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        PermutationReadError(
            const Permutation& _m, std::istream& _is, ptrdiff_t _n) throw() :
            ReadError("Permutation."),
            m(_m), i(0), n(_n),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        PermutationReadError(
            ptrdiff_t _i, const Permutation& _m, std::istream& _is) throw() :
            ReadError("Permutation."),
            m(_m), i(_i), n(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}
        PermutationReadError(
            ptrdiff_t _i, const Permutation& _m, std::istream& _is,
            const std::string& _e, const std::string& _g) throw() :
            ReadError("Permutation."),
            m(_m), i(_i), exp(_e), got(_g), n(_m.size()),
            is(_is), iseof(_is.eof()), isbad(_is.bad()) {}

        PermutationReadError(const PermutationReadError& rhs) :
            m(rhs.m), i(rhs.i), exp(rhs.exp), got(rhs.got), n(rhs.n),
            is(rhs.is), iseof(rhs.iseof), isbad(rhs.isbad) {}
        ~PermutationReadError() throw() {}

        void write(std::ostream& os) const throw()
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
                for(ptrdiff_t k=0;k<i;++k) os<<' '<<m.getValues()[k]<<' ';
                os<<" )\n";
            }
        }
    };
#endif

    inline void Permutation::read(const TMV_Reader& reader)
    {
        std::string exp,got;
        ptrdiff_t temp;
        if (!reader.readCode("P",exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Permutation Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw PermutationReadError(reader.getis(),exp,got);
#endif
        }
        ptrdiff_t n=size();
        if (!reader.readSize(n,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Permutation Read Error: reading size\n";
            exit(1);
#else
            throw PermutationReadError(reader.getis(),exp,got);
#endif
        }
        if (n != size()) resize(n);
        n=size();
        if (!reader.readSimpleSize(n,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Permutation Read Error: reading size\n";
            exit(1);
#else
            throw PermutationReadError(reader.getis(),exp,got);
#endif
        }
        if (n != size()) {
#ifdef NOTHROW
            std::cerr<<"Permutation Read Error: wrong size\n";
            exit(1);
#else
            throw PermutationReadError(*this,reader.getis(),n);
#endif
        }

        if (!reader.readFullSize(temp,exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Permutation Read Error: reading inv\n";
            exit(1);
#else
            throw PermutationReadError(reader.getis(),exp,got);
#endif
        }
        isinv = temp;

        if (!reader.isCompact()) {
#ifdef NOTHROW
            std::cerr<<"NonCompact Read is not supported for Permutation");
            exit(1);
#else
            throw ReadError(
                "Permutation.\n"
                "NonCompact Read is not supported for Permutation.");
#endif
        }

        if (!reader.readStart(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Permutation Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw PermutationReadError(0,*this,reader.getis(),exp,got);
#endif
        }
        if (!reader.readLParen(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Permutation Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw PermutationReadError(0,*this,reader.getis(),exp,got);
#endif
        }
        for(ptrdiff_t i=0;i<itsn;++i) {
            if (i>0 && !reader.readSpace(exp,got)) {
#ifdef NOTHROW
                std::cerr<<"Permutation ReadError: "<<got<<" != "<<exp<<std::endl;
                exit(1);
#else
                throw PermutationReadError(i,*this,reader.getis(),exp,got);
#endif
            }
            if (!reader.readValue(temp)) {
#ifdef NOTHROW
                std::cerr<<"Permutation ReadError: reading value\n";
                exit(1);
#else
                throw PermutationReadError(i,*this,reader.getis());
#endif
            }
            itsmem[i] = temp;
        }
        if (!reader.readRParen(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Permutation Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw PermutationReadError(itsn,*this,reader.getis(),exp,got);
#endif
        }
        if (!reader.readFinal(exp,got)) {
#ifdef NOTHROW
            std::cerr<<"Permutation Read Error: "<<got<<" != "<<exp<<std::endl;
            exit(1);
#else
            throw PermutationReadError(itsn,*this,reader.getis(),exp,got);
#endif
        }
    }

    inline std::istream& operator>>(const TMV_Reader& reader, Permutation& p)
    { p.read(reader); return reader.getis(); }

    inline std::istream& operator>>(std::istream& is, Permutation& m)
    { return is >> IOStyle() >> m; }

    // 
    // Vector::sort 
    // Wait unil here to define the version with Permutation.
    //

    template <class V>
    TMV_INLINE V& BaseVector_Mutable<V>::sort(Permutation& P, ADType ad, CompType comp)
    {
        DoVectorSort(vec(),P,ad,comp);
        return vec();
    }


    //
    // TMV_Text
    //
    
    inline std::string TMV_Text(const Permutation& )
    { return "Permutation"; }

 
} // namespace tmv

#endif
