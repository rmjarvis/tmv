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
// This file defines the TMV DiagMatrix class.
//
// The DiagMatrix class is provided for efficient storage of a diagonal
// matrix.  You can do most of the things that you can do with a 
// regular Matrix, but it will do them more efficiently.
//
// Constructors:
//
//    DiagMatrix<T>(size_t size)
//        Makes a DiagMatrix with column size and row size = size
//        with _uninitialized_ values
//
//    DiagMatrix<T>(size_t size, T x)
//        Makes a DiagMatrix of size n with all values = x
//
//    DiagMatrix<T>(size_t size, T* vv)
//    DiagMatrix<T>(size_t size, const std::vector<T>& vv)
//        Makes a DiagMatrix of size n which copies the values is vv
//
//    DiagMatrix<T>(const Vector<T>& vv)
//        Make a DiagMatrix which copies the elements of vv.
//
//    ConstDiagMatrixView<T>(const Vector<T>& v)
//        Make a constant DiagMatrix view with v as the diagonal.
//        While this view cannon be modified, changing the original v or m
//        will cause corresponding changes in this view.
//
//    DiagMatrixView<T>(Vector<T>& v)
//        Make a mutable DiagMatrix view with v as the diagonal.
//
//
// Access Functions
//
//    size_t colsize() const
//    size_t rowsize() const
//    size_t size() const
//        Return the dimensions of the DiagMatrix
//
//    T& operator()(int i)
//    T operator()(int i) const
//    T& operator()(int i, int j)
//    T operator()(int i, int j) const
//        Return the (i,j) element of the DiagMatrix
//        For the single paramter version, j=i
//
//    VectorView& diag()
//    ConstVectorView& diag() const
//        Return the diagonal of the DiagMatrix as a VectorView
//
//
// Modifying Functions - The same as the regular Matrix counterparts
//
//    DiagMatrix& setZero()
//    DiagMatrix& setAllTo(T x)
//    DiagMatrix<T>& transposeSelf() 
//        (Does nothing.)
//    DiagMatrix& conjugateSelf()
//    DiagMatrix& setToIdentity(x = 1)
//    void Swap(DiagMatrix& m1, DiagMatrix& m2)
//
//
// subDiagMatrix:
//
//    DiagMatrixView subDiagMatrix(int i1, int i2, int istep=1)
//        Returns a Sub-DiagMatrix which extends from i1 to i2 (step istep)
//        which refers to the same physical elements as the original.
//        As usual, i2 is the "one past the end" element.
//
//    DiagMatrixView realPart()
//    DiagMatrixView imagPart()
//        Returns a view to the real/imag elements of a complex DiagMatrix.
//
//    DiagMatrixView view()
//    DiagMatrixView transpose()
//        Returns a view of a DiagMatrix.
//
//    DiagMatrixView conjugate(m)
//    DiagMatrixView adjoint(m)
//        Returns a view to the conjugate of a DiagMatrix.
//
// Functions of DiagMatrices - Same as for regular Matrices:
//
//    m.det() or Det(m)
//    m.logDet() or m.logDet(T* sign) or LogDet(m)
//    m.trace() or Trace(m)
//    m.sumElements() or SumElements(m)
//    m.sumAbsElements() or SumAbsElements(m)
//    m.norm() or m.normF() or Norm(m) or NormF(m)
//    m.normSq() or NormSq(m)
//    m.norm1() or Norm1(m)
//    m.norm2() or Norm2(m)
//    m.normInf() or NormInf(m)
//    m.maxAbsElement() or MaxAbsElements(m)
//    m.maxAbs2Element() or MaxAbs2Elements(m)
//        (Note - for diagonal matrices, 
//        norm1 = norm2 = normInf = maxAbsElement.)
//
//    m.inverse() or Inverse(m)
//    m.invertSelf()
//    m.makeInverse(minv) (takes either Matrix or DiagMatrix argument)
//    m.makeInverseATA(invata) (takes either Matrix or DiagMatrix argument)
//
// I/O: 
//
//    os << d 
//        Writes d to ostream os as a full matrix
//
//    d.writeCompact(os)
//        Writes only the diagonal Vector to os
//
//    is >> d
//        Reads in d in the compact format
//
//

#ifndef TMV_DiagMatrix_H
#define TMV_DiagMatrix_H

#include "tmv/TMV_BaseDiagMatrix.h"
#include "tmv/TMV_BaseTriMatrix.h"
#include "tmv/TMV_Vector.h"
#include "tmv/TMV_TriMatrix.h"
#include "tmv/TMV_Matrix.h"
#include <vector>

namespace tmv {

    template <class T> 
    class GenDiagMatrix : 
        virtual public AssignableToDiagMatrix<T>,
        virtual public AssignableToUpperTriMatrix<T>,
        virtual public AssignableToLowerTriMatrix<T>,
        public BaseMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef T value_type;
        typedef RT real_type;
        typedef CT complex_type;
        typedef GenDiagMatrix<T> type;
        typedef DiagMatrix<T> copy_type;
        typedef ConstDiagMatrixView<T> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstDiagMatrixView<RT> const_realpart_type;
        typedef ConstVectorView<T> const_vec_type;
        typedef DiagMatrixView<T> nonconst_type;

        //
        // Constructors
        //

        inline GenDiagMatrix<T>() {}
        inline GenDiagMatrix(const type&) {}
        virtual inline ~GenDiagMatrix() {}

        //
        // Access Functions
        //

        inline size_t size() const { return diag().size(); }
        inline size_t colsize() const { return size(); }
        inline size_t rowsize() const { return size(); }
        inline DiagType dt() const { return NonUnitDiag; }

        inline T operator()(int i, int j) const 
        {
            TMVAssert(i>=0 && i<int(size()));
            TMVAssert(j>=0 && j<int(size()));
            if (i==j) return diag()(i); 
            else return T(0);
        }

        inline T operator()(int i) const 
        {
            TMVAssert(i>=0 && i<int(size()));
            return diag()(i); 
        }

        inline const_vec_type diag() const { return cdiag(); }

        template <class T2> 
        inline bool isSameAs(const BaseMatrix<T2>& ) const
        { return false; }

        inline bool isSameAs(const GenDiagMatrix<T>& m2) const
        { 
            if (this == &m2) return true;
            else return (diag().isSameAs(m2.diag())); 
        }

        template <class T2>
        TMV_DEPRECATED(bool SameAs(const BaseMatrix<T2>& m2) const);
        TMV_DEPRECATED(bool SameAs(const GenDiagMatrix<T>& m2) const)
        { return isSameAs(m2); }

        inline void assignToM(const MatrixView<RT>& m2) const
        {
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            TMVAssert(isReal(T()));
            m2.diag() = diag(); 
            m2.upperTri().offDiag().setZero();
            m2.lowerTri().offDiag().setZero();
        }
        inline void assignToM(const MatrixView<CT>& m2) const
        {
            TMVAssert(m2.colsize() == size());
            TMVAssert(m2.rowsize() == size());
            m2.diag() = diag(); 
            m2.upperTri().offDiag().setZero();
            m2.lowerTri().offDiag().setZero();
        }

        inline void assignToU(const UpperTriMatrixView<RT>& m2) const
        {
            TMVAssert(m2.size() == size());
            TMVAssert(isReal(T()));
            m2.diag() = diag(); 
            m2.offDiag().setZero();
        }
        inline void assignToU(const UpperTriMatrixView<CT>& m2) const
        {
            TMVAssert(m2.size() == size());
            m2.diag() = diag(); 
            m2.offDiag().setZero();
        }

        inline void assignToL(const LowerTriMatrixView<RT>& m2) const
        {
            TMVAssert(m2.size() == size());
            TMVAssert(isReal(T()));
            m2.diag() = diag(); 
            m2.offDiag().setZero();
        }
        inline void assignToL(const LowerTriMatrixView<CT>& m2) const
        {
            TMVAssert(m2.size() == size());
            m2.diag() = diag(); 
            m2.offDiag().setZero();
        }

        inline void assignToD(const DiagMatrixView<RT>& m2) const
        { 
            TMVAssert(m2.size() == size());
            TMVAssert(isReal(T()));
            m2.diag() = diag(); 
        }
        inline void assignToD(const DiagMatrixView<CT>& m2) const
        {
            TMVAssert(m2.size() == size());
            m2.diag() = diag(); 
        }

        //
        // subDiagMatrix
        //

        inline const_view_type cSubDiagMatrix(int i1, int i2) const
        { return const_view_type(diag().cSubVector(i1,i2)); }

        inline const_view_type subDiagMatrix(int i1, int i2) const
        {
            TMVAssert(diag().hasSubVector(i1,i2,1));
            return cSubDiagMatrix(i1,i2);
        }

        inline const_view_type cSubDiagMatrix(int i1, int i2, int istep) const
        { return const_view_type(diag().cSubVector(i1,i2,istep)); }

        inline const_view_type subDiagMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(diag().hasSubVector(i1,i2,istep));
            return cSubDiagMatrix(i1,i2,istep);
        }

        inline const_realpart_type realPart() const
        { return const_realpart_type(diag().realPart()); }

        inline const_realpart_type imagPart() const
        { return const_realpart_type(diag().imagPart()); }

        inline const_view_type view() const
        { return const_view_type(diag()); }

        inline const_view_type transpose() const
        { return view(); }

        inline const_view_type conjugate() const
        { return const_view_type(diag().conjugate()); }

        inline const_view_type adjoint() const
        { return conjugate(); }

        inline nonconst_type nonConst() const
        { return nonconst_type(diag().nonConst()); }


        TMV_DEPRECATED(const_view_type SubDiagMatrix(int i1, int i2) const)
        { return subDiagMatrix(i1,i2); }
        TMV_DEPRECATED(const_view_type SubDiagMatrix(
                int i1, int i2, int istep) const)
        { return subDiagMatrix(i1,i2,istep); }
        TMV_DEPRECATED(const_realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(const_realpart_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_view_type Adjoint() const)
        { return adjoint(); }
        TMV_DEPRECATED(nonconst_type NonConst() const)
        { return nonConst(); }


        //
        // Functions of DiagMatrix
        //

        T det() const;

        RT logDet(T* sign=0) const;

        inline T trace() const
        { return diag().sumElements(); }

        inline T sumElements() const
        { return diag().sumElements(); }

        inline RT sumAbsElements() const
        { return diag().sumAbsElements(); }

        inline RT norm() const 
        { return normF(); }

        inline RT normF() const 
        { return diag().norm(); }

        inline RT normSq(RT scale = RT(1)) const 
        { return diag().normSq(scale); }

        inline RT norm1() const 
        { return diag().maxAbsElement(); }

        inline RT doNorm2() const 
        { return diag().maxAbsElement(); }

        inline RT doCondition() const 
        { return diag().maxAbsElement()/diag().minAbsElement(); }

        inline RT norm2() const 
        { return doNorm2(); }

        inline RT condition() const 
        { return doCondition(); }

        inline RT normInf() const
        { return diag().maxAbsElement(); }

        inline RT maxAbsElement() const
        { return diag().maxAbsElement(); }

        inline RT maxAbs2Element() const
        { return diag().maxAbs2Element(); }

        inline bool isSingular() const 
        { return diag().minAbsElement() == RT(0); }

        template <class T1> 
        void doMakeInverse(const MatrixView<T1>& minv) const;

        template <class T1> 
        void doMakeInverse(const DiagMatrixView<T1>& minv) const;

        void doMakeInverseATA(const MatrixView<T>& ata) const;
        void doMakeInverseATA(const DiagMatrixView<T>& ata) const;

        inline void makeInverse(const MatrixView<T>& minv) const
        { 
            TMVAssert(minv.colsize() == size());
            TMVAssert(minv.rowsize() == size());
            doMakeInverse(minv);
        }

        template <class T1> 
        inline void makeInverse(const MatrixView<T1>& minv) const
        { 
            TMVAssert(minv.colsize() == size());
            TMVAssert(minv.rowsize() == size());
            doMakeInverse(minv);
        }

        template <class T1> 
        inline void makeInverse(const DiagMatrixView<T1>& minv) const
        { 
            TMVAssert(minv.size() == size());
            doMakeInverse(minv);
        }

        QuotXD<T,T> QInverse() const;
        inline QuotXD<T,T> inverse() const
        { return QInverse(); }

        inline void makeInverseATA(const MatrixView<T>& ata) const
        { 
            TMVAssert(ata.colsize() == size());
            TMVAssert(ata.rowsize() == size());
            doMakeInverseATA(ata);
        }
        inline void makeInverseATA(const DiagMatrixView<T>& ata) const
        { 
            TMVAssert(ata.size() == size());
            doMakeInverseATA(ata);
        }

        template <class T1, IndexStyle I> 
        inline void makeInverse(DiagMatrix<T1,I>& minv) const
        { makeInverse(minv.view()); }

        template <class T1, StorageType S, IndexStyle I> 
        inline void makeInverse(Matrix<T1,S,I>& minv) const
        { makeInverse(minv.view()); }

        template <IndexStyle I> 
        inline void makeInverseATA(DiagMatrix<T,I>& minv) const
        { makeInverseATA(minv.view()); }

        template <StorageType S, IndexStyle I> 
        inline void makeInverseATA(Matrix<T,S,I>& minv) const
        { makeInverseATA(minv.view()); }

        auto_ptr<BaseMatrix<T> > newTranspose() const;
        auto_ptr<BaseMatrix<T> > newConjugate() const;
        auto_ptr<BaseMatrix<T> > newAdjoint() const;
        auto_ptr<BaseMatrix<T> > newInverse() const;
        auto_ptr<BaseMatrix<T> > newView() const;
        auto_ptr<BaseMatrix<T> > newCopy() const;

        typedef QuotXD<T,T> MyQuotXD;
        TMV_DEPRECATED(MyQuotXD Inverse() const)
        { return inverse(); }
        template <class T1>
        TMV_DEPRECATED(void Inverse(const MatrixView<T1>& minv) const);
        template <class T1, StorageType S, IndexStyle I>
        TMV_DEPRECATED(void Inverse(Matrix<T1,S,I>& minv) const);
        template <class T1>
        TMV_DEPRECATED(void Inverse(const DiagMatrixView<T1>& minv) const);
        template <class T1, IndexStyle I> 
        TMV_DEPRECATED(void Inverse(DiagMatrix<T1,I>& minv) const);
        template <StorageType S, IndexStyle I>
        TMV_DEPRECATED(void InverseATA(Matrix<T,S,I>& ata) const);
        template <IndexStyle I>
        TMV_DEPRECATED(void InverseATA(DiagMatrix<T,I>& ata) const);


        // 
        // I/O
        //

        void write(std::ostream& os) const;
        void write(std::ostream& os, RT thresh) const;

        inline void writeCompact(std::ostream& os) const
        { os << "D " << diag() << std::endl; }
        inline void writeCompact(std::ostream& os, RT thresh) const
        { os << "D "; diag().write(os,thresh); os << std::endl; }

        TMV_DEPRECATED(void WriteCompact(std::ostream& os) const)
        { writeCompact(os); }
        TMV_DEPRECATED(void WriteCompact(std::ostream& os, RT thresh) const)
        { writeCompact(os,thresh); }

        // 
        // Arithmetic Helpers
        //

        template <class T1> 
        void doLDivEq(const VectorView<T1>& v) const;

        template <class T1, class T0> 
        void doLDiv(const GenVector<T1>& v1, const VectorView<T0>& v0) const;

        template <class T1> 
        void doLDivEq(const MatrixView<T1>& m) const;

        template <class T1, class T0> 
        void doLDiv(const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const;

        template <class T1> 
        inline void LDivEq(const VectorView<T1>& v) const
        {
            TMVAssert(v.size() == size());
            doLDivEq(v);
        }

        template <class T1, class T0> 
        inline void LDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const
        {
            TMVAssert(v1.size() == size());
            TMVAssert(v0.size() == size());
            doLDiv(v1,v0);
        }

        template <class T1> 
        inline void RDivEq(const VectorView<T1>& v) const
        { LDivEq(v); }

        template <class T1, class T0> 
        inline void RDiv(
            const GenVector<T1>& v1, const VectorView<T0>& v0) const
        { LDiv(v1,v0); }

        template <class T1> 
        inline void LDivEq(const MatrixView<T1>& m) const
        {
            TMVAssert(m.colsize() == size());
            doLDivEq(m);
        }
        template <class T1, class T0> 
        inline void LDiv(
            const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
        {
            TMVAssert(m1.colsize() == size());
            TMVAssert(m0.colsize() == size());
            TMVAssert(m1.rowsize() == m0.rowsize());
            doLDiv(m1,m0);
        }
        template <class T1> 
        inline void RDivEq(const MatrixView<T1>& m) const
        { LDivEq(m.transpose()); }

        template <class T1, class T0> 
        void RDiv(const GenMatrix<T1>& m1, const MatrixView<T0>& m0) const
        { LDiv(m1.transpose(),m0.transpose()); }

        template <class T1> 
        void divEq(const DiagMatrixView<T1>& m) const
        { LDivEq(m.diag()); }

        template <class T1, class T0> 
        void div(
            const GenDiagMatrix<T1>& m1, const DiagMatrixView<T0>& m0) const
        { LDiv(m1.diag(),m0.diag()); }

        // For easier compatibility with regular matrices:
        inline void divideInPlace() const {}
        inline void saveDiv() const {}
        inline void setDiv() const {}
        inline void unsetDiv() const {}
        inline void resetDiv() const {}
        inline bool divIsSet() const { return true; }
        inline void divideUsing(DivType TMV_DEBUGPARAM(dt)) const 
        { TMVAssert(dt == LU); }
        inline bool checkDecomp(std::ostream* fout=0) const { return true; }
        inline bool checkDecomp(
            const BaseMatrix<T>& m2, std::ostream* fout=0) const
        { return true; }

        TMV_DEPRECATED(void DivideInPlace() const) {}
        TMV_DEPRECATED(void SaveDiv() const) {}
        TMV_DEPRECATED(void SetDiv() const) {}
        TMV_DEPRECATED(void UnSetDiv() const) {}
        TMV_DEPRECATED(void ReSetDiv() const) {}
        TMV_DEPRECATED(bool DivIsSet() const) { return divIsSet(); }
        TMV_DEPRECATED(void DivideUsing(DivType dt) const) { divideUsing(dt); }
        TMV_DEPRECATED(bool CheckDecomp(std::ostream* fout=0) const)
        { return true; }
        TMV_DEPRECATED(bool CheckDecomp(
                const BaseMatrix<T>& m2, std::ostream* fout=0) const)
        { return true; }

        inline T cref(int i, int j) const
        { return i==j ? cdiag()(i) : 0; }
        inline const T* cptr() const
        { return cdiag().cptr(); }
        inline int step() const
        { return cdiag().step(); }

    protected :

        virtual ConstVectorView<T> cdiag() const = 0;

    private :

        type& operator=(const type&);

    }; // GenDiagMatrix

    template <class T>
    template <class T2>
    inline bool GenDiagMatrix<T>::SameAs(const BaseMatrix<T2>& m2) const
    { return isSameAs(m2); }

    template <class T>
    template <class T1>
    inline void GenDiagMatrix<T>::Inverse(const MatrixView<T1>& minv) const
    { makeInverse(minv); }

    template <class T>
    template <class T1, StorageType S, IndexStyle I>
    inline void GenDiagMatrix<T>::Inverse(Matrix<T1,S,I>& minv) const
    { makeInverse(minv); }

    template <class T>
    template <class T1>
    inline void GenDiagMatrix<T>::Inverse(const DiagMatrixView<T1>& minv) const
    { makeInverse(minv); }

    template <class T>
    template <class T1, IndexStyle I>
    inline void GenDiagMatrix<T>::Inverse(DiagMatrix<T1,I>& minv) const
    { makeInverse(minv); }

    template <class T>
    template <StorageType S, IndexStyle I>
    inline void GenDiagMatrix<T>::InverseATA(Matrix<T,S,I>& ata) const
    { makeInverseATA(ata); }

    template <class T>
    template <IndexStyle I>
    inline void GenDiagMatrix<T>::InverseATA(DiagMatrix<T,I>& ata) const
    { makeInverseATA(ata); }



    template <class T, IndexStyle I> 
    class ConstDiagMatrixView : public GenDiagMatrix<T>
    {
    public :

        typedef GenDiagMatrix<T> base;
        typedef ConstDiagMatrixView<T,I> type;
        typedef ConstVectorView<T> const_vec_type;

        inline ConstDiagMatrixView(const type& rhs) : itsdiag(rhs.cdiag()) {}

        inline ConstDiagMatrixView(const base& rhs) : itsdiag(rhs.diag()) {}

        explicit inline ConstDiagMatrixView(const GenVector<T>& v) :
            itsdiag(v) {}

        virtual inline ~ConstDiagMatrixView() {}

    protected :

        ConstVectorView<T> itsdiag;
        inline ConstVectorView<T> cdiag() const { return itsdiag; }

    private :

        type& operator=(const type&);

    }; // ConstDiagMatrixView

    template <class T> 
    class ConstDiagMatrixView<T,FortranStyle> : 
        public ConstDiagMatrixView<T,CStyle>
    {
    public :

        typedef TMV_RealType(T) RT;
        typedef GenDiagMatrix<T> base;
        typedef ConstDiagMatrixView<T,FortranStyle> type;
        typedef ConstVectorView<T,FortranStyle> const_vec_type;
        typedef ConstDiagMatrixView<T,FortranStyle> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstDiagMatrixView<RT,FortranStyle> const_realpart_type;
        typedef ConstDiagMatrixView<T,CStyle> c_type;
        typedef DiagMatrixView<T,FortranStyle> nonconst_type;

        inline ConstDiagMatrixView(const type& rhs) : c_type(rhs) {}

        inline ConstDiagMatrixView(const base& rhs) : c_type(rhs) {}

        explicit inline ConstDiagMatrixView(const GenVector<T>& v) :
            c_type(v) {}

        virtual inline ~ConstDiagMatrixView() {}

        //
        // Access Functions
        //

        inline T operator()(int i, int j) const 
        { 
            TMVAssert(i>0 && i<=int(size()));
            TMVAssert(j>0 && j<=int(size()));
            if (i==j) return diag()(i);
            else return T(0);
        }

        inline T operator()(int i) const 
        { 
            TMVAssert(i>0 && i<=int(size()));
            return diag()(i);
        }

        inline const_vec_type diag() const 
        { return const_vec_type(c_type::diag()); }

        //
        // subDiagMatrix
        //

        inline const_view_type subDiagMatrix(int i1, int i2) const
        {
            TMVAssert(diag().hasSubVector(i1,i2,1));
            return base::cSubDiagMatrix(i1-1,i2);
        }

        inline const_view_type subDiagMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(diag().hasSubVector(i1,i2,istep));
            return base::cSubDiagMatrix(i1-1,i2-1+istep,istep);
        }

        inline const_realpart_type realPart() const
        { return const_realpart_type(diag().realPart()); }

        inline const_realpart_type imagPart() const
        { return const_realpart_type(diag().imagPart()); }

        inline const_view_type view() const
        { return const_view_type(diag()); }

        inline const_view_type transpose() const
        { return view(); }

        inline const_view_type conjugate() const
        { return const_view_type(diag().conjugate()); }

        inline const_view_type adjoint() const
        { return conjugate(); }

        inline nonconst_type nonConst() const
        { return nonconst_type(diag().nonConst()); }

        TMV_DEPRECATED(const_view_type SubDiagMatrix(int i1, int i2) const)
        { return subDiagMatrix(i1,i2); }
        TMV_DEPRECATED(const_view_type SubDiagMatrix(
                int i1, int i2, int istep) const)
        { return subDiagMatrix(i1,i2,istep); }
        TMV_DEPRECATED(const_realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(const_realpart_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_view_type Adjoint() const)
        { return adjoint(); }
        TMV_DEPRECATED(nonconst_type NonConst() const)
        { return nonConst(); }

        using c_type::size;

        type& operator=(const type&);

    }; // ConstDiagMatrixView - FortranStyle

    template <class T, IndexStyle I> 
    class DiagMatrixView : public GenDiagMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenDiagMatrix<T> base;
        typedef DiagMatrixView<T,I> type;
        typedef DiagMatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef DiagMatrixView<RT,I> realpart_type;
        typedef VectorView<T,I> vec_type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef TMV_RefType(T) reference;

        //
        // Constructors
        //

        inline DiagMatrixView(const type& rhs) : itsdiag(rhs.diag()) {}

        explicit inline DiagMatrixView(const vec_type& _diag) :
            itsdiag(_diag) {}

        virtual inline ~DiagMatrixView() {} 

        //
        // Op=
        //

        inline const type& operator=(const type& m2) const
        { m2.assignToD(*this); return *this; }

        inline const type& operator=(const type& m2) 
        { m2.assignToD(*this); return *this; }

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { m2.assignToD(*this); return *this; }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { m2.assignToD(*this); return *this; }

        template <class T2> 
        inline const type& operator=(const GenDiagMatrix<T2>& m2) const
        { itsdiag = m2.diag(); return *this; }

        inline const type& operator=(const T& x) const 
        { return setToIdentity(x); }

        inline const type& operator=(const AssignableToDiagMatrix<RT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            m2.assignToD(*this);
            return *this;
        }

        inline const type& operator=(const AssignableToDiagMatrix<CT>& m2) const
        { 
            TMVAssert(size() == m2.size());
            m2.assignToD(*this);
            return *this;
        }


        //
        // Access
        //

        inline reference operator()(int i) const 
        { 
            TMVAssert(i>=0 && i<int(size()));
            return diag()(i); 
        }
        inline reference operator()(int i, int TMV_DEBUGPARAM(j)) const 
        { 
            TMVAssert(i>=0 && i<int(size()));
            TMVAssert(i==j); 
            return diag()(i); 
        }

        inline vec_type diag() const { return itsdiag; }

        //
        // Modifying Functions
        //

        inline const type& setZero() const 
        { diag().setZero(); return *this; }

        inline const type& setAllTo(const T& x) const
        { diag().setAllTo(x); return *this; }

        inline const type& clip(RT thresh) const
        { diag().clip(thresh); return *this; }

        inline const type& transposeSelf() const
        { return *this; }

        inline const type& conjugateSelf() const
        { diag().conjugateSelf(); return *this; }

        const type& invertSelf() const;

        inline const type& setToIdentity(const T& x=T(1)) const 
        { return setAllTo(x); }

        TMV_DEPRECATED(const type& Zero() const)
        { return setZero(); }
        TMV_DEPRECATED(const type& SetAllTo(const T& x) const)
        { return setAllTo(x); }
        TMV_DEPRECATED(const type& Clip(RT thresh) const)
        { return clip(thresh); }
        TMV_DEPRECATED(const type& TransposeSelf() const)
        { return transposeSelf(); }
        TMV_DEPRECATED(const type& ConjugateSelf() const)
        { return conjugateSelf(); }
        TMV_DEPRECATED(const type& InvertSelf() const)
        { return invertSelf(); }
        TMV_DEPRECATED(const type& SetToIdentity(const T& x=T(1)) const)
        { return setToIdentity(x); }


        //
        // subDiagMatrix
        //

        inline view_type cSubDiagMatrix(int i1, int i2) const
        { return view_type(itsdiag.cSubVector(i1,i2)); }

        inline view_type subDiagMatrix(int i1, int i2) const
        {
            TMVAssert(diag().hasSubVector(i1,i2,1));
            return cSubDiagMatrix(i1,i2);
        }

        inline view_type cSubDiagMatrix(int i1, int i2, int istep) const
        { return view_type(itsdiag.cSubVector(i1,i2,istep)); }

        inline view_type subDiagMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(diag().hasSubVector(i1,i2,istep));
            return cSubDiagMatrix(i1,i2,istep);
        }

        inline realpart_type realPart() const
        { return realpart_type(diag().realPart()); }

        inline realpart_type imagPart() const
        { return realpart_type(diag().imagPart()); }

        inline view_type view() const
        { return *this; }

        inline view_type transpose() const
        { return *this; }

        inline view_type conjugate() const
        { return view_type(diag().conjugate()); }

        inline view_type adjoint() const
        { return conjugate(); }

        TMV_DEPRECATED(view_type SubDiagMatrix(int i1, int i2) const)
        { return subDiagMatrix(i1,i2); }
        TMV_DEPRECATED(view_type SubDiagMatrix(
                int i1, int i2, int istep) const)
        { return subDiagMatrix(i1,i2,istep); }
        TMV_DEPRECATED(realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(realpart_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(view_type View() const)
        { return view(); }
        TMV_DEPRECATED(view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(view_type Adjoint() const)
        { return adjoint(); }


        using GenDiagMatrix<T>::size;

    protected:

        VectorView<T> itsdiag;
        inline ConstVectorView<T> cdiag() const { return itsdiag; }

    }; // DiagMatrixView

    template <class T> 
    class DiagMatrixView<T,FortranStyle> : public DiagMatrixView<T,CStyle>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenDiagMatrix<T> base;
        typedef DiagMatrixView<T,FortranStyle> type;
        typedef DiagMatrixView<T,CStyle> c_type;
        typedef DiagMatrixView<T,FortranStyle> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef DiagMatrixView<RT,FortranStyle> realpart_type;
        typedef VectorView<T,FortranStyle> vec_type;

        //
        // Constructors
        //

        inline DiagMatrixView(const type& rhs) : c_type(rhs) {}

        inline DiagMatrixView(const c_type& rhs) : c_type(rhs) {}

        explicit inline DiagMatrixView(const VectorView<T>& _diag) :
            c_type(_diag) {}

        virtual inline ~DiagMatrixView() {} 

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

        inline const type& operator=(const GenDiagMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const GenDiagMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        template <class T2> 
        inline const type& operator=(const GenDiagMatrix<T2>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const T& x) const 
        { c_type::operator=(x); return *this; }

        inline const type& operator=(const AssignableToDiagMatrix<RT>& m2) const
        { c_type::operator=(m2); return *this; }

        inline const type& operator=(const AssignableToDiagMatrix<CT>& m2) const
        { c_type::operator=(m2); return *this; }

        //
        // Access
        //

        inline TMV_RefType(T) operator()(int i) const 
        { 
            TMVAssert(i>0 && i<=int(size()));
            return diag()(i); 
        }
        inline TMV_RefType(T) operator()(int i, int TMV_DEBUGPARAM(j)) const 
        { 
            TMVAssert(i==j); 
            TMVAssert(i>0 && i<=int(size()));
            return diag()(i); 
        }

        inline vec_type diag() const 
        { return c_type::diag(); }

        //
        // Modifying Functions
        //

        inline const type& setZero() const 
        { diag().setZero(); return *this; }

        inline const type& setAllTo(const T& x) const
        { diag().setAllTo(x); return *this; }

        inline const type& clip(RT thresh) const
        { diag().clip(thresh); return *this; }

        inline const type& transposeSelf() const
        { return *this; }

        inline const type& conjugateSelf() const
        { diag().conjugateSelf(); return *this; }

        inline const type& invertSelf() const
        { return c_type::invertSelf(); }

        inline const type& setToIdentity(const T& x=T(1)) const 
        { diag().setAllTo(x); return *this; }

        TMV_DEPRECATED(const type& Zero() const)
        { return setZero(); }
        TMV_DEPRECATED(const type& SetAllTo(const T& x) const)
        { return setAllTo(x); }
        TMV_DEPRECATED(const type& Clip(RT thresh) const)
        { return clip(thresh); }
        TMV_DEPRECATED(const type& TransposeSelf() const)
        { return transposeSelf(); }
        TMV_DEPRECATED(const type& ConjugateSelf() const)
        { return conjugateSelf(); }
        TMV_DEPRECATED(const type& InvertSelf() const)
        { return invertSelf(); }
        TMV_DEPRECATED(const type& SetToIdentity(const T& x=T(1)) const)
        { return setToIdentity(x); }


        //
        // subDiagMatrix
        //

        inline view_type cSubDiagMatrix(int i1, int i2) const
        { return view_type(diag().cSubVector(i1,i2)); }

        inline view_type subDiagMatrix(int i1, int i2) const
        {
            TMVAssert(diag().hasSubVector(i1,i2,1));
            return view_type(diag().subVector(i1,i2)); 
        }

        inline view_type cSubDiagMatrix(int i1, int i2, int istep) const
        { return view_type(diag().cSubVector(i1,i2,istep)); }

        inline view_type subDiagMatrix(int i1, int i2, int istep) const
        {
            TMVAssert(diag().hasSubVector(i1,i2,istep));
            return view_type(diag().subVector(i1,i2,istep)); 
        }

        inline realpart_type realPart() const
        { return realpart_type(diag().realPart()); }

        inline realpart_type imagPart() const
        { return realpart_type(diag().imagPart()); }

        inline view_type view() const
        { return *this; }

        inline view_type transpose() const
        { return *this; }

        inline view_type conjugate() const
        { return view_type(diag().conjugate()); }

        inline view_type adjoint() const
        { return conjugate(); }

        TMV_DEPRECATED(view_type SubDiagMatrix(int i1, int i2) const)
        { return subDiagMatrix(i1,i2); }
        TMV_DEPRECATED(view_type SubDiagMatrix(
                int i1, int i2, int istep) const)
        { return subDiagMatrix(i1,i2,istep); }
        TMV_DEPRECATED(realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(realpart_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(view_type View() const)
        { return view(); }
        TMV_DEPRECATED(view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(view_type Adjoint() const)
        { return adjoint(); }

        using GenDiagMatrix<T>::size;

    }; // FortranStyle DiagMatrixView


    template <class T, IndexStyle I> 
    class DiagMatrix : public GenDiagMatrix<T>
    {
    public:

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;
        typedef GenDiagMatrix<T> base;
        typedef DiagMatrix<T,I> type;
        typedef DiagMatrixView<T,I> view_type;
        typedef view_type transpose_type;
        typedef view_type conjugate_type;
        typedef view_type adjoint_type;
        typedef DiagMatrixView<RT,I> realpart_type;
        typedef VectorView<T,I> vec_type;
        typedef ConstDiagMatrixView<T,I> const_view_type;
        typedef const_view_type const_transpose_type;
        typedef const_view_type const_conjugate_type;
        typedef const_view_type const_adjoint_type;
        typedef ConstDiagMatrixView<RT,I> const_realpart_type;
        typedef ConstVectorView<T,I> const_vec_type;
        typedef T& reference;

        //
        // Constructors
        //

        inline explicit DiagMatrix(size_t _size) : itsdiag(_size) {}

        inline DiagMatrix(size_t _size, const T& x) : itsdiag(_size,x)  {}

        inline DiagMatrix(size_t _size, const T* vv) : itsdiag(_size,vv)  {}

        inline explicit DiagMatrix(const std::vector<T>& vv) : itsdiag(vv) {}

        inline explicit DiagMatrix(const GenVector<T>& rhs) : itsdiag(rhs) {}

        inline explicit DiagMatrix(const GenMatrix<T>& m) : itsdiag(m.diag()) {}

        inline DiagMatrix(const type& rhs) : itsdiag(rhs.diag()) {}

        template <IndexStyle I2> 
        inline DiagMatrix(const DiagMatrix<T,I2>& rhs) : itsdiag(rhs.diag()) {}

        inline DiagMatrix(const GenDiagMatrix<RT>& rhs) : itsdiag(rhs.size()) 
        { rhs.assignToD(view()); }

        inline DiagMatrix(const GenDiagMatrix<CT>& rhs) : itsdiag(rhs.size()) 
        { 
            TMVAssert(isComplex(T()));
            rhs.assignToD(view()); 
        }

        template <class T2> 
        inline DiagMatrix(const GenDiagMatrix<T2>& rhs) : itsdiag(rhs.diag()) {}

        inline DiagMatrix(const AssignableToDiagMatrix<RT>& m2) :
            itsdiag(m2.colsize())
        {
            TMVAssert(m2.colsize() == m2.rowsize());
            m2.assignToD(view()); 
        }

        inline DiagMatrix(const AssignableToDiagMatrix<CT>& m2) :
            itsdiag(m2.colsize())
        {
            TMVAssert(m2.colsize() == m2.rowsize());
            TMVAssert(isComplex(T()));
            m2.assignToD(view()); 
        }

        virtual inline ~DiagMatrix() {}


        //
        // Op=
        //

        inline type& operator=(const type& m2)
        {
            TMVAssert(m2.size() == size());
            m2.assignToD(view());
            return *this; 
        }

        template <IndexStyle I2>
        inline type& operator=(const DiagMatrix<T,I2>& m2)
        {
            TMVAssert(m2.size() == size());
            m2.assignToD(view());
            return *this; 
        }

        inline type& operator=(const GenDiagMatrix<RT>& m2)
        {
            TMVAssert(m2.size() == size());
            m2.assignToD(view()); 
            return *this; 
        }

        inline type& operator=(const GenDiagMatrix<CT>& m2)
        {
            TMVAssert(m2.size() == size());
            TMVAssert(isComplex(T()));
            m2.assignToD(view()); 
            return *this; 
        }

        template <class T2> 
        inline type& operator=(const GenDiagMatrix<T2>& m2)
        { 
            TMVAssert(m2.size() == size());
            view() = m2; 
            return *this; 
        }

        inline type& operator=(const T& x) { view() = x; return *this; }

        inline type& operator=(const AssignableToDiagMatrix<RT>& m2)
        { 
            TMVAssert(m2.size() == size());
            m2.assignToD(view());
            return *this; 
        }

        inline type& operator=(const AssignableToDiagMatrix<CT>& m2)
        { 
            TMVAssert(m2.size() == size());
            TMVAssert(isComplex(T()));
            m2.assignToD(view());
            return *this; 
        }


        //
        // Access
        //

        inline T& operator()(int i) 
        { 
            if (I==CStyle) { 
                TMVAssert(i>=0 && i<int(size()));
                return itsdiag(i); 
            } else { 
                TMVAssert(i>0 && i<=int(size()));
                return itsdiag(i-1); 
            }
        }

        inline T& operator()(int i, int TMV_DEBUGPARAM(j)) 
        { 
            if (I==CStyle) { 
                TMVAssert(i>=0 && i<int(size())); 
                TMVAssert(j>=0 && j<int(size())); 
            } else { 
                TMVAssert(i>0 && i<=int(size())); 
                TMVAssert(j>0 && j<=int(size())); 
            }
            TMVAssert(i==j);
            return operator()(i);
        }

        inline T operator()(int i) const 
        {
            if (I==CStyle) { 
                TMVAssert(i>=0 && i<int(size()));
                return itsdiag(i);
            } else { 
                TMVAssert(i>0 && i<=int(size()));
                return itsdiag(i-1);
            }
        }

        inline T operator()(int i,int j) const 
        {
            if (I==CStyle) { 
                TMVAssert(i>=0 && i<int(size())); 
                TMVAssert(j>=0 && j<int(size())); 
            } else { 
                TMVAssert(i>0 && i<=int(size())); 
                TMVAssert(j>0 && j<=int(size())); 
            }
            if (i==j) return operator()(i);
            else return T(0);
        }

        inline vec_type diag() { return itsdiag.view(); }

        inline const_vec_type diag() const { return itsdiag.view(); }

        //
        // Modifying Functions
        //

        inline type& setZero() 
        { itsdiag.setZero(); return *this; }

        inline type& setAllTo(const T& x) 
        { itsdiag.setAllTo(x); return *this; }

        inline type& clip(RT thresh)
        { diag().clip(thresh); return *this; }

        inline type& transposeSelf() 
        { return *this; }

        inline type& conjugateSelf() 
        { itsdiag.conjugateSelf(); return *this; }

        inline type& invertSelf()
        { view().invertSelf(); return *this; }

        inline type& setToIdentity(const T& x=T(1)) 
        { itsdiag.setAllTo(x); return *this; }

        TMV_DEPRECATED(type& Zero())
        { return setZero(); }
        TMV_DEPRECATED(type& SetAllTo(const T& x))
        { return setAllTo(x); }
        TMV_DEPRECATED(type& Clip(RT thresh))
        { return clip(thresh); }
        TMV_DEPRECATED(type& TransposeSelf())
        { return transposeSelf(); }
        TMV_DEPRECATED(type& ConjugateSelf())
        { return conjugateSelf(); }
        TMV_DEPRECATED(type& InvertSelf())
        { return invertSelf(); }
        TMV_DEPRECATED(type& SetToIdentity(const T& x=T(1)))
        { return setToIdentity(x); }


        //
        // subDiagMatrix
        //

        inline view_type cSubDiagMatrix(int i1, int i2) 
        { return view_type(itsdiag.cSubVector(i1,i2)); }

        inline view_type subDiagMatrix(int i1, int i2) 
        { 
            TMVAssert(diag().hasSubVector(i1,i2,1));
            return view_type(itsdiag.subVector(i1,i2)); 
        }

        inline view_type cSubDiagMatrix(int i1, int i2, int istep) 
        { return view_type(itsdiag.cSubVector(i1,i2,istep)); }

        inline view_type subDiagMatrix(int i1, int i2, int istep) 
        { 
            TMVAssert(diag().hasSubVector(i1,i2,istep));
            return view_type(itsdiag.subVector(i1,i2,istep)); 
        }

        inline realpart_type realPart()
        { return realpart_type(diag().realPart()); }

        inline realpart_type imagPart()
        { return realpart_type(diag().imagPart()); }

        inline const_view_type cSubDiagMatrix(int i1, int i2) const
        { return const_view_type(itsdiag.cSubVector(i1,i2)); }

        inline const_view_type subDiagMatrix(int i1, int i2) const
        {
            TMVAssert(diag().hasSubVector(i1,i2,1));
            return const_view_type(itsdiag.subVector(i1,i2)); 
        }

        inline const_view_type cSubDiagMatrix(int i1, int i2, int istep) const
        { return const_view_type(itsdiag.cSubVector(i1,i2,istep)); }

        inline const_view_type subDiagMatrix(int i1, int i2, int istep) const
        { 
            TMVAssert(diag().hasSubVector(i1,i2,istep));
            return const_view_type(itsdiag.subVector(i1,i2,istep)); 
        }

        inline const_realpart_type realPart() const
        { return const_realpart_type(diag().realPart()); }

        inline const_realpart_type imagPart() const
        { return const_realpart_type(diag().imagPart()); }

        inline const_view_type view() const
        { return const_view_type(itsdiag); }

        inline const_view_type transpose() const
        { return view(); }

        inline const_view_type conjugate() const
        { return const_view_type(itsdiag.conjugate()); }

        inline const_view_type adjoint() const
        { return conjugate(); }

        inline view_type view() 
        { return view_type(itsdiag.view()); }

        inline view_type transpose() 
        { return view(); }

        inline view_type conjugate()
        { return view_type(itsdiag.conjugate()); }

        inline view_type adjoint() 
        { return conjugate(); }

        TMV_DEPRECATED(const_view_type SubDiagMatrix(int i1, int i2) const)
        { return subDiagMatrix(i1,i2); }
        TMV_DEPRECATED(const_view_type SubDiagMatrix(
                int i1, int i2, int istep) const)
        { return subDiagMatrix(i1,i2,istep); }
        TMV_DEPRECATED(const_realpart_type Real() const)
        { return realPart(); }
        TMV_DEPRECATED(const_realpart_type Imag() const)
        { return imagPart(); }
        TMV_DEPRECATED(const_view_type View() const)
        { return view(); }
        TMV_DEPRECATED(const_view_type Transpose() const)
        { return transpose(); }
        TMV_DEPRECATED(const_view_type Conjugate() const)
        { return conjugate(); }
        TMV_DEPRECATED(const_view_type Adjoint() const)
        { return adjoint(); }

        TMV_DEPRECATED(view_type SubDiagMatrix(int i1, int i2))
        { return subDiagMatrix(i1,i2); }
        TMV_DEPRECATED(view_type SubDiagMatrix(int i1, int i2, int istep))
        { return subDiagMatrix(i1,i2,istep); }
        TMV_DEPRECATED(realpart_type Real())
        { return realPart(); }
        TMV_DEPRECATED(realpart_type Imag())
        { return imagPart(); }
        TMV_DEPRECATED(view_type View())
        { return view(); }
        TMV_DEPRECATED(view_type Transpose())
        { return transpose(); }
        TMV_DEPRECATED(view_type Conjugate())
        { return conjugate(); }
        TMV_DEPRECATED(view_type Adjoint())
        { return adjoint(); }

        using base::size;

    protected :

        Vector<T> itsdiag;
        inline ConstVectorView<T> cdiag() const { return itsdiag.view(); }

        template <IndexStyle I2>
        friend void Swap(DiagMatrix<T,I>& m1, DiagMatrix<T,I2>& m2)
        {
            TMVAssert(m1.size() == m2.size());
            Swap(m1.itsdiag,m2.itsdiag);
        }

    }; // DiagMatrix

    //---------------------------------------------------------------------------

    //
    // Special Creators:
    //   DiagMatrixViewOf(v)
    //   DiagMatrixViewOf(T* v, n)
    //

    template <class T> 
    inline ConstDiagMatrixView<T> DiagMatrixViewOf(const GenVector<T>& v)
    { return ConstDiagMatrixView<T>(v); }

    template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> DiagMatrixViewOf(
        const ConstVectorView<T,I>& v)
    { return ConstDiagMatrixView<T,I>(v); }

    template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> DiagMatrixViewOf(const Vector<T,I>& v)
    { return ConstDiagMatrixView<T,I>(v.view()); }

    template <class T, IndexStyle I> 
    inline DiagMatrixView<T,I> DiagMatrixViewOf(const VectorView<T,I>& v)
    { return DiagMatrixView<T,I>(v); }

    template <class T, IndexStyle I> 
    inline DiagMatrixView<T,I> DiagMatrixViewOf(Vector<T,I>& v)
    { return DiagMatrixView<T,I>(v.view()); }

    template <class T> 
    inline ConstDiagMatrixView<T> DiagMatrixViewOf(const T* m, size_t size)
    { return ConstDiagMatrixView<T>(VectorViewOf(m,size)); }

    template <class T> 
    inline DiagMatrixView<T> DiagMatrixViewOf(T* m, size_t size)
    { return DiagMatrixView<T>(VectorViewOf(m,size)); }

    //
    // Swap Matrices
    //

    template <class T> 
    inline void Swap(const DiagMatrixView<T>& m1, const DiagMatrixView<T>& m2)
    { Swap(m1.diag(),m2.diag()); }

    template <class T, IndexStyle I> 
    inline void Swap(const DiagMatrix<T,I>& m1, const DiagMatrixView<T>& m2)
    { Swap(m1.diag(),m2.diag()); }

    template <class T, IndexStyle I> 
    inline void Swap(const DiagMatrixView<T>& m1, const DiagMatrix<T,I>& m2)
    { Swap(m1.diag(),m2.diag()); }

    //
    // Views:
    //

    template <class T> 
    inline ConstDiagMatrixView<T> Transpose(const GenDiagMatrix<T>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> Transpose(const ConstDiagMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> Transpose(const DiagMatrix<T,I>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline DiagMatrixView<T,I> Transpose(const DiagMatrixView<T,I>& m)
    { return m.transpose(); }

    template <class T, IndexStyle I> 
    inline DiagMatrixView<T,I> Transpose(DiagMatrix<T,I>& m)
    { return m.transpose(); }

    template <class T> 
    inline ConstDiagMatrixView<T> Conjugate(const GenDiagMatrix<T>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> Conjugate(const ConstDiagMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> Conjugate(const DiagMatrix<T,I>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline DiagMatrixView<T,I> Conjugate(const DiagMatrixView<T,I>& m)
    { return m.conjugate(); }

    template <class T, IndexStyle I> 
    inline DiagMatrixView<T,I> Conjugate(DiagMatrix<T,I>& m)
    { return m.conjugate(); }

    template <class T> 
    inline ConstDiagMatrixView<T> Adjoint(const GenDiagMatrix<T>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> Adjoint(const ConstDiagMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline ConstDiagMatrixView<T,I> Adjoint(const DiagMatrix<T,I>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline DiagMatrixView<T,I> Adjoint(const DiagMatrixView<T,I>& m)
    { return m.adjoint(); }

    template <class T, IndexStyle I> 
    inline DiagMatrixView<T,I> Adjoint(DiagMatrix<T,I>& m)
    { return m.adjoint(); }

    template <class T> 
    inline QuotXD<T,T> Inverse(const GenDiagMatrix<T>& m)
    { return m.inverse(); }



    //
    // DiagMatrix ==, != DiagMatrix
    //

    template <class T1, class T2> 
    inline bool operator==(
        const GenDiagMatrix<T1>& m1, const GenDiagMatrix<T2>& m2)
    { return m1.diag() == m2.diag(); }

    template <class T1, class T2> 
    inline bool operator!=(
        const GenDiagMatrix<T1>& m1, const GenDiagMatrix<T2>& m2)
    { return !(m1 == m2); }


    //
    // I/O
    //

    template <class T, IndexStyle I> 
    std::istream& operator>>(std::istream& fin, auto_ptr<DiagMatrix<T,I> >& m);

    template <class T> 
    std::istream& operator>>(std::istream& fin, const DiagMatrixView<T>& m);

    template <class T, IndexStyle I> 
    inline std::istream& operator>>(std::istream& fin, DiagMatrix<T,I>& m)
    { return fin >> m.view(); }

} // namespace tmv

#endif
