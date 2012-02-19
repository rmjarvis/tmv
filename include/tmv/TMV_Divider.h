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
// This file defines the TMV Divider class.
//
// There are currently 4 algorithms for doing division (and Inverse and Det)
//
// LU Decomposition
// QR Decomposition (with or without Permutation)
// Singular Value Decomposition (compact or full)
// Cholskey (only for SymMatrix)
//
// To tell a Matrix to use a particular algorithm, use the command:
// m.divideUsing(ALG)
// where ALG is LU, QR, QRP, SV or CH  for the algorithms above.
//
// The default algorithm is LU for square matrices or QR for non-square.
//
// By default, the appropriate Divider class is created the first
// time it is needed (eg. when the statement v = b/m is called).
// However, you can also setup the Divider class beforehand manually
// by calling m.setDiv().
//
// You can also query whether the Divider class is already set up.
// This will only be true, if it was previously set up, _and_ the 
// Matrix hasn't been modified since then.
//
// If you want access to the various Divider functions directly,
// They can be accessed by:
//
// m.lud()
// m.qrd()
// m.svd()
// m.chd()
//
// The one of these that is probably most useful to access is svd(),
// since it is generally a good idea to look for small 
// singular values and zero them out before using SVD for division.
//
// To set to zero all singular value which are less than thresh * 
// the largest singular value use:
//
// m.svd()->setThresh(thresh);
//
// To use only the largest nsv singular values use:
//
// m.svd()->setTop(nsv);
//
// Also, the singular value decomposition can be used for principal
// component analysis of a Matrix.  The principal component vectors
// are the rows of V.  You can access the decomposition using:
//
// m.svd()->getU();
// m.svd()->getS(); // A DiagMatrix
// m.svd()->getV();
//
//


#ifndef TMV_Divider_H
#define TMV_Divider_H

#include "tmv/TMV_BaseMatrix.h"

namespace tmv {

    template <class T> 
    class Divider 
    {

        typedef TMV_RealType(T) RT;
        typedef TMV_ComplexType(T) CT;

    public :

        Divider() {}
        virtual ~Divider() {}

        virtual inline bool isSV() const { return false; }

        virtual T det() const =0;
        virtual RT logDet(T* sign) const =0;
        virtual void makeInverseATA(MatrixView<T> minv) const =0;
        virtual bool isSingular() const =0;
        virtual inline RT norm2() const 
        { TMVAssert(TMV_FALSE); return RT(0); }
        virtual inline RT condition() const 
        { TMVAssert(TMV_FALSE); return RT(0); }

#define DefDivEq(T) \
        virtual void LDivEq(MatrixView<T>) const =0; \
        virtual void RDivEq(MatrixView<T>) const =0; \
        virtual void makeInverse(MatrixView<T> minv) const =0 

        DefDivEq(RT);
        DefDivEq(CT);
#undef DefDivEq

#define DefDiv(T1,T2) \
        virtual void LDiv(const GenMatrix<T1>& b, MatrixView<T2> x) const =0; \
        virtual void RDiv(const GenMatrix<T1>& b, MatrixView<T2> x) const =0 

        DefDiv(RT,RT);
        DefDiv(RT,CT);
        DefDiv(CT,CT);
#undef DefDiv

        virtual bool checkDecomp(
            const BaseMatrix<T>& m, std::ostream* fout) const=0;
    };

    template <class T> 
    inline std::string TMV_Text(const Divider<T>& d)
    { return std::string("Divider<")+TMV_Text(T())+">"; }

} // namespace tmv

#endif
