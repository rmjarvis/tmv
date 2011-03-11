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

#include "tmv/TMV_ScaleU.h"
#include "tmv/TMV_TriMatrix.h"

namespace tmv {

    template <class T, class M>
    static void DoScale(const T x, M& m)
    {
        if (x == T(-1)) {
            InlineScale(Scaling<-1,T>(x),m); 
        } else if (x == T(0)) {
            m.setZero();
        } else if (x != T(1)) {
            InlineScale(Scaling<0,T>(x),m); 
        }
    }

    template <class T, class M>
    static void DoScale(const std::complex<T> x, M& m)
    { 
        if (imag(x) == T(0)) {
            if (real(x) == T(-1)) {
                InlineScale(Scaling<-1,T>(real(x)),m); 
            } else if (real(x) == T(0)) {
                m.setZero();
            } else if (real(x) != T(1)) {
                InlineScale(Scaling<0,T>(real(x)),m); 
            }
        } else  {
            InlineScale(Scaling<0,std::complex<T> >(x),m);
        }
    }

    template <class T>
    void InstScale(const T x, UpperTriMatrixView<T,NonUnitDiag> m)
    {
        if (m.iscm()) {
            UpperTriMatrixView<T,NonUnitDiag,1> mcm = m.cmView();
            DoScale(x,mcm); 
        } else if (m.isrm()) {
            UpperTriMatrixView<T,NonUnitDiag,UNKNOWN,1> mrm = m.rmView();
            DoScale(x,mrm); 
        } else {
            DoScale(x,m); 
        }
    }


#define InstFile "TMV_ScaleU.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


