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

#include "TMV_Blas.h"
#include "tmv/TMV_AddVV.h"
#include "tmv/TMV_MultXV.h"
#include "tmv/TMV_Vector.h"

namespace tmv {

    // Just do the most common cases: v+v, v-v
    // For anything else, do it in two function calls to MultXV.
    template <class V1, class V2, class T3>
    static void DoAddVV(
        const T3 x1, const V1& v1, const T3 x2, const V2& v2,
        VectorView<T3,Unit> v3)
    {
        typedef typename Traits<T3>::real_type RT;
        if (x1 == RT(1)) {
            if (x2 == RT(1)) {
                InlineAddVV(Scaling<1,RT>(),v1,Scaling<1,RT>(),v2,v3);
            } else if (x2 == RT(-1)) {
                InlineAddVV(Scaling<1,RT>(),v1,Scaling<-1,RT>(),v2,v3);
            } else if (x2 == RT(0)) {
                InstMultXV(x1,v1.xView(),v3.xView());
            } else {
                InstMultXV(x1,v1.xView(),v3.xView());
                InstAddMultXV(x2,v2.xView(),v3.xView());
            }
        } else if (x1 == RT(0)) {
            InstMultXV(x2,v2.xView(),v3.xView());
        } else {
            InstMultXV(x1,v1.xView(),v3.xView());
            InstAddMultXV(x2,v2.xView(),v3.xView());
        }
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAddVV(
        const T3 x1, const ConstVectorView<T1,C1>& v1,
        const T3 x2, const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    {
#ifdef BLAS
        InstMultXV(x1,v1,v3);
        InstAddMultXV(x2,v2,v3);
#else
        if (v1.step() == 1 && v2.step() == 1 && v3.step() == 1) 
            DoAddVV(x1,v1.unitView(),x2,v2.unitView(),v3.unitView());
        else if (x1 == T3(0)) 
            InstMultXV(x2,v2,v3);
        else if (x2 == T3(0)) 
            InstMultXV(x1,v1,v3);
        else {
            InstMultXV(x1,v1,v3);
            InstAddMultXV(x2,v2,v3);
        }
#endif
    }

    template <class T1, int C1, class T2, int C2, class T3>
    void InstAliasAddVV(
        const T3 x1, const ConstVectorView<T1,C1>& v1,
        const T3 x2, const ConstVectorView<T2,C2>& v2, VectorView<T3> v3)
    { InlineAliasAddVV(Scaling<0,T3>(x1),v1,Scaling<0,T3>(x2),v2,v3); }

#define InstFile "TMV_AddVV.inst"
#include "TMV_Inst.h"
#undef InstFile

} // namespace tmv


