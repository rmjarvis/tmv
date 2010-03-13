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


#ifndef TMV_Mat_H
#define TMV_Mat_H

#include "tmv/TMV_Matrix.h"
#include "tmv/TMV_SmallMatrix.h"

#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_CopyM.h"
#include "tmv/TMV_SwapM.h"
#include "tmv/TMV_TransposeM.h"
#include "tmv/TMV_PermuteM.h"
#include "tmv/TMV_NormM.h"

#include "tmv/TMV_ScaleM.h"
#include "tmv/TMV_MultXM.h"
#include "tmv/TMV_AddMM.h"
#include "tmv/TMV_Rank1VVM.h"
#include "tmv/TMV_MultMV.h"
#include "tmv/TMV_MultMM.h"

#include "tmv/TMV_Det.h"
#include "tmv/TMV_InvertM.h"

#ifndef TMV_H
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_SumMX.h"
#include "tmv/TMV_OProdVV.h"
#include "tmv/TMV_ProdMV.h"
#include "tmv/TMV_ProdMM.h"
#endif

#include "tmv/TMV_LUDecompose.h"
#if 0
#include "tmv/TMV_LUD.h"
#include "tmv/TMV_QRD.h"
#include "tmv/TMV_QRPD.h"
#include "tmv/TMV_SVD.h"
#include "tmv/TMV_PackedQ.h"
#endif

#endif
