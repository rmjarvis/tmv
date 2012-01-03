

#ifndef TMV_DIAG_H
#define TMV_DIAG_H

#include "tmv/TMV_BaseMatrix_Diag.h"
#include "tmv/TMV_DiagMatrix.h"
#include "tmv/TMV_SmallDiagMatrix.h"
#include "tmv/TMV_DiagMatrixIO.h"
#include "tmv/TMV_AddDD.h"
#include "tmv/TMV_MultXD.h"
#include "tmv/TMV_MultDV.h"
#include "tmv/TMV_MultMD.h"
#include "tmv/TMV_InvertD.h"
#include "tmv/TMV_DivVD.h"

#ifndef TMV_H
// In case TMV_Diag.h is included not as part of TMV.h:
#include "tmv/TMV_MatrixIO.h"
#include "tmv/TMV_ProdMM.h"
#include "tmv/TMV_ProdMV.h"
#include "tmv/TMV_ProdXM.h"
#include "tmv/TMV_SumMM.h"
#include "tmv/TMV_SumMX.h"
#include "tmv/TMV_ElemProdMM.h"
#endif

#endif
