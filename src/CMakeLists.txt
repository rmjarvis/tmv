PROJECT(TMV)

SET(BASIC TMV_Vector.cpp TMV_MultVV.cpp TMV_AddVV.cpp TMV_MultXV.cpp TMV_BaseMatrix.cpp TMV_Matrix.cpp TMV_MultXM.cpp TMV_AddMM.cpp TMV_MultMV.cpp TMV_Rank1_VVM.cpp TMV_MultMM.cpp TMV_MultMM_CCC.cpp TMV_MultMM_CRC.cpp TMV_MultMM_RCC.cpp TMV_Givens.cpp TMV_Householder.cpp TMV_LUD.cpp TMV_LUDecompose.cpp TMV_LUDiv.cpp TMV_LUInverse.cpp TMV_QRD.cpp TMV_QRDecompose.cpp TMV_PackedQ.cpp TMV_QRDiv.cpp TMV_QRInverse.cpp TMV_GetQFromQR.cpp TMV_QRUpdate.cpp TMV_QRDowndate.cpp TMV_QRPD.cpp TMV_QRPDecompose.cpp TMV_SVD.cpp TMV_SVDecompose.cpp TMV_SVDecompose_Bidiag.cpp TMV_SVDecompose_QR.cpp TMV_SVDecompose_DC.cpp TMV_SVDiv.cpp mmgr.cpp)

SET(DIAG TMV_DiagMatrix.cpp TMV_MultDV.cpp TMV_AddDM.cpp TMV_MultDM.cpp)

SET(TRI TMV_TriMatrix.cpp TMV_TriDiv.cpp TMV_TriDiv_V.cpp TMV_TriDiv_M.cpp TMV_TriDiv_L.cpp TMV_TriInverse.cpp TMV_MultUV.cpp TMV_MultXU.cpp TMV_AddUU.cpp TMV_MultUM.cpp TMV_MultUL.cpp TMV_MultUU.cpp TMV_MultDU.cpp)

SET(BAND TMV_BandMatrix.cpp TMV_BandLUD.cpp TMV_BandLUDecompose.cpp TMV_BandLUDiv.cpp TMV_BandLUInverse.cpp TMV_BandTriDiv_V.cpp TMV_BandTriDiv_M.cpp TMV_BandQRD.cpp TMV_BandQRDecompose.cpp TMV_BandQRDiv.cpp TMV_BandQRInverse.cpp TMV_BandSVD.cpp TMV_BandSVDecompose.cpp TMV_MultBV.cpp TMV_MultXB.cpp TMV_AddBB.cpp TMV_MultBM.cpp TMV_MultBB.cpp)

SET(SYM TMV_SymMatrix.cpp TMV_MultSV.cpp TMV_AddSS.cpp TMV_Rank1_VVS.cpp TMV_Rank2_VVS.cpp TMV_MultSM.cpp TMV_RankK_MMS.cpp TMV_Rank2K_MMS.cpp TMV_RankK_ULS.cpp TMV_RankK_LUS.cpp TMV_SymMultMMS.cpp TMV_SymHouseholder.cpp TMV_SymLDLD.cpp TMV_SymLDLDecompose.cpp TMV_SymLDLDiv.cpp TMV_SymLDLInverse.cpp TMV_SymLDLPseudo.cpp TMV_SymSquare.cpp TMV_SymCHD.cpp TMV_SymCHDecompose.cpp TMV_SymCHDiv.cpp TMV_SymCHInverse.cpp TMV_SymSVD.cpp TMV_SymSVInverse.cpp TMV_SymSVDecompose.cpp TMV_IsNaN.cpp TMV_SymSVDecompose_Tridiag.cpp TMV_SymSVDecompose_QR.cpp TMV_SymSVDecompose_DC.cpp)

SET(SYMBAND TMV_SymBandMatrix.cpp TMV_MultsBV.cpp TMV_AddsBsB.cpp TMV_MultsBM.cpp TMV_MultsBB.cpp TMV_SymBandCHD.cpp TMV_SymBandCHDecompose.cpp TMV_SymBandCHDiv.cpp TMV_SymBandCHInverse.cpp TMV_SymBandSVD.cpp TMV_SymBandSVDecompose.cpp)

SET(LIB ${BASIC} ${DIAG} ${TRI})
SET(SBLIB ${BAND} ${SYM} ${SYMBAND})

INCLUDE_DIRECTORIES(../include .)

ADD_LIBRARY(tmv ${LIB})
ADD_LIBRARY(tmv_symband ${SBLIB})
