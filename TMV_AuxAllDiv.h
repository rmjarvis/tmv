#define RT RealType(T)
#define CT ComplexType(T)

#define DefDivEq(T) \
      inline void LDivEq(const MatrixView<T>& m) const \
      { DoLDivEq(m); } \
      inline void RDivEq(const MatrixView<T>& m) const \
      { DoRDivEq(m); }
      DefDivEq(RT);
      DefDivEq(CT);
#undef DefDivEq

#define DefDiv(T1,T2) \
      inline void LDiv(const GenMatrix<T1>& m1, \
	  const MatrixView<T2>& m0) const \
      { DoLDiv(m1,m0); } \
      inline void RDiv(const GenMatrix<T1>& m1, \
	  const MatrixView<T2>& m0) const \
      { DoRDiv(m1,m0); }
      DefDiv(RT,RT);
      DefDiv(RT,CT);
      DefDiv(CT,CT);
#undef DefDiv

#undef RT
#undef CT
