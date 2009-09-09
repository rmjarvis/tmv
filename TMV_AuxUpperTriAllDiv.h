#define RT RealType(T)
#define CT ComplexType(T)

#define DefDivEq(T) \
      inline void LDivEq(const UpperTriMatrixView<T>& m) const \
      { DoLDivEq(m); } \
      inline void RDivEq(const UpperTriMatrixView<T>& m) const \
      { DoRDivEq(m); } \

      DefDivEq(RT);
      DefDivEq(CT);
#undef DefDivEq

#define DefDiv(T1,T2) \
      inline void LDiv(const GenUpperTriMatrix<T1>& m1, \
	  const UpperTriMatrixView<T2>& m0) const \
      { DoLDiv(m1,m0); } \
      inline void RDiv(const GenUpperTriMatrix<T1>& m1, \
	  const UpperTriMatrixView<T2>& m0) const \
      { DoRDiv(m1,m0); } \

      DefDiv(RT,T);
      DefDiv(CT,CT);
#undef DefDiv

#undef RT
#undef CT
