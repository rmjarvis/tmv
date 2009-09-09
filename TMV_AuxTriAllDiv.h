#define RT RealType(T)
#define CT ComplexType(T)

#define DefDivEq(T,C) \
      inline void LDivEq(const TriMatrixView<T,U,C>& m) const \
      { DoLDivEq(m); } \
      inline void RDivEq(const TriMatrixView<T,U,C>& m) const \
      { DoRDivEq(m); } \

      DefDivEq(RT,NonConj);
      DefDivEq(CT,NonConj);
      DefDivEq(CT,Conj);
#undef DefDivEq

#define DefDiv(T1,C1,T2,C2) \
      inline void LDiv(const GenTriMatrix<T1,U,C1>& m1, \
	  const TriMatrixView<T2,U,C2>& m0) const \
      { DoLDiv(m1,Regularize(m0)); } \
      inline void RDiv(const GenTriMatrix<T1,U,C1>& m1, \
	  const TriMatrixView<T2,U,C2>& m0) const \
      { DoRDiv(m1,Regularize(m0)); } \

      DefDiv(RT,NonConj,T,NonConj);
      DefDiv(RT,NonConj,T,Conj);
      DefDiv(CT,NonConj,CT,NonConj);
      DefDiv(CT,NonConj,CT,Conj);
      DefDiv(CT,Conj,CT,NonConj);
      DefDiv(CT,Conj,CT,Conj);
#undef DefDiv

#undef RT
#undef CT
