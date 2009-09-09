
// Need to define the following with #define statements.
// (The given definition is for a regular Matrix*Matrix.  Modify as 
// appropriate for the various other matrices.)
//
// #define PRODMM ProdMM

// x*(x*m*m)
template <class T, class T1, class T2> inline PRODMM<T,T1,T2> operator*(
    T x, const PRODMM<T,T1,T2>& pmm)
{ return PRODMM<T,T1,T2>(x*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    CT x, const PRODMM<T,T,T>& pmm)
{ return PRODMM<CT,T,T>(x*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    CCT x, const PRODMM<T,T,T>& pmm)
{ return PRODMM<CT,T,T>(CT(x)*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    VCT x, const PRODMM<T,T,T>& pmm)
{ return PRODMM<CT,T,T>(CT(x)*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    T x, const PRODMM<CT,T1,T2>& pmm)
{ return PRODMM<CT,T1,T2>(x*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    CCT x, const PRODMM<CT,T1,T2>& pmm)
{ return PRODMM<CT,T1,T2>(CT(x)*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    VCT x, const PRODMM<CT,T1,T2>& pmm)
{ return PRODMM<CT,T1,T2>(CT(x)*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

// (x*m*m)*x
template <class T, class T1, class T2> inline PRODMM<T,T1,T2> operator*(
    const PRODMM<T,T1,T2>& pmm, T x)
{ return PRODMM<T,T1,T2>(x*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    const PRODMM<T,T,T>& pmm, CT x)
{ return PRODMM<CT,T,T>(x*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    const PRODMM<T,T,T>& pmm, CCT x)
{ return PRODMM<CT,T,T>(CT(x)*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T> inline PRODMM<CT,T,T> operator*(
    const PRODMM<T,T,T>& pmm, VCT x)
{ return PRODMM<CT,T,T>(CT(x)*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    const PRODMM<CT,T1,T2>& pmm, T x)
{ return PRODMM<CT,T1,T2>(x*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    const PRODMM<CT,T1,T2>& pmm, CCT x)
{ return PRODMM<CT,T1,T2>(CT(x)*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator*(
    const PRODMM<CT,T1,T2>& pmm, VCT x)
{ return PRODMM<CT,T1,T2>(CT(x)*pmm.GetX(),pmm.GetM1(),pmm.GetM2()); }

// (x*m*m)/x
template <class T, class T1, class T2> inline PRODMM<T,T1,T2> operator/(
    const PRODMM<T,T1,T2>& pmm, T x)
{ return PRODMM<T,T1,T2>(pmm.GetX()/x,pmm.GetM1(),pmm.GetM2()); }

template <class T> inline PRODMM<CT,T,T> operator/(
    const PRODMM<T,T,T>& pmm, CT x)
{ return PRODMM<CT,T,T>(pmm.GetX()/x,pmm.GetM1(),pmm.GetM2()); }

template <class T> inline PRODMM<CT,T,T> operator/(
    const PRODMM<T,T,T>& pmm, CCT x)
{ return PRODMM<CT,T,T>(pmm.GetX()/CT(x),pmm.GetM1(),pmm.GetM2()); }

template <class T> inline PRODMM<CT,T,T> operator/(
    const PRODMM<T,T,T>& pmm, VCT x)
{ return PRODMM<CT,T,T>(pmm.GetX()/CT(x),pmm.GetM1(),pmm.GetM2()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator/(
    const PRODMM<CT,T1,T2>& pmm, T x)
{ return PRODMM<CT,T1,T2>(pmm.GetX()/x,pmm.GetM1(),pmm.GetM2()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator/(
    const PRODMM<CT,T1,T2>& pmm, CCT x)
{ return PRODMM<CT,T1,T2>(pmm.GetX()/CT(x),pmm.GetM1(),pmm.GetM2()); }

template <class T, class T1, class T2> inline PRODMM<CT,T1,T2> operator/(
    const PRODMM<CT,T1,T2>& pmm, VCT x)
{ return PRODMM<CT,T1,T2>(pmm.GetX()/CT(x),pmm.GetM1(),pmm.GetM2()); }

