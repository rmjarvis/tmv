
// Need to define the following with #define statements.
// (The given definition is for a regular Matrix*Vector.  Modify as 
// appropriate for the various other matrices.)
//
// #define PRODMV ProdMV

// x*(x*m*v)
template <class T, class T1, bool tr, class T2> 
inline PRODMV<T,T1,tr,T2> operator*(T x, const PRODMV<T,T1,tr,T2>& pmv)
{ return PRODMV<T,T1,tr,T2>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    CT x, const PRODMV<T,T,tr,T>& pmv)
{ return PRODMV<CT,T,tr,T>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    CCT x, const PRODMV<T,T,tr,T>& pmv)
{ return PRODMV<CT,T,tr,T>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    VCT x, const PRODMV<T,T,tr,T>& pmv)
{ return PRODMV<CT,T,tr,T>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T2,tr,T1> operator*(T x, const PRODMV<CT,T1,tr,T2>& pmv)
{ return PRODMV<CT,T1,tr,T2>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator*(CCT x, const PRODMV<CT,T1,tr,T2>& pmv)
{ return PRODMV<CT,T1,tr,T2>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator*(VCT x, const PRODMV<CT,T1,tr,T2>& pmv)
{ return PRODMV<CT,T1,tr,T2>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

// (x*m*v)*x
template <class T, class T1, bool tr, class T2> 
inline PRODMV<T,T1,tr,T2> operator*(const PRODMV<T,T1,tr,T2>& pmv, T x)
{ return PRODMV<T,T1,tr,T2>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    const PRODMV<T,T,tr,T>& pmv, CT x)
{ return PRODMV<CT,T,tr,T>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    const PRODMV<T,T,tr,T>& pmv, CCT x)
{ return PRODMV<CT,T,tr,T>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator*(
    const PRODMV<T,T,tr,T>& pmv, VCT x)
{ return PRODMV<CT,T,tr,T>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator*(const PRODMV<CT,T1,tr,T2>& pmv, T x)
{ return PRODMV<CT,T1,tr,T2>(x*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator*(const PRODMV<CT,T1,tr,T2>& pmv, CCT x)
{ return PRODMV<CT,T1,tr,T2>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator*(const PRODMV<CT,T1,tr,T2>& pmv, VCT x)
{ return PRODMV<CT,T1,tr,T2>(CT(x)*pmv.GetX(),pmv.GetM(),pmv.GetV()); }

// (x*m*v)/x
template <class T, class T1, bool tr, class T2> 
inline PRODMV<T,T1,tr,T2> operator/(const PRODMV<T,T1,tr,T2>& pmv, T x)
{ return PRODMV<T,T1,tr,T2>(pmv.GetX()/x,pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator/(
    const PRODMV<T,T,tr,T>& pmv, CT x)
{ return PRODMV<CT,T,tr,T>(pmv.GetX()/x,pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator/(
    const PRODMV<T,T,tr,T>& pmv, CCT x)
{ return PRODMV<CT,T,tr,T>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, bool tr> inline PRODMV<CT,T,tr,T> operator/(
    const PRODMV<T,T,tr,T>& pmv, VCT x)
{ return PRODMV<CT,T,tr,T>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator/(const PRODMV<CT,T1,tr,T2>& pmv, T x)
{ return PRODMV<CT,T1,tr,T2>(pmv.GetX()/x,pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator/(const PRODMV<CT,T1,tr,T2>& pmv, CCT x)
{ return PRODMV<CT,T1,tr,T2>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

template <class T, class T1, bool tr, class T2> 
inline PRODMV<CT,T1,tr,T2> operator/(const PRODMV<CT,T1,tr,T2>& pmv, VCT x)
{ return PRODMV<CT,T1,tr,T2>(pmv.GetX()/CT(x),pmv.GetM(),pmv.GetV()); }

