// x*(x*m/m)
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator*(
    T x, const QUOTMM<T,T1,T2>& qmm)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    CT x, const QUOTMM<CT,T,T>& qmm)
{ return QUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    CCT x, const QUOTMM<CT,T,T>& qmm)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    VCT x, const QUOTMM<CT,T,T>& qmm)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    T x, const QUOTMM<CT,T1,T2>& qmm)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    CCT x, const QUOTMM<CT,T1,T2>& qmm)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    VCT x, const QUOTMM<CT,T1,T2>& qmm)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m/m)*x
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator*(
    const QUOTMM<T,T1,T2>& qmm, T x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    const QUOTMM<CT,T,T>& qmm, CT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    const QUOTMM<CT,T,T>& qmm, CCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator*(
    const QUOTMM<CT,T,T>& qmm, VCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    const QUOTMM<CT,T1,T2>& qmm, T x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    const QUOTMM<CT,T1,T2>& qmm, CCT x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator*(
    const QUOTMM<CT,T1,T2>& qmm, VCT x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m/m)/x
template <class T, class T1, class T2> inline QUOTMM<T,T1,T2> operator/(
    const QUOTMM<T,T1,T2>& qmm, T x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator/(
    const QUOTMM<CT,T,T>& qmm, CT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator/(
    const QUOTMM<CT,T,T>& qmm, CCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline QUOTMM<CT,T,T> operator/(
    const QUOTMM<CT,T,T>& qmm, VCT x)
{ return QUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator/(
    const QUOTMM<CT,T1,T2>& qmm, T x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator/(
    const QUOTMM<CT,T1,T2>& qmm, CCT x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline QUOTMM<CT,T1,T2> operator/(
    const QUOTMM<CT,T1,T2>& qmm, VCT x)
{ return QUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

// x*(x*m%m)
template <class T, class T1, class T2> inline RQUOTMM<T,T1,T2> operator*(
    T x, const RQUOTMM<T,T1,T2>& qmm)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    CT x, const RQUOTMM<CT,T,T>& qmm)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    CCT x, const RQUOTMM<CT,T,T>& qmm)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    VCT x, const RQUOTMM<CT,T,T>& qmm)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    T x, const RQUOTMM<CT,T1,T2>& qmm)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    CCT x, const RQUOTMM<CT,T1,T2>& qmm)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    VCT x, const RQUOTMM<CT,T1,T2>& qmm)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m%m)*x
template <class T, class T1, class T2> inline RQUOTMM<T,T1,T2> operator*(
    const RQUOTMM<T,T1,T2>& qmm, T x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    const RQUOTMM<CT,T,T>& qmm, CT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    const RQUOTMM<CT,T,T>& qmm, CCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator*(
    const RQUOTMM<CT,T,T>& qmm, VCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    const RQUOTMM<CT,T1,T2>& qmm, T x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    const RQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator*(
    const RQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()*CT(x),qmm.GetM1(),qmm.GetM2()); }

// (x*m%m)/x
template <class T, class T1, class T2> inline RQUOTMM<T,T1,T2> operator/(
    const RQUOTMM<T,T1,T2>& qmm, T x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator/(
    const RQUOTMM<CT,T,T>& qmm, CT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator/(
    const RQUOTMM<CT,T,T>& qmm, CCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T> inline RQUOTMM<CT,T,T> operator/(
    const RQUOTMM<CT,T,T>& qmm, VCT x)
{ return RQUOTMM<CT,T,T>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator/(
    const RQUOTMM<CT,T1,T2>& qmm, T x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()/x,qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator/(
    const RQUOTMM<CT,T1,T2>& qmm, CCT x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }
    
template <class T, class T1, class T2> inline RQUOTMM<CT,T1,T2> operator/(
    const RQUOTMM<CT,T1,T2>& qmm, VCT x)
{ return RQUOTMM<T,T1,T2>(qmm.GetX()/CT(x),qmm.GetM1(),qmm.GetM2()); }

