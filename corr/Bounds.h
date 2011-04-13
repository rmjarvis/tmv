//---------------------------------------------------------------------------
#ifndef BoundsH
#define BoundsH

#include <complex>
#include "MemDebug.h"

//---------------------------------------------------------------------------

class Position2D 
{

public:
    Position2D() : x(0.),y(0.) {}
    Position2D(const Position2D& rhs) : x(rhs.x), y(rhs.y) {}
    ~Position2D() {}
    Position2D(double xin,double yin) : x(xin), y(yin) {}
    Position2D& operator=(const Position2D& rhs) 
    { x = rhs.x; y = rhs.y; return *this; }

    double GetX() const { return x; }
    double GetY() const { return y; }
    operator std::complex<double>() const { return std::complex<double>(x,y); }

    Position2D& operator+=(const Position2D& p2)
    { x += p2.GetX(); y += p2.GetY(); return *this; }
    Position2D& operator-=(const Position2D& p2)
    { x -= p2.GetX(); y -= p2.GetY(); return *this; }
    Position2D& operator*=(double a)
    { x *= a; y *= a; return *this; }
    Position2D& operator/=(double a)
    { x /= a; y /= a; return *this; }

    Position2D operator+(const Position2D& p2) const
    { return Position2D(x+p2.GetX(),y+p2.GetY()); }
    Position2D operator-(const Position2D& p2) const
    { return Position2D(x-p2.GetX(),y-p2.GetY()); }
    Position2D operator*(double a) const
    { return Position2D(x*a,y*a); }
    Position2D operator/(double a) const
    { return Position2D(x/a,y/a); }

    void Read(std::istream& fin) { fin >> x >> y; }
    void Write(std::ostream& fout) const
    { fout << GetX() << " " << GetY() << " "; }

private:

    double x,y;

}; // Position2D

inline double DistSq(const Position2D& p1, const Position2D& p2)
{ 
    Position2D r = p1-p2;
    return r.GetX()*r.GetX() + r.GetY()*r.GetY(); 
}
inline double Dist(const Position2D& p1, const Position2D& p2)
{ return sqrt(DistSq(p1,p2)); }

inline std::ostream& operator<<(std::ostream& os, const Position2D& pos)
{ pos.Write(os); return os; }

inline std::istream& operator>>(std::istream& os, Position2D& pos)
{ pos.Read(os); return os; }

class Bounds2D 
{
    // Basically just a rectangle.  This is used to keep track of the bounds of
    // catalogs and fields.  You can set values, but generally you just keep
    // including positions of each galaxy or the bounds of each catalog
    // respectively using the += operators

public:
    Bounds2D(double x1, double x2, double y1, double y2):
        defined(1),xmin(x1),xmax(x2),ymin(y1),ymax(y2) {}
    Bounds2D(const Position2D& pos):
        defined(1),xmin(pos.GetX()),xmax(pos.GetX()),
        ymin(pos.GetY()),ymax(pos.GetY()) {}
    Bounds2D(): defined(0),xmin(0.),xmax(0.),ymin(0.),ymax(0.) {}
    ~Bounds2D() {}
    double GetXMin() const {return xmin;}
    double GetXMax() const {return xmax;}
    double GetYMin() const {return ymin;}
    double GetYMax() const {return ymax;}
    bool IsDefined() const {return defined;}
    void operator+=(const Position2D& pos)
        // Expand the bounds to include the given position.
    {
        if(defined) {
            if(pos.GetX() < xmin) xmin = pos.GetX();
            else if (pos.GetX() > xmax) xmax = pos.GetX();
            if(pos.GetY() < ymin) ymin = pos.GetY();
            else if (pos.GetY() > ymax) ymax = pos.GetY();
        }
        else {
            xmin = xmax = pos.GetX();
            ymin = ymax = pos.GetY();
            defined = 1;
        }
    }

    void Write(std::ostream& fout) const
    { fout << xmin << ' ' << xmax << ' ' << ymin << ' ' << ymax << ' '; }
    void Read(std::istream& fin)
    { fin >> xmin >> xmax >> ymin >> ymax; defined = true; }

private:
    bool defined;
    double xmin,xmax,ymin,ymax;

};

inline std::ostream& operator<<(std::ostream& fout, const Bounds2D& b)
{ b.Write(fout); return fout;}

inline std::istream& operator>>(std::istream& fin, Bounds2D& b)
{ b.Read(fin); return fin;}

class Position3D 
{

public:
    Position3D() : x(0.),y(0.),z(0.) {}
    Position3D(const Position3D& rhs) : x(rhs.x),y(rhs.y),z(rhs.z) {}
    ~Position3D() {}
    Position3D(double xin,double yin,double zin) : x(xin),y(yin),z(zin) {}
    Position3D& operator=(const Position3D& rhs) 
    { x = rhs.GetX(); y = rhs.GetY(); z = rhs.GetZ(); return *this; }

    double GetX() const { return(x); }
    double GetY() const { return(y); }
    double GetZ() const { return(z); }

    Position3D& operator+=(const Position3D& p2)
    { x += p2.GetX(); y += p2.GetY(); z += p2.GetZ(); return *this; }
    Position3D& operator-=(const Position3D& p2)
    { x -= p2.GetX(); y -= p2.GetY(); z -= p2.GetZ(); return *this; }
    Position3D& operator*=(double a)
    { x *= a; y *= a; z *= a; return *this; }
    Position3D& operator/=(double a)
    { x /= a; y /= a; z /= a; return *this; }

    Position3D operator+(const Position3D& p2) const
    { return Position3D(x+p2.GetX(),y+p2.GetY(),z+p2.GetZ()); }
    Position3D operator-(const Position3D& p2) const
    { return Position3D(x-p2.GetX(),y-p2.GetY(),z-p2.GetZ()); }
    Position3D operator*(double a) const
    { return Position3D(x*a,y*a,z*a); }
    Position3D operator/(double a) const
    { return Position3D(x/a,y/a,z/a); }

    void Read(std::istream& fin) 
    { fin >> x >> y >> z; }
    void Write(std::ostream& fout) const
    { fout << x << " " << y << " " << z << " "; }

private:
    double x,y,z;

}; // Position3D

inline std::ostream& operator<<(std::ostream& os, const Position3D& pos)
{ pos.Write(os); return os; }

inline std::istream& operator>>(std::istream& os, Position3D& pos)
{ pos.Read(os); return os; }

class Bounds3D 
{

public:
    Bounds3D(double x1, double x2, double y1, double y2, double z1, double z2):
        defined(1),xmin(x1),xmax(x2),ymin(y1),ymax(y2),zmin(z1),zmax(z2) {}
    Bounds3D(const Position3D& pos):
        defined(1),xmin(pos.GetX()),xmax(pos.GetX()),
        ymin(pos.GetY()),ymax(pos.GetY()),
        zmin(pos.GetZ()),zmax(pos.GetZ()) {}
    Bounds3D(): defined(0),xmin(0.),xmax(0.),ymin(0.),ymax(0.),
    zmin(0.),zmax(0.) {}
    ~Bounds3D() {}
    double GetXMin() const {return xmin;}
    double GetXMax() const {return xmax;}
    double GetYMin() const {return ymin;}
    double GetYMax() const {return ymax;}
    double GetZMin() const {return zmin;}
    double GetZMax() const {return zmax;}
    bool IsDefined() const {return defined;}
    void operator+=(const Position3D& pos)
        // Expand the bounds to include the given position.
    {
        if(defined) {
            if(pos.GetX() < xmin) xmin = pos.GetX();
            else if (pos.GetX() > xmax) xmax = pos.GetX();
            if(pos.GetY() < ymin) ymin = pos.GetY();
            else if (pos.GetY() > ymax) ymax = pos.GetY();
            if(pos.GetZ() < zmin) zmin = pos.GetZ();
            else if (pos.GetZ() > zmax) zmax = pos.GetZ();
        }
        else {
            xmin = xmax = pos.GetX();
            ymin = ymax = pos.GetY();
            zmin = zmax = pos.GetZ();
            defined = 1;
        }
    }

    void Write(std::ostream& fout) const
    { fout << xmin << ' ' << xmax << ' ' << ymin << ' ' << ymax << 
        ' ' << zmin << ' ' << zmax << ' ';}
        void Read(std::istream& fin)
        { fin >> xmin >> xmax >> ymin >> ymax >> zmin >> zmax; defined = true; }

private:
        bool defined;
        double xmin,xmax,ymin,ymax,zmin,zmax;

};

inline double DistSq(const Position3D& p1, const Position3D& p2)
{ 
    Position3D r = p1-p2;
    return r.GetX()*r.GetX() + r.GetY()*r.GetY() + r.GetZ()*r.GetZ(); 
}
inline double Dist(const Position3D& p1, const Position3D& p2)
{ return sqrt(DistSq(p1,p2)); }

inline std::ostream& operator<<(std::ostream& fout, const Bounds3D& b)
{ b.Write(fout); return fout;}

inline std::istream& operator>>(std::istream& fin, Bounds3D& b)
{ b.Read(fin); return fin;}

#endif
