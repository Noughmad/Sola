#ifndef PLANETS_H
#define PLANETS_H

#include <QVector>
#include <QDebug>
#include <qmath.h>

static const double invsqrt8 = 1.0/sqrt(8);

class QPaintDevice;
struct Vector
{
    double x;
    double y;
    
    Vector operator-(const Vector& other) const;
    Vector operator+(const Vector& other) const;
    Vector operator*(double factor) const;
    operator QPointF() const;

    void to_polar(double& r, double& phi) const;
};

double operator*(const Vector& one, const Vector& other);
QDebug operator<<(QDebug stream, const Vector& vector);


inline Vector from_polar(double r, double phi)
{
    Vector p;
    
    p.x = r * cos(phi);
    p.y = r * sin(phi);
    return p;
}

inline void Vector::to_polar(double& r, double& phi) const
{
    phi = atan2(y, x);
    if (phi < 0)
    {
        phi += 2*M_PI;
    }
    r = sqrt(x*x+y*y);
}

struct Trajectory
{
    Trajectory();
    Trajectory(int n) : n(n), positions(n+1) {}

    void operator=(const Trajectory& other);
    
    QVector<Vector> positions;
    int n;
};

double distance_sq(const Vector& one, const Vector& other);
Vector from_polar(double r, double phi);

class Planets
{
public:
    Planets(double mu, double phase);
    
    inline Vector planet_one(double t) const;
    inline Vector planet_two(double t) const;
    
    Vector force(const Vector& pos, double t);
    
    double relax_step(Trajectory& trajectory);
    void move_planets();
    
    Trajectory direct_route(int steps);
    Trajectory ellipse(int steps, int circles = 0);
    Trajectory spline(int steps);
    Trajectory orbit(int steps, int circles = 0);
    Trajectory toolbox(int steps, int circles = 0);

    QPair<double, double> burst(const Trajectory& trajectory) const;
    void plot(QPaintDevice* image, const Trajectory& trajectory);
    double adjust_speed(double speed) const;
    
    double dt;
    double mu;
    double phase;
    double M;
    double m;
    double chebishev;
    
    Vector star;
};

Vector Planets::planet_one(double t) const
{
    return from_polar(1, t);
}

Vector Planets::planet_two(double t) const
{
    return from_polar(2, t * invsqrt8 + phase);
}

#endif // PLANETS_H
