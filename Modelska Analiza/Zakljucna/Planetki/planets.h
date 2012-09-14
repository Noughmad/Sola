#ifndef PLANETS_H
#define PLANETS_H

#include <QVector>
#include <QDebug>
#include <qmath.h>

static const double invsqrt8 = 1.0/sqrt(8);

class QImage;
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

    QPair<double, double> burst(const Trajectory& trajectory) const;
    QImage plot(const Trajectory& trajectory);
    
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
