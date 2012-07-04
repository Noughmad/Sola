#ifndef PLANETS_H
#define PLANETS_H

#include <QVector>

struct Vector
{
    double x;
    double y;
    
    Vector operator-(const Vector& other) const;
    Vector operator+(const Vector& other) const;
    Vector operator*(double factor) const;
};

struct Trajectory
{
    QList<Vector> positions;
    int n;
};

double distance_sq(const Vector& one, const Vector& other);
Vector from_polar(double r, double phi);

class Planets
{
public:
    Planets(double r, double mu);
    
    Vector planet_one(double t);
    Vector planet_two(double t);
    
    Vector force(const Vector& pos, double t);
    
    void relax_step(Trajectory& trajectory);
    void move_planets();
    
    double dt;
    double mu;
    double phase;
    double omega;
    double M;
    double m;
    double r;
    
    Vector star;
};

#endif // PLANETS_H
