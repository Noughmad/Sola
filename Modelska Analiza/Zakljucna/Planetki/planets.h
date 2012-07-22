#ifndef PLANETS_H
#define PLANETS_H

#include <QVector>
#include <QDebug>

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

QDebug operator<<(QDebug stream, const Vector& vector);

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
    
    double relax_step(Trajectory& trajectory);
    void move_planets();
    
    Trajectory direct_route(int steps);
    Trajectory ellipse(int steps, int circles);
    
    QImage plot(const Trajectory& trajectory);
    
    double dt;
    double mu;
    double phase;
    double omega;
    double M;
    double m;
    double r;
    double chebishev;
    
    Vector star;
};

#endif // PLANETS_H
