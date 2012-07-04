#include "planets.h"
#include "math.h"

Vector Vector::operator+(const Vector& other) const
{
    Vector v;
    v.x = x + other.x;
    v.y = y + other.y;
    return v;
}

Vector Vector::operator-(const Vector& other) const
{
    Vector v;
    v.x = x - other.x;
    v.y = y - other.y;
}

Vector Vector::operator*(double factor) const
{
    Vector v;
    v.x = x * factor;
    v.y = y * factor;
    return v;
}

Vector gravity(const Vector& object, const Vector& source)
{
    Vector d = source - object;
    double d2 = (d.x) * (d.x) + (d.y) * (d.y);
    double d13 = pow(d2, -1.5);
    return d * d13;
}

Vector from_polar(double r, double phi)
{
    Vector p;
    p.x = r * cos(phi);
    p.y = r * sin(phi);
    return p;
}

Planets::Planets(double r, double mu)
 : r(r)
 , omega(sqrt(1.0/r))
 , mu(mu)
{
}

Vector Planets::planet_one(double t)
{
    return from_polar(1, 2 * omega * t);
}

Vector Planets::planet_two(double t)
{
}

void Planets::relax_step(Trajectory& trajectory)
{
    double change;
    for (int i = 1; i < trajectory.n-1; i += 2)
    {
        const Vector& previous = trajectory.positions.at(i-1);
        Vector& current = trajectory.positions[i];
        const Vector& next = trajectory.positions.at(i);
        const Vector F = force(current, dt * i);
        
        const Vector eps = previous + next - current * 2 - F * 4 * dt * dt;
        current = current + eps * 0.25;
    }
}

Vector Planets::force(const Vector& pos, double t)
{
    return gravity(pos, star) + (gravity(pos, planet_one(t)) + gravity(pos, planet_two(t))) * mu;
}
