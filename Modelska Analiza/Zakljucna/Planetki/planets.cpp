#include "planets.h"
#include "math.h"

#include <QImage>
#include <QPainter>

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
    return v;
}

Vector Vector::operator*(double factor) const
{
    Vector v;
    v.x = x * factor;
    v.y = y * factor;
    return v;
}

QDebug operator<<(QDebug stream, const Vector& vector)
{
    stream << '(' << vector.x << ", " << vector.y << ')';
    return stream;
}

double operator*(const Vector& one, const Vector& other)
{
    return sqrt(one.x * other.x + one.y * other.y);
}

Vector gravity(const Vector& object, const Vector& source)
{
    Q_ASSERT(object.x != source.x || object.y != source.y);
    Vector d = source - object;
    const double sd2 = d * d;
    return d * (1.0 / (sd2 * sd2 * sd2));
}

Vector::operator QPointF() const
{
    return QPointF(400 + 100 * x, 400 - 100 * y);
}

Trajectory::Trajectory() : n(0), positions(0)
{

}

void Trajectory::operator=(const Trajectory& other)
{
    n = other.n;
    positions = other.positions;
}


Planets::Planets(double mu, double phase)
 : mu(mu), phase(phase)
{
}

double Planets::relax_step(Trajectory& trajectory)
{
    double change = 0;
    const double factor = 0.5 * chebishev;
    for (int i = 1; i < trajectory.n; i += 2)
    {
        const Vector& previous = trajectory.positions.at(i-1);
        Vector& current = trajectory.positions[i];
        const Vector& next = trajectory.positions.at(i+1);
        const Vector F = force(current, dt * i);
        
        const Vector eps = previous + next - current * 2 - F * dt * dt;
        change += (eps.x * eps.x) + (eps.y * eps.y);
        current = current + eps * factor;
    }
    for (int i = 2; i < trajectory.n; i += 2)
    {
        const Vector& previous = trajectory.positions.at(i-1);
        Vector& current = trajectory.positions[i];
        const Vector& next = trajectory.positions.at(i+1);
        const Vector F = force(current, dt * i);
        
        const Vector eps = previous + next - current * 2 - F * dt * dt;
        change += (eps.x * eps.x) + (eps.y * eps.y);
        current = current + eps * factor;
    }
    return change;
}

Trajectory Planets::direct_route(int steps)
{
    Trajectory t(steps);
    Vector start = planet_one(0);
    Vector end = planet_two(steps * dt);
    for (int i = 0; i < steps+1; ++i)
    {
        const double u = (double)i/steps;
        t.positions[i] = end * u + start * (1-u);
    }
    return t;
}

Trajectory Planets::ellipse(int steps, int circles)
{
    Trajectory t(steps);

    const double angle = (circles * 2 * M_PI + phase) / steps + dt * invsqrt8;
    const double rad = 1.0 / steps;

    for (int i = 0; i < steps+1; ++i)
    {
        t.positions[i] = from_polar(1.0 + i*rad, i*angle);
    }

    return t;
}

Trajectory Planets::spline(int steps)
{
    Trajectory t(steps);
    int m = steps / 4;

    const double R1 = sqrt(sqrt(2));

    for (int i = 0; i < m; ++i)
    {
        t.positions[i] = planet_one(i * dt * R1);
        t.positions[steps-i] = planet_two((steps-i*R1) * dt);
    }

    double phi1 = m * dt * R1;
    double phi2 = (steps-m*R1) * dt * invsqrt8;

    for (int i = m; i <= (steps-m); ++i)
    {
        const double r = (double)(i - m)/(steps-2*m);
        t.positions[i] = from_polar(1 + r, phi1 * (1-r) + phi2 * r);
    }

    return t;
}

Trajectory Planets::orbit(int steps, int circles)
{
    Trajectory t(steps);

    const double angle = (circles * 2 * M_PI + phase) / steps + dt * invsqrt8;

    t.positions[0] = planet_one(0);
    t.positions[steps] = planet_two(steps * dt);

    double e = 1.0/3.0;
    double r0 = 1+e;
    for (int i = 1; i < steps; ++i)
    {
        t.positions[i] = from_polar(r0 / (1 + e * cos(i*M_PI/steps)), i*angle);
    }
    
    return t;
}

Vector Planets::force(const Vector& pos, double t)
{
    return gravity(pos, star) + (gravity(pos, planet_one(t)) + gravity(pos, planet_two(t))) * mu;
}

void Planets::plot(QPaintDevice* image, const Trajectory& trajectory)
{

}

double Planets::adjust_speed(double dr) const
{
    double speed = dr / dt;
    dr *= 0.6385; // The factor here is magic, I got it by trying when the result is not dependent on number of steps
    double v = sqrt(speed * speed - 2 * mu * (1.0/dr - 1000));
    return v;
}

QPair< double, double > Planets::burst(const Trajectory& trajectory) const
{
    double invt = 1.0/dt;
    const QVector<Vector>& r = trajectory.positions;
    Vector dv1 = r[1] - planet_one(dt);

    int n = trajectory.n;
    Vector dv2 = r[n-1] - planet_two((n-1)*dt);
    return qMakePair(dv1 * dv1, dv2 * dv2);
}
