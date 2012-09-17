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
    return QPointF(400 + 200 * x, 400 - 200 * y);
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
    for (int i = 1; i < trajectory.n; i += 2)
    {
        const Vector& previous = trajectory.positions.at(i-1);
        Vector& current = trajectory.positions[i];
        const Vector& next = trajectory.positions.at(i+1);
        const Vector F = force(current, dt * i);
        
        const Vector eps = previous + next - current * 2 - F * dt * dt;
        change += (eps.x * eps.x) + (eps.y * eps.y);
        current = current + eps * 0.25 * chebishev;
    }
    for (int i = 2; i < trajectory.n; i += 2)
    {
        const Vector& previous = trajectory.positions.at(i-1);
        Vector& current = trajectory.positions[i];
        const Vector& next = trajectory.positions.at(i+1);
        const Vector F = force(current, dt * i);
        
        const Vector eps = previous + next - current * 2 - F * dt * dt;
        change += (eps.x * eps.x) + (eps.y * eps.y);
        current = current + eps * 0.25 * chebishev;
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

QImage Planets::plot(const Trajectory& trajectory)
{
    QImage image(800, 800, QImage::Format_ARGB32);
    image.fill(Qt::black);
    QPainter painter(&image);
    const double T = dt * trajectory.n;
    
    painter.setBrush(Qt::yellow);
    painter.drawEllipse(star, 15, 15);

    painter.setBrush(Qt::blue);
    painter.drawEllipse(planet_one(0), 5, 5);
    
    painter.setBrush(Qt::red);
    painter.drawEllipse(planet_two(trajectory.n*dt), 5, 5);

    painter.setBrush(Qt::NoBrush);
    painter.setPen(Qt::white);
    painter.drawEllipse(star, 200, 200);
    painter.drawEllipse(star, 400, 400);
    
    for (int i = 0; i < trajectory.n; ++i)
    {
        int f = 255 * (double)i / trajectory.n;
        f = qBound(0, f, 255);
        painter.setPen(QColor(qRgb(f, 255, 255-f)));
        painter.drawLine(trajectory.positions[i], trajectory.positions[i+1]);
    }
    painter.end();
    return image;
}

QPair< double, double > Planets::burst(const Trajectory& trajectory) const
{
    double invt = 1.0/dt;
    Vector dv1 = (trajectory.positions[1] - trajectory.positions[0] - planet_one(dt) + planet_one(0)) * invt;
    int n = trajectory.n;
    Vector dv2 = (trajectory.positions[n] - trajectory.positions[n-1] - planet_two(n*dt) + planet_two((n-1)*dt)) * invt;
    return qMakePair(dv1 * dv1, dv2 * dv2);
}
