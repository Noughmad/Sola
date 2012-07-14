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


Vector gravity(const Vector& object, const Vector& source)
{
    Vector d = source - object;
    const double sd2 = sqrt((d.x) * (d.x) + (d.y) * (d.y));
 //   qDebug() << "gravity(): " << d << d2 << d13;
    return d * (1.0 / sd2 / sd2 / sd2);
}

Vector from_polar(double r, double phi)
{
    Vector p;
    p.x = r * cos(phi);
    p.y = r * sin(phi);
    return p;
}

Vector::operator QPointF() const
{
    return QPointF(400 + 200 * x, 400 - 200 * y);
}

Planets::Planets(double r, double mu)
 : r(r)
 , omega(sqrt(1.0/r))
 , mu(mu)
{
}

Vector Planets::planet_one(double t)
{
    return from_polar(r, omega * t);
}

Vector Planets::planet_two(double t)
{
    return from_polar(2*r, omega * t / sqrt(8) + phase);
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
    for (int i = trajectory.n-2; i > 1; i -= 2)
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
    Trajectory t;
    t.n = steps;
    Vector start = planet_one(0);
    Vector end = planet_two(steps * dt);
    qDebug() << start << end;
    for (int i = 0; i < steps+1; ++i)
    {
        const double u = (double)i/steps;
        t.positions << (end * u + start * (1-u));
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
    painter.drawEllipse(star, 10, 10);
    
    for (int i = 0; i < trajectory.n; ++i)
    {
        int f = 255 * (double)i / trajectory.n;
        f = qBound(0, f, 255);
        painter.setPen(QColor(qRgb(f, 255, 255-f)));
        painter.drawLine(trajectory.positions[i], trajectory.positions[i+1]);
        
        painter.setBrush(QColor(qRgb(0, 0, 255 - f/2)));
        painter.drawEllipse(planet_one(i*dt), 5, 5);
        
        painter.setBrush(QColor(qRgb(255 - f/2, 0, 0)));
        painter.drawEllipse(planet_two(i*dt), 5, 5);
    }
    painter.end();
    return image;
}
