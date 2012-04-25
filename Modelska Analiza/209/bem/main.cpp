#include <iostream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <gsl/gsl_odeiv2.h>

#include <QtGui/QPolygonF>
#include <QtGui/qtransform.h>
#include <QtCore/qmath.h>
#include <QtCore/QPair>
#include <QtCore/qline.h>
#include <QtCore/QFile>
#include <QtCore/QDebug>
#include <QtCore/qtextstream.h>

inline QString fileName(const QString& name)
{
    return "../../g_" + name + ".dat";
}

inline bool compare(double p1, double p2)
{
    return (qAbs(p1 - p2) <= 0.001 * qMin(qAbs(p1), qAbs(p2)));
}

inline double arctan(double y, double x)
{
    double a = qAtan2(y, x);
    if (x > 0)
    {
        return a;
    }
    else if (y <= 0)
    {
        return a + M_PI;
    }
    else if (y > 0)
    {
        return a - M_PI;
    }
    return 0;
}

typedef QPair<double,double> K;

struct Sistem
{
    gsl_matrix* A;
    gsl_matrix* B;
    gsl_vector* b;
    QPolygonF delitev;
    gsl_vector* resitev;
    int N;
    
    Sistem& resi();
    Sistem& shrani(const QString& filename);
    
    Sistem& el_polje(const QString& filename);
    Sistem& hitrostno_polje(const QString& filename);
    Sistem& tokovnice(const QString& filename);
    
    double epot(double x, double y) const;
    
    QPair<double, double> hitrost(double x, double y);
    
    void preveri_hitrosti();
    Sistem& preveri();
    double kot(int i);
    
    Sistem& naboj(const QString& filename);
    
    Sistem& matrika(const QString& filename);
    
    double kapaciteta();
    
    void test();
};

void save(const QPolygonF& poly, const QString& name)
{
    QFile file(name);
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    
    foreach (const QPointF& point, poly)
    {
        stream << point.x() << " " << point.y() << endl;
    }
    
    file.close();
}

QPolygonF trak(int n)
{
    QPolygonF poly(n+1);
    for (int i = 0; i <= n; ++i)
    {
        poly[i] = QPointF(i*1.0/n, 0);
    }
    return poly;
}

double n_y(double x)
{
    return 1.457112*sqrt(x) - 0.624424*x - 1.727016 * x * x + 1.384087 * x * x * x - 0.489769 * x * x * x * x;
}

QPolygonF naca(int n, double t, double kot)
{
    QPolygonF poly(n+1);
    
    double F = 2.75;
    
    double k = F / n;
    double m = F / 3.0;
    double l = (double)2*(3.0-F) / n;
    
    for (int i = 0; i <= n; ++i)
    {
        double x = 1 - cos(M_PI * (i-n/2) * 1.0 / n);
        if (i < n/2)
        {
            poly[i] = QPointF(2*x-1, t / 50.0 * n_y(x));
        }
        else
        {
            poly[i] = QPointF(2*x-1, -t / 50.0 * n_y(x));
        }
    }
    
    return QTransform().rotate(kot).map(poly);
}

inline QPointF center(const QPolygonF& poly, int i)
{
    return (poly[i] + poly[i+1])/2;
}

double u(double x, double y, double ksi1, double ksi2)
{
    double x1 = x-ksi1;
    double x2 = x-ksi2;
    return 0.25 * M_1_PI * (2*x1 - 2*x2 + 2*y*(arctan(x1, y) - arctan(x2, y)) + x1 * log(x1*x1 + y*y) - x2 * log(x2*x2 + y*y));
}

inline double u(const QPointF& xy, const QPair<double, double>& center)
{
    return u(xy.x(), xy.y(), center.first, center.second);
}

void rotate(double angle, double* x, double* y)
{    
    double rr = (*x) * (*x) + (*y) * (*y);
    
    double xn = (*x) * cos(angle) + (*y) * sin(angle);
    double yn = (*y) * cos(angle) - (*x) * sin(angle);
    
    Q_ASSERT(qFuzzyCompare(xn*xn + yn*yn, rr));
    
    *x = xn;
    *y = yn;
}

void rotate(double angle, K* k)
{
    rotate(angle, &(k->first), &(k->second));
}

double angle(const QPointF& p1, const QPointF& p2)
{
    return -M_PI / 180.0 * QLineF(p1, p2).angle();
}

void transform(const QPointF& p1, const QPointF& p2, double* x, double* y)
{
    rotate(angle(p1, p2), x, y);
}

double Sistem::kot(int i)
{
    return angle(delitev[i], delitev[i+1]);
}

K v(const QPointF& ap1, const QPointF& ap2, double x, double y)
{
    QPointF c = (ap1+ap2)/2;
    
    x -= c.x();
    y -= c.y();
    
    double ksi1 = ap1.x() - c.x();
    double ksi2 = ap2.x() - c.x();
    double eta1 = ap1.y() - c.y();
    double eta2 = ap2.y() - c.y();
    
    double kot = angle(ap1, ap2);
    
    rotate(kot, &x, &y);
    rotate(kot, &ksi1, &eta1);
    rotate(kot, &ksi2, &eta2);
    
    Q_ASSERT(qFuzzyIsNull(eta1));
    Q_ASSERT(qFuzzyIsNull(eta2));
    
    
    double x1 = x - ksi1;
    double x2 = x - ksi2;
    
    K v;
    v.first = 0.25 * M_1_PI * log( (x1*x1 + y*y) / (x2*x2 + y*y) );
    v.second = (y == 0) ? -0.5 : 0.5 * M_1_PI * (arctan(x1, y) - arctan(x2, y));
    
    rotate(-kot, &v);
    return v;
}

Sistem elektroda(const QPolygonF& delitev)
{
    int n = delitev.size() - 1;
    gsl_matrix* matrika = gsl_matrix_alloc(n, n);
    gsl_vector* rhs = gsl_vector_alloc(n);
    
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            gsl_matrix_set(matrika, i, j, u(center(delitev, j), qMakePair(delitev[i].x(), delitev[i+1].x()) ));
        }
        gsl_vector_set(rhs, i, 1);
    }
        
    Sistem s;
    s.A = matrika;
    s.b = rhs;
    s.resitev = 0;
    s.delitev = delitev;
    s.N = n;
    return s;
}

Sistem obtekanje(const QPolygonF& delitev)
{
    int n = delitev.size() - 1;
    Sistem s;
    s.A = gsl_matrix_alloc(n, n);
    s.B = gsl_matrix_alloc(n, n);
    s.b = gsl_vector_alloc(n);
    s.delitev = delitev;
    s.resitev = 0;
    s.N = n;
    
    K u_inf;
    u_inf.first = 1;
    u_inf.second = 0;
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (i == j)
            {
                gsl_matrix_set(s.A, i, j, -0.5);
            }
            else
            {
                QPointF p = center(delitev, i);
                K e = v(delitev[j], delitev[j+1], p.x(), p.y());
                rotate(s.kot(i), &e);
                gsl_matrix_set(s.A, i, j, e.second);
                gsl_matrix_set(s.B, i, j, -e.first);
            }
        }
        K u = u_inf;
        rotate(-s.kot(i), &u);
        gsl_vector_set(s.b, i, u.second);
    }
    
    return s;
}

Sistem& Sistem::resi()
{
    Q_ASSERT(N);
    resitev = gsl_vector_alloc(b->size);
    gsl_matrix* T = gsl_matrix_alloc(N, N);
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            gsl_matrix_set(T, i, j, gsl_matrix_get(A, i, j));
        }
    }
    gsl_linalg_HH_solve(T, b, resitev);
    
    return *this;
}

void napisi(QTextStream& stream, Sistem& sistem, int i, int j)
{
    stream << sistem.delitev[i].x() << " " << sistem.delitev[i].y() << " ";
    stream << gsl_vector_get(sistem.resitev, j);
    
    if (i < sistem.N)
    {
        QPointF p = center(sistem.delitev, i);
        K h = sistem.hitrost(p.x(), p.y());
       // rotate(-sistem.kot(i), &h);
        stream << " " << h.first << " " << h.second;
    }
    
    stream << endl;
}

Sistem& Sistem::shrani(const QString& filename)
{
    if (!filename.isEmpty())
    {
        int n = delitev.size() - 1;
        QFile file(filename);
        file.open(QIODevice::WriteOnly);
        QTextStream stream(&file);
        napisi(stream, *this, 0, 0);
        for (int i = 1; i < n; ++i)
        {
            napisi(stream, *this, i, i-1);
            napisi(stream, *this, i, i);
        }
        napisi(stream, *this, n, n-1);        
        
        file.close();
    }
    return *this;
}

QPolygonF zukovski(int n, double A, double B)
{
    double r = sqrt((A-1)*(A-1) + B*B);
    
    QPolygonF poly(n+1);
    for (int i = 0; i <= n; ++i)
    {
        double x = A + r * cos(i * 2.0 * M_PI / n);
        double y = B + r * sin(i * 2.0 * M_PI / n);
        
        double t = x*x + y*y;
        QPointF p(x + x/t, y - y/t);
        poly[i] = p/2;
    }
    return poly;
}

QPolygonF elipsoid(int n, double b)
{
    double phi = 2 * M_PI / n;
    QPolygonF poly(n+1);
    for (int i = 0; i <= n; ++i)
    {
        poly[i] = QPointF(cos(i*phi), b * sin(i*phi));
    }
    return poly;
}

double Sistem::epot(double x, double y) const
{
    double pot = 0;
    for (int i = 0; i < N; ++i)
    {
        const QPointF p = center(delitev, i);
        double rr = (p.x() - x) * (p.x() - x) + (p.y() - y) * (p.y() - y);
        pot += gsl_vector_get(resitev, i) * 0.25 * M_1_PI * log(rr);
    }
    return pot;
}

Sistem& Sistem::el_polje(const QString& filename)
{
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    const double L = 2;
    
    for (double y = -L; y <= L; y += 0.02)
    {
        for (double x = -L; x <= L; x += 0.02)
        {
            stream << x << " " << y << " " << epot(x, y) << endl;
        }
        stream << endl;
    }
    
    file.close();
    return *this;
}

K Sistem::hitrost(double x, double y)
{
    double vx = 1, vy = 0;
    Q_ASSERT(N > 0);
    for (int i = 0; i < N; ++i)
    {
        K t;
        if (QPointF(x,y) == center(delitev, i))
        {
            t = qMakePair(0.0, -0.5);
            rotate(-kot(i), &t);
        }
        else
        {
            t = v(delitev[i], delitev[i+1], x, y);
        }
        
        vx += gsl_vector_get(resitev, i) * t.first;
        vy += gsl_vector_get(resitev, i) * t.second;
    }
    
    return qMakePair(vx, vy);
}


Sistem& Sistem::hitrostno_polje(const QString& filename)
{
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    const double Lx = 2;
    const double Ly = 1;
    
    for (double y = -Ly; y <= Ly; y += Ly/15.5)
    {
        for (double x = -Lx; x <= Lx; x += Lx/10.0)
        {
            K v = hitrost(x, y);
            stream << x << " " << y << " " << v.first << " " << v.second << endl;
        }
        stream << endl;
    }
    
    file.close();
    return *this;
}

int odvod_hitrost(double t, const double y[], double dydt[], void* params)
{
    Sistem* s = (Sistem*)params;
    
    QPair< double, double > p = s->hitrost(y[0], y[1]);
    dydt[0] = p.first;
    dydt[1] = p.second;
    
    return GSL_SUCCESS;
}

Sistem& Sistem::tokovnice(const QString& filename)
{
    
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    const double Lx = 2;
    const double Ly = 1;
    
    gsl_odeiv2_system sys = {odvod_hitrost, 0, 2, this};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-2, 1e-3, 1e-3);
    
    for (double y = -Ly; y <= Ly; y += Ly / 10.0)
    {
        double v[2] = {-Lx, y};
        double t = 0;
        double dt = 1e-2;
        
        while (v[0] < Lx && t < 10)
        {
            gsl_odeiv2_driver_apply(driver, &t, t+dt, v);
            stream << v[0] << " " << v[1] << endl;
        }
        stream << endl;
    }
    
    file.close();
    return *this;
}

void Sistem::preveri_hitrosti()
{
    for (int i = 0; i < N; ++i)
    {
        QPointF p = center(delitev, i);
        QPair< double, double > a = hitrost(p.x(), p.y());
        /*
        
        double phi = -M_PI / 180.0 * QLineF(delitev[i], delitev[i+1]).angle();
        qDebug() << a.first * cos(phi) + a.second * sin(phi) << a.second * cos(phi) - a.first * sin(phi);
        */
        qDebug() << a;
    }
}

void test()
{
    QPointF one(0,0);
    QPointF two(-1,-1);
    
    QFile file("g_test.dat");
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    
    double L = 3.0;
    
    for (double x = -L; x < L; x += 0.1)
    {
        for (double y = -L; y < L; y += 0.1)
        {
            K h = v(one, two, x, y);
            stream << x << " " << y << " " << h.first << " " << h.second << endl;
        }
    }
    file.close();
}

Sistem& Sistem::matrika(const QString& filename)
{
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            stream << i << " " << j << " " << gsl_matrix_get(A, i, j) << " " << gsl_matrix_get(B, i, j) << endl;
        }
        stream << endl;
    }
    
    file.close();
    return *this;
}

Sistem& Sistem::preveri()
{
    // Najprej preverimo resite sistema
    gsl_blas_dgemv(CblasNoTrans, 1, A, resitev, -1, b);
    
    double r;
    gsl_blas_ddot(b, b, &r);
    qDebug() << "Residual =" << r;
    
    // Potem pa se hitrosti po povrsini
    
    for (int i = 0; i < N; ++i)
    {
        const QPointF p = center(delitev, i);
        K h = hitrost(p.x(), p.y());
        rotate(kot(i), &h);
        qDebug() << "Tok [" << i << "] =" << h.first << h.second;
    }
    
    return *this;
}

void vse_trak()
{
    
    QString naboj = fileName("trak_naboj_%1");
    QString potencial = fileName("trak_potencial_%1");
    
    
    QFile file(fileName("kapaciteta"));
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    
    for (int i = 100; i <= 2000; i += 100)
    {
        Sistem s = elektroda(trak(i)).resi().el_polje(potencial.arg(i)).naboj(naboj.arg(i));
        stream << i << " " << s.kapaciteta() << endl;
    }
    
    file.close();
}

Sistem& Sistem::naboj(const QString& filename)
{
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    
    for (int i = 0; i < N; ++i)
    {
        stream << center(delitev, i).x() << " " << gsl_vector_get(resitev, i) << endl;
    }
    stream << endl;
    
    file.close();
    return *this;
}

void naredi_obtekanje(const QPolygonF& delitev, const QString& ime)
{
    QString telo = fileName("%1_telo");
    QString polje = fileName("%1_polje");
    QString tok = fileName("%1_tok");
    
    obtekanje(delitev).resi().shrani(telo.arg(ime)).hitrostno_polje(polje.arg(ime)).tokovnice(tok.arg(ime));
    qDebug() << "Naredil" << ime;
}

double Sistem::kapaciteta()
{
    double C = 0;
    for (int i = 0; i < N; ++i)
    {
        C += gsl_vector_get(resitev, i);
    }
    return C * 0.5 / N;
}

double analiticna(double b, double x, double y)
{
    return (1+b)*y / sqrt(y*y + b*b*b*b*x*x);
}

void tangencialna(int n)
{
    double b = 0.4;
    Sistem s = obtekanje(elipsoid(n, b));
    s.resi();
    
    QFile file(fileName(QString("tangencialna_%1").arg(n)));
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    
    for (int i = 0; i < n; ++i)
    {
        const QPointF p = center(s.delitev, i);
        K v = s.hitrost(p.x(), p.y());
        rotate(s.kot(i), &v);
        
        stream << p.x() << " " << p.y() << " " << -v.first << " " << v.second << " " << analiticna(b, p.x(), p.y()) << endl;
    }
    
    file.close();
}

int main(int argc, char **argv) {
    vse_trak();
    
    /*
    naredi_obtekanje(elipsoid(100, 0.2), "elipsoid");
    naredi_obtekanje(naca(100, 15, 0), "naca");
    naredi_obtekanje(naca(100, 15, -20), "naca-r");
    naredi_obtekanje(zukovski(100, -0.2, 0.1), "zukovski");
    */
    
    return 0;
}
