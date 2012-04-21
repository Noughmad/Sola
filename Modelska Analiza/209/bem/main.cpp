#include <iostream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

#include <QtGui/QPolygonF>
#include <QtGui/qtransform.h>
#include <QtCore/qmath.h>
#include <QtCore/QPair>
#include <QtCore/qline.h>
#include <QtCore/QFile>
#include <QtCore/QDebug>
#include <QtCore/qtextstream.h>

struct Sistem
{
    gsl_matrix* A;
    gsl_vector* b;
    QPolygonF delitev;
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
    poly[0] = QPointF(0,0);
    
    double F = 2.75;
    
    double k = F / n;
    double m = F / 3.0;
    double l = (double)2*(3.0-F) / n;
    
    for (int i ; i <= n; ++i)
    {
        double x = 1 - cos(M_PI * (i-n/2) * 1.0 / n);
        if (i < n/2)
        {
            poly[i] = QPointF(x, t / 50.0 * n_y(x));
        }
        else
        {
            poly[i] = QPointF(x, -t / 50.0 * n_y(x));
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
    return 0.25 * M_1_PI * (2*x1 - 2*x2 + 2*y*(atan2(x1, y) - atan2(x2, y)) + x1 * log(x1*x1 + y*y) - x2 * log(x2*x2 + y*y));
}

inline double u(const QPointF& xy, const QPair<double, double>& center)
{
    return u(xy.x(), xy.y(), center.first, center.second);
}

void transform(const QPointF& p1, const QPointF& p2, double* x, double* y)
{
    double angle = -M_PI / 180.0 * QLineF(p1, p2).angle();
    
    double xn = (*x) * cos(angle) + (*y) * sin(angle);
    double yn = (*y) * cos(angle) - (*x) * sin(angle);
    *x = xn;
    *y = yn;
}

double v_T(const QPointF& ap1, const QPointF& ap2, double x, double y)
{
    QPointF c = (ap1+ap2)/2;
    QTransform t;
    t.translate(-c.x(), -c.y());
    QPointF p1 = t.map(ap1);
    QPointF p2 = t.map(ap2);
    
    x -= c.x();
    y -= c.y();
    
    
    double ksi1 = p1.x();
    double ksi2 = p2.x();
    double eta1 = p1.y();
    double eta2 = p2.y();
    
    transform(p1, p2, &x, &y);
    transform(p1, p2, &ksi1, &eta1);
    transform(p1, p2, &ksi2, &eta2);
    
    Q_ASSERT(qFuzzyIsNull(eta1));
    Q_ASSERT(qFuzzyIsNull(eta2));
    
    return 0.5 * M_1_PI * (atan2(x-ksi1, y) - atan2(x-ksi1, y));
}

double v_II(const QPointF& p1, const QPointF& p2, double x, double y)
{
    double ksi1 = p1.x();
    double ksi2 = p2.x();
    double eta1 = p1.y();
    double eta2 = p2.y();
    
    transform(p1, p2, &x, &y);
    transform(p1, p2, &ksi1, &eta1);
    transform(p1, p2, &ksi2, &eta2);
    
    Q_ASSERT(qFuzzyIsNull(eta1));
    Q_ASSERT(qFuzzyIsNull(eta2));
    
    double x1 = x - ksi1;
    double x2 = x - ksi2;
    
    return 0.25 * M_1_PI * log( (x1*x1 + y*y) / (x2*x2 + y*y) );
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
            gsl_matrix_set(matrika, i, j, u(center(delitev, i), qMakePair(delitev[i].x(), delitev[i+1].x()) ));
        }
        gsl_vector_set(rhs, i, 1);
    }
    
    Sistem s;
    s.A = matrika;
    s.b = rhs;
    s.delitev = delitev;
    return s;
}

Sistem obtekanje(const QPolygonF& delitev)
{
    int n = delitev.size() - 1;
    Sistem s;
    s.A = gsl_matrix_alloc(n, n);
    s.b = gsl_vector_alloc(n);
    s.delitev = delitev;
    
    for (int i = 0; i < n; ++i)
    {
        gsl_matrix_set(s.A, i, i, -0.5);
        for (int j = 0; j < i; ++j)
        {
            QPointF p = center(delitev, j);
            double e = v_T(delitev[i], delitev[i+1], p.x(), p.y());
            gsl_matrix_set(s.A, i, j, e);
            gsl_matrix_set(s.A, j, i, e);
        }
        gsl_vector_set(s.b, i, -sin(M_PI / 180.0 * QLineF(delitev[i], delitev[i+1]).angle() ));
    }
    return s;
}

void resi(const Sistem& sistem, const QString& filename)
{
    gsl_vector* x = gsl_vector_alloc(sistem.b->size);
    gsl_linalg_HH_solve(sistem.A, sistem.b, x);
    
    if (!filename.isEmpty())
    {
        int n = sistem.delitev.size() - 1;
        QFile file(filename);
        file.open(QIODevice::WriteOnly);
        QTextStream stream(&file);
        stream << sistem.delitev[0].x() << " " << sistem.delitev[0].y() << " " << gsl_vector_get(x, 0) << endl;
        for (int i = 1; i < n; ++i)
        {
            stream << sistem.delitev[i].x() << " " << sistem.delitev[i].y() << " " << gsl_vector_get(x, i-1) << endl;
            stream << sistem.delitev[i].x() << " " << sistem.delitev[i].y() << " " << gsl_vector_get(x, i) << endl;
        }
        stream << sistem.delitev[n].x() << " " << sistem.delitev[n].y() << " " << gsl_vector_get(x, n-1) << endl;
        
        file.close();
    }
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

int main(int argc, char **argv) {
  //  resi(obtekanje(naca(120, 15, -10)), "g_naca_15.dat");
 //   save(zukovski(120, -0.1, 0.1), "g_slika_zukovski.dat");
 //   resi(obtekanje(zukovski(120, -0.1, 0.025)), "g_zukovski_0.5_0.dat");
    resi(obtekanje(elipsoid(500, 0.3)), "g_elipsoid.dat");
    return 0;
}
