#include <iostream>

#include "delitev.h"
#include "QtCore/qmath.h"

#include <cholmod.h>
#include <QtGui/QImage>
#include <QtGui/QPainter>

using namespace std;

bool operator==(const Trikotnik& t1, const Trikotnik& t2)
{
    return (t1.i == t2.i) && (t1.j == t2.j) && (t1.k == t2.k);
}

quint16 qHash(const Trikotnik& t)
{
    return (t.i << 10) | (t.j << 5) | t.k;
}

Tocka::operator QPointF()
{
    return QPointF(400 + 400 * x, 400 * y);
}


Delitev::Delitev() : cc(0)
{

}

Delitev::~Delitev()
{

}

int Delitev::indeks(int i, int j) const
{
    Q_ASSERT(tocke[i].noter);
    Q_ASSERT(tocke[j].noter);
    Q_ASSERT(not_indeksi.contains(i));
    Q_ASSERT(not_indeksi.contains(j));
    
    int ni = not_indeksi[i];
    int nj = not_indeksi[j];
    return ni * not_indeksi.size() + nj;
}

double Delitev::ploscina(const Trikotnik& t)
{
    if (ploscine.contains(t))
    {
        return ploscine.value(t);
    }

    double s = 0;
    s += tocke[t.i].x * y(t.j,t.k);
    s += tocke[t.k].x * y(t.i,t.j);
    s += tocke[t.j].x * y(t.k,t.i);
    s = fabs(s)/2;
    ploscine.insert(t, s);
    return s;
}

double Delitev::dd(int i, int j) const
{
    return x(i,j)*x(i,j) + y(i,j)*y(i,j);
}

double Delitev::x(int i, int j) const
{
    return tocke[i].x - tocke[j].x;
}

double Delitev::y(int i, int j) const
{
    return tocke[i].y - tocke[j].y;
}

int Delitev::dodaj_tocko(double x, double y)
{
    const int i = tocke.size();
    
    Tocka t;
    t.x = x;
    t.y = y;
    if (qFuzzyIsNull(y) || (x*x+y*y) > 0.999)
    {
        t.noter = false;
        t.z = 0;

    }
    else
    {
        t.noter = true;
        not_indeksi.insert(i,notranje.size());
        notranje << i;
    }
    tocke << t;
        
    return i;
}

int Delitev::dodaj_trikotnik(int i, int j, int k)
{
    Trikotnik t = {i,j,k};
    trikotniki << t;
    return trikotniki.size() - 1;
}

double Delitev::pretok()
{
    double phi = 0;
    foreach (const Trikotnik& t, trikotniki)
    {
        phi += (tocke[t.i].z + tocke[t.j].z + tocke[t.k].z) * ploscina(t) / 3;
    }
    return phi;
}

cholmod_dense* Delitev::desne_strani()
{
    int n = notranje.size();
    cholmod_dense* b = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, cc);
    double* bx = (double*)b->x;
    foreach (const Trikotnik& t, trikotniki)
    {
        const double s = ploscina(t)/3;
        if (tocke[t.i].noter)
        {
            bx[not_indeksi[t.i]] += s;
        }
        if (tocke[t.j].noter)
        {
            bx[not_indeksi[t.j]] += s;
        }
        if (tocke[t.k].noter)
        {
            bx[not_indeksi[t.k]] += s;
        }
    }
    return b;
}

cholmod_sparse* Delitev::matrika()
{
    int n = notranje.size();
    int N = 7*n;
    cholmod_dense* D = cholmod_allocate_dense(n, n, n, CHOLMOD_REAL, cc);
    // cholmod_sparse* A = cholmod_allocate_sparse(n, n, N, 1, 1, 0, CHOLMOD_REAL, c);
    
    double* Dx = (double*)D->x;
    for (int i = 0; i < n*n; ++i)
    {
        Dx[i] = 0;
    }
        
    foreach (const Trikotnik& t, trikotniki)
    {
        const int& i = t.i;
        const int& j = t.j;
        const int& k = t.k;
        
        bool ni = tocke[i].noter;
        bool nj = tocke[j].noter;
        bool nk = tocke[k].noter;
        
        double s = 0.25 / ploscina(t);
        
        if (ni)
        {
            Dx[indeks(t.i, t.i)] += dd(t.j, t.k) * s;
        }
        
        if (nj)
        {
            Dx[indeks(t.j, t.j)] += dd(t.k, t.i) * s;
        }
        
        if (nk)
        {
            Dx[indeks(t.k, t.k)] += dd(t.i, t.j) * s;
        }
        
        if (ni && nj)
        {
            const double ts = ( x(j,k)*x(k,i) + y(j,k)*y(k,i) ) * s;
            Dx[indeks(i,j)] += ts;
            Dx[indeks(j,i)] += ts;
        }
        
        if (nj && nk)
        {
            const double ts = ( x(j,i)*x(i,k) + y(j,i)*y(i,k) ) * s;
            Dx[indeks(j,k)] += ts;
            Dx[indeks(k,j)] += ts;
        }
        
        if (ni && nk)
        {
            const double ts = ( x(i,j)*x(j,k) + y(i,j)*y(j,k) ) * s;
            Dx[indeks(k,i)] += ts;
            Dx[indeks(i,k)] += ts;
        }
    }
    
    cholmod_sparse* A = cholmod_dense_to_sparse(D, 1, cc);
    cholmod_free_dense(&D, cc);
    return A;
}

void Delitev::resi_poisson()
{
    cholmod_common Common;
    cc = &Common;
    
    cholmod_start(cc);
    
    cholmod_sparse* A = matrika();
    cholmod_dense* b = desne_strani();
    
    cholmod_factor* L = cholmod_analyze(A, cc);
    cholmod_factorize(A, L, cc);
    cholmod_dense* x = cholmod_solve(CHOLMOD_A, L, b, cc);
    double* xx = (double*)x->x;
    
    int n = x->nrow;
    Q_ASSERT(n == notranje.size());
    
    for (int i = 0; i < n; ++i)
    {
        tocke[notranje[i]].z = xx[i];
    }
}

void Delitev::narisi(const QString& file)
{
    QImage image(801, 401, QImage::Format_ARGB32);
    image.fill(Qt::white);
    QPainter painter;
    painter.begin(&image);
        
    foreach (const Trikotnik& t, trikotniki)
    {
        painter.setPen(QPen(Qt::black, 1));
        painter.drawLine(tocke[t.i], tocke[t.j]);
        painter.drawLine(tocke[t.j], tocke[t.k]);
        painter.drawLine(tocke[t.k], tocke[t.i]);
        
        painter.setPen(QPen(Qt::red, 2));
        painter.drawPoint(QPointF(tocke[t.i] + tocke[t.j] + tocke[t.k])/3);
    }
    painter.end();
    image.save(file);
}
