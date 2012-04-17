#include <iostream>

#include "delitev.h"
#include "QtCore/qmath.h"

#include <cholmod.h>
#include <QtGui/QImage>
#include <QtGui/QPainter>
#include <QtCore/QMap>

#include <QtCore/QDebug>
#include <QFile>




using namespace std;

bool operator==(const Trikotnik& t1, const Trikotnik& t2)
{
    return (t1.i == t2.i) && (t1.j == t2.j) && (t1.k == t2.k);
}

quint16 qHash(const Trikotnik& t)
{
    return (t.i << 10) | (t.j << 5) | t.k;
}

Tocka::operator QPointF() const
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

int Delitev::dodaj_tocko(double x, double y, bool noter)
{
    const int i = tocke.size();
    
    Tocka t;
    t.x = x;
    t.y = y;
    if (!noter)
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
    int ti, tj, tk;
    if (i < qMin(j,k))
    {
        ti = i;
        tj = qMin(j,k);
        tk = qMax(j,k);
    }
    else if (j < qMin(i,k))
    {
        ti = j;
        tj = qMin(i,k);
        tk = qMax(i,k);
    }
    else 
    {
        ti = k;
        tj = qMin(i,j);
        tk = qMax(i,j);
    }
    Trikotnik t = {ti,tj,tk};
    
    double s = 0;
    s += tocke[t.i].x * y(t.j,t.k);
    s += tocke[t.k].x * y(t.i,t.j);
    s += tocke[t.j].x * y(t.k,t.i);
    t.ploscina = fabs(s)/2;
    
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
    
    for (int i = 0; i < n; ++i)
    {
        bx[i] = 0;
    }
    
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

Matrika Delitev::matrika()
{
    int n = notranje.size();
    
    Matrika elementi;
        
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
            elementi[t.i][t.i] += dd(t.j, t.k) * s;
        }
        
        if (nj)
        {
            elementi[t.j][t.j] += dd(t.k, t.i) * s;
        }
        
        if (nk)
        {
            elementi[t.k][t.k] += dd(t.i, t.j) * s;
        }
        
        if (ni && nj)
        {
            elementi[i][j] += ( x(j,k)*x(k,i) + y(j,k)*y(k,i) ) * s;
        }
        
        if (nj && nk)
        {
            elementi[j][k] += ( x(j,i)*x(i,k) + y(j,i)*y(i,k) ) * s;
        }
        
        if (ni && nk)
        {
            elementi[i][k] += ( x(i,j)*x(j,k) + y(i,j)*y(j,k) ) * s;
        }
    }
        
    Q_ASSERT(elementi.size() == n);
        
    return elementi;
}

cholmod_sparse* Delitev::sparse(const Matrika& matrika, bool symmetric)
{    
    Matrika elementi(matrika);
        
    if (!symmetric)
    {
        // Insert the other half of values
        for (Matrika::const_iterator git = matrika.constBegin(); git != matrika.constEnd(); ++git)
        {
            Vrstica::const_iterator end = git.value().constEnd();
            Vrstica::const_iterator it = git.value().constBegin();
            for (; it != end; ++it)
            {
                if (it.key() != git.key())
                {
                    elementi[it.key()][git.key()] = it.value();
                }
            }
        }
    }
    
    Matrika::const_iterator gend = elementi.constEnd();
    Matrika::const_iterator git = elementi.constBegin();
    
    int n = notranje.size();
    int m = 0;
    foreach (const Vrstica& vrstica, elementi)
    {
        m += vrstica.size();
    }
    
    cholmod_sparse* A = cholmod_allocate_sparse(n, n, m, 1, 1, symmetric ? -1 : 0, CHOLMOD_REAL, cc);
    
    double* Ax = (double*)A->x;
    int* Ap = (int*)A->p;
    int* Ai = (int*)A->i;
    
    int t = 0;
    int i = 0;
    
    
    for (git = elementi.constBegin(); git != gend; ++git)
    {
        Ap[i] = t;
        Vrstica::const_iterator end = git.value().constEnd();
        Vrstica::const_iterator it = git.value().constBegin();
        for (; it != end; ++it)
        {          
            Ai[t] = not_indeksi[it.key()];
            Ax[t] = it.value();
            ++t;
        }
        ++i;
    }
    Ap[n] = t;
    Q_ASSERT(t == m);
    
    return A;
}

void Delitev::resi_poisson(bool risi)
{
    cholmod_common Common;
    cc = &Common;
    
    cholmod_start(cc);
    Common.dtype = CHOLMOD_DOUBLE;
        
    cholmod_sparse* A = sparse(matrika(), true);
    cholmod_dense* b = desne_strani();
    
   // cholmod_print_sparse(A, "A", cc);
        
    cholmod_factor* L = cholmod_analyze(A, cc);
    cholmod_factorize(A, L, cc);
    cholmod_dense* x = cholmod_solve(CHOLMOD_A, L, b, cc);
    
    cholmod_free_sparse(&A, cc);
    cholmod_free_dense(&b, cc);
    
    double* xx = (double*)x->x;
    
    int n = x->nrow;
    Q_ASSERT(n == notranje.size());
    
    for (int i = 0; i < n; ++i)
    {
        tocke[notranje[i]].z = xx[i];
    }
    
    if (risi)
    {
        plot(x, 0, 0);
    }
    
    cholmod_free_dense(&x, cc);
    
    cholmod_finish(cc);
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
    
    painter.setPen(QPen(Qt::green, 2));
    foreach (const Tocka& t, tocke)
    {
        if (t.noter)
        {
            continue;
        }
        painter.drawPoint(t);
    }
    
    painter.end();
    image.save(file);
}

int Delitev::st_notranjih() const
{
    return notranje.size();
}

int Delitev::st_tock() const
{
    return tocke.size();
}

int Delitev::st_trikotnikov() const
{
    return trikotniki.size();
}

Matrika Delitev::masa()
{
    Matrika elementi;
    foreach (const Trikotnik& t, trikotniki)
    {
        const int& i = t.i;
        const int& j = t.j;
        const int& k = t.k;
        
        bool ni = tocke[i].noter;
        bool nj = tocke[j].noter;
        bool nk = tocke[k].noter;
        
        double s = ploscina(t) / 12.0;
        
        if (ni)
        {
            elementi[t.i][t.i] += 2 * s;
        }
        
        if (nj)
        {
            elementi[t.j][t.j] += 2 * s;
        }
        
        if (nk)
        {
            elementi[t.k][t.k] += 2 * s;
        }
        
        if (ni && nj)
        {
            elementi[i][j] += s;
        }
        
        if (nj && nk)
        {
            elementi[j][k] += s;
        }
        
        if (ni && nk)
        {
            elementi[i][k] += s;
        }
    }
    return elementi;
}

void Delitev::shrani(const QString& file)
{
    QFile f(file);
    f.open(QIODevice::WriteOnly);
    QTextStream stream(&f);
    foreach (const Tocka& t, tocke)
    {
        stream << t.x << " " << t.y << " " << t.z << endl;
    }
    f.close();
}


