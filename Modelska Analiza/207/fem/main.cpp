#include <iostream>
#include <fstream>
#include <stdio.h>
#include <ctime>
#include <cassert>

#include <QtCore/qmath.h>
#include <QtCore/qglobal.h>

#include <QtCore/QList>
#include <QtCore/QQueue>
#include <QtCore/QHash>
#include <QtCore/QPair>

#include <QtGui/QImage>
#include <QtGui/QPainter>

#include <SuiteSparseQR.hpp>
#include <cholmod.h>

#include <gsl/gsl_matrix_char.h>
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sys.h>

using namespace std;

void shrani_matriko(gsl_matrix* m, int n, const char* filename)
{
    ofstream out;
    out.open(filename);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            out << gsl_matrix_get(m, i, j) << " ";
        }
        out << std::endl;
    }
    out.close();
}

struct Tocka
{
    double x,y;
};

struct Trikotnik
{
    int i,j,k;
};

quint16 qHash(const Trikotnik& t)
{
    return (t.i << 8) ^ (t.j << 4) ^ t.k;
}

bool operator==(const Trikotnik& t1, const Trikotnik& t2)
{
    return (t1.i == t2.i) && (t1.j == t2.j) && (t1.k == t2.k);
}

class Delitev
{
public:
    
    
    Delitev() : povezave(0)
    {
        
    }
    
    void obrni()
    {
        QList<Tocka> novi;
        while (!tocke.isEmpty())
        {
            novi << tocke.takeLast();
        }
        tocke = novi;
    }
    
    ~Delitev()
    {
        gsl_matrix_char_free(povezave);
    }
    
    double pretok(const cholmod_dense* v)
    {
        double phi;
        double* vx = (double*)v->x;
        for (int i = 0; i < n; ++i)
        {
            if (!noter(i))
            {
                continue;
            }
            for (int j = 0; j < n; ++j)
            {
                for (int k = 0; k < j; ++k)
                {
                    if (povezano(i,j,k))
                    {
                        phi += vx[i] * ploscina(i,j,k) / 3;
                    }
                }
            }
        }
        return phi;
    }
    
    void init()
    {
        n = tocke.size();
        cout << "Delitev ima " << n << " tock" << endl;
        povezave = gsl_matrix_char_alloc(n,n);
        gsl_matrix_char_set_zero(povezave);
    }
    
    int dodaj_tocko(const Tocka& tocka)
    {
        tocke << tocka;
        return tocke.size() - 1;
    }
    
    int dodaj_tocko(double x, double y)
    {
        Tocka t = {x, y};
        return dodaj_tocko(t);
    }
    
    void konec_vrstice()
    {
        vrstice << tocke.size();
    }
    
    int dodaj_tocko_noter(double x, double y, double h)
    {
        
    }
    
    void vez(int i, int j)
    {
        gsl_matrix_char_set(povezave, i, j, 1);
        gsl_matrix_char_set(povezave, j, i, 1);
    }
    
    double ploscina(int i, int j, int k)
    {
        Trikotnik t = {i,j,k};
        if (ploscine.contains(t))
        {
            return ploscine.value(t);
        }
        
        double s = 0;
        s += tocke[i].x * y(j,k);
        s += tocke[k].x * y(i,j);
        s += tocke[j].x * y(k,i);
        s = fabs(s)/2;
        ploscine.insert(t, s);
        return s;
    }
    
    double dd(int i, int j)
    {
        return x(i,j)*x(i,j) + y(i,j)*y(i,j);
    }
    
    double skupna_ploscina()
    {
        double S = 0;
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < i; ++j)
            {
                if (!gsl_matrix_char_get(povezave, i, j))
                {
                    continue;
                }
                for (int k = 0; k < j; ++k)
                {
                    if (povezano(i,j,k))
                    {
                        S += ploscina(i,j,k);
                    }
                }
            }
        }
        return S;
    }
    
    double x(int i, int j)
    {
        return tocke[i].x - tocke[j].x;
    }
    
    double y(int i, int j)
    {
        return tocke[i].y - tocke[j].y;
    }
    
    double element(int i, int j)
    {
        double a = 0;
        int st = 0;
        if (i == j)
        {
            for (int k = 0; k < n; ++k)
            {
                if (!gsl_matrix_char_get(povezave, i, k))
                {
                    continue;
                }
                for (int l = 0; l < k; ++l)
                {
                    if (povezano(i,k,l))
                    {
                        a += dd(k,l) / 4.0 / ploscina(i,k,l);
                        ++st;
                    }
                }
            }
            if (st > 6)
            {
                cout << "Tocka z vec kot 6 trikotniki " << i << endl;
            }
        }
        else if (gsl_matrix_char_get(povezave, i, j) && noter(i))
        {
            for (int k = 0; k < n; ++k)
            {
                if (gsl_matrix_char_get(povezave, i, k) && gsl_matrix_char_get(povezave, j, k))
                {
                    a += ( x(j,k)*x(k,i) + y(j,k)*y(k,i) ) / 4.0 / ploscina(i,j,k);
                    ++st;
                }
            }
            if (st > 2)
            {
                cout << "Povezava z vec kot 2 trikotnikoma " << i << ", " << j << endl; 
            }
        }
        return a;
    }
    
    cholmod_sparse* matrika(cholmod_common* c)
    {
        int N = 7*n;
        cholmod_sparse* A = cholmod_allocate_sparse(n, n, N, 1, 1, 0, CHOLMOD_REAL, c);
                
        
        double* Ax = (double*)A->x;
        int* Ap = (int*)A->p;
        int* Ai = (int*)A->i;   
        
        int t = 0;
        for (int i = 0; i < n; ++i)
        {
            Ap[i] = t;
            for (int j = 0; j < n; ++j)
            {
                const double E = element(i,j);
                if (!qFuzzyIsNull(E))
                {
                    assert(t < N);
                    Ax[t] = E;
                    Ai[t] = j;
                    ++t;
                }
            }
        }
        Ap[n] = t;
        
        
        Ax = (double*)A->x;
        A->stype = 0;
        return A;
    }
    
    virtual bool noter(int i)
    {
        const Tocka& t = tocke[i];
        if (qFuzzyIsNull(t.y))
        {
            return false;
        }
        double rr = t.x*t.x + t.y*t.y;
        return rr < 0.9999;
    }
    
    bool povezano(int i, int j, int k)
    {
        if (i == j || j == k || k == i)
        {
            return false;
        }
        return gsl_matrix_char_get(povezave, i, j) && gsl_matrix_char_get(povezave, j, k) && gsl_matrix_char_get(povezave, k, i);
    }
    
    cholmod_dense* desne(cholmod_common* c)
    {
        cholmod_dense* b = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, c);
        double* bx = (double*)b->x;
        
        for (int i = 0; i < n; ++i)
        {
            bx[i] = 0;
            if (!noter(i))
            {
                continue;
            }
            for (int j = 0; j < n; ++j)
            {
                if (!gsl_matrix_char_get(povezave, i, j))
                {
                    continue;
                }
                for (int k = 0; k < j; ++k)
                {
                    if (povezano(i,j,k))
                    {
                        bx[i] += ploscina(i, j, k) / 3;
                    }
                }
            }
        }
        
        return b;
    }
    
    
    double resitev(const QString& filename = QString())
    {
        cholmod_common cc;
        cholmod_start(&cc);
        cc.dtype = CHOLMOD_DOUBLE;
        cc.nmethods = CHOLMOD_MAXMETHODS;
        cc.print = 5;
        
        cholmod_dense* b = desne(&cc);
        cholmod_sparse* A = matrika(&cc);
        cholmod_sparse* T = cholmod_copy_sparse(A, &cc);
        double rnorm, one [2] = {1.0/n,0}, minusone [2] = {-1,0} ;

        
        cout << "Norme: |A| = " << cholmod_norm_sparse(A, 1, &cc)/n << ", |b| = " << cholmod_norm_dense(b, 1, &cc)/3 << endl;
        
        
        cholmod_factor* L = cholmod_analyze(A, &cc);
        cholmod_factorize(A, L, &cc);
        cholmod_dense* x = cholmod_solve(CHOLMOD_A, L, b, &cc);
        double* xx = (double*)x->x;
        
        if (!filename.isEmpty())
        {
            FILE *f = fopen(filename.toLatin1(), "wt");
            for (int i = 0; i < n; ++i)
            {
                if (vrstice.contains(i))
                {
                    fprintf(f, "\n");
                }
                const Tocka& t = tocke[i];
                fprintf(f, "%g %g %g\n", t.x, t.y, xx[i]);
            }
            fclose(f);
        }
        
        int nt = 0;
        double phi = 0;
        
        double* bx = (double*)b->x;
        for (int i = 0; i < n; ++i)
        {
            phi += bx[i] * xx[i];
        }
        
        phi /= n;
        
        double S = skupna_ploscina();
        
        cholmod_dense* Residual = cholmod_copy_dense (b, &cc) ;
        cholmod_sdmult (A, 0, minusone, one, x, Residual, &cc) ;
        
     //   cholmod_print_dense(Residual, "Residual", &cc);

        
        cout << "Ploscina je " << S << ", pretok je " << phi << ", norma x je " << cholmod_norm_dense(x, 1, &cc)/n << endl;
                
        cholmod_finish(&cc);
        
        return phi;
    }
    
    void odrezi(double h)
    {
        for (int i = 0; i < tocke.size(); ++i)
        {
            const Tocka& t = tocke[i];
            double r2 = t.x * t.x + t.y * t.y;
            if (r2 > 1+1.8*h)
            {
                tocke.removeAt(i);
                --i;
                for (int j = 0; j < vrstice.size(); ++j)
                {
                    if (vrstice[j] > i)
                    {
                        vrstice[j] = vrstice[j] - 1;
                    }
                }
                
                vrstice << i+1;
            }
        }
    }
    
    virtual void popravi(double h)
    {
        for (int i = 0; i < tocke.size(); ++i)
        {
            Tocka& t = tocke[i];
            double r2 = t.x * t.x + t.y * t.y;
            if (r2 > 1)
            {
                double s = 1.0 / sqrt(r2);
                t.x *= s;
                t.y *= s;
            }
        }
    }
    
    void samo_povezi(double h)
    {
        for (int i = 0; i < n; ++i)
        {
            int si = 0;
            for (int j = 0; j < n; ++j)
            {
                if (gsl_matrix_char_get(povezave, i, j))
                {
                    ++si;
                }
            }
            if (si >= 6)
            {
                continue;
            }
            for (int j = i+1; j < n; ++j)
            {
                if (dd(i,j) < h)
                {
                    vez(i,j);
                }
            }
        }
    }
    
    virtual void narisi(const QString& filename)
    {
        QImage image(801, 401, QImage::Format_ARGB32);
        image.fill(Qt::white);
        QPainter painter;
        painter.begin(&image);
        
        for (int i = 0; i < n; ++i)
        {
            const Tocka& t = tocke[i];
            QPointF p = QPointF(400 + 400 * t.x, 400 * t.y);
            
            for (int j = 0; j < i; ++j)
            {
                if (gsl_matrix_char_get(povezave, i, j))
                {
                    QPointF q = QPointF(400 + 400 * tocke[j].x, 400 * tocke[j].y);
                    painter.drawLine(p, q);
                }
            }
        }
        
        for (int i = 0; i < n; ++i)
        {
            const Tocka& t = tocke[i];
            QPointF p = QPointF(400 + 400 * t.x, 400 * t.y);
            
            if (!noter(i))
            {
                painter.setBrush(Qt::red);
            }
            painter.drawEllipse(p, 2, 2);
            painter.setBrush(Qt::black);
        }
        
        painter.end();
        image.save(filename);
    }
    
    int stevilo()
    {
        return n;
    }
    
protected:
    QList<Tocka> tocke;
    gsl_matrix_char* povezave;
    int n;
    QQueue<int> vrstice;
    QHash<Trikotnik, double> ploscine;
};

class DelitevBatman : public Delitev
{
public:
    virtual bool noter(int i)
    {
        Tocka& t = tocke[i];
        if (t.x < 0.001 || t.x > 0.999 || t.y < 0.001 || t.y > 0.999)
        {
            return false;
        }
        if (t.x < 1.0/3)
        {
            return t.y < 0.999-t.x;
        }
        else if (t.x < 2.0/3)
        {
            return t.y < 0.666;
        }
        else
        {
            return t.y < t.x;
        }
    }
    
    virtual void narisi(const QString& filename)
    {
        QImage image(401, 401, QImage::Format_ARGB32);
        image.fill(Qt::white);
        QPainter painter;
        painter.begin(&image);
        
        for (int i = 0; i < n; ++i)
        {
            const Tocka& t = tocke[i];
            QPointF p = QPointF(400 * t.x, 400 * t.y);
            
            for (int j = 0; j < i; ++j)
            {
                if (gsl_matrix_char_get(povezave, i, j))
                {
                    QPointF q = QPointF(400 * tocke[j].x, 400 * tocke[j].y);
                    painter.drawLine(p, q);
                }
            }
        }
        for (int i = 0; i < n; ++i)
        {
            const Tocka& t = tocke[i];
            QPointF p = QPointF(400 * t.x, 400 * t.y);
            
            if (!noter(i))
            {
                painter.setBrush(Qt::red);
            }
            painter.drawEllipse(p, 2, 2);
            painter.setBrush(Qt::black);
        }
        
        painter.end();
        image.save(filename);
    }
    
};

Delitev srediscna(int n)
{
    double h = 1.0/n;
    Delitev d;
    d.dodaj_tocko(0,0);
    for (int i = 0; i < n; ++i)
    {
        double r = (i+1)*h;
        int nf = 3*(i+1);
        for (int i = 0; i <= nf; ++i)
        {
            double f = i * M_PI / nf;
            double x = r*cos(f);
            double y = r*sin(f);
            d.dodaj_tocko(x,y);
        }
        d.konec_vrstice();
    }
    
    d.init();
    d.popravi(h);
    d.samo_povezi(h*h*2.1);
    return d;
}

Delitev heksagonalna(int n)
{
    double h = 1.0/n;
    double v = h * sqrt(3) / 2;
    Delitev d;
    for (int i = 0; i*v < 1.2; ++i) // Vrstice
    {
        d.dodaj_tocko(i*h, 0);
        for (int j = 1; j < i+1; ++j)
        {
            d.dodaj_tocko((i-0.5*j)*h, j*v);
        }
        for (int j = 1; j < i+1; ++j)
        {
            d.dodaj_tocko((i*0.5-j)*h, i*v);
        }
        for (int j = i; j > 0; --j)
        {
            d.dodaj_tocko((-i+0.5*(j-1))*h, (j-1)*v);
        }
        d.konec_vrstice();
    }
    
    d.odrezi(h);

    d.init();
    d.samo_povezi(h*h*1.2);
    d.popravi(h);
    return d;
}

DelitevBatman batman(int n)
{
    double h = 1.0/n;
    double v = h * sqrt(3) / 2;
    int k = n/3;
    
    DelitevBatman d;
    double y = 0;
    
    QList<int> lineStarters;
    
    for (int i = 0; i < k; ++i)
    {
        lineStarters << d.dodaj_tocko(i*h, 0);
        for (int j = 1; j < n-i+1; ++j)
        {
            d.dodaj_tocko(i*h, j*h);
        }
        d.konec_vrstice();
    }
    for (int i = k; i < 2*k; ++i)
    {
        lineStarters << d.dodaj_tocko(i*h, 0);
        for (int j = 1; j < 2*k+1; ++j)
        {
            d.dodaj_tocko(i*h, j*h);
        }
        d.konec_vrstice();
    }
    for (int i = 2*k; i < n+1; ++i)
    {
        lineStarters << d.dodaj_tocko(i*h, 0);
        for (int j = 1; j < i+1; ++j)
        {
            d.dodaj_tocko(i*h, j*h);
        }
        d.konec_vrstice();
    }
    
    d.init();
    d.samo_povezi(h*h*1.01);
    
    const int m = 2*k;
    for (int i = 0; i < n/2; ++i)
    {
        const int s1 = lineStarters[i];
        const int s2 = lineStarters[i+1];
        for (int j = 0; j < m; ++j)
        {
            d.vez(s1+j, s2+j+1);
        }
    }
    for (int i = n/2; i < n; ++i)
    {
        const int s1 = lineStarters[i];
        const int s2 = lineStarters[i+1];
        for (int j = 0; j < m; ++j)
        {
            d.vez(s1+j+1, s2+j);
        }
    }
    
    for (int i = 0; i < k; ++i)
    {
        const int s1 = lineStarters[i];
        const int s2 = lineStarters[i+1];
        
        const int z1 = lineStarters[n-i];
        const int z2 = lineStarters[n-i-1];
        for (int j = 0; j < k-i; ++j)
        {
            d.vez(s1 + 2*k + j + 1, s2 + 2*k + j);
            d.vez(z1 + 2*k + j + 1, z2 + 2*k + j);
        }
    }
        
    return d;
}

void srediscna_vse()
{
    FILE* f = fopen("g_konv_srediscna.dat", "wt");
    for (int i = 10; i < 80; i += 3)
    {
        Delitev d = srediscna(i);
        int n = d.stevilo();
        clock_t start = clock();
        double r = d.resitev();
        fprintf(f, "%d %d %g %g\n", i, n, r, (double)(clock() - start)/(CLOCKS_PER_SEC));
    }
    fclose(f);
}

void heksa_vse()
{
    FILE* f = fopen("g_konv_hex.dat", "wt");
    for (int i = 10; i < 80; i += 3)
    {
        Delitev d = heksagonalna(i);
        int n = d.stevilo();
        clock_t start = clock();
        double r = d.resitev();
        fprintf(f, "%d %d %g %g\n", i, n, r, (double)(clock() - start)/(CLOCKS_PER_SEC));
    }
    fclose(f);
}

void batman_vse()
{
    FILE* f = fopen("g_konv_batman.dat", "wt");
    for (int i = 12; i < 121; i += 6)
    {
        DelitevBatman d = batman(i);
        int n = d.stevilo();
        clock_t start = clock();
        double r = d.resitev();
        fprintf(f, "%d %d %g %g\n", i, n, r, (double)(clock() - start)/(CLOCKS_PER_SEC));
    }
    fclose(f);
}


int main(int argc, char **argv) {
    srediscna_vse();
    heksa_vse();
    batman_vse();
  /* 
  heksagonalna(10).narisi("g_povezave_hex.png");
  srediscna(10).narisi("g_povezave_sred.png");
  batman(18).narisi("g_batman_18.png");
  batman(30).narisi("g_batman_30.png");
**/
    return 0;
}
