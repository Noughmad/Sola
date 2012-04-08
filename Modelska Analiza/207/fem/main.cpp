#include <iostream>
#include <fstream>
#include <stdio.h>

#include <QtCore/qmath.h>
#include <QtCore/qglobal.h>

#include <QtCore/QList>
#include <QtCore/QPair>

#include <QtGui/QImage>
#include <QtGui/QPainter>

#include <SuiteSparseQR.hpp>

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

/**
 * Tri tocke s ciklicnimi indeksi
 **/
class Trikotnik
{
public:
    Trikotnik(Tocka* a, Tocka* b, Tocka* c)
    {
        tocke[0] = 0;
        tocke[1] = b;
        tocke[2] = c;
    }
    
    Tocka& operator[](int i)
    {
        return *(tocke[i % 3]);
    }
    
    const Tocka& operator[](int i) const
    {
        return *(tocke[i % 3]);
    }
    
    double ploscina()
    {
    }
        
private:
    Tocka* tocke[3];
};

class Delitev
{
public:
    Delitev() : povezave(0)
    {
        
    }
    
    ~Delitev()
    {
        gsl_matrix_char_free(povezave);
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
        double s = 0;
        s += tocke[i].x * ( tocke[j].y - tocke[k].y );
        s += tocke[k].x * ( tocke[i].y - tocke[j].y );
        s += tocke[j].x * ( tocke[k].y - tocke[i].y );
        if (gsl_isnan(s) || gsl_isinf(s) || s == 0)
        {
            printf("%d, %d, %d => S = %g\n", i, j, k, s);
        }
        return fabs(s)/2;
    }
    
    double dd(int i, int j)
    {
        return pow(tocke[i].x - tocke[j].x, 2) + pow(tocke[i].y - tocke[j].y, 2);
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
        if (i == j)
        {
            for (int k = 0; k < n; ++k)
            {
                for (int l = 0; l < k; ++l)
                {
                    if (povezano(i,k,l))
                    {
                        a += dd(k,l) / 4 / ploscina(i,k,l);
                    }
                }
            }
        }
        else if (gsl_matrix_char_get(povezave, i, j) && noter(i) && noter(j))
        {
            for (int k = 0; k < n; ++k)
            {
                if (gsl_matrix_char_get(povezave, i, k) && gsl_matrix_char_get(povezave, j, k))
                {
                    a += ( x(j,k)*x(k,i) + y(j,k)*y(k,i) ) / 4 / ploscina(i,j,k);
                }
            }
        }
        return a;
    }
    
    gsl_matrix* matrika()
    {
        gsl_matrix* m = gsl_matrix_alloc(n,n);
        for (int i = 0; i < n; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                gsl_matrix_set(m, i, j, element(i,j));
            }
        }
        return m;
    }
    
    bool noter(int i)
    {
        const Tocka& t = tocke[i];
        if (qFuzzyIsNull(t.y))
        {
            return false;
        }
        double rr = t.x*t.x + t.y*t.y;
        return rr < 0.99;
    }
    
    bool povezano(int i, int j, int k)
    {
        if (i == j || j == k || k == i)
        {
            return false;
        }
        return gsl_matrix_char_get(povezave, i, j) && gsl_matrix_char_get(povezave, j, k) && gsl_matrix_char_get(povezave, k, i);
    }
    
    gsl_vector* desne()
    {
        gsl_vector* b = gsl_vector_alloc(n);
        gsl_vector_set_zero(b);
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
                        gsl_vector_set(b, i, gsl_vector_get(b, i) + ploscina(i, j, k));
                    }
                }
            }
        }
        return b;
    }
    
    QList<int> sosede(int i)
    {
        QList<int> list;
        for (int j = 0; j < n; ++j)
        {
            if (gsl_matrix_char_get(povezave, i, j))
            {
                list << j;
            }
        }
    }
    
    gsl_vector* resitev(const char* filename)
    {
        gsl_vector* v = desne();
        gsl_matrix* m = matrika();
        
        gsl_vector_fprintf(stdout, v, "%g");
        
        shrani_matriko(m, n, "g_matrika_A.dat");

        int s;
        gsl_permutation * p = gsl_permutation_alloc (n);
        gsl_linalg_LU_decomp (m, p, &s);
        gsl_linalg_LU_svx (m, p, v);
       
        FILE *f = fopen(filename, "wt");
        for (int i = 0; i < n; ++i)
        {
            const Tocka& t = tocke[i];
            fprintf(f, "%g %g %g\n", t.x, t.y, gsl_vector_get(v, i));
        }
        fclose(f);
        
        return v;
    }
    
    void popravi_krog(double h)
    {
        for (int i = 0; i < tocke.size(); ++i)
        {
            Tocka& t = tocke[i];
            double r2 = t.x * t.x + t.y * t.y;
            if (r2 > 1+1.5*h)
            {
                tocke.removeAt(i);
                --i;
            }
            else if (r2 > 1)
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
    
    void narisi(const QString& filename)
    {
        QImage image(800, 400, QImage::Format_ARGB32);
        image.fill(Qt::white);
        QPainter painter;
        painter.begin(&image);
        
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
            
            for (int j = 0; j < i; ++j)
            {
                if (gsl_matrix_char_get(povezave, i, j))
                {
                    QPointF q = QPointF(400 + 400 * tocke[j].x, 400 * tocke[j].y);
                    painter.drawLine(p, q);
                }
            }
        }
        
        painter.end();
        image.save(filename);
    }

private:
    QList<Tocka> tocke;
    gsl_matrix_char* povezave;
    int n;
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
    }
    
    d.init();
    d.popravi_krog(h);
    d.samo_povezi(h*h*2.05);
    d.narisi("g_povezave_20.png");
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
            d.dodaj_tocko((i*0.5-j)*h, i*v);
            d.dodaj_tocko((-i+0.5*(j-1))*h, (j-1)*v);
        }
    }
    
    d.popravi_krog(h);

    d.init();
    d.samo_povezi(h*h*1.2);
    d.narisi("g_povezave_hex.png");
    return d;
}

int main(int argc, char **argv) {
    /*
    srediscna(20).resitev("g_srediscna_20.dat");
//    srediscna(10).resitev("g_srediscna_10.dat");
*/
    heksagonalna(20).resitev("g_hex_20.dat");
    return 0;
}
