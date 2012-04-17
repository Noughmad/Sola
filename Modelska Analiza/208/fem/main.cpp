#include <iostream>
#include <fstream>
#include <ctime>

#include <QtCore/qmath.h>
#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QDebug>

#include "delitev.h"
#include "galerkin.h"

using namespace std;

Delitev* srediscna(int n)
{
    double h = 1.0/n;
    Delitev* d = new Delitev();
    QList<int> zadnja;
    zadnja << d->dodaj_tocko(0,0, false);
    
    for (int i = 0; i < n; ++i)
    {
        QList<int> nova;
        double r = (i+1)*h;
        int nf = 3*(i+1);
        for (int j = 0; j <= nf; ++j)
        {
            double f = j * M_PI / nf;
            double x = r*cos(f);
            double y = r*sin(f);
            nova << d->dodaj_tocko(x, y, i < n-1 && j != 0 && j != nf);
        }
        for (int j = 0; j < i+1; ++j)
        {
            d->dodaj_trikotnik(nova[j], nova[j+1], zadnja[j]);
        }
        for (int j = 1; j < i+1; ++j)
        {
            d->dodaj_trikotnik(nova[j], zadnja[j-1], zadnja[j]);
        }
        for (int j = i+1; j < 2*i+2; ++j)
        {
            d->dodaj_trikotnik(nova[j], nova[j+1], zadnja[j-1]);
        }
        for (int j = i+2; j < 2*i+2; ++j)
        {
            d->dodaj_trikotnik(nova[j], zadnja[j-2], zadnja[j-1]);
        }
        for (int j = 2*i+2; j < 3*i+3; ++j)
        {
            d->dodaj_trikotnik(nova[j], nova[j+1], zadnja[j-2]);
        }
        for (int j = 2*i+3; j < 3*i+3; ++j)
        {
            d->dodaj_trikotnik(nova[j], zadnja[j-3], zadnja[j-2]);
        }
        zadnja = nova;
    }
    
    d->rect = QRectF(-1, 0, 2, 1);
    d->name = "srediscna";
    
    return d;
}

double xc(double h)
{
    return sqrt(1 - h*h);
}

int meja(int i, int n)
{
    if (i < n/3)
    {
        return n-i;
    }
    else if (i < 2*n/3)
    {
        return 2*n/3;
    }
    else
    {
        return i;
    }
}

Delitev* batman(int n)
{
    Delitev* d = new Delitev();
    double h = 1.0/n;
    
    QList<int> zadnja;
    for (int i = 0; i <= n; ++i)
    {
        QList<int> nova;
        nova << d->dodaj_tocko(i*h, 0, false);
        int mj = meja(i,n);
        for (int j = 1; j < mj; ++j)
        {
            nova << d->dodaj_tocko(i*h, j*h, i > 0 && i < n);
        }
        nova << d->dodaj_tocko(i*h, mj*h, false);
        
        int m = 2*n/3;
        if (i > 0 && i > n/2)
        {
            for (int j = 0; j < m; ++j)
            {
                d->dodaj_trikotnik(nova[j], nova[j+1], zadnja[j]);
                d->dodaj_trikotnik(nova[j+1], zadnja[j], zadnja[j+1]);
            }
        }
        else if (i > 0)
        {
            for (int j = 0; j < m; ++j)
            {
                d->dodaj_trikotnik(nova[j], nova[j+1], zadnja[j+1]);
                d->dodaj_trikotnik(nova[j], zadnja[j], zadnja[j+1]);
            }            
        }
        
        if (i > 0 && i <= n/3)
        {
            for (int j = m; j < mj; ++j)
            {
                d->dodaj_trikotnik(nova[j], nova[j+1], zadnja[j+1]);
            }
            for (int j = m; j <= mj; ++j)
            {
                d->dodaj_trikotnik(nova[j], zadnja[j], zadnja[j+1]);
            }
        }
        else if (i > 0 && i > 2*n/3)
        {
            for (int j = m; j < mj; ++j)
            {
                d->dodaj_trikotnik(nova[j], nova[j+1], zadnja[j]);
            }
            for (int j = m; j < mj-1; ++j)
            {
                d->dodaj_trikotnik(nova[j+1], zadnja[j], zadnja[j+1]);
            }
        }
        zadnja = nova;
    }
    
    Q_ASSERT(d->st_trikotnikov() == 2*7*n*n/9);
    
    d->rect = QRectF(0, 0, 1, 1);
    d->name = "batman";
    
    return d;
}

void vse_n(Generator g, const QString& filename)
{
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    
    QString vrstica = "%1 %2 %3 %4\n";
    
    const int NV = 8;
    
    for (int n = 10; n < 70; n += 3)
    {
        clock_t start = clock();
        
        Delitev* d = g(n);
        double vrednosti[NV];
        d->resi_nihanje(NV, false, vrednosti);
        
        double t = (double)(clock() - start)/CLOCKS_PER_SEC;
        stream << vrstica.arg(n)
                         .arg(d->st_notranjih())
                         .arg(vrednosti[0])
                         .arg(t);
        cout << "Resil za n = " << n << ", tock = " << d->st_notranjih() << ", lambda_0 = " << vrednosti[0] << endl;
    }
    file.close();
}

void vse_p(Generator g, const QString& filename)
{
    QFile file(filename);
    file.open(QIODevice::WriteOnly);
    QTextStream stream(&file);
    
    QString vrstica = "%1 %2 %3 %4\n";
    
    for (int n = 9; n < 300; n += 9)
    {
        clock_t start = clock();
        Delitev* d = g(n);
        d->resi_poisson(false);
        double phi = d->pretok();
        double t = (double)(clock() - start)/CLOCKS_PER_SEC;
        
        stream << vrstica.arg(n)
                         .arg(d->st_notranjih())
                         .arg(phi)
                         .arg(t);
        cout << "Resil za n = " << n << ", tock = " << d->st_notranjih() << ", pretok = " << phi << endl;
    }
    file.close();
}

Delitev* galerkin(int m, int n)
{
    Galerkin* g = new Galerkin(m, n);
    return g;
}

int main(int argc, char **argv) {
    double lv[5];
    galerkin(1, 10)->resi_nihanje(5, true, lv);
    return 0;
}
