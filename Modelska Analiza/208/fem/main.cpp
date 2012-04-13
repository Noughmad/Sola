#include <iostream>
#include <fstream>
#include <ctime>

#include <QtCore/qmath.h>
#include "delitev.h"

using namespace std;



Delitev srediscna(int n)
{
    double h = 1.0/n;
    Delitev d;
    QList<int> zadnja;
    zadnja << d.dodaj_tocko(0,0);
    
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
            nova << d.dodaj_tocko(x,y);
        }
        for (int j = 0; j < i+1; ++j)
        {
            d.dodaj_trikotnik(nova[j], nova[j+1], zadnja[j]);
        }
        for (int j = 1; j < i+1; ++j)
        {
            d.dodaj_trikotnik(nova[j], zadnja[j-1], zadnja[j]);
        }
        for (int j = i+1; j < 2*i+2; ++j)
        {
            d.dodaj_trikotnik(nova[j], nova[j+1], zadnja[j-1]);
        }
        for (int j = i+2; j < 2*i+2; ++j)
        {
            d.dodaj_trikotnik(nova[j], zadnja[j-2], zadnja[j-1]);
        }
        for (int j = 2*i+2; j < 3*i+3; ++j)
        {
            d.dodaj_trikotnik(nova[j], nova[j+1], zadnja[j-2]);
        }
        for (int j = 2*i+3; j < 3*i+3; ++j)
        {
            d.dodaj_trikotnik(nova[j], zadnja[j-3], zadnja[j-2]);
        }
        zadnja = nova;
    }
    return d;
}

int main(int argc, char **argv) {
    Delitev d = srediscna(8);
    d.narisi("g_test.png");
    d.resi_poisson();
  //  cout << d.pretok();
    return 0;
}
