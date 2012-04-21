
#include "galerkin.h"

#include <QtCore/QDebug>
#include <QtCore/qmath.h>
#include <QFile>

double g(double m, double k, double r, double phi)
{
    return pow(r, m+k) * (1-r) * sin(m*phi);
}

Galerkin::Galerkin(int m, int n) : Delitev(), 
m(m),
n(n)
{
    for (int i = 0; i < n; ++i)
    {
        not_indeksi[i] = i;
    }
}

Galerkin::~Galerkin()
{

}

Matrika Galerkin::masa()
{
    Matrika mat;
    for (int i = 0; i < n; ++i)
    {
        for (int j = i; j < n; ++j)
        {
            mat[i][j] = 1.0/(2*m+i+j+2) - 2.0/(2*m+i+j+3) + 1.0/(2*m+i+j+4);
        }
    }
    return mat;
}


Matrika Galerkin::matrika()
{
    Matrika mat;
    for (int k = 0; k < n; ++k)
    {
        for (int l = k; l < n; ++l)
        {
            mat[k][l] = 1.0*k*l/(2*m+k+l) - 1.0*(2*k*l+k+l)/(2*m+k+l+1) + 1.0*(k+1)*(l+1)/(2*m+k+l+2);
        }
    }
    return mat;
}

void Galerkin::plot(cholmod_dense* a, int d, double k)
{
    qDebug() << "Resil galerkin: " << d << k;
    
    QFile file(QString("g_galerkin_%1_%2.dat").arg(m).arg(d));
    file.open(QIODevice::WriteOnly);
    
    QTextStream stream(&file);
    double* ax = (double*)a->x;
    
    QString vector = "(";
    for(int i = 0; i < n; ++i)
    {
        vector += QString::number(ax[i], 'g', 3) += ',';
    }
    vector.chop(1);
    vector += ")";
    qDebug() << m << d << ": => " << vector;
    
    for (double r = 0; r <= 1; r += 0.01)
    {
        for (double phi = 0; phi < M_PI; phi += M_PI*0.01)
        {
            double t = 0;
            for (int i = 0; i < n; ++i)
            {
                t += ax[i] * g(m, d, r, phi);
            }
            stream << t << " ";
        }
        stream << endl;
    }
    file.close();
}

int Galerkin::st_notranjih() const
{
    return n;
}
