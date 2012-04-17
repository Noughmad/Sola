
#include "galerkin.h"
#include <QDebug>

Galerkin::Galerkin(int m, int n) : Delitev(), 
m(m),
n(n)
{

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
    qDebug() << "Masa" << n << mat.size();
    return mat;
}


Matrika Galerkin::matrika()
{

    Matrika mat;
    for (int k = 0; k < n; ++k)
    {
        for (int l = k; l < n; ++l)
        {
            mat[k][l] = 1.0*k*l/(2*m+k+l) - 1.0*(k+l)/(2*m+k+l+1) + 1.0*(k+l+k*l)/(2*m+k+l+2);
        }
    }
    qDebug() << "Matrika" << n << mat.size();
    qDebug() << mat;
    return mat;
}

void Galerkin::plot(cholmod_dense* a, int d, double k)
{
    qDebug() << "Resil galerkin: " << d << k;
}

int Galerkin::st_notranjih() const
{
    return n;
}

