
#include "galerkin.h"
#include <QDebug>

#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_vector_complex_double.h>
#include <gsl/gsl_eigen.h>

#include <gsl/gsl_sort_vector.h>

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
    return mat;
}


Matrika Galerkin::matrika()
{
    Matrika mat;
    for (int k = 0; k < n; ++k)
    {
        for (int l = k; l < n; ++l)
        {
            mat[k][l] = 1.0*k*l/(2*m+k+l) - 1.0*(k+l)/(2*m+k+l+1) + 1.0*(k+1)*(l+1)/(2*m+k+l+2);
        }
    }
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

// #if not defined WITH_ARPACK
void Galerkin::resi_nihanje(int stevilo, bool risi, double lastne_vrednosti[])
{
   gsl_matrix* A = gsl(matrika());
   gsl_matrix* B = gsl(masa());
   
   gsl_vector* lambda = gsl_vector_alloc(n);
   gsl_matrix* evec = gsl_matrix_alloc(n, n);
   
    gsl_eigen_gensymmv_workspace* w = gsl_eigen_gensymmv_alloc(n);
    gsl_eigen_gensymmv(A, B, lambda, evec, w);
    gsl_eigen_gensymmv_sort(lambda, evec, GSL_EIGEN_SORT_ABS_ASC);

    for (int i = 0; i < stevilo; ++i)
    {
        lastne_vrednosti[i] = gsl_vector_get(lambda, i);
        qDebug() << gsl_vector_get(lambda, i);
    }
    
    gsl_eigen_gensymmv_free(w);
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(evec);
    gsl_vector_free(lambda);
}

gsl_matrix* Galerkin::gsl(const Matrika& mat)
{
    gsl_matrix* g = gsl_matrix_alloc(n,n);
    for (int i = 0; i < n; ++i)
    {
        gsl_matrix_set(g, i, i, mat[i][i]);
        for (int j = i; j < n; ++j)
        {
            gsl_matrix_set(g, i, j, mat[i][j]);
            gsl_matrix_set(g, j, i, mat[i][j]);
        }
    }
    return g;
}

// #endif