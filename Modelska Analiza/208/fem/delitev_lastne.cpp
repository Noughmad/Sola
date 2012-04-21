#include "delitev.h"

#include <QtCore/qmath.h>

#include <QtCore/QDebug>

#ifndef WITH_ARPACK
#  include <gsl/gsl_matrix_double.h>
#  include <gsl/gsl_vector_complex_double.h>
#  include <gsl/gsl_eigen.h>
#endif


void Delitev::resi_nihanje(int stevilo, bool risi, double lastne_vrednosti[])
{
#ifdef WITH_ARPACK
    cholmod_common Common;
    cc = &Common;
    cholmod_start(cc);
    cc->dtype = CHOLMOD_DOUBLE;
    
    cholmod_sparse* A = sparse(matrika(), true);
    cholmod_sparse* B = sparse(masa(), true);
    
    int n = st_notranjih();
    double* EigVal = new double[stevilo];
    double* EigVec = new double[stevilo * n];
    
    int nconv = lastne_gen(EigVal, EigVec, n, A->nzmax, A, B->nzmax, B, stevilo, "SM");
    
    qDebug() << "Izracunal" << nconv << "lastnih nihanj";
    for (int i = 0; i < nconv; ++i)
    {
        qDebug() << "k [" << i << "] =" << sqrt(EigVal[i]);
        lastne_vrednosti[i] = sqrt(EigVal[i]);
    }
    
    if (risi)
    {
        for (int vec = 0; vec < nconv; ++vec)
        {
            cholmod_dense* b = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, cc);
            double* bx = (double*)b->x;
            for (int i = 0; i < n; ++i)
            {
                bx[i] = EigVec[vec*n+i];
            }
            plot(b, vec+1, sqrt(EigVal[vec]));
        }
    }
    
    cholmod_finish(cc);
    
#else // WITH_ARPACK
    
    int n = st_notranjih();
    gsl_matrix* A = gsl(matrika());
    gsl_matrix* B = gsl(masa());
    
    // print_matrix(A);
    // print_matrix(B);
    
    gsl_vector* lambda = gsl_vector_alloc(n);
    gsl_matrix* evec = gsl_matrix_alloc(n, n);
   
    gsl_eigen_gensymmv_workspace* w = gsl_eigen_gensymmv_alloc(n);
    gsl_eigen_gensymmv(A, B, lambda, evec, w);
    gsl_eigen_gensymmv_sort(lambda, evec, GSL_EIGEN_SORT_ABS_ASC);

    for (int i = 0; i < stevilo; ++i)
    {
        lastne_vrednosti[i] = gsl_vector_get(lambda, i);
        qDebug() << sqrt(gsl_vector_get(lambda, i));
    }
    
    if (risi)
    {
        for (int vec = 0; vec < stevilo; ++vec)
        {
            cholmod_common Common;
            cc = &Common;
            cholmod_start(cc);
            
            cholmod_dense* b = cholmod_allocate_dense(n, 1, n, CHOLMOD_REAL, cc);
            double* bx = (double*)b->x;
            for (int i = 0; i < n; ++i)
            {
                bx[i] = gsl_matrix_get(evec, vec, i);
            }
            
            plot(b, vec+1, sqrt(gsl_vector_get(lambda, vec)));
            cholmod_finish(cc);
        }
    }
    
    gsl_eigen_gensymmv_free(w);
    gsl_matrix_free(A);
    gsl_matrix_free(B);
    gsl_matrix_free(evec);
    gsl_vector_free(lambda);
    
#endif // WITH_ARPACK
}