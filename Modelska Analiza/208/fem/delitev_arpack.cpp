#include "delitev.h"

#if defined WITH_ARPACK
#include "areig.h"
#endif


int lastne_arpack(double EigVal[], double EigVec[], int n, int nnz, double A[],
          int irow[], int pcol[], char uplo, int nev,
          char* which)
{
#if defined WITH_ARPACK
    return AREig<double>(EigVal, EigVec, n, nnz, A, irow, pcol, uplo, nev, which);
#else
    return 0;
#endif
}


int lastne(double EigVal[], double EigVec[], int n, int nnz, cholmod_sparse* A, int nev, char* which)
{
    return lastne_arpack(EigVal, EigVec, n, nnz, A->x, A->i, A->p, (A->stype == 1) ? 'U' : 'L', nev, which);
}
int lastne_gen(double EigVal[], double EigVec[], int n, int nnzA, cholmod_sparse* A, int nnzB, cholmod_sparse* B, int nev, char* which)
{
#if defined WITH_ARPACK
    return AREig<double>(EigVal, EigVec, n, nnzA, (double*)A->x, (int*)A->i, (int*)A->p, 
                                         nnzB, (double*)B->x, (int*)B->i, (int*)B->p, (A->stype == 1) ? 'U' : 'L', nev, which);
#else
    return 0;
#endif
}

