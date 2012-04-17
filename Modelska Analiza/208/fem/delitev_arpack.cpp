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
