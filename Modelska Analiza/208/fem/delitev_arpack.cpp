#include "delitev.h"

#include "areig.h"

int lastne_arpack(double EigVal[], double EigVec[], int n, int nnz, double A[],
          int irow[], int pcol[], char uplo, int nev,
          char* which = "LM")
{
    return AREig<double>(EigVal, EigVec, n, nnz, A, irow, pcol, uplo, nev, which);
}