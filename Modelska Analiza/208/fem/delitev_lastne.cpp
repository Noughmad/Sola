#include "delitev.h"

#include <QtCore/qmath.h>

#include <QtCore/QDebug>

void Delitev::resi_nihanje(int stevilo, bool risi, double lastne_vrednosti[])
{
    cholmod_common Common;
    cc = &Common;
    cholmod_start(cc);
    cc->dtype = CHOLMOD_DOUBLE;
    
    cholmod_sparse* A = sparse(matrika(), false);
    cholmod_sparse* B = sparse(masa(), true);
    
    cholmod_factor* L = cholmod_analyze(B, cc);
    cholmod_factorize(B, L, cc);
    Q_ASSERT(L);
    
    cholmod_sparse* Yt = cholmod_spsolve(CHOLMOD_L, L, A, cc);
    Q_ASSERT(Yt);
    
    cholmod_sparse* Y = cholmod_transpose(Yt, 1, cc);
    cholmod_sparse* C = cholmod_spsolve(CHOLMOD_L, L, Y, cc);
    Q_ASSERT(C);
        
    // Matrika C mora biti simetricna, ceprav cholmod_spsolve tega ne uposteva
    
    int n = st_notranjih();
    int nzmax = (C->nzmax + n + 1) / 2;
    
    int t = 0;
    int c = 0;
    
    int* Cp = (int*)C->p;
    int* Ci = (int*)C->i;
    double* Cx = (double*)C->x;
    
    cholmod_sparse* R = cholmod_allocate_sparse(n, n, nzmax, 1, 1, 1, CHOLMOD_REAL, cc);
    
    int* Rp = (int*)R->p;
    int* Ri = (int*)R->i;
    double* Rx = (double*)R->x;
    
    Rp[0] = 0;
    Ri[0] = 0;
    Rx[0] = Cx[0];
    ++t;
    
    for (int col = 1; col < n; ++col)
    {
        Rp[col] = t;
        for (int j = Cp[col]; j < Cp[col+1]; ++j)
        {
            int row = Ci[j];
            if (row <= col)
            {
                Ri[t] = row;
                Rx[t] = Cx[j];
                ++t;
            }
        }
    }
    Rp[n] = t;
    
    int nconv = stevilo;
    double* EigVal = new double[nconv];
    double* EigVec = new double[nconv * n];
    
    nconv = lastne_arpack(EigVal, EigVec, n, 
                            t, Rx, Ri, Rp, 
                            'U', nconv, "SM");
    
    qDebug() << "Izracunal" << nconv << "lastnih nihanj";
    for (int i = 0; i < nconv; ++i)
    {
        qDebug() << "lambda [" << i << "] =" << EigVal[i];
        lastne_vrednosti[i] = EigVal[i];
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
            
            cholmod_dense* a = cholmod_solve(CHOLMOD_Lt, L, b, cc);
            double* ax = (double*)a->x;
            
            plot(a, vec+1, sqrt(EigVal[vec]));
        }
    }
    
    cholmod_finish(cc);
}


int Delitev::arpack(const Matrika& elementi, int** i, int** p, double** x)
{
    int nnz = 0;
    int n = notranje.size();
    int t = 0;
    
    Matrika::const_iterator gend = elementi.constEnd();
    Matrika::const_iterator git = elementi.constBegin();
    
    for (git = elementi.constBegin(); git != gend; ++git)
    {
        nnz += git.value().size();
    }
    
    qDebug() << n << nnz;
    
    int* Ai = new int[nnz];
    int* Ap = new int[n+1];
    double* Ax = new double[nnz];
    
    int col = 0;
    for (git = elementi.constBegin(); git != gend; ++git)
    {
        Ap[col] = t;
        Vrstica::const_iterator end = git.value().constEnd();
        Vrstica::const_iterator it = git.value().constBegin();
        for (; it != end; ++it)
        {          
            Ai[t] = not_indeksi[it.key()];
            Ax[t] = it.value();
            ++t;
        }
        ++col;
    }
    
    Ap[n] = t;
    
    for (int j = 0; j < n+1; ++j)
    {
        qDebug() << "Ap [" << j << "] = " << Ap[j];
    }
    
    for (int j = 0; j < nnz; ++j)
    {
        qDebug() << "Ai [" << j << "] = " << Ai[j] << ", Ax [" << j << "] = " << Ax[j];
    }
    
    
    
    *i = Ai;
    *p = Ap;
    *x = Ax;
    return nnz;
}