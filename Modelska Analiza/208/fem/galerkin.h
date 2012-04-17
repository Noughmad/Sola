
#ifndef GALERKIN_H
#define GALERKIN_H

#include "delitev.h"

#include <gsl/gsl_matrix.h>

class Galerkin : public Delitev
{

public:
    /**
     * @param m Indeks kotnega dela resitve
     * @param n Stevilo k-jev
     **/
    Galerkin(int m, int n);
    virtual ~Galerkin();
    
    virtual Matrika matrika();
    virtual Matrika masa();
    
    virtual void plot(cholmod_dense* a, int d, double k);
    virtual int st_notranjih() const;
    
#ifndef WITH_ARPACK
    virtual void resi_nihanje(int stevilo, bool risi, double lastne_vrednosti[]);
    gsl_matrix* gsl(const Matrika& mat);
#endif
    
private:
    int m;
    int n;
};

#endif // GALERKIN_H
