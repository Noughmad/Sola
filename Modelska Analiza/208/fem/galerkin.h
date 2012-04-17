
#ifndef GALERKIN_H
#define GALERKIN_H

#include "delitev.h"


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
    
private:
    int m;
    int n;
};

#endif // GALERKIN_H
