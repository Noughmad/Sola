#ifndef DELITEV_H
#define DELITEV_H

#include <QtCore/QList>
#include <QtCore/QHash>

#include <cholmod.h>

class QPointF;
struct Tocka
{
    double x, y, z;
    bool noter;
    
    operator QPointF() const;
};

struct Trikotnik
{
    int i, j, k;
    double ploscina;
};

quint16 qHash(const Trikotnik& t);
bool operator==(const Trikotnik& t1, const Trikotnik& t2);

class Delitev
{
public:
    Delitev();
    virtual ~Delitev();
    
    inline double ploscina(const Trikotnik& t) const
    {
        return t.ploscina;
    }
    
    double dd(int i, int j) const;
    double pretok();
    
    int dodaj_tocko(double x, double y, bool noter);
    int dodaj_trikotnik(int i, int j, int k);
    
    cholmod_sparse* matrika();
    cholmod_dense* desne_strani();
    cholmod_sparse* masa();
    
    double x(int i, int j) const;
    double y(int i, int j) const;
    
    void resi_poisson();
    void resi_nihanje();
    
    void narisi(const QString& file);
    
    int st_notranjih() const;
    int st_tock() const;
    
private:
    int indeks(int i, int j) const;
        
private:
    QList<Tocka> tocke;
    QList<int> notranje;
    QHash<int,int> not_indeksi;
    QList<Trikotnik> trikotniki;
    QHash<Trikotnik, double> ploscine;
    cholmod_common* cc;
};

typedef Delitev (*Generator)(int);

#endif // DELITEV_H

struct Tocka;
