#ifndef DELITEV_H
#define DELITEV_H

#include <QtCore/QList>
#include <QtCore/QHash>
#include <QtCore/QMap>
#include <QtCore/QRectF>

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

typedef QMap<int, double> Vrstica;
typedef QMap<int, Vrstica> Matrika;

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
    
    virtual Matrika matrika();
    cholmod_dense* desne_strani();
    virtual Matrika masa();
    
    cholmod_sparse* sparse(const Matrika& elementi, bool symmetric);
    int arpack(const Matrika& m, int** i, int** p, double** x);
    
    double x(int i, int j) const;
    double y(int i, int j) const;
    
    void resi_poisson(bool risi);
    virtual void resi_nihanje(int stevilo, bool risi, double lastne_vrednosti[]);
    
    void narisi(const QString& file);
    void shrani(const QString& file);
    
    virtual int st_notranjih() const;
    int st_tock() const;
    int st_trikotnikov() const;
    
    virtual void plot(cholmod_dense* a, int d, double k);
    
    QRectF rect;
    const char* name;
    
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

typedef Delitev* (*Generator)(int);

int lastne_arpack(double EigVal[], double EigVec[], int n, int nnz, double A[],
          int irow[], int pcol[], char uplo, int nev,
          char* which = "LM");

#endif // DELITEV_H

struct Tocka;
