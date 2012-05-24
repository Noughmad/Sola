#include <iostream>
#include <fstream>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

#include <math.h>

class Matrix
{
public:
    Matrix(size_t n1, size_t n2);
    ~Matrix();
    
    inline double get(int i, int j);
    inline void set(int i, int j, double v);
    
    void save(const char* filename);
        
private:
    gsl_matrix* m;
};

Matrix::Matrix(size_t n1, size_t n2) : m(gsl_matrix_alloc(n1, n2))
{

}

Matrix::~Matrix()
{
    gsl_matrix_free(m);
}

double Matrix::get(int i, int j)
{
    return gsl_matrix_get(m, i, j);
}

void Matrix::set(int i, int j, double v)
{
    gsl_matrix_set(m, i, j, v);
}

struct NsWorkspace
{
    NsWorkspace(int n);
    ~NsWorkspace();
    
    void zacetni_pogoj();
    inline void korak()
    {
        izracunaj_zeta();
        izracunaj_psi();
        izracunaj_v();
    }
    
    void izracunaj_zeta();
    void izracunaj_psi();
    void izracunaj_v();
    double korak_psi();
    
    inline void shrani(const char* filename)
    {
        Psi.save(filename);
    }
    
    inline void popravi_psi(int i, int j, double eps, double* vsota)
    {
        Psi.set(i, j, Psi.get(i,j) + eps * 0.25 * omega);
        *vsota += eps * eps;
    }
    
    int N;
    double k;
    double R;
    double h;
    double omega;
    
    Matrix Psi, Zeta, vx, vy;
};

NsWorkspace::NsWorkspace(int n) : 
N(n),
Psi(n,n),
Zeta(n,n),
vx(n,n),
vy(n,n)
{
    h = 1.0 / (N-1);
    k = h / 3.0;
    omega = 2.0/(1+M_PI/N);
    zacetni_pogoj();
}

NsWorkspace::~NsWorkspace()
{
    
}

void NsWorkspace::zacetni_pogoj()
{
    /*
     * Najenostavnejsi mozni zacetni pogoj: Vse je nic
     * Je tudi fizikalno smiselen, saj lahko recemo, da ob casu T=0
     * zacnemo premikati zgornjo plosco
     */
    
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            Zeta.set(i,j,0);
            Psi.set(i, j, 0);
            vx.set(i, j, 0);
            vy.set(i, j, 0);
        }
    }
    
    // Zgonja plosca
    
    for (int i = 0; i < N; ++i)
    {
        vx.set(i, 0, 1);
    }
    
    /*
     * Na zacetku predpostavimo nekaj enostavnega, a se vedno dovolj smiselnega
     * 
     * * Vrtinec na sredini opisemo kot vrh hitrostnega potenciala Psi
     * * Hitrost izracunamo iz psi
     * * Vrtincnost lahko izrazimo neposredno iz hitrosti
     */
    
    /*
    
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            Psi.set(i, j, sin(M_PI * i / N) * sin(M_PI * j / N));
        }
    }
    
    izracunaj_v();
    
    for (int i = 1; i < N-1; ++i)
    {
        for (int j = 1; j < N-1; ++j )
        {
            double udy = vx.get(i, j+1) - vx.get(i, j-1);
            double vdx = vy.get(i+1, j) - vy.get(i-1, j);
            Zeta.set(i, j, 0.5 / h * (udy - vdx));
        }
    }
    
    */
}


void NsWorkspace::izracunaj_zeta()
{
    // zeta = - (rot v) = dvx/dy - dvy/dx
    
    // Robni pogoji: Enacba (10) v navodilih
    // NOTE: Mogoce bo treba robne pogoje upostevati pred racunom na sredini
    for (int i = 0; i < N; ++i)
    {
        Zeta.set(i, 0, 2.0/h/h * (Psi.get(i, 1) - h));
        Zeta.set(i, N-1, 2.0/h/h * Psi.get(i, N-2));
    }
    for (int j = 0; j < N; ++j)
    {
        Zeta.set(0, j, 2.0/h/h * Psi.get(1, j));
        Zeta.set(N-1, j, 2.0/h/h * Psi.get(N-2, j));
    }
    
    double uzx, vzy, nbz;
    for (int i = 1; i < N-1; ++i)
    {
        for (int j = 1; j < N-1; ++j)
        {
            uzx = (Zeta.get(i+1,j) * vx.get(i+1,j) - Zeta.get(i-1,j) * vx.get(i-1,j)) / 2;
            vzy = (Zeta.get(i,j+1) * vy.get(i,j+1) - Zeta.get(i,j-1) * vy.get(i,j-1)) / 2;
            nbz = ( -4 * Zeta.get(i,j) + Zeta.get(i-1,j) + Zeta.get(i+1,j) + Zeta.get(i,j-1) + Zeta.get(i,j+1) ) / h;
        
            Zeta.set(i, j, Zeta.get(i,j) - k / h * (uzx + vzy - nbz / R ));
            if (isnan(Zeta.get(i,j)))
            {
                std::cout << " NaN: " << i << ", " << j << std::endl;
                exit(-4);
            }
        }
    }
    
}

void NsWorkspace::izracunaj_v()
{
    for (int i = 1; i < N-1; ++i)
    {
        for (int j = 1; j < N-1; ++j)
        {
            vx.set(i, j, (Psi.get(i+1,j) - Psi.get(i-1,j)) / 2 / h);
            vy.set(i, j, (Psi.get(i,j+1) - Psi.get(i,j-1)) / 2 / h);
        }
    }
    
    /*
     * Za hitrost ni treba pazit na robne pogoje
     * Hitrost na robovih je vedno enaka nic (ali pa hitrosti stene)
     * vsekakor se pa s casom ne spreminja
     */
}

double NsWorkspace::korak_psi()
{   
    double popravek = 0;
    
    // Lihe tocke v notranjosti
    for (int i = 1; i < N-1; ++i)
    {
        for (int j = 1+(i%2); j < N-1; j += 2)
        {
            double eps = Psi.get(i+1, j) + Psi.get(i-1, j) + Psi.get(i, j+1) + Psi.get(i, j-1) - 4*Psi.get(i, j) - h*h*Zeta.get(i,j);
            popravi_psi(i, j, eps, &popravek);
        }
    }
    
    // Sode tocke v notranjosti
    for (int i = 1; i < N-1; ++i)
    {
        for (int j = 2-(i%2); j < N-1; j += 2)
        {
            double eps = Psi.get(i+1, j) + Psi.get(i-1, j) + Psi.get(i, j+1) + Psi.get(i, j-1) - 4*Psi.get(i, j) - h*h*Zeta.get(i,j);
            popravi_psi(i, j, eps, &popravek);
        }
    }
    
    /*
     * Robni pogoj zahteva, da ni hitrosti v steno
     * To je ekvivalentno robnemu pogoju prve vrste, torej Psi=0
     * 
     * Podobno kot pri racunu Zeta se robne tocke s casom ne spreminjajo. 
     */
    return popravek;
}

void NsWorkspace::izracunaj_psi()
{
    double eps;
    do 
    {
        eps = korak_psi();
    }
    while (eps > 1e-7);
}

void Matrix::save(const char* filename)
{
    std::ofstream stream(filename);
    for (int i = 0; i < m->size1; ++i)
    {
        for (int j = 0; j < m->size2; ++j)
        {
            stream << get(i,j) << " ";
        }
        stream << std::endl;
    }
    stream.close();
}


double postopek(int N)
{
    NsWorkspace workspace(N);
    workspace.R = 50;
    
    const int IterationsPerSave = 100;
    const int Saves = 10;
    
    char buf[64];
    for (int i = 0; i < IterationsPerSave * Saves; ++i)
    {
        if (i % IterationsPerSave == 0)
        {
            sprintf(buf, "g_tok_50_%d.dat", i / IterationsPerSave);
            workspace.shrani(buf);
            std::cout << "Shranil " << buf << std::endl;
            sprintf(buf, "g_zeta_50_%d.dat", i / IterationsPerSave);
            workspace.Zeta.save(buf);
        }
        
        workspace.korak();
    }
}

int main(int argc, char **argv) {
    postopek(40);
    return 0;
}
