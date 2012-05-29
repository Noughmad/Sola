#include <iostream>
#include <fstream>
#include <stdlib.h>

#include <math.h>
#include <limits>
#include <assert.h>

#define CHECK_VALID(matrix)     \
if (!matrix.isValid()) { std::cout << "Neveljavna matrika " << #matrix << " v koraku " << i << std::endl; exit(5); }

class Matrix
{
public:
    Matrix(int n);
    ~Matrix();
    
    inline const double& get(int i, int j) const
    {
        return data[i*n+j];
    }
    inline void set(int i, int j, double v)
    {
        data[i*n+j] = v;
    }
    
    double interpolate(double x, double y) const;
    
    void save(const char* filename);
    
    Matrix operator*(int f);
    
    inline int size() const
    {
        return n;
    }
    
    bool isValid();
        
private:
    int n;
    double* data;
};

Matrix::Matrix(int n) : n(n), data(new double[n*n])
{

}

Matrix::~Matrix()
{
    delete[] data;
}

double Matrix::interpolate(double x, double y) const
{
    int lx = floor(x);
    int ux = ceil(x);
    int ly = floor(y);
    int uy = ceil(y);
    
    double fQ11 = get(lx, ly);
    double fQ21 = get(ux, ly);
    double fQ12 = get(lx, uy);
    double fQ22 = get(ux, uy);

    //if point exactly found on a node do not interpolate
    if ((lx == ux) && (ly == uy))  
        return fQ11;
    
    double x1 = lx;
    double x2 = ux;
    double y1 = ly;
    double y2 = uy;

    //if xcoord lies exactly on an xAxis node do linear interpolation
    if (lx == ux) 
            return fQ11 + (fQ12 - fQ11) * (y - y1) / (y2 - y1);
    //if ycoord lies exactly on an xAxis node do linear interpolation
    if (ly == uy) 
            return fQ11 + (fQ22 - fQ11) * (x - x1) / (x2 - x1);

    double fxy = fQ11 * (x2 - x) * (y2 - y);
    fxy = fxy + fQ21 * (x - x1) * (y2 - y);
    fxy = fxy + fQ12 * (x2 - x) * (y - y1);
    fxy = fxy + fQ22 * (x - x1) * (y - y1);
    fxy = fxy / ((x2 - x1) * (y2 - y1));

    return fxy;
}

bool Matrix::isValid()
{
    for (int i = 0; i < n*n; ++i)
    {
        if (isnan(data[i]) || isinf(data[i]))
        {
            std::cout << "Nan v elementu " << i << std::endl;
            return false;
        }
    }
    return true;
}


struct NsWorkspace
{
    NsWorkspace(int n, double R);
    ~NsWorkspace();
    
    void zacetni_pogoj();
    inline void korak(int i)
    {
        izracunaj_zeta();
        
        CHECK_VALID(Zeta)
        
        izracunaj_psi();
        CHECK_VALID(Psi)
        
        izracunaj_v();
        assert(vx.isValid());
        assert(vy.isValid());
        
        /*
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                assert(vx.get(i,j) >= -1.1);
                assert(vy.get(i,j) >= -1.1);
                assert(vx.get(i,j) <= 1.1);
                assert(vy.get(i,j) <= 1.1);
            }
        }
        */
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
    
    Matrix Psi, Zeta, vx, vy, tmp;
};

NsWorkspace::NsWorkspace(int n, double R) : 
N(n),
Psi(n),
Zeta(n),
vx(n),
vy(n),
tmp(n),
R(R)
{
    h = 1.0 / (double)(N-1);
    k = h / 1000.0;
    omega = 2.0/(1+M_PI/N);
    std::cout << "Omega = " << omega << ", h = " << h << ", k = " << k << std::endl;
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
    double uzx, vzy, nbz;
    for (int i = 1; i < N-1; ++i)
    {
        for (int j = 1; j < N-1; ++j)
        {
            vzy = (Zeta.get(i+1,j) * vy.get(i+1,j) - Zeta.get(i-1,j) * vy.get(i-1,j)) / 2;
            uzx = (Zeta.get(i,j+1) * vx.get(i,j+1) - Zeta.get(i,j-1) * vx.get(i,j-1)) / 2;
            nbz = ( -4 * Zeta.get(i,j) + Zeta.get(i-1,j) + Zeta.get(i+1,j) + Zeta.get(i,j-1) + Zeta.get(i,j+1) ) / h;
        
            tmp.set(i, j, Zeta.get(i,j) - k / h * (uzx + vzy - nbz / R ));
        }
    }
    
    for (int i = 1; i < N-1; ++i)
    {
        for (int j = 1; j < N-1; ++j)
        {
            Zeta.set(i, j, tmp.get(i, j));
        }
    }
    
    // Robni pogoji: Enacba (10) v navodilih
    for (int i = 0; i < N; ++i)
    {
        Zeta.set(i, 0, 2.0/h/h * (Psi.get(i, 1) - h * vx.get(i, 0)));
        Zeta.set(i, N-1, 2.0/h/h * Psi.get(i, N-2) - h * vx.get(i, N-1));
    }
    for (int j = 0; j < N; ++j)
    {
        Zeta.set(0, j, 2.0/h/h * Psi.get(1, j));
        Zeta.set(N-1, j, 2.0/h/h * Psi.get(N-2, j));
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
    
    double t = omega;
    omega = 1;
    
    for (int i = 0; i < 5; ++i)
    {
        korak_psi();
    }
    omega = t;
    
    double lastEps = std::numeric_limits<double>::max();
    int i = 0;
    do 
    {
        ++i;
        eps = korak_psi();
        
      //  std::cout << "Popravek " << i << " : " << eps << std::endl;
    }
    while (eps > 1e-6 && i < 1000);
}

void Matrix::save(const char* filename)
{
    std::ofstream stream(filename);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            stream << get(i,j) << " ";
        }
        stream << std::endl;
    }
    stream.close();
}

void interpoliraj(const Matrix& mala, Matrix& velika)
{
    int n = mala.size();
    int m = velika.size();
    
    double f = (double)(n-1)/(m-1); // f < 1
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < m; ++j)
        {
            velika.set(i, j, mala.interpolate(f*i, f*j));
        }
    }
}

NsWorkspace* razdeli_mrezo(const NsWorkspace* stari, int n)
{
    NsWorkspace* novi = new NsWorkspace(n, stari->R);
    
    interpoliraj(stari->Psi, novi->Psi);
    interpoliraj(stari->Zeta, novi->Zeta);
    interpoliraj(stari->vx, novi->vx);
    interpoliraj(stari->vy, novi->vy);
    
    delete stari;
    return novi;
}

double postopek(int N, double R)
{
    NsWorkspace* workspace = new NsWorkspace(10, R);
    
    const int IterationsPerSave = 1e5;
    const int Saves = 10;
    
    const int Iterations = 1e4;
    
    char buf[64];

    for (int n = workspace->N; n < N; n = 2*n-2)
    {
        for (int i = 0; i < Iterations; ++i)
        {   
            workspace->korak(i);
        }
        
        sprintf(buf, "g_po_%d.dat", n);
        workspace->shrani(buf);
        
        if (n < N/2)
        {
            std::cout << "Povecujem mrezo na " << 2*n-2 << std::endl;
            workspace = razdeli_mrezo(workspace, 2*n-2);
            sprintf(buf, "g_prej_%d.dat", 2*n-2);
            workspace->shrani(buf);
        }
        else
        {
            break;
        }
    }
    
    workspace = razdeli_mrezo(workspace, N);

    for (int i = 0; i < IterationsPerSave * Saves; ++i)
    {
        if (i % IterationsPerSave == 0)
        {
            sprintf(buf, "g_tok_50_%d.dat", i / IterationsPerSave);
            workspace->shrani(buf);
            std::cout << "Shranil " << buf << std::endl;
            sprintf(buf, "g_zeta_50_%d.dat", i / IterationsPerSave);
            workspace->Zeta.save(buf);
            sprintf(buf, "g_vx_50_%d.dat", i / IterationsPerSave);
            workspace->vx.save(buf);
            sprintf(buf, "g_vy_50_%d.dat", i / IterationsPerSave);
            workspace->vy.save(buf);
        }
        
        workspace->korak(i);
    }
    
    delete workspace;
}

int main(int argc, char **argv) {
    postopek(100, 2000);
    return 0;
}
