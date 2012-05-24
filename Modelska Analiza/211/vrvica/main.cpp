#include <iostream>
#include <fstream>
#include <stdio.h>

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_roots.h>

#include <QtGui/QImage>
#include <QtGui/QPainter>
#include <QtGui/QBitmap>
#include <QtGui/QApplication>
#include <QtCore/QDebug>

inline double sq(double x)
{
    return x*x;
}

const int VelikostSlike = 200;
QString ImageFileNameFormat;

int StClenov;
int StClenovMinusDva;
double k = 1e-5;
double h;

double C; // C = k*k/h/h
double D; // D = h*h/k/k

const int IterationsPerFrame = 1e5 / 25;

double *Phi, *Phi1, *pn;
gsl_vector* F;

QImage* image;
QPainter* painter;

gsl_vector *d, *u, *l, *rhs;

gsl_root_fdfsolver* solver;
gsl_function_fdf function;

double fn_f(double x, void* params)
{
    double* p = (double*)params;
    return x - p[1] + h*gsl_vector_get(F,0) * cos(x);
}

double fn_df(double x, void* params)
{
    return 1 - h*gsl_vector_get(F,0) * sin(x);
}

void fn_fdf (double x, void *params, double *y, double *dy)
{
    double* p = (double*)params;
    const double b = h * gsl_vector_get(F,0);
    *y = x - p[1] + b*cos(x);
    *dy = 1 - b*sin(x);
}

void polozaj_repa(double& x, double& y)
{
    x = 0;
    y = 0;
    for (int i = 0; i < StClenov; ++i)
    {
        x += cos(Phi[i]);
        y += sin(Phi[i]);
    }
    
    x *= h;
    y *= h;
}

struct energy_t
{
    double potencialna;
    double kineticna;
    double rotacijska;
};

std::ostream& operator<<(std::ostream& stream, const energy_t& energy)
{
    stream << energy.potencialna << " " << energy.kineticna << " " << energy.rotacijska;
    return stream;
}

energy_t energija()
{   
    double x = cos(Phi[0]);
    double y = sin(Phi[0]);
    double xk = cos(Phi1[0]);
    double yk = sin(Phi1[0]);
    
    energy_t E;
    E.potencialna = -y;
    E.kineticna = (x-xk) * (x-xk) + (y-yk) * (y-yk);
    E.rotacijska = (Phi[0] - Phi1[0]) * (Phi[0] - Phi1[0]);
        
    for (int i = 1; i < StClenov; ++i)
    {
        x += cos(Phi[i-1]) + cos(Phi[i]);
        y += sin(Phi[i-1]) + sin(Phi[i]);
        xk += cos(Phi1[i-1]) + cos(Phi1[i]);
        yk += sin(Phi1[i-1]) + sin(Phi1[i]);

        E.potencialna -= y;
        E.kineticna += (x-xk) * (x-xk) + (y-yk) * (y-yk);
        E.rotacijska += (Phi[i] - Phi1[i]) * (Phi[i] - Phi1[i]);
    }
    
    E.potencialna *= 0.5 * h * h;
    E.kineticna *= 0.125 * D * h;
    E.rotacijska *= 1.0/12.0 * h * D;

    return E;
}

void izracunaj_silo()
{    
    gsl_vector_set(d, 0, -1);
    gsl_vector_set(u, 0, 1);
    gsl_vector_set(rhs, 0, -h*sin(Phi[0]));
    
    for (int i = 1; i < StClenovMinusDva; ++i)
    {
        gsl_vector_set(d, i, -2.0 - 0.25*sq(Phi[i+1]-Phi[i-1]));
        gsl_vector_set(l, i-1, 1);
        gsl_vector_set(u, i, 1);
        gsl_vector_set(rhs, i, -D*sq(Phi[i] - Phi1[i]));
    }
    
    const int i = StClenovMinusDva;
    gsl_vector_set(d, i, -2.0 - 0.25*sq(Phi[i+1]-Phi[i-1]));
    gsl_vector_set(l, i-1, 1);
    gsl_vector_set(rhs, i, -D*sq(Phi[i] - Phi1[i]));
    
    gsl_linalg_solve_tridiag(d, u, l, rhs, F);
}

void izracunaj_kot()
{
    double* pm = Phi1;
    double* p = Phi;
    
    for (int i = 1; i < StClenovMinusDva; ++i)
    {
        const double O1 = (gsl_vector_get(F, i+1)-gsl_vector_get(F, i-1))*(p[i+1]-p[i-1]) * 0.5;
        const double O2 = (p[i-1] - 2*p[i] + p[i+1])*gsl_vector_get(F, i);
        pn[i] = 2*p[i] - pm[i] + C * (O1 + O2);
    }

    const int n = StClenovMinusDva;
    const double O1 = (0 -gsl_vector_get(F, n-1))*(p[n+1]-p[n-1]) * 0.5;
    const double O2 = (p[n-1] - 2*p[n] + p[n+1])*gsl_vector_get(F, n);
    pn[n] = 2*p[n] - pm[n] + C * (O1 + O2);
    
    // Robni pogoj na koncu
    pn[StClenov-1] = 2*pn[StClenov-2] - pn[StClenov-3];
    
    
    // Robni pogoj na zacetku

    // Privzamem da se vrednost ne bo veliko spreminjala
    /*
    function.params = pn;
    gsl_root_fdfsolver_set (solver, &function, Phi[0]);
    
    int status;
    double x, x0;
    do
    {
           gsl_root_fdfsolver_iterate (solver);
           x0 = x;
           x = gsl_root_fdfsolver_root (solver);
           status = gsl_root_test_delta (x, x0, 0, 1e-3);
    }
    while (status == GSL_CONTINUE);  
    
    pn[0] = gsl_root_fdfsolver_root(solver);
    
    */
    
    pn[0] = pn[2] + 2*h/gsl_vector_get(F,1) * cos(pn[1]);
        
    Phi1 = p;
    Phi = pn;
    pn = pm;
}

void zacetni_pogoj(int n, double phi_0)
{
    StClenov = n;
    StClenovMinusDva = n-2;
    h = 1.0 / n;
    C = k*k/h/h;
    D = h*h/k/k;
    
    Phi = new double[StClenov];
    Phi1 = new double[StClenov];
    
    double x = 0;
    
    double c = 0.7;
    double a = 1;
    double bb = -c * cosh((x-a)/c);
    for (int i = 0; i < StClenov; ++i)
    {
        Phi[i] = -atan(sinh((x-a)/c));
        Phi1[i] = Phi[i];
        
        x += h*cos(Phi[i]);
    }
    
    pn = new double[StClenov];
    F = gsl_vector_alloc(StClenov-1);
    
    d = gsl_vector_alloc(StClenov-1);
    u = gsl_vector_alloc(StClenov-2);
    l = gsl_vector_alloc(StClenov-2);
    rhs = gsl_vector_alloc(StClenov-1);
    
    image = new QImage(2*VelikostSlike, VelikostSlike, QImage::Format_RGB16);
    painter = new QPainter(image);
    painter->setRenderHint(QPainter::Antialiasing, true);
    
    solver = gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);
    function.fdf = fn_fdf;
    function.df = fn_df;
    function.f = fn_f;
    function.params = 0;
}

void shrani_sliko(int frame)
{
    image->fill(QColor(Qt::black));
    
    QPointF lastPos = QPointF(VelikostSlike, 0);
    QPointF currentPos;
    
    const double h = 1.0/StClenov * VelikostSlike;
    
    QPen pen;
    pen.setWidth(5);
    for (int i = 0; i < StClenov-1; ++i)
    {
        if (gsl_vector_get(F, i) < 0)
        {
            exit(27);
        }
        pen.setColor(QColor(255, 150 * gsl_vector_get(F, i), 0));
        painter->setPen(pen);
        
        currentPos = lastPos + h*QPointF(cos(Phi[i]), sin(Phi[i]));
        painter->drawLine(lastPos, currentPos);
        lastPos = currentPos;
    }
    
    pen.setColor(QColor(255, 0, 0));
    painter->setPen(pen);
        
    currentPos = lastPos + h*QPointF(cos(Phi[StClenov-1]), sin(Phi[StClenov-1]));
    painter->drawLine(lastPos, currentPos);
    
    image->save(ImageFileNameFormat.arg(frame, 4, 10, QLatin1Char('0')));
    
    std::cout << "Shranil sliko " << frame << std::endl;
}

void korak()
{
    izracunaj_silo();
    izracunaj_kot();
}

void postopek(int n, double phi, int slike)
{
    zacetni_pogoj(n, phi);
    double t = 0;
    int iteration = 0;
    while (iteration < slike * IterationsPerFrame)
    {
        izracunaj_silo();

        if (iteration % IterationsPerFrame == 0)
        {
            shrani_sliko(iteration/IterationsPerFrame+1);
        }
        
        
        if (gsl_vector_get(F, 0) < 0)
        {
            std::cout << "Negativna sila " << iteration << std::endl;
            exit(4);
        }
        
        izracunaj_kot();
        
        
        ++iteration;
    }
}

void natancnost()
{
    int iteracije = 10000;
    const double phi = 1.0;
    
    int delitve[] = {10, 20, 30, 50, 75, 100, 200, 300, 500, 750, 1000, 0};
    
    std::ofstream f("g_natancnost_delitve.dat");
    
    for (int i = 0; delitve[i]; ++i)
    {
        zacetni_pogoj(delitve[i], phi);
        
        int iter = 0;
        while (iter < iteracije)
        {
            izracunaj_silo();
            izracunaj_kot();
            ++iter;
        }
        
        double x, y;
        polozaj_repa(x, y);
        f << delitve[i] << " " << x << " " << y << " " << energija() << std::endl;
    }
    
    f.close();
}

void zacetek_energije()
{
    int iteracije = 20;
    int IterationsPerEnergy = 1;
    const double phi = 1.0;
    
    int delitve[] = {10, 30, 100, 300, 1000, 0};
    char buf[64];
    
    for (int i = 0; delitve[i]; ++i)
    {
        sprintf(buf, "g_energija_zacetek_%d.dat", delitve[i]);
        std::ofstream f(buf);
        zacetni_pogoj(delitve[i], phi);
        
        f << "# Should be: E_0 = " << -sin(phi) / 2 << std::endl;
        f.precision(12);
        
        int iter = 0;
        while (iter < iteracije)
        {
            if (iter % IterationsPerEnergy == 0)
            {
                f << iter << " " << energija() << std::endl;
            }
            
            izracunaj_silo();
            izracunaj_kot();
            
            ++iter;
        }
    }
}
void ohranitev_energije()
{
    int iteracije = 1e6;
    int IterationsPerEnergy = 1e4;
    
    std::cout << iteracije << " " << IterationsPerEnergy << std::endl;
    
    const double phi = 1.0;
    
    int delitve[] = {10, 30, 100, 300, 1000, 0};
    char buf[64];
    
    for (int i = 0; delitve[i]; ++i)
    {
        sprintf(buf, "g_energija_%d.dat", delitve[i]);
        std::ofstream f(buf);
        zacetni_pogoj(delitve[i], phi);
        
        f << "# Should be: E_0 = " << -sin(phi) / 2 << std::endl;
        f.precision(12);
                
        for (int j = 0; j < 10; ++j)
        {
            izracunaj_silo();
            izracunaj_kot();
        }
        
        energy_t energy = energija();
        double E = energy.kineticna + energy.potencialna + energy.rotacijska;
        
        for (int iter = 0; iter < iteracije; ++iter)
        {
            if (iter % IterationsPerEnergy == 0)
            {
                f << iter+10 << " " << energija() << " " << E << std::endl;
            }
            
            izracunaj_silo();
            izracunaj_kot();
        }
    }
}

void ohranitev_korak()
{
    int ks[] = {1, 4, 10, 40, 100, 0};
    
    
    const double phi = 1.0;
    char buf[64];
    
    const double T = 10.0;
    
    
    for (int i = 0; ks[i]; ++i)
    {
        k = ks[i] * 1e-5;
        int IterationsPerEnergy = 0.01 / k;
        zacetni_pogoj(100, phi);
        
        double t = 0;
        for (int j = 0; j < 10; ++j)
        {
            izracunaj_silo();
            izracunaj_kot();
            t += k;
        }
        
        sprintf(buf, "g_energija_korak_%d.dat", ks[i]);
        std::ofstream f(buf);
        
        energy_t erg = energija();
        double E = erg.kineticna + erg.potencialna + erg.rotacijska;
        
        int iter = 0;
        while (t < T)
        {
            if (iter % IterationsPerEnergy == 0)
            {
                f << t << " " << energija() << " " << E << std::endl;
            }
            
            izracunaj_silo();
            izracunaj_kot();
            
            ++iter;
            t += k;
        }
        
        f.close();
    }
    
}

int main(int argc, char **argv) {
    
    QApplication app(argc, argv);
    
    
    int cleni = argc > 1 ? atoi(argv[1]) : 100;
    double kot = argc > 2 ? atof(argv[2]) : 0.9;
    int slike = argc > 3 ? atoi(argv[3]) : 1000;
    ImageFileNameFormat = argc > 4 ? argv[4] : "g_frame_%1.jpg";
    
    postopek(cleni, kot, slike);
    
    
    return 0;
}
