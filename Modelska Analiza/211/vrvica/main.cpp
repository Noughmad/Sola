#include <iostream>
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
double k = 6e-6;
double h;

double C; // C = k*k/h/h
double D; // D = h*h/k/k

const int IterationsPerFrame = 1e4;

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
    
    pn[0] = pn[2] + h/gsl_vector_get(F,1) * cos(pn[1]);
        
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
    for (int i = 0; i < StClenov; ++i)
    {
        Phi[i] = phi_0;
        Phi1[i] = phi_0;
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
        pen.setColor(QColor(255, 200 * gsl_vector_get(F, i), 0));
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

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    int cleni = argc > 1 ? atoi(argv[1]) : 100;
    double kot = argc > 2 ? atof(argv[2]) : 0.9;
    int slike = argc > 3 ? atoi(argv[3]) : 1000;
    ImageFileNameFormat = argc > 4 ? argv[4] : "g_frame_%1.jpg";
    
    postopek(cleni, kot, slike);
    return 0;
}
