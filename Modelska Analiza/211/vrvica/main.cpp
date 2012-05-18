#include <iostream>

#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_roots.h>

#include <QtGui/QImage>
#include <QtGui/QPainter>
#include <QtGui/QBitmap>
#include <QtGui/QApplication>

inline double sq(double x)
{
    return x*x;
}

const int VelikostSlike = 200;
const QString ImageFileNameFormat = "image_%1.jpg";

int StClenov;
int StClenovMinusDva;
double k = 0.001;
double h;

double C; // C = k*k/h/h
double D; // D = h*h/k/k

const int IterationsPerFrame = 1e4;

double *Phi, *Phi1, *pn;
gsl_vector* F;

QImage* image;
QBitmap* bitmap;
QPainter* painter;
QPainter* imagePainter;

gsl_vector *d, *u, *l, *rhs;

gsl_root_fsolver* solver;
gsl_function function;

double fn(double x, void* params)
{
    return x - Phi[1] - h*gsl_vector_get(F,0) * cos(x);
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
    gsl_root_fsolver_set (solver, &function, 0.5 * Phi[0], 2.0 * Phi[0]);
    
    int status;
    do
    {
           gsl_root_fsolver_iterate (solver);
           double r = gsl_root_fsolver_root (solver);
           double x_lo = gsl_root_fsolver_x_lower (solver);
           double x_hi = gsl_root_fsolver_x_upper (solver);
           status = gsl_root_test_interval (x_lo, x_hi, 0, 1e-5);
    }
    while (status == GSL_CONTINUE);  
    
    pn[0] = gsl_root_fsolver_root(solver);
        
    Phi1 = Phi;
    Phi = pn;
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
    
    bitmap = new QBitmap(2*VelikostSlike, VelikostSlike);
    image = new QImage(2*VelikostSlike, VelikostSlike, QImage::Format_RGB16);
    painter = new QPainter(bitmap);
    QPen pen;
    pen.setColor(Qt::color1);
    pen.setWidth(5);
    painter->setPen(pen);
    
    imagePainter = new QPainter(image);
    imagePainter->setPen(Qt::white);
    imagePainter->setBackgroundMode(Qt::TransparentMode);
    
    solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
    function.function = fn;
    function.params = 0;
}

void shrani_sliko(int frame)
{
    bitmap->clear();
    QPointF lastPos = QPointF(VelikostSlike, 0);
    QPointF currentPos;
    
    const double h = 1.0/StClenov * VelikostSlike;
    
    for (int i = 0; i < StClenov; ++i)
    {
        currentPos = lastPos + h*QPointF(cos(Phi[i]), sin(Phi[i]));
        painter->drawLine(lastPos, currentPos);
        lastPos = currentPos;
    }
    
    image->fill(QColor(Qt::black));
    imagePainter->drawPixmap(0,0,*bitmap);
    image->save(ImageFileNameFormat.arg(frame, 4, 10, QLatin1Char('0')));
    
    std::cout << "Shranil sliko " << frame << std::endl;
}

void korak()
{
    izracunaj_silo();
    izracunaj_kot();
}

void postopek(int n, double phi, double T)
{
    zacetni_pogoj(n, phi);
    double t = 0;
    int iteration = 0;
    while (t < T)
    {
        if (iteration % IterationsPerFrame == 0)
        {
            shrani_sliko(iteration/IterationsPerFrame+1);
        }
        
        izracunaj_silo();
        izracunaj_kot();
        
        
        ++iteration;
        t += k;
    }
}

int main(int argc, char **argv) {
    QApplication app(argc, argv);
    postopek(100, 0.5, 10000.0);
    return 0;
}
