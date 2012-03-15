#include <iostream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

double k(double x, double l, double e)
{
    if (x == 0)
    {
        // TODO
    }
    return 2/x - l*(l+1)/x/x + e;
}

int numerov(double* y, double* y1, double* y2, double* k, double* k1, double* k2, double h)
{
    double yy = (2*(1-5*h*h* (*k1)/12)* (*y1) - (1+h*h* (*k2)/12)* (*y2)) / (1+h*h* (*k)/12);
    
    *y2 = *y1;
    *y1 = *y;
    *y = yy;
    *k2 = *k1;
    *k1 = *k;
}

int neskoncnost(const gsl_vector* y, void* params, gsl_vector* f)
{
    double y0, y1, y2;
    double k0, k1, k2;
    y2 = 0;
    y1 = gsl_vector_get(y, 0);
    double t = 0;
    double h = 1e-6;
    
    while (t < 20)
    {
        k0 = k(t, 0, 0);
        numerov(&y0, &y1, &y2, &k0, &k1, &k2, h);
        t += h;
    }
    gsl_vector_set(f, 0, y0);
    return GSL_SUCCESS;
}

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    return 0;
}
