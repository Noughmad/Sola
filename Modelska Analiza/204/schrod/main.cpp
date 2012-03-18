#include <iostream>
#include <assert.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

inline double k(double x, double l, double e)
{
    return 2.0/x - l*(l+1)/x/x + e;
}

inline int numerov(double& y, double& y1, double& y2, double& k, double& k1, double& k2, double h)
{
    const double yy = (2*(1-5*h*k1)*y - (1+h*k2)*y1) / (1+h*k);
    
    y2 = y1;
    y1 = y;
    y = yy;
    
    k2 = k1;
    k1 = k;
    
    return GSL_SUCCESS;
}

int neskoncnost(const gsl_vector* y, void* params, gsl_vector* f)
{
    double h = 1e-4;
    double y0, y1, y2;
    double k0, k1, k2;
    y1 = 0;
    y0 = gsl_vector_get(y, 0) * h;
    double e = gsl_vector_get(y, 1);
    if (gsl_isnan(y0) || gsl_isnan(e))
    {
        gsl_vector_set(f, 0, y0);
        gsl_vector_set(f, 1, e);
        return GSL_EDOM;
    }
    double t = 2*h;
    double l = *(double*)params;
    double integral = 0;
    
    printf("Zacenjamo integracijo, y1 = %g, e = %g\n", y0, e);
    
   // FILE* file = fopen("g_testfile.dat", "wt");
    
    k2 = 0; // To je ubistvu k(x=0)
    k1 = k(h, l, e);
   // fprintf(file, "%g, %g, %g, %g, %g, %g, %g\n", t, y0, y1, y2, k0, k1, k2);
    double Phi = 0;
    while (t < 10)
    {
        k0 = k(t, l, e);
        numerov(y0, y1, y2, k0, k1, k2, h*h/12.0);
        integral += h*(fabs(y0*y0) + 4*fabs(y1*y1) + fabs(y2*y2))/6;
   //     fprintf(file, "%g, %g, %g, %g, %g, %g, %g\n", t, y0, y1, y2, k0, k1, k2);
        t += h;
    }
    
   // fclose(file);
    
    double r1 = y0/y1;
    double r2 = y1/y2;
    
    if (r1 < 1 && r2 < 1 && r1 > -1 && r2 > -1)
    {
        gsl_vector_set(f, 0, r1 - r2);
        integral += h*y0*y0 / (1-r1*r1);
    }
    else
    {
        gsl_vector_set(f, 0, y0);
    }
    
    // NOTE: To je manj pravilna formula, ker mora vrednost eksponentno padat

    
    gsl_vector_set(f, 1, integral-1);
    printf("Iteracija: %g, %g -> %g, %g, %g, %g\n", gsl_vector_get(y,0), e, y0, r1-1, r2-1, integral - 1);
    return GSL_SUCCESS;
}

int obestrani(const gsl_vector* y, void* params, gsl_vector* values)
{
    double h = 1e-4;
    double e = gsl_vector_get(y, 0);
    const double l = *(double*)params;
    
    double y2 = 0;
    double y1 = -h * gsl_vector_get(y, 1);
    double y0 = 0;
    
    double k2 = 0;
    double k1 = k(h, 0, e);
    double k0;
    
    double tm = 30;
    double z2 = gsl_vector_get(y, 2);
    double kk = fabs(k(20, l, e));
    double z1 = z2 * (1 + h*sqrt(kk));
    double z0 = 0;
    
    double integral = z2*z2 / (h*kk);
    double t1;
    
    FILE* f = fopen("g_test_obestrani.dat", "wt");
    
    for (t1 = 2*h; t1 < 5 + 1.5*h; t1 += h)
    {
        k0 = k(t1, l, e);
        numerov(y0, y1, y2, k0, k1, k2, h*h/12);
        integral += y0 * y0 * h;
        fprintf(f, "%g %g %g\n", t1, y0, k0);
    }
    
    k2 = k(tm, l, e);
    k1 = k(tm-h, l, e);
    fprintf(f, "\n");
    
    for (double t2 = tm-2*h; t2 > 5 - 0.5*h; t2 -= h)
    {
        k0 = k(t2, l, e);
        numerov(z0, z1, z2, k0, k1, k2, h*h/12);
        if (t2 > t1)
        {
            integral += z0 * z0 * h;
            fprintf(f, "%g %g %g\n", t2, z0, k0);
        }
    }
    fclose(f);
    integral -= (z0*z0+z1*z1+y0*y0+y1*y1)/2 * h;
    
    gsl_vector_set(values, 0, z0 + z1 - y0 - y1);
    gsl_vector_set(values, 1, z0 - z1 + y0 - y1);
    gsl_vector_set(values, 2, integral - 1);
    
    
    return GSL_SUCCESS;
}

void vodik()
{
    double l = 1;
    gsl_vector* v = gsl_vector_alloc(3);
    gsl_vector_set(v, 0, -0.3); // energija
    gsl_vector_set(v, 1, 1); // odvod v 0
    gsl_vector_set(v, 2, exp(-20)); // vrednost v 20
    
    gsl_multiroot_function function = {obestrani, 3, &l};
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 3);
    gsl_multiroot_fsolver_set(solver, &function, v);
    
    int status;
    int iter = 0;
    do
    {
        status = gsl_multiroot_fsolver_iterate(solver);
        status = gsl_multiroot_test_residual(solver->f, 1e-10);
        ++iter;
    }
    while (status == GSL_CONTINUE);
    printf("Koncal po %d iteracijah, status = %s\n", iter, gsl_strerror(status));
    
    gsl_vector_fprintf(stdout, solver->x, "%f");
    gsl_vector_fprintf(stdout, solver->f, "%f");
}

int main(int argc, char **argv) {
    vodik();
    std::cout << "Hello, world!" << std::endl;
    return 0;
}
