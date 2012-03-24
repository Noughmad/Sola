#include <iostream>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_min.h>
#include <math.h>
#include <ctime>

int N = 90; // Stevilo tock na enem robu, vseh je NxN
gsl_matrix* U;
double h = 1.0/N;

inline double omega(int korak, double alpha = 1)
{
    if (korak < 3)
    {
        return 1;
    }
    else
    {
        return 2.0/(1+alpha*M_PI/N);
    }
}

inline int meja(int i)
{
    if (i < N/3)
    {
        return N-i;
    }
    else if (i < 2*N/3)
    {
        return 2*N/3;
    }
    else
    {
        return i;
    }
}

double korak_cev(double rho, double omega)
{
    double popravek = 0;
    
    // Popravimo sode tocke, razen tistih na robu
    
    for (int i = 1; i < N/2; ++i)
    {
        const int mj = meja(i);
        for (int j = 1+(i%2); j < mj; j += 2)
        {
            double eps = gsl_matrix_get(U, i+1, j) + gsl_matrix_get(U, i-1, j) + gsl_matrix_get(U, i, j+1) + gsl_matrix_get(U, i, j-1) - 4*gsl_matrix_get(U, i, j) - h*h*rho;
            
            gsl_matrix_set(U, i, j, gsl_matrix_get(U, i,j) + eps * 0.25 * omega);
            popravek += eps * eps;
        }
    }
    
    // Zagotovimo se robni pogoj na odprtem robu (sredini cevi)
    
    int ti = N/2;
    int tmj = meja(ti);
    for (int j = 1+((N/2)%2); j < tmj; j += 2)
    {
        double eps = 2*gsl_matrix_get(U, ti-1, j) + gsl_matrix_get(U, ti, j+1) + gsl_matrix_get(U, ti, j-1) - 4*gsl_matrix_get(U, ti, j) - h*h*rho;
        
        gsl_matrix_set(U, ti, j, gsl_matrix_get(U, ti,j) + eps * 0.25 * omega);
        popravek += eps * eps;
    }

    // In se lihe tocke
    
    for (int i = 1; i < N/2; ++i)
    {
        const int mj = meja(i);
        for (int j = 2-(i%2); j < mj; j += 2)
        {
            double eps = gsl_matrix_get(U, i+1, j) + gsl_matrix_get(U, i-1, j) + gsl_matrix_get(U, i, j+1) + gsl_matrix_get(U, i, j-1) - 4*gsl_matrix_get(U, i, j) - h*h*rho;
            
            gsl_matrix_set(U, i, j, gsl_matrix_get(U, i,j) + eps * 0.25 * omega);
            popravek += eps * eps;
        }
    }
    
    // In se lihe tocke na sredini cevi
    
    for (int j = 2-((N/2)%2); j < tmj; j += 2)
    {
        double eps = 2*gsl_matrix_get(U, ti-1, j) + gsl_matrix_get(U, ti, j+1) + gsl_matrix_get(U, ti, j-1) - 4*gsl_matrix_get(U, ti, j) - h*h*rho;
        
        gsl_matrix_set(U, ti, j, gsl_matrix_get(U, ti,j) + eps * 0.25 * omega);
        popravek += eps * eps;
    }
    
    return popravek;
}

double valj(double rho, double omega)
{
    double popravek = 0;
    
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            
        }
    }
}

int cev(double omg, bool print)
{
    U = gsl_matrix_alloc(N/2+1,N+1);
    gsl_matrix_set_zero(U);
    int korak = 0;
    double eps;
    do
    {
        eps = korak_cev(1.0, omg);
        ++korak;
//        std::cout << "Popravek pri koraku " << korak << " je bil " << eps << std::endl;
    }
    while (eps > 1e-8);
    
    if (print)
    {
        FILE* f = fopen("g_test_cev.dat", "wt");
        for (int j  = 0; j < N; ++j)
        {
            for (int i = 0; i < N/2; ++i)
            {
                fprintf(f, "%g\t", gsl_matrix_get(U, i, j));
            }
            for (int i = N/2-1; i > 0; --i)
            {
                fprintf(f, "%g\t", gsl_matrix_get(U, i, j));
            }
            fprintf(f, "%g\n", gsl_matrix_get(U, 0, j));
        }
        fclose(f);
    }
    
    gsl_matrix_free(U);
    U = 0;
    return korak;
}

double trajanje(double x, void* param)
{
    return cev(x, false);
}

void najdi_omega(FILE* f)
{
    int status;
    int iter = 0, max_iter = 100;
    gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (gsl_min_fminimizer_brent);
    gsl_function F;

    F.function = trajanje;
    F.params = 0;

    gsl_min_fminimizer_set (s, &F, 2.0/(1+M_PI/N), 1.0, 1.999);
    
    double m, a, b;
    do
    {
        iter++;
        status = gsl_min_fminimizer_iterate (s);

        m = gsl_min_fminimizer_x_minimum (s);
        a = gsl_min_fminimizer_x_lower (s);
        b = gsl_min_fminimizer_x_upper (s);

        status 
            = gsl_min_test_interval (a, b, 1e-6, 0.0);

        if (status == GSL_SUCCESS)
            printf ("Converged:\n");

        printf ("%5d [%.7f, %.7f] "
                "%.7f %+.7f %.7f\n",
                iter, a, b,
                m, m - 1.0, b - a);
    }
    while (status == GSL_CONTINUE);
    
    fprintf(f, "%d, %g, %g\n", N, m, gsl_min_fminimizer_f_minimum(s));
}

void hitrost()
{
    FILE* f = fopen("g_hitrost_omega.dat", "wt");
    int Ns[] = {30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 0};
    for (int i = 0; Ns[i] != 0; ++i)
    {
        N = Ns[i];
        najdi_omega(f);
    }
    fclose(f);
}

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    N = 90;
    cev(1.8, true);
    return 0;
}
