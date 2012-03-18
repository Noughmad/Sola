#include <iostream>
#include <assert.h>
#include <stdio.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

double* R;
double* Phi;
double* Rint;
double* Ryint;
int N;
double step;

inline double k(double x, double l, double e)
{
    return 2.0/x - l*(l+1)/x/x + e;
}

inline double kHe(int n, double Z, double e)
{
  const double k = 2.0 * Z / (step*n) + 2*Phi[n] + e;
  // std::cout << "K_He(" << n << ", " << Z << ", " << e <<") = " << k << std::endl;
  return k;
}

inline double kPot(int n)
{
  return - R[n] * R[n] / n / step;
}

inline void numerov_phi(size_t n, double k, double& k1, double& k2)
{
    const double h = step*step/12;
    Phi[n] = (2*(1-5*h*k1)*Phi[n-1] - (1+h*k2)*Phi[n-2]) / (1+h*k);
    
    k2 = k1;
    k1 = k;
}

void plot_data(const char* filename)
{
  FILE* f = fopen(filename, "wt");
  
  for (int i = 0; i < N; ++i)
  {
    fprintf(f, "%g, %g, %g\n", i*step, R[i], Phi[i]);
  }
  
  fclose(f);
}

double integral_RR(size_t from, size_t to)
{
  // TODO: use more advanced formula
  double i = (R[from] + R[to-1])/2;
  for (int j = from+1; j < to-1; ++j)
  {
    i += R[j]*R[j];
  }
  return i*step;
}

double integral_RR_x(size_t from, size_t to)
{
  // TODO: use more advanced formula
  double i = (R[from] + R[to-1])/2;
  for (int j = from+1; j < to-1; ++j)
  {
    i += R[j]*R[j]/j;
  }
  return i;
}

int potencial(const gsl_vector* y, void* params, gsl_vector* v)
{
  Phi[0] = -0.25;
  Phi[1] = Phi[0] + step * gsl_vector_get(y, 0);
  double k2 = 0;
  double k1 = kPot(1);
  double k0;
  
  for (int i = 2; i < N; ++i)
  {
    k0 = kPot(i);
    numerov_phi(i, k0, k1, k2);
  }
  
  gsl_vector_set(v, 0, Phi[N-1]);
}


inline void R_to_Phi()
{
  /*
   * Namesto integriranja je mogoce boljse resevat enacbo
  Rint[0] = 0;
  Ryint[N-1] = 0;
  for (int i = 1; i < N-1; ++i)
  {
    // TODO: Higher-order formula
    Rint[i] = Rint[i-1] + R[i] * R[i] * step;
    Ryint[N-i-1] = Ryint[N-i] + R[N-i-1] * R[N-i-1] / (N-i-1);
  }
  
  R[N-1] = R[N-2] + R[N-1] * R[N-1];
  Ryint[0] = Ryint[1];
  
  Phi[0] = Ryint[0];
  for (int i = 1; i < N; ++i)
  {
    Phi[i] = -Rint[i]/step/i - Ryint[i];
  }
  */
  
  gsl_multiroot_function function = {potencial, 1, 0};
  gsl_vector* v = gsl_vector_alloc(1);
  gsl_vector_set(v, 0, 1);
  gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 1);
  
  gsl_multiroot_fsolver_set(solver, &function, v);
  int status;
  int iter = 0;
  do
  {
      status = gsl_multiroot_fsolver_iterate(solver);
      status = gsl_multiroot_test_residual(solver->f, 1e-8);
      ++iter;
  }
  while (status == GSL_CONTINUE);
  
  for (int i = 1; i < N; ++i)
  {
    Phi[i] = Phi[i] / (i*step);
  }
}

inline void numerov(size_t n, double k, double& k1, double& k2, bool reverse)
{
    const double h = step*step/12;
    if (reverse)
    {      
      assert(n < N-2);
      assert(n >= 0);
      R[n] = (2*(1-5*h*k1)*R[n+1] - (1+h*k2)*R[n+2]) / (1+h*k);
    }
    else
    {
      assert(n > 1);
      assert(n <= N-1);
      R[n] = (2*(1-5*h*k1)*R[n-1] - (1+h*k2)*R[n-2]) / (1+h*k);
    }
    
    k2 = k1;
    k1 = k;
}

int obestrani(const gsl_vector* y, void* params, gsl_vector* values)
{
    const int M = N/5;
    
    double h = step;
    double e = gsl_vector_get(y, 0);
    const double l = *(double*)params;
    
    R[0] = 0;
    R[1] = -h * gsl_vector_get(y, 1);
    
    double k2 = 0;
    double k1 = k(h, 0, e);
    double k0;
    
    double tm = (N-1)*h;
    
    for (int i = 2; i < M + 2; ++i)
    {
        k0 = k(i*h, l, e);
        numerov(i, k0, k1, k2, false);
        if (y->size != 3)
        {
          std::cout << y->size << " " << i << " " << M << std::endl;
        }
    }
    
    double y0 = R[M];
    double y1 = R[M+1];
    
    k2 = k(tm, l, e);
    k1 = k(tm-h, l, e);
    
    std::cout << y->size << std::endl;

    R[N-1] = gsl_vector_get(y, 2);
    R[N-2] = R[N-1] * (1 + h*sqrt(fabs(k1)));
    
    for (int i = N-3; i >= M; --i)
    {
        k0 = k(i*step, l, e);
        numerov(i, k0, k1, k2, true);
    }
    
    gsl_vector_set(values, 0, R[M+1] + R[M] - y0 - y1);
    gsl_vector_set(values, 1, R[M+1] - R[M] + y0 - y1);
    gsl_vector_set(values, 2, integral_RR(0,N) - 1);
    
    
    return GSL_SUCCESS;
}


int obestrani_helij(const gsl_vector* y, void* params, gsl_vector* values)
{
    const int M = N/2;
    
    
    double h = step;
    double e = gsl_vector_get(y, 0);
    
    
    for (int i = 0; i < 3; ++i)
    {
      if (gsl_isnan(gsl_vector_get(y, i)))
      {
        plot_data("g_fail_data.dat");
        exit(i+1);
      }
    }
    
    if (e > 0)
    {
      for (int i = 0; i < 3; ++i)
      {
        gsl_vector_set(values, i, exp(e));
      }
      return GSL_SUCCESS;
    }
    
    const double Z = *(double*)params;
        
    R[0] = 0;
    R[1] = gsl_vector_get(y, 1);
    
    double k2 = 0;
    double k1 = k(h, 0, e);
    double k0;
    
    double tm = (N-1)*h;
    
    for (int i = 2; i < M + 2; ++i)
    {
        k0 = kHe(i, Z, e);
        numerov(i, k0, k1, k2, false);
    }
    
    double y0 = R[M];
    double y1 = R[M+1];
    double y_test = R[M-1];
    
  //  std::cout << "Integracija naprej: " << M << " => " << y0 << ", " << y1 << ", " << R[M-1] << std::endl;
    
    k2 = kHe(N-1, Z, e);
    k1 = kHe(N-2, Z, e);
    
    R[N-1] = gsl_vector_get(y, 2);
    R[N-2] = R[N-1] * (1 + h*sqrt(fabs(k1)));
    
  //  std::cout << "k1 = " << k1 << ", R[N-1] = " << R[N-1] << std::endl;
    
    for (int i = N-3; i >= M; --i)
    {
        k0 = kHe(i, Z, e);
        numerov(i, k0, k1, k2, true);
    }
    
    assert(y_test = R[M-1]);
    
 //   std::cout << "Integracija nazaj: " << M << " => " << R[M+1] << ", " << R[M] << ", " << R[M-1] << std::endl;
    
    gsl_vector_set(values, 0, R[M+1] + R[M] - y0 - y1);
    gsl_vector_set(values, 1, 1e3 * (R[M+1] - R[M] + y0 - y1));
    
    double I = integral_RR(0,N);
    if (gsl_isinf(I))
    {
      I = 1e200;
    }
    gsl_vector_set(values, 2, I - 1);
    
   printf("%g, %g, %g --> %g, %g, %g\n", e, R[1], R[N-1], gsl_vector_get(values, 0), gsl_vector_get(values, 1), gsl_vector_get(values, 2));
   std::cout << std::flush;
    
    return GSL_SUCCESS;
}

void vodik()
{
    N = 30 * 1e3;
    step = 1e-3;
    R = new double[N];
    Phi = new double[N];
    
    double l = 0;
    gsl_vector* v = gsl_vector_alloc(3);
    gsl_vector_set(v, 0, -1); // energija
    gsl_vector_set(v, 1, 1); // odvod v 0
    gsl_vector_set(v, 2, exp(-N*step)); // vrednost v 20
    
    
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
    // printf("Koncal po %d iteracijah, status = %s\n", iter, gsl_strerror(status));
    
    delete[] R;
    delete[] Phi;
}

void helij()
{
    N = 100 * 1e3;
    step = 1e-3;
    R = new double[N];
    Phi = new double[N];
    Rint = new double[N];
    Ryint = new double[N];
    
    assert(R);
    assert(Phi);
    
    
    double Z = 2;
    
    double ZZ = Z - 5.0/16.0;
    for (int i = 0; i < N; ++i)
    {
        // Zacetni priblizek za R(x)
        R[i] = 2*sqrt(ZZ) * i*step * ZZ * exp(-ZZ*i*step);
    }
    
    R_to_Phi();
    plot_data("g_initial.dat");
    
    std::cout << "Zacetni priblizek ima integral " << integral_RR(0, N) << std::endl;
    
    gsl_vector* v = gsl_vector_alloc(3);
    gsl_vector_set(v, 0, -0.1); // energija
    gsl_vector_set(v, 1, R[1]); // odvod v 0
    gsl_vector_set(v, 2, 0.01); // vrednost v 20
    
    
    gsl_multiroot_function function = {obestrani_helij, 3, &Z};
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 3);
    
    
    int r = 0;
    while (r < 100)
    {
        R_to_Phi();
        
        std::cout << "Izracunal Phi(x), zdaj racunam R" << std::endl;

        gsl_multiroot_fsolver_set(solver, &function, v);
        int status;
        int iter = 0;
        do
        {
            status = gsl_multiroot_fsolver_iterate(solver);
            status = gsl_multiroot_test_residual(solver->f, 1e-4);
            ++iter;
        }
        while (status == GSL_CONTINUE);
        
        plot_data("g_med_iteracijami.dat");
        
        std::cout << "Izracunal R(x), zdaj racunam Phi" << std::endl;
        
        for (int i = 0; i < 3; ++i)
        {
          gsl_vector_set(v, i, gsl_vector_get(solver->x, i));
        }
        
        ++r;
    }
    
    
    delete[] R;
    delete[] Phi;
    delete[] Rint;
    delete[] Ryint;
}

int main(int argc, char **argv) {
    vodik();
    std::cout << "Hello, world!" << std::endl;
    helij();
    return 0;
}
