#include <iostream>
#include <map>
#include <assert.h>
#include <stdio.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_odeiv2.h>

double* R;
double* Phi;
double* Rint;
double* Ryint;
int N;
double step;

double energ_He(double Z)
{
  double e = R[1]*R[1]/step/step;
  for (int i = 1; i < N-1; ++i)
  {
    e += pow((R[i+1]-R[i])/step, 2) - 2*Z/i/step * R[i] * R[i] - Phi[i]*R[i]*R[i]/i/step;
  }
  return 2 * e * step;
}

inline double k(double x, double l, double e)
{
    double ret = 2.0/x + e;
    if (fabs(l) > 0.1 && fabs(l+1) > 0.1)
    {
      ret -= l*(l+1)/x/x;
    }
    return ret;
}

inline double kHe(int n, double Z, double e)
{
  const double k = 2.0 * (Z + Phi[n])/(step*n) + e;
  // std::cout << "K_He(" << n << ", " << Z << ", " << e <<") = " << k << std::endl;
  return k;
}

inline int odvod_pot(int n, const double y[], double f[], void* params)
{
  f[0] = y[1];
  
  if (n == 0)
  {
    f[1] = 0;
  }
  else
  {
    f[1] = R[n] * R[n] / step / n;
  }
  
  return GSL_SUCCESS;
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
  double yy[2] = {gsl_vector_get(y, 0), gsl_vector_get(y, 1)};
  Phi[0] = yy[0];
  Phi[1] = Phi[0] + step * yy[1];
  
  for (int i = 2; i < N; ++i)
  {
      double f[2];
      
      yy[0] = Phi[i-2];
      yy[1] = (Phi[i-1] - Phi[i-2])/step;
      
      odvod_pot(i-2, yy, f, 0);
      yy[0] += f[0] * step/6;
      yy[1] += f[1] * step/6;
      odvod_pot(i-1, yy, f, 0);
      yy[0] += f[0] * 4*step/6;
      yy[1] += f[1] * 4*step/6;
      odvod_pot(i, yy, f, 0);
      yy[0] += f[0] * step/6;
      yy[1] += f[1] * step/6;
      
      Phi[i] = yy[0];
      
      if (gsl_isnan(Phi[i]) || gsl_isinf(Phi[i]))
      {
        std::cout << "Nan v potencialu: " << i << ", " << Phi[i-1] << std::endl;
        exit(33);
      }
  }
    
  gsl_vector_set(v, 0, Phi[N-1]);
  gsl_vector_set(v, 1, Phi[N-2]);
  
  return GSL_SUCCESS;
}


inline void R_to_Phi()
{
  /*
  gsl_multiroot_function function = {potencial, 2, 0};
  gsl_vector* v = gsl_vector_alloc(2);
  gsl_vector_set(v, 0, Phi[0]);
  gsl_vector_set(v, 1, (Phi[1]-Phi[0])/step);
  gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 2);
  
  gsl_multiroot_fsolver_set(solver, &function, v);
  int status;
  int iter = 0;
  do
  {
      status = gsl_multiroot_fsolver_iterate(solver);
      if (status != GSL_SUCCESS)
      {
          break;
      }
      status = gsl_multiroot_test_residual(solver->f, 1e-12);
      ++iter;
  }
  while (status == GSL_CONTINUE);
  
  std::cout << "Izracunal Phi(x) po " << iter << " iteracijah, status = " << gsl_strerror(status) << std::endl; 
  */
  
  Rint[0] = 0;
  Ryint[0] = 0;
  for (int i = 1; i < N-1; ++i)
  {
    Rint[i] = Rint[i-1] + step * R[i] * R[i];
    Ryint[i] = Ryint[i-1] + R[N-i-1] * R[N-i-1] / (N-i-1);
  }
  Rint[N-1] = Rint[N-2] + step * R[N-1] * R[N-1];
  Ryint[N-1] = Ryint[N-2];
  
  for (int i = 0; i < N; ++i)
  {
    Phi[i] = -Rint[i] - i*step*Ryint[N-i-1];
    if (gsl_isnan(Phi[i]) || gsl_isinf(Phi[i]))
    {
      std::cout << "Potencial je nan " << i << ", " << Phi[i-1] << std::endl;
    }
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
    const int M = N/2;
    
    double h = step;
    double e = gsl_vector_get(y, 0);
    const double l = *(double*)params;
    
    R[0] = 0;
    R[1] = h * gsl_vector_get(y, 1);
    
    std::cout << "e = " << e << ", R(x1) = " << R[1] << std::endl;
    
    for (int i = 0; i < 3; ++i)
    {
      if (gsl_isnan(gsl_vector_get(y, i)))
      {
        return GSL_EDOM;
      }
    }
    
    if (e > 0)
    {
      for (int i = 0; i < 3; ++i)
      {
        std::cout << "Positive energy, returning a big number" << std::endl;
        gsl_vector_set(values, i, exp(10*e));
      }
      return GSL_SUCCESS;
    }
    
    double k2 = 0;
    double k1 = k(h, l, e);
    double k0;
    
    double tm = (N-1)*h;
    
    for (int i = 2; i < M + 2; ++i)
    {
        k0 = k(i*h, l, e);
        numerov(i, k0, k1, k2, false);
    }
    
    double y0 = R[M];
    double y1 = R[M+1];
    
    k2 = k(tm, l, e);
    k1 = k(tm-h, l, e);
    
    R[N-1] = gsl_vector_get(y, 2);
    R[N-2] = R[N-1] * (1 + h*sqrt(fabs(k1)));
    
    for (int i = N-3; i >= M; --i)
    {
        k0 = k(i*step, l, e);
        numerov(i, k0, k1, k2, true);
    }
    gsl_vector_set(values, 0, R[M+1] + R[M] - y0 - y1);
    gsl_vector_set(values, 1, 1e4 * (R[M+1] - R[M] + y0 - y1));
    
    double I = integral_RR(0,N);
    if (gsl_isinf(I))
    {
      I = 1e200;
    }
    gsl_vector_set(values, 2, I - 1);
    
    gsl_vector_fprintf(stdout, values, "%g");
    
    
    return GSL_SUCCESS;
}


void normiraj()
{
  const double I = 1.0/sqrt(integral_RR(0, N));
  
  for (int i = 0; i < N; ++i)
  {
    R[i] *= I;
  }
}

double enastran_vodik(double e, double l)
{
  R[0] = 0;
  R[1] = step;
  
  double k2 = 0;
  double k1 = k(step, l, e);
  double k0;
    
  for (int i = 2; i < N; ++i)
  {
      k0 = k(i*step, l, e);
      numerov(i, k0, k1, k2, false);
  }
  
  return R[N-1];
}

typedef double (*bis_fun)(double, double);

double bisekcija_res(double plus, double minus, double l, bis_fun fun)
{
  double center = (plus+minus)/2;
  if (fabs(plus-minus) < 1e-15)
  {
    return center;
  }
  if (fun(center, l) > 0)
  {
    return bisekcija_res(center, minus, l, fun);
  }
  else
  {
    return bisekcija_res(plus, center, l, fun);
  }
}

double bisekcija(double plus, double minus, double l, bis_fun fun)
{
  double ret;
  if (enastran_vodik(plus, l) < 0)
  {
    ret = bisekcija_res(minus, plus, l, fun);
  }
  else
  {
    ret = bisekcija_res(plus, minus, l, fun);
  }
  normiraj();
  
  std::cout << "Bisekcija uspela, energija je " << ret << std::endl;
  
  return ret;
}

double bisekcija_vodik(double plus, double minus, double l)
{
    return bisekcija(plus, minus, l, enastran_vodik);
}


void en_vodik(double l)
{
  
}

double enastran_helij(double e, double Z)
{
  R[0] = 0;
  R[1] = step;
  
  double k2 = 0;
  double k1 = kHe(1, Z, e);
  double k0;
    
  for (int i = 2; i < N; ++i)
  {
      k0 = kHe(i, Z, e);
      numerov(i, k0, k1, k2, false);
      
      if (gsl_isnan(R[i]) || gsl_isinf(R[i]))
      {
        std::cout << i << " is nan: " << R[i-1] << std::endl;
        exit(15);
      }
  }
  
  return R[N-1];
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
          std::cout << "Domain error, something is nan" << std::endl;
          return GSL_EDOM;
      }
    }
    
    if (e > 0)
    {
      for (int i = 0; i < 3; ++i)
      {
        std::cout << "Positive energy, returning a big number: " << e << std::endl;
        gsl_vector_set(values, i, exp(100*e));
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
    
    double I = integral_RR(0,N) + R[N-1]/sqrt(fabs(k0));
    if (gsl_isinf(I))
    {
      I = 1e200;
    }
    gsl_vector_set(values, 2, I - 1);
    
   printf("%g, %g, %g --> %g, %g, %g\n", e, R[1], R[N-1], gsl_vector_get(values, 0), gsl_vector_get(values, 1), gsl_vector_get(values, 2));
   std::cout << std::flush;
    
    return GSL_SUCCESS;
}

void vodik(double E, double L)
{
    int n = E < -0.5 ? 30 : E < -0.2 ? 50 : 100;
    N = n * 1e3;
    step = 1e-3;
    R = new double[N];
    Phi = new double[N];

    double l = L;
    gsl_vector* v = gsl_vector_alloc(3);
    gsl_vector_set(v, 0, E); // energija
    gsl_vector_set(v, 1, 0.1); // odvod v 0
    gsl_vector_set(v, 2, exp(-30)); // vrednost v 20
    
    
    gsl_multiroot_function function = {obestrani, 3, &l};
    gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_dnewton, 3);
    gsl_multiroot_fsolver_set(solver, &function, v);
    
    int status;
    int iter = 0;
    do
    {
        status = gsl_multiroot_fsolver_iterate(solver);
        if (status != GSL_SUCCESS)
        {
            break;
        }
        status = gsl_multiroot_test_residual(solver->f, 1e-8);
        
        gsl_vector_fprintf(stdout, solver->f, "%g");
        
        ++iter;
    }
    while (status == GSL_CONTINUE);
    printf("Koncal po %d iteracijah, status = %s\n", iter, gsl_strerror(status));
    
    char buf[32];
    sprintf(buf, "g_vodik_%g_%g.dat", fabs(E), L);
    plot_data(buf);
    
    delete[] R;
    delete[] Phi;
}

bool preveri_helij(double eps)
{
    double check = 0;
    double o;
    for (int i = 1; i < N-1; ++i)
    {
        o = (Phi[i-1] + Phi[i+1] - 2*Phi[i])/step/step - R[i]*R[i]/i/step;
        check += o*o;
    }
    std::cout << " Preveri_helij: povprecno odstopanje je " << check*step << std::endl;
    return check*step > eps;
}



void helij(double plus, double minus, int n, const char* filename)
{
    N = n * 1e4;
    step = 1e-4;
    R = new double[N];
    Phi = new double[N];
    Rint = new double[N];
    Ryint = new double[N];
    
    assert(R);
    assert(Phi);
    
    
    double Z = 2;
    
    double ZZ = Z - 6.0/16.0;
    for (int i = 0; i < N; ++i)
    {
        // Zacetni priblizek za R(x)
        R[i] = 2*sqrt(ZZ) * i*step * ZZ * exp(-ZZ*i*step);
    }
    
    int r = 0;
    double e;
    do
    {
        R_to_Phi();
        plot_data("g_helij_med.dat");
        e = bisekcija(plus, minus, Z, enastran_helij);
        normiraj();
        ++r;
    }
    while (r < 30 && preveri_helij(5e-7));
    
    std::cout << "Resil helij po " << r << " iteracijah, energija je " << e;
    
    delete[] R;
    delete[] Phi;
    delete[] Rint;
    delete[] Ryint;
}

void spekter()
{
  FILE* f = fopen("g_vodik_spekter.dat", "wt");
  for (double energ = 0; energ > -1.1; energ -= 0.01)
  {
    fprintf(f, "%g, %g", energ, enastran_vodik(energ, 0));
  }
  fclose(f);
}


int main(int argc, char **argv) {

  /*
  step = 1e-3;
  N = 50 * 1e3;
  
  R = new double[N];
  Phi = new double[N];
  
  // spekter();
  
  std::cout << bisekcija_vodik(-1.1, -0.9, 0) << std::endl;
  plot_data("g_bisekcija_vodik_1_0.dat");
  
  std::cout << bisekcija_vodik(-0.24, -0.26, 0) << std::endl;
  plot_data("g_bisekcija_vodik_2_0.dat");
  std::cout << bisekcija_vodik(-0.24, -0.26, 1) << std::endl;
  plot_data("g_bisekcija_vodik_2_1.dat");
  
  std::cout << bisekcija_vodik(-0.12, -0.1, 0) << std::endl;
  plot_data("g_bisekcija_vodik_3_0.dat");
  std::cout << bisekcija_vodik(-0.12, -0.1, 1) << std::endl;
  plot_data("g_bisekcija_vodik_3_1.dat");
  std::cout << bisekcija_vodik(-0.12, -0.1, 2) << std::endl;
  plot_data("g_bisekcija_vodik_3_2.dat");
  
  delete[] R;
  delete[] Phi;
  
  N = 100 * 1e3;
  R = new double[N];
  Phi = new double[N];
  
  std::cout << bisekcija_vodik(-0.06, -0.07, 0) << std::endl;
  plot_data("g_bisekcija_vodik_4_0.dat");
  std::cout << bisekcija_vodik(-0.06, -0.07, 1) << std::endl;
  plot_data("g_bisekcija_vodik_4_1.dat");
  std::cout << bisekcija_vodik(-0.06, -0.07, 2) << std::endl;
  plot_data("g_bisekcija_vodik_4_2.dat");
  std::cout << bisekcija_vodik(-0.06, -0.07, 3) << std::endl;
  plot_data("g_bisekcija_vodik_4_3.dat");
  
  std::cout << bisekcija_vodik(-0.042, -0.038, 0) << std::endl;
  plot_data("g_bisekcija_vodik_5_0.dat");
  std::cout << bisekcija_vodik(-0.042, -0.038, 1) << std::endl;
  plot_data("g_bisekcija_vodik_5_1.dat");
  std::cout << bisekcija_vodik(-0.042, -0.038, 2) << std::endl;
  plot_data("g_bisekcija_vodik_5_2.dat");
  std::cout << bisekcija_vodik(-0.042, -0.038, 3) << std::endl;
  plot_data("g_bisekcija_vodik_5_3.dat");
  std::cout << bisekcija_vodik(-0.042, -0.038, 4) << std::endl;
  plot_data("g_bisekcija_vodik_5_4.dat");
  
  delete[] R;
  delete[] Phi;
  */
  
  // Ti dve zacetni stanji najdeta ustrezne resitve
  helij(-0.1, -0.3, 20, "g_helij_2.dat");
  helij(-0.1, -2, 20, "g_helij_2.dat");
  
  // TODO: Najdi se ostale mozne energije
}
