#include <iostream>
#include <cassert>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_blas.h>

#define STEP_START 0.1
#define STEP 0.05
#define FFOR(var) for(double var = STEP_START; var < 0.5; var += STEP)

typedef int (*odvod)(double, const double[], double[], void*);

typedef struct {
   gsl_odeiv2_system system;
   gsl_vector* initial_state; // zacetna (x,y)
   gsl_vector* final_state; // koncna (x,y)
} m_system;

int orbita_n;

int f_vrv(double t, const double y[], double f[], void* params)
{
    // x, y, alpha, F
    const double beta = *(double*)params;
    f[0] = cos(y[3]);
    f[1] = sin(y[3]);
    
    f[2] = -beta * y[0] * f[0] + f[1];
    f[3] = (beta * y[0] * f[1] + f[0]) / y[2];
    
    return GSL_SUCCESS;
}

int f_kristal(double t, const double y[], double f[], void* params)
{
    if (y[0] < -0.5 || y[0] > 0.5 || y[1] < -0.5 || y[1] > 0.5)
    {
        return GSL_EDOM;
    }
    
    const double D = 1.0 / (sin(M_PI*y[0]) * sin(M_PI*y[0]) + sin(M_PI*y[1]) * sin(M_PI*y[1]));
    f[0] = y[2];
    f[1] = y[3];
    
    f[2] = -M_PI_2 * sin(2 * M_PI * y[0]) * D;
    f[3] = -M_PI_2 * sin(2 * M_PI * y[1]) * D;
    
 //   std::cout << "t, x, y = " << t << ", " << y[0] << ", "<< y[1] << ", " << y[2] << ", " << y[3] << std::endl;
    
    return GSL_SUCCESS;
}


void print_orbit(const double v[])
{
  std::cout << "Printing orbit " << orbita_n << std::endl;
  ++orbita_n;
  char buf[32];
  sprintf(buf, "g_orbita_%d.dat", orbita_n);
  FILE* f = fopen(buf, "wt");
  
  fprintf(f, "## %g, %g, %g, $g\n", v[0], v[1], v[2], v[3]);
  
  double y[4] = {v[1], v[1], v[2], v[3]};
  const double t1 = v[0];
  const double step = t1 / 1000;
  
  gsl_odeiv2_system sys = {f_kristal, 0, 4, 0};
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-8, 1e-8);
  gsl_odeiv2_driver_set_hmin(driver, 1e-7);
  
  double t = 0;
  fprintf(f, "%g, %g, %g\n", t, y[0], y[1]);
  while (t < t1)
  {
    if (gsl_odeiv2_driver_apply(driver, &t, t+step, y) != GSL_SUCCESS)
    {
      fprintf(f, " == ERROR == \n");
      --orbita_n;
      break;
    }
    fprintf(f, "%g, %g, %g\n", t, y[0], y[1]);
  }
  
  printf("Resitev: %g, %g, %g, %g\n", t1, v[1], v[2], v[3]);
  
  fclose(f);
}

int konec(const gsl_vector* state, void* params, gsl_vector* values)
{
    m_system m = *(m_system*)params;    
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&m.system, gsl_odeiv2_step_rk4, 1e-6, 1e-12, 1e-12);
    
    double t = 0.0;
    int n = state->size + m.initial_state->size;
    double* v = new double[n];
    int i = 0;
    for (; i < m.initial_state->size; ++i)
    {
      v[i] = gsl_vector_get(m.initial_state, i);
    }
    for (; i < n; ++i)
    {
      v[i] = gsl_vector_get(state, i - m.initial_state->size);
    }
 //   std::cout << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;

    int status;
    const double t1 = 1.0;
    while (t < t1)
    {
      status = gsl_odeiv2_driver_apply(driver, &t, t1, v);
      if (status != GSL_SUCCESS)
      {
        std::cout << "Error:" << gsl_strerror(status) << std::endl;
        break;
      }
    }
  //  std::cout << status << GSL_SUCCESS << std::endl;
    if (status == GSL_SUCCESS)
    {
 //     std::cout << values->size << std::endl;
  //    std::cout << v[0] << " " << v[1] << std::endl;
      // gsl_vector_set_basis(values, m.final_state->size);
      for (i = 0; i < m.final_state->size; ++i)
      {
        gsl_vector_set(values, i, v[i] - gsl_vector_get(m.final_state, i));
      }
    }
    delete[] v;
    return status;
}

void vrv(double beta, double y0, double alpha, double F, const char* filename)
{
  gsl_vector* i = gsl_vector_alloc(2);
  gsl_vector_set(i, 0, 0);
  gsl_vector_set(i, 1, 0);
  gsl_vector* f = gsl_vector_alloc(2);
  gsl_vector_set(f, 0, 0);
  gsl_vector_set(f, 1, y0);
  gsl_odeiv2_system s = {f_vrv, 0, 4, (void*)&beta};
  m_system m = {s, i, f};
  
  gsl_vector* Fa = gsl_vector_alloc(2);
  gsl_vector_set(Fa, 0, F); // F
  gsl_vector_set(Fa, 1, alpha); // alpha
  
  gsl_multiroot_function function = {konec, 2, &m};
  gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 2);
  gsl_multiroot_fsolver_set(solver, &function, Fa);
  
  int iter = 0;
  int status;
  int fstatus;
  do
  {
      iter++;
      fstatus = gsl_multiroot_fsolver_iterate(solver);

      if (fstatus)   /* check if solver is stuck */
        break;

      status = gsl_multiroot_test_residual(solver->f, 1e-12);
   }
  while (status == GSL_CONTINUE && iter < 10000);
  
  FILE* file = fopen(filename, "wt");
  gsl_vector* x = gsl_multiroot_fsolver_root(solver);
  double y[] = {0, 0, gsl_vector_get(x, 0), gsl_vector_get(x, 1) };
  double t = 0.0;
  fprintf(file, "%g\t%g\t%g\n", t, y[0], y[1]);
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&s, gsl_odeiv2_step_rk4, 1e-6, 1e-8, 1e-8);
  while (t < 1)
  {
    gsl_odeiv2_driver_apply(driver, &t, t+0.01, y);
    fprintf(file, "%g\t%g\t%g\n", t, y[0], y[1]);
  }
  fclose(file);
  
  gsl_vector_fprintf(stdout, solver->x, "%f");
  
  gsl_multiroot_fsolver_free(solver);
}

bool cez_diagonalo(double x1, double y1, double x2, double y2)
{
  if (x1 < 0 || x2 < 0 || y1 < 0 || y2 < 0)
  {
    return false;
  }
  return (x1 > y1 && x2 < y2);
}

int orbita(const gsl_vector* state, void* params, gsl_vector* values)
{
    int n = state->size;
    double t1 = gsl_vector_get(state, 0);
    
    if (t1 < 0.2)
    {
    //  std::cout << "Cas je manjsi od 0.5" << std::endl;
      for(int j = 0; j < n; ++j)
      {
        gsl_vector_set(values, j, (1-t1));
      }
      return GSL_SUCCESS;
    }
    
    double test = gsl_vector_get(state, 2) * gsl_vector_get(state, 3);
    if (test > -1e-3)
    {
   //   std::cout << "Oba odvoda kazeta v isto smer " << gsl_vector_get(state, 2) << " " << gsl_vector_get(state, 3) << std::endl;
      for(int j = 0; j < n; ++j)
      {
        gsl_vector_set(values, j, 2 + test);
      }
      return GSL_SUCCESS;
    }
    
    gsl_odeiv2_system* system = (gsl_odeiv2_system*)params;
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(system, gsl_odeiv2_step_rk4, 1e-6, 1e-8, 1e-8);
    gsl_odeiv2_driver_set_hmin(driver, 1e-7);
    double t = 0.0;
    
    double* v = new double[n];
    int i = 1;
    for (; i < n; ++i)
    {
      v[i] = gsl_vector_get(state, i);
    }
    v[0] = gsl_vector_get(state, 1);
    
  //  std::cout << t1 << ";  " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;

    int status;
    
    status = gsl_odeiv2_driver_apply(driver, &t, t1, v);
  //  std::cout << t << " " << status << " " << GSL_SUCCESS << std::endl;
    if (status == GSL_SUCCESS || status == GSL_EDOM)
    {
  //    std::cout << state->size << " " << values->size << std::endl;
  //    std::cout << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;
      gsl_vector_set(values, 0, v[0] - gsl_vector_get(state, 1));
      for (i = 1; i < n; ++i)
      {
        gsl_vector_set(values, i, v[i] - gsl_vector_get(state, i));
      }
   //   gsl_vector_fprintf(stdout, values, "%f");
      status = GSL_SUCCESS;
    }
    delete[] v;
    gsl_odeiv2_driver_free(driver);
    return status;
}

void kristal(double x0, double v0, double u0, double t, const char* filename)
{
  gsl_odeiv2_system sys = {f_kristal, 0, 4, 0};
  gsl_vector* v = gsl_vector_alloc(4);
  gsl_vector_set(v, 0, t); 
  gsl_vector_set(v, 1, x0); 
  gsl_vector_set(v, 2, -u0); 
  gsl_vector_set(v, 3, v0); 
  
//  std::cout << "kristal zacetek " << x0 << std::endl;
  
  gsl_multiroot_function function = {orbita, 4, &sys};
  gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 4);
  gsl_multiroot_fsolver_set(solver, &function, v);
  
 // std::cout << "kristal nastavil sistem " << x0 << std::endl;
  
  int iter = 0;
  int status;
  int fstatus;
  
//  FILE* out = fopen(filename, "wt");

  do
  {
      ++iter;
    //  std::cout << "Starting iteration " << iter << std::endl << std::flush;
    //  std::cout << "  with u = " << gsl_vector_get(solver->x, 2) << ", v = " << gsl_vector_get(solver->x, 3) << std::endl;
      status = gsl_multiroot_fsolver_iterate(solver);
      
     // std::cout << std::endl << "Vector dx: " << std::endl;
     // gsl_vector_fprintf(stdout, solver->dx, "%f");
     // std::cout << std::endl;

      if (status)   /* check if solver is stuck */
        break;

      status = gsl_multiroot_test_residual(solver->f, 1e-6);
   //   fprintf(out, "%d %g %g %g %g\n", iter, gsl_vector_get(solver->x, 0), gsl_vector_get(solver->x, 1), gsl_vector_get(solver->f, 0), gsl_vector_get(solver->f, 1));
   }
  while (status == GSL_CONTINUE && iter < 1000);
  if (status == GSL_SUCCESS)
  {
    print_orbit(solver->x->data);
  }
  gsl_multiroot_fsolver_free(solver);
}

void vsi_kristali()
{
  printf("# x0 = y0, u0, v0, t, x, u, v, delta");
  orbita_n = 0;
  char buf[64];
  for( double x = 0.01; x < 0.45; x += 0.01)
  {
    FFOR(u)
    {
      FFOR(v)
      {
        sprintf(buf, "g_kristal_%g_%g_%g", x, u, v);
        kristal(x, u, v, 10*x, buf);
        std::cout << "Narejeno za kristal " << x << ", " << u << ", " << v << std::endl;
      }
    }
  }
  
  orbita_n = 1000;
  for( double x = 0.01; x < 0.45; x += 0.01)
  {
    FFOR(u)
    {
      FFOR(v)
      {
        sprintf(buf, "g_kristal_%g_%g_%g", x, u, v);
        kristal(x, u, v, 50*x, buf);
        std::cout << "Narejeno za kristal " << x << ", " << u << ", " << v << std::endl;
      }
    }
  }
}

void vse_vrvi()
{
  char buf[64];
  for (double y = 0.5; y <= 0.5; y += 0.1)
  {
    for (double beta = 0; beta < 10; beta += 1)
    {
      sprintf(buf, "g_vrv_%g_%g.dat", y, beta);
      vrv(beta, y, 0, 2, buf);
      sprintf(buf, "g_vrv_visenje_%g_%g.dat", y, beta);
      vrv(beta, y, -(M_PI_2-0.2), 0.5, buf);
      std::cout << "Narejeno za " << y << ", " << beta << std::endl << std::flush;
    }
  }
}

int main(int argc, char **argv) {
  
  /*
    vrv(1.0, "g_vrv_1.0.dat");
    vrv(2.0, "g_vrv_2.0.dat");
    vrv(10.0, "g_vrv_10.0.dat");
    vrv(0.1, "g_vrv_0.1.dat");
    */
  
   vse_vrvi();
   vsi_kristali();
  
  
  return 0;
}
