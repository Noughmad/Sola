#include <iostream>
#include <cassert>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>

typedef int (*odvod)(double, const double[], double[], void*);

typedef struct {
   gsl_odeiv2_system system;
   gsl_vector* initial_state; // zacetna (x,y)
   gsl_vector* final_state; // koncna (x,y)
} m_system;

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
    const double D = 1.0 / (sin(M_PI*y[0]) * sin(M_PI*y[0]) + sin(M_PI*y[0]) * sin(M_PI*y[0]));
    f[0] = y[2];
    f[1] = y[3];
    
    f[2] = M_PI_2 * sin(M_2_PI * y[0]) * D;
    f[2] = M_PI_2 * sin(M_2_PI * y[1]) * D;
    
    return GSL_SUCCESS;
}

int konec(const gsl_vector* state, void* params, gsl_vector* values)
{
    m_system m = *(m_system*)params;    
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&m.system, gsl_odeiv2_step_rk4, 1e-6, 1e-8, 1e-8);
    
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
    std::cout << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << std::endl;

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
    std::cout << status << GSL_SUCCESS << std::endl;
    if (status == GSL_SUCCESS)
    {
      std::cout << values->size << std::endl;
      std::cout << v[0] << " " << v[1] << std::endl;
      // gsl_vector_set_basis(values, m.final_state->size);
      for (i = 0; i < m.final_state->size; ++i)
      {
        gsl_vector_set(values, i, v[i] - gsl_vector_get(m.final_state, i));
      }
    }
    delete[] v;
    return status;
}

void vrv(double beta, const char* filename)
{
  gsl_vector* i = gsl_vector_alloc(2);
  gsl_vector_set(i, 0, 0);
  gsl_vector_set(i, 1, 0);
  gsl_vector* f = gsl_vector_alloc(2);
  gsl_vector_set(f, 0, 0);
  gsl_vector_set(f, 1, -0.5);
  gsl_odeiv2_system s = {f_vrv, 0, 4, (void*)&beta};
  m_system m = {s, i, f};
  
  gsl_vector* Fa = gsl_vector_alloc(2);
  gsl_vector_set(Fa, 0, 0.4); // F
  gsl_vector_set(Fa, 1, 0.3); // alpha
  
  gsl_multiroot_function function = {konec, 2, &m};
  gsl_multiroot_fsolver* solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrid, 2);
  gsl_multiroot_fsolver_set(solver, &function, Fa);
  
  int iter = 0;
  int status;
  int fstatus;
  
  FILE* out = fopen(filename, "wt");

  do
  {
      iter++;
      fstatus = gsl_multiroot_fsolver_iterate(solver);

      if (fstatus)   /* check if solver is stuck */
        break;

      status = gsl_multiroot_test_residual(solver->f, 1e-7);
      fprintf(out, "%d %g %g %g %g\n", iter, gsl_vector_get(solver->x, 0), gsl_vector_get(solver->x, 1), gsl_vector_get(solver->f, 0), gsl_vector_get(solver->f, 1));
   }
  while (status == GSL_CONTINUE && iter < 1000);
  fclose(out);
}

void kristal()
{
  
}

int main(int argc, char **argv) {
    vrv(1.0, "g_vrv_1.0.dat");
    vrv(2.0, "g_vrv_2.0.dat");
    vrv(10.0, "g_vrv_10.0.dat");
    vrv(0.1, "g_vrv_0.1.dat");
    return 0;
}
