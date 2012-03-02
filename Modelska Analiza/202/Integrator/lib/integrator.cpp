#include "integrator.h"

#include <iostream>
#include <iomanip>
#include <math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

State createState(double x, double y, double u, double v)
{
  State s(4);
  s[0] = x;
  s[1] = y;
  s[2] = u;
  s[3] = v;
  return s;
}


Integrator::~Integrator()
{

}

Solution GslIntegrator::integrate(YDot dot, const State& initialState, const Interval& interval)
{
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;
  const size_t dim = initialState.size();

  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, dim);
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (0.0, 1e-14);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (dim);

  gsl_odeiv2_system sys = {dot, 0, dim, 0};

  double t = 0.0, t1 = 30.0;
  double h = 1e-6;
  double* y = (double*)malloc(sizeof(double) * initialState.size());
  for (int i = 0; i < initialState.size(); ++i)
  {
    y[i] = initialState[i];
  }
  
  Solution solution;
  solution[t] = initialState;

  while (t < t1)
  {
      int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);

      if (status != GSL_SUCCESS)
	  break;
      
      State s(dim);
      for (int i = 0; i < dim; ++i)
      {
	s[i] = y[i];
      }
      
      solution[t] = s;
  }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  delete[] y;
  
  return solution;
}

double radius(const State& state)
{
  return sqrt(state[0] * state[0] + state[1] * state[1]);
}

double energy(const State& state)
{
  // TODO: Normalization
  // both parts should come with different constant factors
  return 0.5*(state[2] * state[2] + state[3] * state[3]) - 1.0/radius(state);
}

double momentum(const State& state)
{
  // TODO: Calculate
  return 0;
}



