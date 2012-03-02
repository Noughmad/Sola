#include <iostream>
#include <iomanip>
#include <math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "lib/integrator.h"

int func(double t, const double y[], double f[], void* params)
{
    const double r = sqrt(y[0] * y[0] + y[1] * y[1]);
    const double rr = 1.0 / (r*r*r);
    f[0] = y[2];
    f[1] = y[3];
    f[2] = -y[0] * rr;
    f[3] = -y[1] * rr;

    return GSL_SUCCESS;
}

double z(double t)
{
  return 2 * t - 20.0;
}

int zvezda(double t, const double y[], double f[], void* params)
{
  const double r = sqrt(y[0] * y[0] + y[1] * y[1]);
  const double rr = 1.0 / (r*r*r);
  f[0] = y[2];
  f[1] = y[3];
  f[2] = -y[0] * rr;
  f[3] = -y[1] * rr;
  const double dx = y[0] - z(t);
  const double dy = y[1] - 1.5; 
  
  const double R = sqrt(dx * dx + dy * dy);
  const double RR = 1.0 / (R*R*R);
  
  f[2] -= dx * RR;
  f[3] -= dy * RR;
  
  return GSL_SUCCESS;
}

State zacetni(double phi, double v)
{
  State y(4);
  phi = phi - 10;
  y[0] = cos(phi);
  y[1] = sin(phi);
  y[2] = -v * sin(phi);
  y[3] = v * cos(phi);
  
  return y;
}

int main(int argc, char **argv) {
    Integrator* integrator = new GslIntegrator();
    Interval i = qMakePair(0.0, 30.0);
    Solution sol = integrator->integrate(func, zacetni(0.0, 1.0), i);
    Solution::ConstIterator it = sol.constBegin();
    Solution::ConstIterator end = sol.constEnd();
    for(; it != end; ++it)
    {
      std::cout << it.key() << "\t" << it.value()[0] << "\t" << it.value()[1] << std::endl;
    }
    return 0;
}
