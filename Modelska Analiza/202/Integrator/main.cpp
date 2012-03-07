#include <iostream>
#include <iomanip>
#include <math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "lib/integrator.h"

#include <QtCore/QDebug>

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

double energija2(double t, const State& state)
{
  Q_ASSERT(state.size() == 4);
  State s(4);
  s[0] = state[0] - z(t);
  s[1] = state[1] + 1.5;
  s[2] = state[2] - 2;
  s[3] = state[3];
  return energy(s);
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
  const double dy = y[1] + 1.5; 
  
  const double R = sqrt(dx * dx + dy * dy);
  const double RR = 1.0 / (R*R*R);
  
  f[2] -= dx * RR;
  f[3] -= dy * RR;
  
  return GSL_SUCCESS;
}

State zacetni(double phi, double v)
{
  State y(4);
  phi = phi - 10 - M_PI_2;
  y[0] = cos(phi);
  y[1] = sin(phi);
  y[2] = -v * sin(phi);
  y[3] = v * cos(phi);
  
  return y;
}


int main(int argc, char **argv) {
    Integrator* integrator = new GslIntegrator();
    Interval i = qMakePair(0.0, 50.0);
    
    /*
    double eps[] = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16 };
    State s(4);
    s[0] = 1.0;
    s[1] = 0.0;
    s[2] = 0.0;
    s[3] = 0.1;
    for (int j = 0; j < 11; ++j)
    {
      Solution sol = integrator->integrate(func, s, i, eps[j]);
      qDebug() << eps[j] << " " << odstopanje(sol, energy) << " " << odstopanje(sol, momentum)
	<< " " << obhodni_cas(sol) << " " << povratek(sol) << " " << sol.size();
     // saveToFile(QString("test-%1").arg(j+6), sol, false);
    }
    
    return 0;
    
    qDebug() << "Naredil primer brez zvezde";
    
    */
    
    for (double phi = 0; phi < 2*M_PI; phi += 0.02)
    {
      const State state = (integrator->integrate(zvezda, zacetni(phi, 1.0), i, 1e-14).constEnd()-1).value();
      qDebug() << phi << " " << energy(state) << " "<< energija2(i.second, state);
      
      // Testing only
      // break;
    }
    
    return 0;
}
