#include <iostream>
#include <iomanip>
#include <math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

int func(double t, const double y[], double f[], void* params)
{
    const double r = sqrt(y[0] * y[0] + y[1] * y[1]);
    const double rr = 1.0 / (r*r*r);
    f[0] = y[2];
    f[1] = y[3];
    f[2] = -y[0] * rr;
    f[3] = -y[1] * rr;
    
    // std::cout << "Evaluating " << y[0] << " " << y[1] << " " << y[2] << " " << y[3]  << std::endl;

    return GSL_SUCCESS;
}

double z(double t)
{
  return 2 * t - 10;
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

void zacetni(double y[], double phi, double v)
{
  y[0] = cos(phi);
  y[1] = sin(phi);
  y[2] = -v * sin(phi);
  y[3] = v * cos(phi);
}

void integrate()
{
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;
    const size_t dim = 4;

    gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, dim);
    gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (0.0, 1e-12);
    gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (dim);

    gsl_odeiv2_system sys = {zvezda, 0, dim, 0};

    double t = 0.0, t1 = 10.0;
    double h = 1e-6;
    double y[4];
    zacetni(y, 0.1, 0.8);

    while (t < t1)
    {
        int status = gsl_odeiv2_evolve_apply (e, c, s,
                                              &sys,
                                              &t, t1,
                                              &h, y);

        if (status != GSL_SUCCESS)
            break;

        std::cout << t << "\t" << y[0] << "\t" << y[1] << std::endl;
    }

    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);
}

void integriraj()
{
    gsl_odeiv2_system sys = {func, 0, 4, 0};
    gsl_odeiv2_driver * d =  gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk4, 1e-6, 1e-12, 0.0);

    int i;
    double t = 0.0;
    double t1 = 100.0;
    double y[4] = { 1.0, 0.0, 0.0, 1.0 };
    for (i = 1; i <= 100; i++)
    {
        double ti = i * t1 / 100.0;
        int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

        if (status != GSL_SUCCESS)
        {
            std::cout << "error, return value = " << status << std::endl;
            break;
        }

        std::cout << t << " " << y[0] << " " << y[1] << std::endl;
    }

    gsl_odeiv2_driver_free (d);
}

int main(int argc, char **argv) {
    integrate();
    return 0;
}
