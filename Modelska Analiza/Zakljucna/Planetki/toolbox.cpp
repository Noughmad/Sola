#include "planets.h"

#include "KEP/src/lambert_problem.h"

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

using namespace kep_toolbox;

array3D vectorToArray(const Vector& v)
{
    array3D a;
    a[0] = v.x;
    a[1] = v.y;
    a[2] = 0;
    return a;
}

int odvod(double t, const double y[], double dydt[], void* params)
{
    dydt[0] = y[2];
    dydt[1] = y[3];

    const double r = sqrt(y[0] * y[0] + y[1] * y[1]);
    const double rr = 1.0 / (r*r*r);
    dydt[2] = -y[0] * rr;
    dydt[3] = -y[1] * rr;

    return GSL_SUCCESS;
}

Trajectory Planets::toolbox(int steps, int circles)
{
    Trajectory tr(steps);
    Vector r1 = planet_one(0);
    lambert_problem p(vectorToArray(r1), vectorToArray(planet_two(steps*dt)), steps*dt);

    tr.positions[0] = r1;
    tr.positions[steps] = planet_two(steps*dt);

    int n = p.get_Nmax();
    int m = 0;
    double v = std::numeric_limits<double>::max();
    for (int i = 0; i < 2*n+1; ++i)
    {
        array3D one = p.get_v1()[i];
        array3D two = p.get_v2()[i];
        const double test = sqrt(one[0]*one[0]+one[1]*one[1]) + sqrt(two[0]*two[0] + two[1]*two[1]);
        if (test < v)
        {
            v = test;
            m = i;
        }
    }

    array3D v1 = p.get_v1()[m];
    double y[4] = {r1.x, r1.y, v1[0], v1[1]};
    gsl_odeiv2_system sys = {odvod, 0, 4, 0};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-3, 1e-8, 1e-8);

    double t = 0;
    int i = 1;
    while (i < steps)
    {
        gsl_odeiv2_driver_apply(driver, &t, t+dt, y);
        tr.positions[i].x = y[0];
        tr.positions[i].y = y[1];
        ++i;
    }

    return tr;
}
