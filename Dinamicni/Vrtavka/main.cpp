#include <iostream>
#include <fstream>
#include <math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#define l(x) y[x]
#define n(x) y[x+3]

using namespace std;

struct top_params
{
    double D;
    double a;
    double L;
    double mg;
};

inline double interpolate(double t1, double t2, double x1, double x2, double t)
{
  return ((t-t1)*x2 + (t2-t)*x1) / (t2-t1);
}

int odvod(double t, const double y[], double dy[], void* params)
{
    top_params& p = *(top_params*)params;
    
    const double R = 1.0-1.0/p.D;
    
    dy[0] = y[1] * y[2] * R + p.mg * p.a * y[4];
    dy[1] = -y[2] * y[0] * R + p.mg * (p.L * y[5] - p.a * y[3]);
    dy[2] = -p.mg * p.L * y[4];
    
    dy[3] = y[1] * y[5] - y[2] * y[4] / p.D;
    dy[4] = y[2] * y[3] / p.D - y[0] * y[5];
    dy[5] = y[0] * y[4] - y[1] * y[3];
    
    return GSL_SUCCESS;
}

class TopWorkspace
{
public:
    TopWorkspace(double D, double a, double L);
    ~TopWorkspace();
    void poincare();
    inline void apply()
    {
        for (int i = 0; i < 6; ++i)
        {
          last[i] = y[i];
        }
        gsl_odeiv2_driver_apply(driver, &t, t+t1, y);
    }
    
private:
    gsl_odeiv2_system system;
    gsl_odeiv2_driver* driver;
    top_params params;
    double y[6];
    double t;
    double t1;
    double last[6];
    ofstream output;
};

TopWorkspace::TopWorkspace(double D, double a, double L)
{
    params = {D, a, L, 1.0};
    system = {odvod, 0, 6, &params};
    driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rk4, 1e-6, 1e-8, 1e-8);
    
    t1 = 1e-4;
    
    output.open("g_vrtavka.dat");
    
    // TODO: Zacetni pogoji
    y[0] = 0.0;
    y[1] = 0.0;
    y[2] = 1.0;
    
    y[3] = 0.9;
    y[4] = 0.0;
    y[5] = sqrt(1.0 - y[3] * y[3] - y[4] * y[4]);
    cout << y[5] << endl;
}

TopWorkspace::~TopWorkspace()
{
    gsl_odeiv2_driver_free(driver);
    output.close();
}

void TopWorkspace::poincare()
{
    last[3] = y[3];
    cout << "Zacnem en korak preslikave" << endl;
    cout << t << " " << y[3] << " " << y[4] << " " << y[5] << endl;

    while (last[3] * y[3] > 0 || y[3] > last[3])
    {
        apply();
    }
    
    double tS = t - t1 * fabs(y[3]) / (fabs(y[3]) + fabs(last[3]));
    
    for (int i = 0; i < 6; ++i)
    {
        output << interpolate(t-t1, t, last[i], y[i], tS) << " ";
    }
    output << endl;
}

int main(int argc, char **argv) {
    TopWorkspace w(2.0, 1.0, 0.1);
    for (int i = 0; i < 100; ++i)
    {
      w.poincare();
    }
    return 0;
}
