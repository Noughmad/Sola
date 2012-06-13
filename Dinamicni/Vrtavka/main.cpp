#include <iostream>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

struct top_params
{
    double D;
    double a;
    double L;
    double mg;
};

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
        last = y[3];
        gsl_odeiv2_driver_apply(driver, &t, t1, y);
    }
    
private:
    gsl_odeiv2_system system;
    gsl_odeiv2_driver* driver;
    top_params params;
    double y[];
    double t;
    double t1;
    double last;
};

TopWorkspace::TopWorkspace(double D, double a, double L)
{
    params = {D, a, L, 1.0};
    system = {odvod, 0, 6, &params};
    driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rk4, 1e-6, 1e-8, 1e-8);
    
    t1 = 1e-4;
    
    // TODO: Zacetni pogoji
}

TopWorkspace::~TopWorkspace()
{
    gsl_odeiv2_driver_free(driver);
}

void TopWorkspace::poincare()
{
    last = y[3];
    while (last * y[3] > 0)
    {
        apply();
    }
    apply();
    last = y[3];
    while (last * y[3] > 0)
    {
        apply();;
    }    
}

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    return 0;
}
