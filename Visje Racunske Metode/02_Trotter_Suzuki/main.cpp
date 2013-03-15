#include <vector>
#include <iostream>

#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

using namespace std;

const double lambda = 0.5;

typedef vector<double> shema;

class X
{
public:
    double energija();
    double korak_T(double t);
    double korak_V(double t);
    double korak(double t, shema s);
    
    double q1;
    double q2;
    double p1;
    double p2;
};

double X::energija()
{
    return (p1*p1 + p2*p2 + q1*q1 + q2*q2) / 2 + lambda * q1*q1*q2*q2;
}

double X::korak_T(double t)
{
    q1 += p1 * t;
    q2 += p2 * t;
}

double X::korak_V(double t)
{
    p1 -= q1 * (1 + 2 * lambda * q2 * q2) * t;
    p2 -= q2 * (1 + 2 * lambda * q1 * q1) * t;
}

double X::korak(double t, shema s)
{
    const size_t n = s.size();
    for (size_t i = 0; i < n; ++i)
    {
        if (i % 2)
        {
            korak_T(t * s[i]);
        }
        else
        {
            korak_V(t * s[i]);
        }
    }
}

shema s2()
{
    shema s;
    s.push_back(0.5);
    s.push_back(1.0);
    s.push_back(0.5);
    
    return s;
}

shema s4()
{
    const double x0 = -pow(2, 1.0/3.0) / (2 - pow(2, 1.0/3.0));
    const double x1 = 1 / (2 - pow(2, 1.0/3.0));

    shema s;
    s.push_back(x1/2);
    s.push_back(x1);
    s.push_back((x1+x0)/2);
    s.push_back(x0);
    s.push_back((x1+x0)/2);
    s.push_back(x1);
    s.push_back(x1/2);
    
    return s;
}

double energija(const double y[])
{
    return (y[0]*y[0] + y[1]*y[1] + y[2]*y[2] + y[3]*y[3]) / 2 + lambda * y[0]*y[0]*y[1]*y[1];
}

int odvod(double t, const double y[], double dydt[], void* param)
{
    dydt[0] = y[2];
    dydt[1] = y[3];
    
    dydt[2] = -y[0] * (1 + 2 * lambda * y[1] * y[1]);
    dydt[3] = -y[1] * (1 + 2 * lambda * y[0] * y[0]);
    
    return GSL_SUCCESS;
}

void vse_metode(double T, int N, int n)
{
    const double p0 = 0.5;
    X x2 = {1, 0, 0, p0};
    X x4 = {1, 0, 0, p0};
        
    double y[4] = {1, 0, 0, p0};
    double t = 0;
    double step = T / (N * n);
    
    shema sm2 = s2();
    shema sm4 = s4();
    
    gsl_odeiv2_system sys = {odvod, 0, 4, 0};
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, T/N, 1e-3, 1e-3);
    
    for (int i = 0; i < N; ++i)
    {
        gsl_odeiv2_driver_apply_fixed_step(driver, &t, step, n, y);
        
        for (int j = 0; j < n; ++j)
        {
            x2.korak(step, sm2);
            x4.korak(step, sm4);
        }
        
        cout << t << " " 
             << y[0] << " " << y[1] << " " << energija(y) << " " 
             << x2.q1 << " " << x2.q2 << " " << x2.energija() << " " 
             << x4.q1 << " " << x4.q2 << " " << x4.energija() << endl;
    }
}

void ekviparticija(double T, double N)
{
    X x = {1, 0, 0, 2};
    shema s = s4();
    
    double pp1 = 0;
    double pp2 = 0;
    
    const int Interval = 1;
    
    double step = T / N;
    for (int i = 0; i < N; ++i)
    {
        x.korak(step, s);
        pp1 = i*pp1/(i+1) + x.p1 * x.p1 / (i+1);
        pp2 = i*pp2/(i+1) + x.p2 * x.p2 / (i+1);
        
        if ((i % Interval) == 0)
        {
            cout << i * step << " " << pp1 << " " << pp2 << endl;
        }
    }
}

int main(int argc, char** argv)
{
    // vse_metode(atof(argv[1]), atoi(argv[2]), atoi(argv[3]));
    ekviparticija(atof(argv[1]), atoi(argv[2]));
    return 0;
}