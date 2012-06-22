#include <iostream>
#include <fstream>
#include <math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#define l(x) y[x-1]
#define n(x) y[x+2]

using namespace std;


struct top_params
{
    double D;
    double a;
    double L;
    double mg;
};

inline double energy(const top_params& top, double y[])
{
    return 0.5 * (y[0]*y[0] + y[1]*y[1] + top.D*y[2]*y[2]) + top.mg * (top.L * n(1) + top.a * n(2));
}

void random_state(double y[], const top_params& top, double emax);

double lyapunov(const top_params& top, const double initial[]);

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

    void setInitial(double l1, double l2, double l3, double n1, double n2, double n3);
    void setInitial(const double initial[]);

    double operator-(const TopWorkspace& other);

    void poincare();
    inline void apply()
    {
        for (int i = 0; i < 6; ++i)
        {
          last[i] = y[i];
        }
        gsl_odeiv2_driver_apply(driver, &t, t+t1, y);
    }

    inline double energy()
    {
	return ::energy(params, y);
    }

public:
    double y[6];
    double t;

private:
    gsl_odeiv2_system system;
    gsl_odeiv2_driver* driver;
    top_params params;
    double t1;
    double last[6];
    ofstream output;
};

TopWorkspace::TopWorkspace(double D, double a, double L)
{
    params = {D, a, L, 1.0};
    system = {odvod, 0, 6, &params};
    driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rk4, 1e-4, 1e-4, 1e-4);

    t1 = 1e-3;

    output.open("g_vrtavka.dat");
}

TopWorkspace::~TopWorkspace()
{
    gsl_odeiv2_driver_free(driver);
    output.close();
}

void TopWorkspace::setInitial(double l1, double l2, double l3, double n1, double n2, double n3)
{
  l(1) = l1;
  l(2) = l2;
  l(3) = l3;

  n(1) = n1;
  n(2) = n2;
  n(3) = n3;
}

void TopWorkspace::setInitial(const double initial[])
{
  for (int i = 0; i < 6; ++i)
  {
    y[i] = initial[i];
  }
}

double TopWorkspace::operator-(const TopWorkspace& other)
{
  double d = 0;
  for (int i = 0; i < 6; ++i)
  {
      d += (y[i] - other.y[i]) * (y[i] - other.y[i]);
  }
  return d;
}


void TopWorkspace::poincare()
{
    last[3] = y[3];
    apply();
    apply();
    double tt = t;

    while (last[3] * y[3] > 0 || y[3] > last[3])
    {
        apply();
    }

    double tS = t - t1 * fabs(y[3]) / (fabs(y[3]) + fabs(last[3]));
    
 //   cout << "Poincare: perioda je " << tS - tt << endl; 

    /*
    output << t;
    for (int i = 0; i < 6; ++i)
    {
        output << " " << interpolate(t-t1, t, last[i], y[i], tS);
    }
    output << endl;
    */
    for (int i = 0; i < 6; ++i)
    {
        y[i] = interpolate(t-t1, t, last[i], y[i], tS);
    }
    t = tS;
}

double lyapunov(const top_params& top, const double initial[])
{
  double factor = 1;
  const double step = 1e-12;
  const double limit = 1e-9;

  const int N = 100;

  double z[6];
  for (int i = 0; i < 6; ++i)
  {
    z[i] = initial[i];
  }

  TopWorkspace top1(top.D, top.a, top.L);
  top1.setInitial(z);

  for (int i = 0; i < 6; ++i)
  {
    z[i] += (double)rand() / RAND_MAX * step;
  }

  TopWorkspace top2(top.D, top.a, top.L);
  top2.setInitial(z);
  double d0 = top1 - top2;

  for (int i = 0; i < N; ++i)
  {
      double d = top1 - top2;
      
      if (d > limit)
      {
          for (int i = 0; i < 6; ++i)
          {
              top2.y[i] = (1.0-step) * top1.y[i] + step * top2.y[i];
          }
          factor /= step;
      }
      top1.poincare();
      top2.poincare();
  }
  
  double d_mid = (top1 - top2) * factor;
  factor = 1;
  
  for (int i = 0; i < N; ++i)
  {
      double d = top1 - top2;
      
      if (d > limit)
      {
          cout << "Manjsam razliko" << endl;
          for (int i = 0; i < 6; ++i)
          {
              top2.y[i] = (1.0-step) * top1.y[i] + step * top2.y[i];
          }
          factor /= step;
      }
      top1.poincare();
      top2.poincare();
  }

  double d = (top1 - top2) * factor;
  double l1 = log(d_mid / d0) / N;
  double l2 = log(d / d_mid) / N;
  return l2;
}

void random_state(double y[], const top_params& top, double emax)
{
    do 
    {
	for (int i = 1; i < 4; ++i)
	{
	    l(i) = 10 * (2 * (double)rand() / RAND_MAX - 1);
	    n(i) = (double)rand() / RAND_MAX;
	}
    }
    while (energy(top, y) > emax);
}

double chaos_part(const top_params& top, double emax)
{
    int N = 100;
    double y[6];
    double lambda_limit = 1e-3;
    int C = 0;
    
    for (int i = 0; i < N; ++i)
    {
	random_state(y, top, emax);
	double l = lyapunov(top, y);
	cout << "Eksponent je " << l << endl;
	if (l > lambda_limit)
	{
	    C++;
	}
    }
    
    return (double)C/N;
}

int main(int argc, char **argv) {
    int seed = time(0);
    srand(seed);
    
    top_params params = {2.0, 1.0, 0.5, 1.0};
    cout << "Delez kaosa je " << chaos_part(params, 10) << endl;
    top_params params2 = {2.0, 1.0, 0.0, 1.0};
    cout << "Delez kaosa je " << chaos_part(params2, 10) << endl;
    
    
    return 0;
}
