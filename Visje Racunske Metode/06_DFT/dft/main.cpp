#include <iostream>
#include <fstream>

#include <gsl/gsl_roots.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

using namespace std;

const int N = 2000;
const double Rmax = 20.0;
const double h = Rmax / N;

double r(int i)
{
  return h * i + 1e-12;
}

class DFT
{
public:
  DFT(int n);
  ~DFT();

  double getEnergy();
  void getElectronDensity();
  void getPotential();
  
  void printElectronDensity(const char* file);
  void printPotential(const char* file);

private:
  gsl_root_fsolver *rootSolver;
  gsl_odeiv2_driver *poissonSolver;

public:
  const int n;
  double* u;
  double* U;
};

double numerov_end(double E, void* params);
double numerov_step(double up, double uc, double k2p, double k2c, double k2n, double h2);
int xdot_poisson(double t, const double y[], double dydt[], void* params);

double numerov_end(double E, void* params)
{
  DFT* dft = static_cast<DFT*>(params);

  double a, b;
  a = 0;
  b = 1;

  for (int i = 2; i < N; ++i)
  {
    a = numerov_step(a, b, 2 * (E - dft->U[i-2]), 2 * (E - dft->U[i-1]), 2 * (E - dft->U[i]), h*h);
    std::swap(a, b);
  }

  return b;
}

double numerov_step(double up, double uc, double k2p, double k2c, double k2n, double h2)
{
  return (2 * (1 - 5 * h2 * k2c / 12) * uc - (1.0 + h2 * k2p / 12) * up) / (1.0 + h2 * k2n / 12.0);
}

DFT::DFT(int n)
: n(n)
, rootSolver(gsl_root_fsolver_alloc(gsl_root_fsolver_bisection))
{
  u = new double[N];
  U = new double[N];
}

DFT::~DFT()
{
  gsl_root_fsolver_free(rootSolver);
  delete[] u;
  delete[] U;
}

double DFT::getEnergy()
{
  gsl_function f = {&numerov_end, this};
  gsl_root_fsolver_set(rootSolver, &f, -2, -0);

  int status, x_lo, x_hi;
  int iter = 0;
  
  
    x_lo = gsl_root_fsolver_x_lower(rootSolver);
    x_hi = gsl_root_fsolver_x_upper(rootSolver);
      cout << x_hi << ", " << x_lo << endl;

  double r, lr;
  lr = -1;
  
  do
  {
    ++iter;
    status = gsl_root_fsolver_iterate(rootSolver);
    cout << iter << ", " << status << endl;
    
    r = gsl_root_fsolver_root(rootSolver);
    status = gsl_root_test_delta(lr, r, 1e-12, 0);
    lr = r;
  }
  while (status == GSL_CONTINUE);
  
  cout << "Found energy after " << iter << " with status = " << status << endl;
  cout << x_hi << ", " << x_lo << endl;

  return gsl_root_fsolver_root(rootSolver);
}

void DFT::getElectronDensity()
{
  double E = getEnergy();
  cout << "Energy is " << E << endl;
  
  double sum = 0;

  u[0] = 0;
  u[1] = 1;

  sum += u[1] * u[1];

  for (int i = 2; i < N; ++i)
  {
    u[i] = numerov_step(u[i-2], u[i-1], 2*(E-U[i-2]), 2*(E-U[i-1]), 2*(E-U[i]), h*h);
    sum += u[i] * u[i];
  }

  sum = 1.0 / sqrt(sum);
  for (int i = 0; i < N; ++i)
  {
    u[i] *= sum;
  }
}

void DFT::getPotential()
{
  double Uprime = 1;
  U[0] = 0;

  for (int i = 1; i < N; ++i)
  {
    U[i] = U[i-1] + Uprime * h;
    Uprime += (- u[i] * u[i]) / r(i);
  }

  double adjust = (1 - U[N-1]) / (N-1);
  for (int i = 1; i < N; ++i)
  {
    U[i] += adjust * i;
  }

  // Transform U to V = N*U/r
  for (int i = 1; i < N; ++i)
  {
    U[i] *= n / r(i);
  }
}

void DFT::printElectronDensity(const char* file)
{
  ofstream out(file);

  for (int i = 0; i < N; ++i)
  {
    out << r(i) << ", " << u[i] << endl;
  }

  out.close();
}

void DFT::printPotential(const char* file)
{
  ofstream out(file);

  for (int i = 0; i < N; ++i)
  {
    out << r(i) << ", " << U[i] << endl;
  }

  out.close();
}

void test_vodik()
{
  DFT dft(1);

  for (int i = 0; i < N; ++i)
  {
    dft.U[i] = -1.0 / r(i);
  }

  dft.getElectronDensity();
  dft.printElectronDensity("g_test_vodik.dat");

  dft.getPotential();
  dft.printPotential("g_test_vodik_potencial.dat");
}

int main(int argc, char **argv) {
    test_vodik();
    return 0;
}
