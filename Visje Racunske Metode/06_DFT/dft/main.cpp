#include <iostream>

#include <gsl/gsl_roots.h>
#include <gsl/gsl_odeiv2.h>

const int N = 200;
const double Rmax = 5.0;
const double h = Rmax / N;

class DFT
{
public:
  DFT();
  ~DFT();

private:
  gsl_root_fsolver *rootSolver;
  gsl_odeiv2_driver *poissonSolver;

public:
  double u[];
  double U[];
};

double numerov_end(double E, void* params);
double numerov_step(double up, double uc, double k2p, double k2c, double k2n, double h2);

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
  return (2 * (1 + 5 * h2 * k2c / 12) * uc - (1.0 + h2 * k2p / 12) * up) / (1.0 + h2 * k2n / 12.0);
}

DFT::DFT()
: rootSolver(gsl_root_fsolver_alloc(gsl_root_fsolver_brent))
{
  gsl_function f = {&numerov_end, this};
  gsl_root_fsolver_set(rootSolver, &f, -2, 0);
}

DFT::~DFT()
{
  gsl_root_fsolver_free(rootSolver);
  gsl_odeiv2_driver_free(poissonSolver);
}


int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    return 0;
}
