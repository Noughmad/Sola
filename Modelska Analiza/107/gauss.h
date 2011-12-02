#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "common.h"

typedef double (*gauss_rng_f)(const gsl_rng*, const double);

double* gauss(gauss_rng_f f)
{
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, time(0));
  double* d = create();
  int i;
  for (i=0; i<N; ++i)
  {
    d[i] = f(r, 1);
  }
  gsl_rng_free(r);
  return d;
}

double* gauss_bm()
{
  return gauss(gsl_ran_gaussian);
}

double* gauss_zig()
{
  return gauss(gsl_ran_gaussian_ziggurat);
}

double* gauss_ratio()
{
  return gauss(gsl_ran_gaussian_ratio_method);
}

void gauss_time(gauss_rng_f f)
{
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, time(0));
  double t;
  int i;
  for (i=0; i<N; ++i)
  {
    t = f(r, 1);
  }
  gsl_rng_free(r);
}


int* bin_gauss(double* v)
{
  int* i = zeros();
  int j;
  for (j=0; j<N; ++j)
  {
    if (v[j] > 10)
    {
      v[j] = 10;
    }
    else if (v[j] < -10)
    {
      v[j] = -10;
    }
    int b = v[j] * B / 20 + B/2;
    i[b]++;
  }
  return i;
}

