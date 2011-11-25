#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

#include "common.h"

double* gauss()
{
  gsl_rng* r = gsl_rng_alloc(gsl_rng_mt19937);
  double* d = create();
  int i;
  for (i=0; i<N; ++i)
  {
    d[i] = gsl_ran_ugaussian(r);
  }
  gsl_rng_free(r);
  return d;
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

