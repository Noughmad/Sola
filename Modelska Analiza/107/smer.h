#include "common.h"
#include "math.h"
#include <complex.h>

#include <gsl/gsl_sf_legendre.h>

typedef struct _smer
{
  double theta;
  double phi;
} smer;

smer*  smeri();

smer* smeri()
{
  smer* s = malloc(N * sizeof(smer));
  int j;
  double x,y,z;
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  for (j=0; j<N; ++j)
  {
    gsl_ran_dir_3d(r, &x, &y, &z);
    s[j].phi = atan2(y,x);
    s[j].theta = acos(z);
  }
  return s;
}

typedef double (*s_test)(smer*);

double moment(smer* s, int l, int m)
{
  complex avg = 0;
  int j;
  for (j=0; j<N; ++j)
  {
    avg += gsl_sf_legendre_sphPlm(l,m,cos(s[j].theta)) * exp( I*m*s[j].phi );
  }
  return avg / N;
}

