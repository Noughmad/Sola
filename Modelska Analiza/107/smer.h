#include "common.h"
#include "math.h"
#include <complex.h>

#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_poly.h>

typedef struct _smer
{
  double theta;
  double phi;
} smer;

smer*  smeri_gsl();

smer* smeri_gsl()
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
  gsl_rng_free(r);
  return s;
}

smer* smeri_jaz()
{
  smer* s = malloc(N*sizeof(smer));
  int j;  
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  for (j=0; j<N; ++j)
  {
    s[j].phi = 2 * M_PI * gsl_rng_uniform(r);
    s[j].theta = acos ( 2*gsl_rng_uniform(r) -1 );
  }
  gsl_rng_free(r);
  return s;
}

smer* dipol_jaz()
{
  smer* s = malloc(N*sizeof(smer));
  int j;  
  double A, x0;
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  double x, y, z, g;
  for (j=0; j<N; ++j)
  {
    s[j].phi = 2 * M_PI * gsl_rng_uniform(r);
    int roots = gsl_poly_solve_cubic(0, -3, 3*A*(gsl_rng_uniform(r)-x0), &x, &y, &z);
    if (roots == 1)
    {
      g = x;
    }
    else
    {
      printf("More than one real root\n");
      g = x;
    }
    s[j].theta = acos ( x );
  }
  gsl_rng_free(r);
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

