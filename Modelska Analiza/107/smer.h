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

typedef smer* (*gen_s)(void);

smer* smeri_gsl();

smer* smeri_gsl()
{
  smer* s = malloc(N * sizeof(smer));
  int j;
  double x,y,z;
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(0));
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
  gsl_rng_set(r, time(0));
  for (j=0; j<N; ++j)
  {
    s[j].phi = (2 * gsl_rng_uniform(r) - 1) * M_PI;
    s[j].theta = acos ( 2*gsl_rng_uniform_pos(r) -1 );
  }
  gsl_rng_free(r);
  return s;
}

smer* dipol_jaz()
{
  smer* s = malloc(N*sizeof(smer));
  int j;  
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(0));
  double x, y, z, g;
  for (j=0; j<N; ++j)
  {
    s[j].phi = (2 * gsl_rng_uniform(r) - 1) * M_PI;
    int roots = gsl_poly_solve_cubic(0, -3, 2.0*(2*gsl_rng_uniform_pos(r)-1), &x, &y, &z);
    if (roots < 3)
    {
      printf("WARNING: not enough real roots\n");
      g = x;
    }
    else
    {
      g = y;
    }
    s[j].theta = acos ( g );
  }
  gsl_rng_free(r);
  return s;
}


typedef double (*s_test)(smer*);

double moment(smer* s, int l, int m)
{
  double avg = 0;
  int j;
  for (j=0; j<N; ++j)
  {
    double R = gsl_sf_legendre_sphPlm(l,abs(m),cos(s[j].theta));
    if ( m < 0 )
    {
      avg += R * sin ( m*s[j].phi );
    }
    else
    {
      avg += R * cos ( m*s[j].phi );
    }
  }
  return avg / N;
}

double avgphi(smer* s)
{
  double avg = 0;
  int j;
  for (j=0; j<N; ++j)
  {
    avg += s[j].phi;
  }
  
  return avg / N;
}

double avgc2t(smer* s)
{
  double avg = 0;
  int j;
  for (j=0; j<N; ++j)
  {
    double u = cos(s[j].theta);
    avg += u*u;
  }
  
  return avg / N;
}

void test_smeri()
{
  smer* data;
  double t;
  gen_s generators[] = { smeri_gsl, smeri_jaz, dipol_jaz, 0 };
  int j;
  for (j=0; generators[j]; ++j)
  {
    data = generators[j]();
    printf("%d, l=1: %g\t%g\t%g\n", j, moment(data, 1, -1), moment(data, 1, 0), moment(data, 1, 1));
    printf("%d, l=2: %g\n", j, moment(data, 2, 0));
    printf("%d, <phi>: %g\n", j, avgphi(data));
    printf("%d, <c2t>: %g\n", j, avgc2t(data));
  }
}

