
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "common.h"
#include "gauss.h"
#include "smer.h"
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_statistics_double.h>


double constant(double x)
{
  return AVG;
}

double constant_c(double x)
{
  return AVG*x/B;
}

double chisq(gsl_histogram* values, gsl_histogram* expected);
double chisq2d(gsl_histogram* values, gsl_histogram* expected);
double kolsmir(gsl_histogram* values, gsl_histogram* expected);
double korelacija(gsl_histogram* value, gsl_histogram* expected);

double int2double(int i);

double* builtin();
double* kalkulator();
double* from_file(const char* name);
double* devrandom();
double* devurandom();
double* mersenne();


double chisq ( gsl_histogram* values, gsl_histogram* expected )
{
  gsl_histogram*t  = gsl_histogram_clone(values);
  gsl_histogram_sub(t, expected);
  int j;
  double chi = 0;
  for (j=0; j<B; ++j)
  {
    chi += pow(gsl_histogram_get(t, j), 2)/gsl_histogram_get(expected, j);
  }
  
  return gsl_cdf_chisq_Q(chi, B-1);
}

double kolsmir ( gsl_histogram* values, gsl_histogram* expected )
{  
  double D = 0;
  double cv = 0, ce = 0;
  int j;
  for (j=0; j<B; ++j)
  {
    cv += gsl_histogram_get(values, j);
    ce += gsl_histogram_get(expected, j);
    D = fmax(D, fabs(cv-ce));
  }
  D /= B;
  return D;
}

int* zeros()
{
  int* z = malloc(B * sizeof(int));
  int j;
  for (j = 0; j < B; ++j)
  {
    z[j] = 0;
  }
  return z;
}

double* create()
{
  return (double*)malloc(N * sizeof(double));
}

double* builtin()
{
  double* d = create();
  srand(time(0));
  int i;
  for(i=0; i<N; ++i)
  {
    d[i] = int2double(rand());
  }
  return d;
}

double* kalkulator()
{
  double x = time(0) / RAND_MAX;
  double* d = create();
  int i;
  for(i=0; i<N; ++i)
  {
    x = normalized( pow(x + M_PI,5) );
    d[i] = x;
  }
  return d;
}

double int2double(int i)
{
  double m = RAND_MAX;
  double M = 1.0 / (m + 1);
  return i * M;
}

double* devurandom()
{
  return from_file("/dev/urandom");
}

double* from_file(const char* name)
{
  FILE *f = fopen(name, "rt");
  double* d = create();
  int i;
  int t;
  for(i=0; i<N; ++i)
  {
    fread(&t, sizeof(int), 1, f);
    d[i] = int2double(abs(t));
  }
  fclose(f);
  return d;
}

double* mersenne()
{
  gsl_rng* g = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(g, time(0));
  double* d = create();
  int i;
  for (i=0; i<N; ++i)
  {
    d[i] = gsl_rng_uniform(g);
  }
  gsl_rng_free(g);
  return d;
}

double normalized(double x)
{
  return x - floor(x);
}

int main(int argc, char** argv)
{
  gen g_uni[] = {mersenne, kalkulator, devurandom, builtin, 0};
  test t[] = {chisq, kolsmir, 0};
  task* prob = task_new(0, g_uni, t, 0 );
  
  gsl_histogram* h = gsl_histogram_alloc(B);
  gsl_histogram_set_ranges_uniform(h, 0, 1);
  gsl_histogram_shift(h, AVG);
  prob->density = h;
  
  prob->min = 0;
  prob->max = 1;
  task_perform(prob);
  task_free(prob);
  
  gen g_gauss[] = {gauss, 0};
  prob = task_new(2, g_gauss, t, 0);
  prob->min = -5;
  prob->max = 5;
  
  h = gsl_histogram_alloc(B);
  gsl_histogram_set_ranges_uniform(h, -5, 5);
  int j;
  for (j=0; j<B; ++j)
  {
    double min, max;
    gsl_histogram_get_range(h, j, &min, &max);
    double x = (max + min)/2;
    
    gsl_histogram_accumulate(h, x, N*gsl_ran_ugaussian_pdf(x)/100*(prob->max-prob->min));
  }
  prob->density = h;
  task_perform(prob);
  task_free(prob);

  return 0;
}