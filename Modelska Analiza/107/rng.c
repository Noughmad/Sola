
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "common.h"
#include "gauss.h"
#include "smer.h"
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort_double.h>

double constant(double x)
{
  return AVG;
}

double constant_c(double x)
{
  return AVG*x/B;
}

double chisq(gsl_histogram* values, gsl_histogram* expected);
double kolsmir(gsl_histogram* values, gsl_histogram* expected);
double korelacija(gsl_histogram* value, gsl_histogram* expected);

double int2double(int i);

double* builtin();
double* kalkulator();
double* from_file(const char* name);
double* devrandom();
double* devurandom();
double* mersenne();

double kolmogorov_konst(double* values)
{
  double m = 0;
  gsl_sort(values, 1, N);
  int i;
  for (i = 0; i < N; ++i)
  {
    m = fmax( m, fabs( i*1.0/N - values[i] ) );
  }
  return m * sqrt(N);
}

double kolmogorov_gauss(double* values)
{
  double m = 0, cv=0, ce=0;
  gsl_sort(values, 1, N);
  int i;
  for (i = 0; i < N; ++i)
  {
    m = fmax( m, fabs( i*1.0/N - 0.5 - 0.5*gsl_sf_erf(values[i] / sqrt(2.0)) ) );
  }
  return m * sqrt(N);
}

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
  return D / sqrt(1.0 * B);
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

void builtin_time()
{
  srand(time(0));
  double t;
  int i;
  for(i=0; i<N; ++i)
  {
    t = int2double(rand());
  }
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

void kalkulator_time()
{
  double x = time(0) / RAND_MAX;
  int i;
  for (i=0; i<N; ++i)
  {
    x = normalized( pow(x+M_PI, 5) );
  }
}

double int2double(int i)
{
  double m = RAND_MAX;
  double M = 1.0 / (m + 1);
  return i * M;
}

void devurandom_time()
{
  FILE *f = fopen("/dev/urandom", "rb");
  int i;
  int t;
  for(i=0; i<N; ++i)
  {
    fread(&t, sizeof(int), 1, f);
  }
  fclose(f);
}

double* devurandom()
{
  return from_file("/dev/urandom");
}

double* devrandom()
{
  return from_file("/dev/random");
}

double* randorg()
{
  return from_file("/home/miha/Prejeto/RandomNumbers");
}

double* from_file(const char* name)
{
  FILE *f = fopen(name, "rb");
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

void mersenne_time()
{
  gsl_rng* g = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(g, time(0));
  double t;
  int i;
  for (i=0; i<N; ++i)
  {
    t = gsl_rng_uniform(g);
  }
  gsl_rng_free(g);
}

double normalized(double x)
{
  return x - floor(x);
}

void test_kolmogorov( gen* generators_k, gen* generators_g )
{
  int i;
  for (i=0; generators_k[i]; ++i)
  {
    double* d = generators_k[i]();
    printf("konstantni %d: %lf\n", i, kolmogorov_konst(d));
    free(d);
  }
  for (i=0; generators_g[i]; ++i)
  {
    double* d = generators_g[i]();
    printf("Gaussov %d: %lf\n", i, kolmogorov_gauss(d));
    free(d);
  }
}

void test_normal()
{
  gen g_uni[] = {builtin, devurandom, kalkulator, mersenne, 
#if N <= (1 << 12)
  randorg,
#endif
  0};
  test t[] = {chisq, kolsmir, 0};
  task* prob = task_new(0, g_uni, t, 0 );
  
  gsl_histogram* h = gsl_histogram_alloc(B);
  gsl_histogram_set_ranges_uniform(h, 0, 1);
  gsl_histogram_shift(h, AVG);
  prob->density = h;
  
  gsl_histogram2d* hd = gsl_histogram2d_alloc(B2, B2);
  gsl_histogram2d_set_ranges_uniform(hd, 0, 1, 0, 1);
  gsl_histogram2d_shift(hd, AVG2/2);
  prob->density2d = hd;
  
  prob->min = 0;
  prob->max = 1;
  task_perform(prob);
  task_free(prob);
  
  gen g_gauss[] = {gauss_bm, gauss_zig, gauss_ratio, 0};
  prob = task_new(2, g_gauss, t, 0);
  prob->max = 3;
  prob->min = -prob->max;
  
  h = gsl_histogram_alloc(B);
  gsl_histogram_set_ranges_uniform(h, prob->min, prob->max);
  int j;
  for (j=0; j<B; ++j)
  {
    double min, max;
    gsl_histogram_get_range(h, j, &min, &max);
    double x = (max + min)/2;
    
    gsl_histogram_accumulate(h, x, N*gsl_ran_ugaussian_pdf(x) / B * (prob->max - prob->min));
  }
  
  hd = gsl_histogram2d_alloc(B2, B2);
  gsl_histogram2d_set_ranges_uniform(hd, prob->min, prob->max, prob->min, prob->max);
  int k;
  for (j=0; j < B2; ++j)
    for (k=0; k<B2; ++k)
    {
      double min, max;
      gsl_histogram2d_get_xrange(hd, j, &min, &max);
      double x = (max + min)/2;
      
      gsl_histogram2d_get_yrange(hd, k, &min, &max);
      double y = (max + min)/2;
      
      gsl_histogram2d_accumulate(hd, x, y, gsl_ran_ugaussian_pdf(x) * gsl_ran_ugaussian_pdf(y));
    }
    
  gsl_histogram2d_scale( hd, N/2/gsl_histogram2d_sum(hd) );
  prob->density = h;
  prob->density2d = hd;
  task_perform(prob);
  task_free(prob);
}

int main(int argc, char** argv)
{
 // test_normal();
 // test_smeri();
 // mersenne_time();
 // devurandom_time();
 // builtin_time();
 // kalkulator_time();
 // gauss_time(gsl_ran_gaussian_ratio_method);
 gen gk[] = {builtin, devurandom, kalkulator, mersenne, 
#if N <= (1 << 12)
  randorg,
#endif
0};

  gen gg[] = { gauss_bm, gauss_zig, gauss_ratio, 0};
  test_kolmogorov(gk, gg);
  return 0;
}