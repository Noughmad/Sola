#ifndef VREDNOSTI_H
#define VREDNOSTI_H
#include <stdlib.h>
#include <stdio.h>

#include <gsl/gsl_histogram.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <math.h>

#define N (1 << 24)
#define B (1 << 12)
#define B2 (1 << 8)
#define AVG (N/B)
#define AVG2 (N/B2/B2)

#define D_UNIFORM 1
#define D_VECTOR 2
#define D_GAUSS 3

double* create();
int* zeros();
double normalized(double x);

double chisq2d(gsl_histogram2d* values, gsl_histogram2d* expected)
{
  double chi;
  int i,j;
  for (i=0; i<B2; ++i)
    for (j=0; j<B2; ++j)
    {
      chi += pow(gsl_histogram2d_get(values, i, j) - gsl_histogram2d_get(expected, i,j), 2) / gsl_histogram2d_get(expected, i, j);
    }
  return gsl_cdf_chisq_Q(chi, B2*B2-1);
}


typedef double* (*gen)(void);
typedef double (*test)(gsl_histogram*,gsl_histogram*);

typedef struct _task
{
  int id;
  gen* generators;
  test* tests;
  gsl_histogram* density;
  gsl_histogram2d* density2d;
  double min, max;
} task;

task* task_new(int id, gen* generators, test* tests, gsl_histogram* dens)
{
  task* t = malloc(sizeof(task));
  t->id = id;
  t->generators = generators;
  int i;

  t->tests = tests;
  t->density = dens;
  
  return t;
}

void task_free(task* t)
{
  free(t);
  gsl_histogram_free(t->density);
  gsl_histogram2d_free(t->density2d);
}


double kol_p(double x)
{
  double p = 1, t;
  int j;
  for (j=1; j<100; ++j)
  {
    t = 2*exp(-2*j*j*x*x);
    if (j%2)
      p -= t;
    else
      p += t;
  }
  return p;
}

void task_perform(task* t)
{
  int i,j;
  for (i=0; t->generators[i]; ++i)
  {
    double* v = t->generators[i]();
    double r;
    
    printf("%d-corr: %lf\n", t->id, gsl_stats_lag1_autocorrelation(v, 0, N));
    
    gsl_histogram* h = gsl_histogram_alloc(B);
    gsl_histogram_set_ranges_uniform(h, t->min, t->max);
    
    gsl_histogram2d* hd = gsl_histogram2d_alloc(B2, B2);
    gsl_histogram2d_set_ranges_uniform(hd, t->min, t->max, t->min, t->max);
    for (j=0; j<N; ++j)
    {
      gsl_histogram_increment(h, v[j]);
      if ( j % 2 )
	gsl_histogram2d_increment(hd, r, v[j]);
      else
	r = v[j];
    }
        
    for (j=0; t->tests[j]; ++j)
    {
      printf("%d-%d-%d: %lf\n", t->id, i, j, t->tests[j](h, t->density));
    }
    
    printf("%d-%d-chi2: %lf\n", t->id, i, chisq2d(hd, t->density2d ));
    free(v);
    
    gsl_histogram_free(h);
    gsl_histogram2d_free(hd);
  }
}


#endif // VREDNOSTI_H