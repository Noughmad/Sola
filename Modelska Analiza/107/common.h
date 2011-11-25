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

#define N 1e6
#define B 1e2
#define B2 30
#define AVG (N/B)
#define AVG2 (N/B2/B2)

#define D_UNIFORM 1
#define D_VECTOR 2
#define D_GAUSS 3

double* create();
int* zeros();
double normalized(double x);

typedef double* (*gen)(void);
typedef double (*test)(gsl_histogram*,gsl_histogram*);

typedef struct _task
{
  int id;
  gen* generators;
  test* tests;
  gsl_histogram* density;
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
}

void task_perform(task* t)
{
  int i,j;
  for (i=0; t->generators[i]; ++i)
  {
    double* v = t->generators[i]();
    
    printf("%d-corr: %lf\n", t->id, gsl_stats_lag1_autocorrelation(v, 0, N));
    
    gsl_histogram* h = gsl_histogram_alloc(B);
    gsl_histogram_set_ranges_uniform(h, t->min, t->max);
    for (j=0; j<N; ++j)
    {
      gsl_histogram_increment(h, v[j]);
    }
        
    for (j=0; t->tests[j]; ++j)
    {
      printf("%d-%d-%d: %lf\n", t->id, i, j, t->tests[j](h, t->density));
    }
    free(v);
  }
}


#endif // VREDNOSTI_H