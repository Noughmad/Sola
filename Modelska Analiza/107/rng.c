
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_cdf.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#define N 1e7
#define B 1e2
#define AVG (N/B)
#define AVG2 (N/B/B)

int* zeros();
double* create();
double normalized(double x);

double chi1(int n, double* values);
double chi2(int n, double* values);
double ks(int n, double* values);

int* bin(int n, double* values);
double int2double(int i);

double* builtin();
double* kalkulator();
double* devrandom();
double* devurandom();
double* mersenne();

int* bin(int n, double* values)
{
  int* a = zeros();
  int b, j = 0;
  for (j; j < N; ++j)
  {
    b = values[j] * B;
    a[b]++;
  }
  free(values);
  return a;
}

double chi1(int n, double* values)
{
  int* b = bin(n, values);
  double chi = 0;
  int i = 0;
  for (i; i < B; ++i)
  {
    chi += pow( (b[i] - AVG)/sqrt(AVG), 2 );
  }
  free(b);
  return gsl_cdf_chisq_Q(chi, B-1);
}

double chi2(int n, double* values)
{
  int** b = malloc(B * sizeof(int*));
  int i;
  int k,l;
  for (i=0; i<B; ++i)
  {
    b[i] = zeros();
  }
  for (i=0; i<n; i+=2)
  {
    k = values[i] * B;
    l = values[i+1] * B;
    b[k][l]++;
  }
  
  double chi;
  
  for(k=0;k<B;++i)
  {
    for(l=0;l<B;++l)
    {
      chi += pow( (b[k][l] - AVG2) / sqrt(AVG2), 2 );
    }
  }
  
  return gsl_cdf_chisq_Q(chi, B*B-1);
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
  FILE *f = fopen("/dev/urandom", "rb");
  double* d = create();
  int i;
  double t;
  for(i=0; i<N; ++i)
  {
    fscanf(f, "%lf", &t);
    d[i] = normalized(t);
  }
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
  printf("Vgrajen: %lf\n", chi1(N, builtin()));
  printf("Kalkulator: %lf\n", chi1(N, kalkulator()));
  printf("Mersenne: %lf\n", chi1(N, mersenne()));
 // printf("Urandom: %lf\n", chi1(N, devurandom()));
  return 0;
}