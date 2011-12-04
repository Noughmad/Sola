#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <math.h>

using namespace std;

struct Polozaj 
{
  Polozaj(gsl_rng* rng) 
  {
    r = pow(gsl_rng_uniform_pos(rng), 1.0/3.0);
    phi = gsl_rng_uniform(rng) * 2 * M_PI;
    theta = acos( 2*gsl_rng_uniform_pos(rng) - 1 );
  }
  
  double r;
  double phi;
  double theta;
};

struct Foton
{
  Foton(gsl_rng* rng)
  {
    r = pow(gsl_rng_uniform_pos(rng), 1.0/3.0);
    u = 2*gsl_rng_uniform_pos(rng) - 1;
  }
  
  double r;
  double u;
};

bool je_znotraj(Polozaj p)
{
  return p.r * sin(p.theta) > cos(p.phi);
}

void telo()
{
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(0));
  int Z = 1e5;
  int N = 0;
  double X = 0;
  for (int i = 0; i < Z; ++i) 
  {
    Polozaj p(r);
    if (je_znotraj(p))
    {
      ++N;
      X += p.r * cos(p.phi);
    }
  }
  
  gsl_rng_free(r);
  
  double p = 1.0 * N / Z;
  double sigma = sqrt(p*(1-p)/Z);
  
  printf("Enakomerna gostota: m = %g +/- %g\n", p, sigma);
  printf("Enakomerna gostota: x* = %g\n", X/N);
}

double rand_exp(gsl_rng* r)
{
  return -log(1-gsl_rng_uniform(r));
}

void sevanje()
{
  int B = 100;
  
  gsl_histogram* h = gsl_histogram_alloc(B);
  gsl_histogram_set_ranges_uniform(h, 0, 2);
  
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(0));
  int Z = 1e7;
  double l;
  
  for (int i = 0; i < Z; ++i)
  {
    Foton f(r);
    l = sqrt(1 - f.r * f.r + f.r * f.r * f.u * f.u) - f.r * f.u;
    gsl_histogram_increment(h, l);
  }
  
  FILE* hf = fopen("sevanje_histogram.dat", "wt");
  gsl_histogram_fprintf(hf, h, "%g", "%g");
  fclose(hf);
  
  gsl_histogram_pdf* pdf = gsl_histogram_pdf_alloc(B);
  gsl_histogram_pdf_init(pdf, h);
  
  int m = 100;
  double mu[100];
  int N[100];
  
  for (int j = 0; j < m; ++j)
  {
    mu[j] = 5.0 * j / m;
    N[j] = 0;
  }
  
  for (int i = 0; i < Z; ++i)
  {
    double l = gsl_histogram_pdf_sample(pdf, gsl_rng_uniform(r));
    double y = rand_exp(r) / l;
    int k = fmin(y * m / 5.0 + 1, m);
    
    for (int j = 0; j < k; ++j)
    {
      ++N[j];
    }
  }
  
  FILE* f = fopen("sevanje2.dat", "wt");
  
  for (int j = 0; j < m; ++j)
  {
    fprintf(f, "%g\t%g\n", mu[j], 1.0 * N[j] / Z);
  }
  
  fclose(f);

}

void reflektor()
{
  const int P = 100; // stevilo pregrad
  const int M = 1000; // najvecje stevilo odbojev
  const int D = 1e7; // stevilo delcev
  
  double pregrade[P];
  gsl_histogram* histogrami[P];
  for (int i = 0; i < P; ++i) 
  {
    pregrade[i] = 0.2 * i;
    histogrami[i] = gsl_histogram_alloc(M);
    gsl_histogram_set_ranges_uniform(histogrami[i], 0, M);
  }
  
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(0));
  for (int i = 0; i < D; ++i)
  {
    int dosegel = 0;
    double x = rand_exp(r); // Zacetni premik do prvega odboja

    for (int j = 0; j < M; ++j)
    {
      const int trenutna = floor(5 * x);
      for (int k = dosegel; k < trenutna && k < P; ++k)
      {
	gsl_histogram_increment(histogrami[k], j);
      }
      dosegel = max(dosegel, trenutna);
      if (trenutna >= P)
      {
	// Delec je zletel ven
	break;
      }
      if (x < 0)
      {
	break;
      }
      
      int smer = (gsl_rng_uniform(r) < 0.5) ? 1 : -1;
      x += smer * rand_exp(r);
    }
  }
  
  for (int i = 0; i < P; ++i)
  {
    char name[64];
    sprintf(name, "reflektor_%d.dat", i);
    FILE* f = fopen(name, "wt");
    gsl_histogram_fprintf(f, histogrami[i], "%g", "%g");
    fclose(f);
  }
}

int main(int argc, char **argv) 
{
  // telo();
  // sevanje();
  reflektor();
  return 0;
}
