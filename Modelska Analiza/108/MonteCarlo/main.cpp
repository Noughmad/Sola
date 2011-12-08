#include <iostream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram.h>
#include <math.h>

const double V = 4.0*M_PI/3;


using namespace std;

struct Polozaj 
{
  Polozaj(gsl_rng* rng) 
  {
    r = pow(gsl_rng_uniform_pos(rng), 1.0/3.0);
    phi = gsl_rng_uniform(rng) * 2 * M_PI;
    theta = acos( 2*gsl_rng_uniform_pos(rng) - 1 );
  }
  
  double gostota();
  
  double r;
  double phi;
  double theta;
  static const char* name;
  static double mass;
};


struct Polozaj_3
{
  Polozaj_3(gsl_rng* rng) 
  {
    r = pow(gsl_rng_uniform(rng), 1.0/6.0);
    phi = gsl_rng_uniform(rng) * 2 * M_PI;
    theta = acos( 2*gsl_rng_uniform_pos(rng) - 1 );
  }
  
  double gostota();
  
  double r;
  double phi;
  double theta;
  static const char* name;
  static double mass;
};

const char* Polozaj::name = "Enakomerna";
const char* Polozaj_3::name = "Spremenljiva";

double Polozaj::mass = V;
double Polozaj_3::mass = V/2;


double Polozaj::gostota()
{
  return 1;
}

double Polozaj_3::gostota()
{
  return r * r * r;
}


struct Foton
{
  Foton(gsl_rng* rng)
  {
    r = pow(gsl_rng_uniform(rng), 1.0/3.0);
    u = 2*gsl_rng_uniform_pos(rng) - 1;
  }
  
  double r;
  double u;
};

template <class P>
bool je_znotraj(P p)
{
  return p.r * sin(p.theta) > -cos(p.phi);
}

template <class P>
void telo()
{
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(0));
  int Z = 1e8;
  int N = 0;
  double X = 0;
  double sumsq = 0;
  for (int i = 0; i < Z; ++i) 
  {
    P p(r);
    if (je_znotraj(p))
    {
      ++N;
      const double v = p.r * cos(p.phi);
      X += v;
      sumsq += v * v;
    }
  }
  
  gsl_rng_free(r);
  
  double p = 1.0 * N / Z;
  double sigma = sqrt(p*(1-p)/Z);
  double m = p * P::mass;
  
  double xa = X/N;
  double var = sumsq / N - xa * xa;
  
  printf("%s gostota: m = %g +/- %g\n", P::name, m, sigma*P::mass);
  printf("%s gostota: x* = %g +/ %g\n", P::name, xa, sqrt(var/N));
}

double rand_exp(gsl_rng* r)
{
  return -log(1-gsl_rng_uniform(r));
}

void sevanje()
{
  const int B = 100;
  
  gsl_histogram* h = gsl_histogram_alloc(B);
  gsl_histogram_set_ranges_uniform(h, 0, 2);
  
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(0));
  const int Z = 1e8;
  double l;
  
  for (int i = 0; i < Z; ++i)
  {
    Foton f(r);
    l = sqrt(1 - f.r * f.r + f.r * f.r * f.u * f.u) - f.r * f.u;
    gsl_histogram_increment(h, l);
  }
  
  FILE* hf = fopen("sevanje_histogram.dat", "wt");
  gsl_histogram_scale(h, 0.5*B/Z);
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

void reflektor(bool izotropno)
{
  const int P = 100; // stevilo pregrad
  const int M = 100; // najvecje stevilo odbojev
  const int D = 1e7; // stevilo delcev
  
  double pregrade[P];
  gsl_histogram* histogrami[P];
  gsl_histogram* kot_o[P];
  gsl_histogram* kot_p[P];
  for (int i = 0; i < P; ++i) 
  {
    pregrade[i] = 0.2 * i;
    histogrami[i] = gsl_histogram_alloc(M);
    gsl_histogram_set_ranges_uniform(histogrami[i], -0.5, M-0.5);
    
    if (izotropno)
    {
      kot_o[i] = gsl_histogram_alloc(100);
      gsl_histogram_set_ranges_uniform(kot_o[i], -1, 0);
      kot_p[i] = gsl_histogram_alloc(100);
      gsl_histogram_set_ranges_uniform(kot_p[i], 0, 1);
    }
  }
  
  gsl_histogram* prepustnost = gsl_histogram_alloc(P);
  gsl_histogram_set_ranges_uniform(prepustnost, -0.1, 0.2 * P-0.1);
  
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(0));
  for (int i = 0; i < D; ++i)
  {
    int dosegel = 0;
    double smer = 1;

    double x = rand_exp(r); // Zacetni premik do prvega odboja

    for (int j = 0; true; ++j)
    {
      const int trenutna = floor(5 * x);
      for (int k = dosegel; k < trenutna && k < P; ++k)
      {
	if (izotropno)
	  gsl_histogram_increment(kot_p[k], smer);
	gsl_histogram_increment(histogrami[k], j);
      }
      dosegel = max(dosegel, min(trenutna, P) );
      if (trenutna >= P)
      {
	// Delec je zletel ven
	break;
      }
      if (x < 0)
      {
	if (izotropno)
	for (int k = dosegel; k < P; ++k)
	{
	  gsl_histogram_increment(kot_o[k], smer);
	}
	break;
      }
      
      if (izotropno)
      {
	smer = 2 * gsl_rng_uniform_pos(r) - 1;
      } 
      else 
      {
	smer = (gsl_rng_uniform(r) < 0.5) ? 1 : -1;
      }
      x += smer * rand_exp(r);
    }
    
    for (int k = 0; k < dosegel; ++k )
    {
      gsl_histogram_increment(prepustnost, k * 0.2);
    }
  }
  
  gsl_histogram_scale(prepustnost, 1.0/D);
  
  char rname[64];
  sprintf(rname, "odbojnost_%s.dat", izotropno ? "izo" : "nn");
  FILE* o = fopen(rname, "wt");
  
  char mean_name[64];
  sprintf(mean_name, "mean_%s.dat", izotropno ? "izo" : "nn");
  FILE* mf = fopen(mean_name, "wt");
  
  gsl_histogram_fprintf(o, prepustnost, "%g", "%g");
  fclose(o);
  
  char mname[64];
  sprintf(mname, "matrix_%s.dat", izotropno ? "izo" : "nn");
  FILE* mat = fopen(mname, "wt");
  
  for (int i = 0; i < P; ++i)
  {
    char name[64];
    sprintf(name, "reflektor_%s_%d.dat", izotropno ? "izo" : "nn", i);
    FILE* f = fopen(name, "wt");
    gsl_histogram_scale(histogrami[i], 1.0/D);
    gsl_histogram_fprintf(f, histogrami[i], "%g", "%g");
    fclose(f);
    
    if (izotropno)
    {
      sprintf(name, "kot_o_%d.dat", i);
      f = fopen(name, "wt");
      gsl_histogram_scale(kot_o[i], 1.0/gsl_histogram_sum(kot_o[i]));
      gsl_histogram_fprintf(f, kot_o[i], "%g", "%g");
      fclose(f);
      
      sprintf(name, "kot_p_%d.dat", i);
      f = fopen(name, "wt");
      gsl_histogram_scale(kot_p[i], 1.0/gsl_histogram_sum(kot_p[i]));
      gsl_histogram_fprintf(f, kot_p[i], "%g", "%g");
      fclose(f);
    }
    
    for (int j = 0; j < M; ++j)
    {
      fprintf(mat, "%g\t", gsl_histogram_get(histogrami[i], j));
    }
    fprintf(mat, "\n");
    
    fprintf(mf, "%g\t%g\n", pregrade[i], gsl_histogram_mean(histogrami[i]));
  }
  
  fclose(mat);
  fclose(mf);
}

int main(int argc, char **argv) 
{
  /*
 telo<Polozaj>();
 telo<Polozaj_3>();
 sevanje();
 */
  reflektor(false);
  reflektor(true);
  return 0;
}
