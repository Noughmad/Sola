#include <iostream>

#include <gsl/gsl_rng.h>
#include <math.h>

#define NV 15
#define N_SAMPLES 1e3
#define SAMPLE_INTERVAL 1e2

#define NS 256
#define NS_AND (NS-1)

const double m = 1;
const double k = 0.1;

double T = 60;

struct Veriga
{
  Veriga();
  double clenki[NV];
};

struct Spini
{
  bool s[NS][NS];
};

double J = 1;
double H = 0;

double energija(Spini spini)
{
  double e = 0;
  for (int i = 0; i < NS; ++i)
  {
    for (int j = 0; j < NS; ++j)
    {
      e += (spini.s[i][j] == spini.s[(i+1) & NS_AND][j]) ? -J : J;
      e += (spini.s[i][j] == spini.s[i][(j+1) & NS_AND]) ? -J : J;
      e += spini.s[i,j] ? -H : H;
    }
  }
  return e;
}

Spini korak(Spini spini, gsl_rng* rng)
{
  int i = gsl_rng_uniform(rng) * NS;
  int j = gsl_rng_uniform(rng) * NS;
  
  Spini t = spini;
  t.s[i][j] = !t.s[i][j];
  
  double dE = energija(t) - energija(spini);
  if ( dE < 0 || exp(- dE/T) > gsl_rng_uniform(rng) )
  {
    return t;
  }
  else
  {
    return spini;
  }
}

Veriga::Veriga()
{
  for(int i=0; i<NV; ++i)
  {
    clenki[i] = 0;
  }
}

double energija(Veriga v);
Veriga korak(Veriga& v, double xi);

double energija(Veriga v)
{
  double e = 0;
  for (int i=0; i<NV; ++i)
  {
    e -= m * v.clenki[i];
  }
  for (int i=0; i<NV-1; ++i)
  {
    e += k * (v.clenki[i] - v.clenki[i+1]) * (v.clenki[i] - v.clenki[i+1]);
  }
  return e;
}

Veriga korak(Veriga& v, gsl_rng* r)
{
  int n = 5*pow(gsl_rng_uniform_pos(r),5); // Stevilo clenkov, ki jih premaknemo naenkrat
  
  int k = gsl_rng_get(r) % (NV-n);
  int s = (gsl_rng_get(r) % 2) ? 1 : -1;
  
  Veriga t = v;
  for (int i = k; i < k+n; ++i)
  {
    int tmp = t.clenki[i] + s;
    t.clenki[i] = std::min(std::max(tmp, 0), 19);
  }
  double dE = energija(t) - energija(v);
  if ( dE < 0 || exp(-T * dE) > gsl_rng_uniform(r) )
  {
    return t;
  }
  else
  {
    return v;
  }
}

template <class Sistem>
void postopek()
{
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(0));
  double avgE = 0;
  
  Sistem veriga;
  for (int i = 0; i < 100; ++i)
  {
    // Na zacetku, pred samplanjem, naredimo nekaj ponovitev
    veriga = korak(veriga, r);
  }
  
  for (int i=0; i < N_SAMPLES; ++i)
  {
    for (int j = 0; j < 100; ++j)
    {
      // Na zacetku, pred samplanjem, naredimo nekaj ponovitev
      veriga = korak(veriga, r);
    }
  avgE += energija(veriga); 
  }
  
  avgE /= N_SAMPLES;

  std::cout << "Povprecna energija pri T=" << T << " je E=" << avgE << std::endl;
}

int main(int argc, char **argv) {
   // verizica();
    postopek<Veriga>();
    return 0;
}
