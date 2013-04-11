#include <iostream>
#include <fstream>

#include <gsl/gsl_rng.h>
#include <math.h>

#define NV 15
#define N_SAMPLES 1e3
 
#define NS 64
#define NS_AND (NS-1)

#define N NS*NS

#define SAMPLE_INTERVAL 1e3
#define SAMPLE_WAIT 5e3

// #define UREJENO
#define COOLING

using namespace std;

const double m = 1;
double k = 0.1;

double T = 10;

struct Veriga
{
  Veriga();
  double clenki[NV];
};

struct Spini
{ 
  Spini()
  {
#ifndef UREJENO
    srand(time(0));
    int r;
    int t = 0;

#endif
    
    for (int i = 0; i < NS; ++i)
      for (int j = 0; j < NS; ++j)
      {
#ifndef UREJENO
	if (t == 0)
	{
	  r = rand();
	}
	s[i][j] = (r & ( 1 << t )) ? 1 : -1;
	t = (t+1) & 31;
#else
	s[i][j] = -1;
#endif
      }
  }
  
  
  short s[NS][NS];
  
  inline int sp(int i, int j) const
  {
    return s[i & NS_AND][j & NS_AND] ? 1 : -1;
  }
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
      e -= spini.sp(i,j) * ( J * ( spini.sp(i+1,j) + spini.sp(i,j+1) ) + H );
    }
  }
  return e;
}

double energ_sosedi(Spini spini, int i, int j)
{
  return -spini.sp(i,j) * ( J * ( spini.sp(i+1,j) + spini.sp(i,j+1) + spini.sp(i-1,j) + spini.sp(i,j-1) ) + H );
}

Spini korak(Spini spini, gsl_rng* rng)
{
  int i = gsl_rng_uniform(rng) * NS;
  int j = gsl_rng_uniform(rng) * NS;
  
  double e_prej = 0, e_pol = 0;
  
  Spini t = spini;
  t.s[i][j] = !t.s[i][j];
  
  double dE = 2 * energ_sosedi(t,i,j);
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

double korelacija(Spini s1, Spini s2)
{
  double kor = 0;
  for (int i = 0; i < NS; ++i)
    for (int j = 0; j < NS; ++j)
    {
      kor += (s1.s[i][j] == s2.s[i][j]);
    }
  return 2 * kor / NS / NS - 1;
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
  if ( dE < 0 || exp(-dE / T) > gsl_rng_uniform(r) )
  {
    return t;
  }
  else
  {
    return v;
  }
}

template <class Sistem>
Sistem postopek(const Sistem& s, ostream& stream)
{
  gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
  gsl_rng_set(r, time(0));
  double avgE = 0;
  double avgE2 = 0;
  
  Sistem veriga = s;
      
  for (int i=0; i < N_SAMPLES; ++i)
  {
    for (int j = 0; j < SAMPLE_INTERVAL; ++j)
    {
      // Na zacetku, pred samplanjem, naredimo nekaj ponovitev
      veriga = korak(veriga, r);
    }
    double e = energija(veriga);
  avgE += e;
  avgE2 += e*e;
  }
  
  avgE /= N_SAMPLES;
  avgE2 /= N_SAMPLES;
  
  stream << T << " " << avgE << " " << (avgE2 - avgE * avgE) / N / T / T << endl;
  
  return veriga;
}

double magnetizacija(Spini s)
{
  double m = 0;
  for (int i = 0; i < NS; ++i)
    for (int j = 0; j < NS; ++j)
      m += s.sp(i,j);
  return m;
}

void histereza()
{
  ofstream o;
  o.open("../../sp-meritve-gor.dat");
  double temp[] = {1, 2, 3, 5, 10};
  Spini s;
  for (int i = 0; i < 30; ++i)
  {
      T = 1.5 + i * 0.05;
      s = postopek<Spini>(s, o);
  }
  o.close();
  o.open("../../sp-meritve-dol.dat");

  for (int i = 0; i < 30; ++i)
  {
      T = 3 - i * 0.05;
      s = postopek<Spini>(s, o);
  }
  o.close();
}

void siroko()
{
  ofstream o;
  o.open("../../sp-siroko.dat");
  
  Spini s;
  for (int i = 0; i < 30; ++i)
  {
      T = 1 + i * 0.3;
      s = postopek<Spini>(s, o);
  }
  o.close();
}

void polje()
{
  ofstream o;
  char buf[64];
  sprintf(buf, "../../sp-polje-%g.dat", H);
  o.open(buf);
  Spini s;
  for (int i = 0; i < 20; ++i)
  {
    T = 2 + i * 0.1;
    s = postopek(s, o);
  }
  
  o.close();
}

void ver()
{
  Veriga v;
  ofstream o;
  char buf[64];
  sprintf(buf, "../../veriga-%g.dat", k);
  o.open(buf);
  for (int i = 0; i < 100; ++i)
  {
    T = i * 0.2;
    v = postopek(v, o);
  }
  o.close();
}

int main(int argc, char **argv) {
  double kji[] = {0.1, 0.3, 1, 3, 10};
  for (int i = 0; i < 5; ++i)
  {
    k = kji[i];
    ver();
  }
  return 0;
}
