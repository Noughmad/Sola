#include <iostream>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sort.h>
#include <math.h>

#include <vector>
#include <fstream>

#include <algorithm>
#include <numeric>
#include <iterator>

#include <cstring>

#define forto(i, n) for(int i = 0; i < n; ++i)
#define foreach(it, vector) for (auto it = vector.begin(); it != vector.end(); ++it)

#define G_OUTPUT 0
#define G_STAT 1

double dt = 0.01;
double beta = 1.0;

double StStat = 1e4;

typedef int (*korak)(int N, gsl_rng* r);
typedef void (*korak_zl)(int *Z, int *L, gsl_rng* r);
typedef std::vector<double> data;
typedef double (*test)(int n, data D);

double kolmogorov(int n, data d1, data d2)
{
    std::sort(d1.begin(), d1.end());
    std::sort(d2.begin(), d2.begin());
    
    int d = 0;
    int i=0,j=0;
    while (i != n || j != n)
    {
        if (d1[i] >= d2[j])
        {
            ++j;
        }
        if (d1[i] <= d2[j])
        {
            ++i;
        }
        d = std::max(d, abs(i - j));
    }
    
    return d/sqrt(n);
}

double smrt_exp(int N);
double smrt_rs(int N);

int korak_exp(int N, gsl_rng* r)
{
    return N - gsl_ran_poisson(r, N * beta * dt);
}

int korak_rs(int N, gsl_rng* r)
{
    return N - gsl_ran_poisson(r, N * 5 * beta * dt) + gsl_ran_poisson(r, N * 4 * beta * dt);
}

double smrt(int N, korak k, gsl_rng* r)
{
    double T = 0;
    while (N > 0)
    {
        N = k(N,r);
        T += dt;
    }
    return T;
}

double smrt_zl(int Z, int L, korak_zl k, gsl_rng* r)
{
    double T = 0;
    while (Z > 0 && L > 0)
    {
        k(&Z, &L, r);
        T += dt;
    }
    return T;
}

void korak_zajlis(int *Z, int *L, gsl_rng* r)
{
    int z = *Z;
    int l = *L;
    int dz = gsl_ran_poisson(r, 5*beta*z*dt) - gsl_ran_poisson(r, 4*beta*z*dt + beta*z*l/50.0*dt);
    int dl = gsl_ran_poisson(r, 4*beta*l*dt + beta*z*l/200.0*dt) - gsl_ran_poisson(r, 5*beta*l*dt);
    *Z += dz;
    *L += dl;
}

data statistika(int N, korak k)
{
    data stat(StStat);
    gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, time(0));
    forto(i,StStat)
    {
        stat[i] = smrt(N,k,r);
    }
    gsl_rng_free(r);
    return stat;
}

data stat_zl(int Z, int L, korak_zl k)
{
    data stat(StStat);
    gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r, time(0));
    forto(i,StStat)
    {
        stat[i] = smrt_zl(Z, L, k, r);
    }
    gsl_rng_free(r);
    return stat;
}

gsl_histogram* histogram(data d, int max)
{
    // double max = *(std::max_element(d.begin(), d.end()));
    // double min = *(std::min_element(d.begin(), d.end()));
		
		int bins = max / dt;
		
		std::cout << "Making histogram with " << bins << " bins" << std::endl;
		
    gsl_histogram* h = gsl_histogram_alloc(bins);

    gsl_histogram_set_ranges_uniform(h, dt/2, max+dt/2);
    
    foreach (it, d)
    {
        gsl_histogram_increment(h, *it);
    }
    
    gsl_histogram_scale(h, 1.0 / dt / StStat);
    
    return h;
}

void zajlis(int Z, int L, korak_zl k, std::ostream& out)
{
    data d = stat_zl(Z, L, k);
#if G_OUTPUT
    gsl_histogram* h = histogram(d, 60);
    int B = gsl_histogram_bins(h);
    forto(i,B)
    {
        double min, max;
        gsl_histogram_get_range(h, i, &min, &max);
        out << (min+max)/2 << '\t' << gsl_histogram_get(h,i) << std::endl;
    }
#endif

#if G_STAT    
    std::cout << "Rezultat:\t" << gsl_stats_mean(d.data(), 1, StStat) << " & " << gsl_stats_sd(d.data(), 1, StStat) << std::endl;
#endif
}

void izumrtje(int N, korak k, test* t, std::ostream& out)
{
    data d = statistika(N, k);
#if G_OUTPUT
    gsl_histogram* h = histogram(d, 15);
    int B = gsl_histogram_bins(h);
    forto(i,B)
    {
        double min, max;
        gsl_histogram_get_range(h, i, &min, &max);
        out << (min+max)/2 << '\t' << gsl_histogram_get(h,i) << std::endl;
    }
#endif

#if G_STAT
		std::cout << "Rezultat:\t" << gsl_stats_mean(d.data(), 1, StStat) << " & " << gsl_stats_sd(d.data(), 1, StStat) << std::endl;
#endif
}

int main(int argc, char **argv) {
		
		if (argc < 6)
		{
				std::cout << "Usage: populacije <outfile> <dt> <StStat> izumrtje <rs|exp> N" << std::endl;
				std::cout << "   or: populacije <outfile> <dt> <StStat> zajlis Z L" << std::endl;
				return 0;
		}
		
		forto(i,7)
				std::cout << argv[i] << " ";
		std::cout << std::endl;
		
		std::ofstream o;
#if G_OUTPUT
		o.open(argv[1]);
#endif
		dt = atof(argv[2]);
		StStat = atoi(argv[3]);
		char* op = argv[4];
		if (strcmp(op, "izumrtje") == 0)
		{
				if (strcmp(argv[5], "rs") == 0)
				{
						izumrtje(atoi(argv[6]), korak_rs, 0, o);
				}
				else
				{
						izumrtje(atoi(argv[6]), korak_exp, 0, o);
				}
		}
		else if (strcmp(op, "zajlis") == 0)
		{
				zajlis(atoi(argv[5]), atoi(argv[6]), korak_zajlis, o);
		}
		o.close();
    return 0;
}
