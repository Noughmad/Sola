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

#define forto(i, n) for(int i = 0; i < n; ++i)
#define foreach(it, vector) for (auto it = vector.begin(); it != vector.end(); ++it)

const double dt = 0.01;
const double beta = 1.0;

const double StStat = 1e3;

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
    std::cout << T << std::endl;
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
    std::vector<double> stat(StStat);
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

gsl_histogram* histogram(data d)
{
    int max = *(std::max_element(d.begin(), d.end()));
    int min = *(std::min_element(d.begin(), d.end()));
    gsl_histogram* h = gsl_histogram_alloc(max - min + 1);
    gsl_histogram_set_ranges_uniform(h, min, max);
    
    foreach (it, d)
    {
        gsl_histogram_increment(h, *it);
    }
    return h;
}

void zajlis(int Z, int L, korak_zl k)
{
    data d = stat_zl(Z, L, k);
    std::ofstream o;
    char buffer[64];
    sprintf(buffer, "../../zajlis-%d-%d.dat", Z, L);
    o.open(buffer);
    gsl_histogram* h = histogram(d);
    int B = gsl_histogram_bins(h);
    forto(i,B)
    {
        double min, max;
        gsl_histogram_get_range(h, i, &min, &max);
        o << (min+max)/2 << '\t' << gsl_histogram_get(h,i) << std::endl;
    }
    o.close();
}

void izumrtje(int N, korak k, test* t)
{
    data d = statistika(N, k);
    std::ofstream o;
    char buffer[64];
    sprintf(buffer, "../../izumrtje-%d.dat", N);
    o.open(buffer);
    gsl_histogram* h = histogram(d);
    int B = gsl_histogram_bins(h);
    forto(i,B)
    {
        double min, max;
        gsl_histogram_get_range(h, i, &min, &max);
        o << (min+max)/2 << '\t' << gsl_histogram_get(h,i) << std::endl;
    }
    o.close();
}

int main(int argc, char **argv) {
    izumrtje(25, korak_exp, 0);
    izumrtje(250, korak_exp, 0);
    zajlis(200, 30, korak_zajlis);
    return 0;
}
