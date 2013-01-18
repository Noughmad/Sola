#include "system.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_histogram.h>

double simulation(int L, int tau, double beta, double H)
{
    System sys(L);
    sys.beta = beta;
    sys.h = 0;
    sys.randomState();

    sys.metropolisSteps(L*L*5);

    double work = 0;
    const double hdot = H/tau;

    for (int t = 0; t < tau; ++t)
    {
        sys.h += hdot;
        sys.metropolisSteps(2*L*L);
        work -= hdot * sys.magnetization();
    }

    work /= (L*L);

    std::cout << "Work: " << work << std::endl;

    return work;
}

void complete(int L, int tau, double beta, double H, int N)
{
    char buf[64];
    sprintf(buf, "../data/g_work_%d_%g_%g.dat", tau, beta, H);
    FILE* f = fopen(buf, "wt");

    for (int i = 0; i < N; ++i)
    {
        fprintf(f, "%g\n", simulation(L, tau, beta, H));
    }

    fclose(f);
}

void timeDependence(int L, double beta, double H, int N)
{
    int TT[] = {15, 25, 50, 100, 150, 0};
    for (int i = 0; TT[i]; ++i)
    {
        complete(L, TT[i], beta, H, N);
        std::cout << "Completed for time " << TT[i] << std::endl;
    }
}

void repeatFromArticle()
{
    const int L = 256;
    const int N = 1000000;
    const double betaC = log(1 + sqrt(2)) / 2;

    timeDependence(L, 0.2, 1.5, N);
    timeDependence(L, betaC, 0.1/betaC, N);
    timeDependence(L, betaC, 0.3/betaC, N);
    timeDependence(L, 0.7, 0.1/0.7, N);
}

int main(int argc, char **argv) {
    srand(time(0));

    repeatFromArticle();

    return 0;
}
