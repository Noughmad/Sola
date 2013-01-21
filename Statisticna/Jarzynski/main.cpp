#include "system.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

using namespace std;

void testEnergyChanges()
{
    System sys(32);
    sys.setBeta(0.2);
    sys.setH(0);
    for (int i = 0; i < 64; ++i)
    {
        for (int j = 0; j < 32; ++j)
        {
            sys.value(i, j) = 1;
        }
    }

    for (double h = 0; h < 1; h += 0.1)
    {
        sys.setH(h);
        std::cout << h << ", " << sys.energyChange(3, 3) << std::endl;
    }

    sys.setBeta(0.4);

    for (double h = 0; h < 1; h += 0.1)
    {
        sys.setH(h);
        std::cout << h << ", " << sys.energyChange(3, 3) << std::endl;
    }
}

double simulation(int L, int tau, double beta, double H)
{
    System sys(L);
    sys.setBeta(beta);
    sys.setH(0);
    sys.randomState();

    int EqSteps;
    if (beta < 0.3)
    {
        EqSteps = 10;
    }
    else if (beta > 0.6)
    {
        EqSteps = 80;
    }
    else
    {
        EqSteps = 300;
    }

    sys.metropolisSteps(L*L*EqSteps);

    double work = 0;
    const double hdot = H/tau;
    const int n = L*L / 16;

    for (int t = 0; t < tau; ++t)
    {
        sys.setH(sys.h() + hdot);
        sys.metropolisSteps(n);
        work -= hdot * sys.magnetization();
    }
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
    int TT[] = {25, 50, 100, 0};
    for (int i = 0; TT[i]; ++i)
    {
        complete(L, TT[i], beta, H, N);
        std::cout << "Completed for time " << TT[i] << std::endl;
    }
}

void repeatFromArticle()
{
    const int L = 32;
    const int N = 1e5;
    const double betaC = log(1 + sqrt(2)) / 2;

    timeDependence(L, 0.2, 1.5, N);
    timeDependence(L, betaC, 0.1/betaC, N);
    timeDependence(L, betaC, 0.3/betaC, N);
    timeDependence(L, 0.7, 0.1/0.7, N);
}

int main(int argc, char **argv) {
    srand(time(0));

    if (argc < 3)
    {
        testEnergyChanges();
    }
    else
    {
        int L = 32;
        int N = 1e4;

        double beta;
        if (strcmp(argv[1], "c") == 0)
        {
            beta = log(1 + sqrt(2)) / 2;
        }
        else
        {
            beta = atof(argv[1]);
        }

        double betaH = atof(argv[2]);
        timeDependence(L, beta, betaH/beta, N);
    }
    return 0;
}
