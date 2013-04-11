#include <iostream>
#include <cmath>
#include <string.h>
#include "veriga.h"

using namespace std;

inline double sqr(double x)
{
    return x*x;
}

void run_sim(Veriga& h, int InitialSteps, int AverageSteps, int MeasureInterval)
{
    h.setup();

    for (int i = 0; i < InitialSteps; ++i)
    {
        h.step();
    }
    
    int N = h.size();

    double J[N];
    double T[N];

    double J2[N];
    double T2[N];

    for (int j = 0; j < N; ++j)
    {
        J[j] = 0;
        T[j] = 0;
    }

    for (int i = 0; i < AverageSteps; ++i)
    {
        h.step();
        if ((i % MeasureInterval) == 0)
        {
            const int d = i / MeasureInterval;
            const double f = 1.0 / (d+1);
            
            #pragma omp parallel for
            for (int j = 0; j < N; ++j)
            {
                T[j] = T[j] * d * f + h.y[N+j] * h.y[N+j] * f;
                T2[j] = T2[j] * d * f + pow(h.y[N+j], 4) * f;
            }
            
            #pragma omp parallel for
            for (int j = 0; j < N; ++j)
            {
                J[j] = J[j] * d * f + h.flux(j) * f;
                J2[j] = J2[j] * d * f + sqr(h.flux(j)) * f;
            }
        }
    }

    const double M = sqrt((double)AverageSteps / MeasureInterval);

    for (int j = 0; j < N; ++j)
    {
        std::cout << J[j] << ", " << sqrt(J2[j] - J[j] * J[j]) / M << ", " << T[j] << ", " << sqrt(T2[j] - T[j] * T[j]) / M << ", " << h.y[j] << ", " << h.y[N+j] << endl;
    }
}

int main(int argc, char **argv) {
    
    string chainType(argv[1]);
    int numAtoms = atoi(argv[2]);
    int InitialSteps = atof(argv[3]);
    const int AverageSteps = atof(argv[4]);
    const int MeasureInterval = atof(argv[5]);

    double lambda = atof(argv[6]);

    if (chainType == "maxwell")
    {
        Maxwell m(numAtoms, 1, 2);
        m.K = 1;
        m.Q = 1;
        m.lambda = lambda;
        m.resetInterval = 100;
        m.h = 1e-2;
        run_sim(m, InitialSteps, AverageSteps, MeasureInterval);
    }
    else if (chainType == "hoover")
    {
        Hoover h(numAtoms, 1, 2);
        h.K = 1;
        h.Q = 1;
        h.lambda = lambda;
        h.invTau = 1;
        h.h = 1e-2;
        run_sim(h, InitialSteps, AverageSteps, MeasureInterval);
    }
}
