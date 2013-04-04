#include <iostream>
#include <cmath>
#include "veriga.h"

using namespace std;

const int InitialSteps = 1e8;
const int AverageSteps = 1e7;
const int MeasureInterval = 10;

const int N = 30;

inline double sqr(double x)
{
    return x*x;
}

typedef Maxwell V;

int main(int argc, char **argv) {
    
    /*
    Hoover h(N, 1, 2);
    */
    
    V h(N, 6, 3);
    
    // h.invTau = 10;
    h.resetInterval = 10;
    h.h = 1e-2;
    
    h.K = 1;
    h.Q = 1;
    h.lambda = 100;
        
    h.setup();
    
    for (int i = 0; i < InitialSteps; ++i)
    {
        h.step();
    }
    
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
            
            J[0] = J[0] * d * f + h.y[N] * h.Vprime(h.y[0], 0, h.y[1]) * f;
            J2[0] = J2[0] * d * f + sqr(h.y[N] * h.Vprime(h.y[0], 0, h.y[1])) * f;
            
#pragma omp parallel for
            for (int j = 1; j < N-1; ++j)
            {
                J[j] = J[j] * d * f + h.y[N+j] * h.Vprime(h.y[j], h.y[j-1], h.y[j+1]) * f;
                J2[j] = J2[j] * d * f + sqr(h.y[N+j] * h.Vprime(h.y[j], h.y[j-1], h.y[j+1])) * f;
            }
            J[N-1] = J[N-1] * d * f + h.y[2*N-1] * h.Vprime(h.y[N-1], h.y[N-2], 0) * f;
            J2[N-1] = J2[N-1] * d * f + sqr(h.y[2*N-1] * h.Vprime(h.y[N-1], h.y[N-2], 0)) * f;
        }
    }
    
    const double M = sqrt((double)AverageSteps / MeasureInterval);
        
    for (int j = 0; j < N; ++j)
    {
        std::cout << J[j] << ", " << sqrt(J2[j] - J[j] * J[j]) / M << ", " << T[j] << ", " << sqrt(T2[j] - T[j] * T[j]) / M << ", " << h.y[j] << ", " << h.y[N+j] << endl;
    }
        
    return 0;
}
