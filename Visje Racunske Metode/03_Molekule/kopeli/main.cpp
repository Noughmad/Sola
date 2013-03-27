#include <iostream>
#include "veriga.h"

using namespace std;

const int InitialSteps = 1e7;
const int AverageSteps = 1e7;
const int MeasureInterval = AverageSteps / 1;

const int N = 100;

int main(int argc, char **argv) {
    
    Hoover h(N, 1, 2);
    h.K = 1;
    h.Q = 1;
    h.invTau = 1;
    h.lambda = 1;
        
    h.setup();
    
    for (int i = 0; i < InitialSteps; ++i)
    {
        h.step();
    }
    
    double J[N];
    double T[N];
    
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
            for (int j = 0; j < N; ++j)
            {
                T[j] = T[j] * d * f + h.y[N+j] * h.y[N+j] * f;
            }

            J[0] = J[0] * d * f + h.y[N] * h.Vprime(h.y[0], 0, h.y[1]) * f;
            for (int j = 0; j < N-1; ++j)
            {
                J[j] = J[j] * d * f + h.y[N+j] * h.Vprime(h.y[j], h.y[j-1], h.y[j+1]) * f;
            }
            J[N-1] = J[N-1] * d * f + h.y[2*N-1] * h.Vprime(h.y[N-1], h.y[N-2], 0) * f;
        }
    }
        
    for (int j = 0; j < N; ++j)
    {
        std::cout << J[j] << ", " << T[j] << ", " << h.y[j] << ", " << h.y[N+j] << endl;
    }
        
    return 0;
}
