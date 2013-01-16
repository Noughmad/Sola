#include <iostream>
#include "system.h"

void simulation(int L, int tau, double beta)
{
    System sys(L);
    sys.beta = beta;

    for (int i = 0; i < L*L*5; ++i)
    {
        sys.metropolis();
    }

    for (int t = 0; t < tau; ++t)
    {
        for (int i = 0; i < L*L; ++i)
        {
            sys.metropolis();
        }
    }
}

int main(int argc, char **argv) {
    simulation(1024, 100, 1);

    std::cout << "Hello, world!" << std::endl;
    return 0;
}
