#include <iostream>
#include "system.h"

void simulation(int L, int tau, double beta, double H)
{
    System sys(L);
    sys.beta = beta;
    sys.h = 0;
    sys.randomState();

    for (int i = 0; i < L*L*5; ++i)
    {
        sys.metropolis();
    }

    double work = 0;
    for (int t = 0; t < tau; ++t)
    {
        sys.h += H/tau;
        for (int i = 0; i < L*L; ++i)
        {
            sys.metropolis();
        }
        std::cout << sys.magnetization() << std::endl;
        work -= H/tau * sys.magnetization();
    }

    std::cout << "Work: " << work << std::endl;
}

int main(int argc, char **argv) {
    simulation(256, 100, 0.02, 0.3);

    return 0;
}
