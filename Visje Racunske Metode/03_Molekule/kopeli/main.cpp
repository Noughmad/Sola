#include <iostream>
#include "veriga.h"

using namespace std;

int main(int argc, char **argv) {
    
    Hoover h(10, 1, 2);
    h.K = 1;
    h.Q = 1;
    h.invTau = 1;
        
    h.setup();
    
    for (int i = 0; i < 100000; ++i)
    {
        h.step();
    }
    
    std::cout << h.t << endl;
    
    for (int i = 0; i < 22; ++i)
    {
        std::cout << h.y[i] << ", ";
    }
    
    std::cout << std::endl;
    return 0;
}
