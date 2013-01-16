#include <iostream>
#include "system.h"

int main(int argc, char **argv) {
    System sys(128);
    sys.beta = 0.2;

    sys.metropolis();

    std::cout << "Hello, world!" << std::endl;
    return 0;
}
