#include <iostream>
#include <cassert>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>



int konec_vrvi(const gsl_vector* state, void* params, gsl_vector* values)
{
    assert(state->size == 2);
    values->size = 2;
    
    return GSL_SUCCESS;
}

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    return 0;
}
