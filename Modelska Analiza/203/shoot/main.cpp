#include <iostream>
#include <cassert>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_multiroots.h>

typedef int (*odvod)(double, const double[], double[], void*);

typedef struct {
   gsl_odeiv2_system system; 
   gsl_vector* final_state;
} m_system;

int f_vrv(double t, const double y[], double f[], void* params)
{
    
}

int konec(const gsl_vector* state, void* params, gsl_vector* values)
{
    m_system m = *(m_system*)params;
    assert(state->size == m.system.dimension);
    
    gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&m.system, gsl_odeiv2_step_rk4, 1e-6, 1e-8, 1e-8);
    
    double t = 0.0;
    gsl_vector* v = gsl_vector_alloc(state->size);
    for (int i = 0; i < state->size; ++i)
    {
        gsl_vector_set(v, i, gsl_vector_get(state, i));
    }
    int status = gsl_odeiv2_driver_apply(driver, &t, 1.0, v->data);
    if (status == GSL_SUCCESS)
    {
        gsl_vector_sub(v, m.final_state);
        *values = *v;
    }
    return status;
}

void vrv

int main(int argc, char **argv) {
    std::cout << "Hello, world!" << std::endl;
    return 0;
}
