/*
    <one line to give the program's name and a brief idea of what it does.>
    Copyright (C) 2013  Miha Cancula <email>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/


#include "veriga.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

Veriga::Veriga(int N) : N(N)
{
    y = new double[2*N];
}

Veriga::~Veriga()
{
    delete[] y;
}

size_t Veriga::size() const
{
    return N;
}

double Veriga::Uprime(double x)
{
    return (Q + 4 * lambda * x * x) * x;
}

double Veriga::Vprime(double x)
{
    return K * x;
}

double Veriga::Vprime(double current, double previous, double next)
{
    return K * (2 * current - previous - next);
}

int odvod(double t, const double y[], double dydt[], void* params)
{
    Hoover* h = static_cast<Hoover*>(params);
    int n = h->size();
    
    for (int i = 0; i < n; ++i)
    {
        // Prvih n so q_i
        dydt[i] = y[n+i];
    }
    
    // p_1 ima zraven clen z zeta_L = y[2*n]
    dydt[n] = -h->Vprime(y[0], 0, y[1]) - h->Uprime(y[0]) - y[2*n]*y[n];
    
    for (int i = 1; i < n-1; ++i)
    {
        dydt[n+i] = -h->Vprime(y[i], y[i-1], y[i+1]) - h->Uprime(y[i]);
    }
    
    // p_n ima zraven clen z zeta_R = y[2*n+1]
    dydt[2*n-1] = -h->Vprime(y[n-1], y[n-2], 0) - h->Uprime(y[n-1]) - y[2*n+1]*y[2*n-1];
    
    dydt[2*n] = h->invTau * (y[n]*y[n] - h->T_L);
    dydt[2*n+1] = h->invTau * (y[2*n-1]*y[2*n-1] - h->T_R);
    
    return GSL_SUCCESS;
}

Hoover::Hoover(int N, double T_L, double T_R)
: Veriga(N+1)
, T_L(T_L)
, T_R(T_R)
{
    
}

Hoover::~Hoover()
{
    gsl_odeiv2_driver_free(driver);
    delete sys;
}

void Hoover::setup()
{
    t = 0;
    sys = new gsl_odeiv2_system {odvod, 0, 2*N, static_cast<void*>(this)};
    driver = gsl_odeiv2_driver_alloc_y_new(sys, gsl_odeiv2_step_rk4, 1e-2, 1e-5, 0);
    
    gsl_rng* r = gsl_rng_alloc(gsl_rng_default);
    
    int n = size();
    for (int i = 0; i < n; ++i)
    {
        y[i] = 0;
        y[n+i] = gsl_rng_uniform(r);
    }
    y[2*n] = 0;
    y[2*n+1] = 0;
    
    gsl_rng_free(r);
}

void Hoover::step()
{
    gsl_odeiv2_driver_apply(driver, &t, t + 1, y);
}

size_t Hoover::size() const
{
    return N-1;
}

Maxwell::Maxwell(int N, double T_L, double T_R): Veriga(N)
{
    rng = gsl_rng_alloc(gsl_rng_default);
}

Maxwell::~Maxwell()
{
    gsl_rng_free(rng);
}


void Maxwell::setup()
{
    const int n = size();
    rng = gsl_rng_alloc(gsl_rng_default);
    for (int i = 0; i < n; ++i)
    {
        y[i] = 0;
        y[n+i] = gsl_ran_chisq(rng, 3); // TODO: Normalization
    }
    
    stepNumber = 0;
}

void Maxwell::step()
{
    // TODO: This is S2. Use S4 instead
    stepKinetic(0.5 * h);
    stepPotential(1 * h);
    stepKinetic(0.5 * h);
    
    if ((stepNumber % resetInterval) == 0)
    {
        
    }
}

void Maxwell::stepKinetic(double x)
{
    const int n = size();
    for (int i = 0; i < n; ++i)
    {
        y[i] += x * y[i+n];
    }
}

void Maxwell::stepPotential(double x)
{
    const int n = size();
    y[n] -= x * (Vprime(y[n], 0, y[n+1] + Uprime(y[n]));
    for (int i = 0; i < n; ++i)
    {
        y[n+i] -= x * (Vprime(y[n+i], y[n+i-1], y[n+i+1])  + Uprime(y[n]));
    }
    y[n] -= x * (Vprime(y[2*n-1], y[2*n-2], 0 + Uprime(y[2*n-1]));
}


