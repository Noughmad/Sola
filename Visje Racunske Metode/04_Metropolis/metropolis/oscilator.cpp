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


#include "oscilator.h"
#include <stdlib.h>
#include <boost/iterator/iterator_concepts.hpp>
#include <gsl/gsl_randist.h>

#include <iostream>

using namespace std;

Oscilator::Oscilator()
{
    rng = gsl_rng_alloc(gsl_rng_default);
    // setRandomPath();
    
    for (int j = 0; j < M; ++j)
    {
        mPath[j] = gsl_ran_gaussian(rng, 1);
    }
}

Oscilator::~Oscilator()
{
    gsl_rng_free(rng);
}

bool Oscilator::step()
{
    int i = rand() % M;
    int prev = (i+M-1) % M;
    int next = (i+1) % M;

    State oldState = mPath[i];
    State newState = oldState + randomState();
    
    double dE = logMatrixElement(mPath[prev], newState) + logMatrixElement(newState, mPath[next]) - logMatrixElement(mPath[prev], oldState) - logMatrixElement(oldState, mPath[next]);
     
    if (dE <= 0 || gsl_rng_uniform(rng) < exp(- M * dE))
    {
        mPath[i] = newState;
        return true;
    }
    return false;
}

int Oscilator::manySteps(int steps)
{
    int accepted = 0;
    for (int i = 0; i < steps; ++i)
    {
        if (step())
        {
            ++accepted;
        }
    }
    
    return accepted;
}

double Oscilator::energy()
{
    double e = 0;
    for (int j = 0; j < M; ++j)
    {
        int next = (j+1) % M;
        e += energy(mPath[j], mPath[next]);
    }
    return e;
}

double Oscilator::energy(Oscilator::State one, Oscilator::State two)
{
    return M * 0.5 / beta / beta * (one - two) * (one - two) + 1.0 / M * potential(one);
}

double Oscilator::logMatrixElement(Oscilator::State one, Oscilator::State two)
{
    return  M * 0.5 / beta * (one - two) * (one - two) + beta / M * potential(one);
}
    
Oscilator::State Oscilator::randomState()
{
    return gsl_ran_gaussian(rng, epsilon);
}

double Oscilator::potential(State state)
{
    return state * state * (0.5 + lambda * state * state);
}

void Oscilator::measure()
{
    
    E.insertValue(0, energy());
    
    double x = 0;
    for (int j = 0; j < M; ++j)
    {
        x += pow(mPath[j], 2);
    }
    X.insertValue(0, x / M);
}

