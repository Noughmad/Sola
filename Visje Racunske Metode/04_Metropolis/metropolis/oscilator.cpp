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
        mPath[j] = 0;
    }
}

Oscilator::~Oscilator()
{
    gsl_rng_free(rng);
}

void Oscilator::setRandomPath()
{
    for (int i = 0; i < M; ++i)
    {
        mPath[i] = randomState();
    }
}

bool Oscilator::step()
{
    int i = rand() % M;
    int prev = (i+M-1) % M;
    int next = (i+1) % M;

    State oldState = mPath[i];
    State newState = oldState + epsilon * randomState();
    
    double dE = energy(mPath[prev], newState) + energy(newState, mPath[next]) - energy(mPath[prev], oldState) - energy(mPath[prev], oldState);
     
    if (dE <= 0 || (double)rand() / RAND_MAX < exp(-beta * dE))
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
    
    cout << "Steps accepted: " << accepted << " out of " << steps << endl;
    return accepted;
}


double Oscilator::energy(Oscilator::State one, Oscilator::State two)
{
    return + M / 2.0 / beta / beta * (one - two) * (one - two) + 1.0 / M * potential(one);
}


double Oscilator::matrixElement(Oscilator::State one, Oscilator::State two)
{
    return exp(-beta * energy(one, two));
}
    
Oscilator::State Oscilator::randomState()
{
    return gsl_ran_gaussian(rng, 1);
}

double Oscilator::potential(State state)
{
    return state * state * (0.5 + lambda * state * state);
}

void Oscilator::measure()
{
    double e = M / 2 / beta;
    for (int j = 0; j < M; ++j)
    {
        int next = (j+1) % M;
        e += energy(mPath[j], mPath[next]);
    }
    E.insertValue(0, e);
}

