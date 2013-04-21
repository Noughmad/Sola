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


#ifndef OSCILATOR_H
#define OSCILATOR_H

#include <NUtils/observable.h>
#include <NUtils/vector.h>
#include <gsl/gsl_rng.h>

template <int M>
class Oscilator
{
    
    typedef ::Observable<int, double, double> Observable;
    typedef double State;
    typedef ::Vector<M, State> Path;

public:
    Oscilator();
    virtual ~Oscilator();
    
    double energy(State one, State two);
    double energy();
    double logMatrixElement(State one, State two);
    double logMatrixElementTotal(int start, int end);
    
    double potentialEnergy();
    
    double lambda;
    double beta;
    double epsilon;

    Observable E;
    Observable X;
    Observable PotE;
    Path mPath;
    
    void setRandomPath();
    bool step();
    int manySteps(int steps);
    void measure();
    
    State randomState();
    double potential(State state);
    
    gsl_rng* rng;
};


#include <stdlib.h>
#include <gsl/gsl_randist.h>
#include <iostream>

using namespace std;

template <int M>
Oscilator<M>::Oscilator()
{
    rng = gsl_rng_alloc(gsl_rng_default);
    // setRandomPath();
    
    for (int j = 0; j < M; ++j)
    {
        mPath[j] = /* gsl_ran_gaussian(rng, 2) */ 0;
    }
}

template <int M> Oscilator<M>::~Oscilator()
{
    gsl_rng_free(rng);
}

template <int M> bool Oscilator<M>::step()
{
    int r = rand() % 20;
    if (r == 0)
    {
        // Move ALL the points
        
        Path oldPath = mPath;
        double oldEnergy = logMatrixElement(0, M);
        
        double delta = randomState();
        for (int j = 0; j < M; ++j)
        {
            mPath[j] += delta;
        }
        
        double dE = logMatrixElement(0, M) - oldEnergy;
        if (dE <= 0 || gsl_rng_uniform(rng) < exp(-dE))
        {
            return true;
        }
        else
        {
            mPath = oldPath;
            return false;
        }
    }
    else if (r == 1)
    {
        // Move some points
        
        Path oldPath = mPath;
        
        int start = rand() % M;
        int end = rand() % M;
        
        if (end <= start)
        {
            end += M;
        }
        double oldE = logMatrixElementTotal(start, end);
        
        double delta = randomState();
        for (int j = start; j < end; ++j)
        {
            mPath[j % M] += delta;
        }
        
        double dE = logMatrixElementTotal(start, end) - oldE;
        if (dE <= 0 || gsl_rng_uniform(rng) < exp(-dE))
        {
            return true;
        }
        else
        {
            mPath = oldPath;
            return false;
        }
    }
    else 
    {
        // Move just one point
        int i = rand() % M;
        int prev = (i+M-1) % M;
        int next = (i+1) % M;
        
        State oldState = mPath[i];
        State newState = oldState + randomState();
        
        double dE = logMatrixElement(mPath[prev], newState) + logMatrixElement(newState, mPath[next]) - logMatrixElement(mPath[prev], oldState) - logMatrixElement(oldState, mPath[next]);
        
        if (dE <= 0 || gsl_rng_uniform(rng) < exp(-dE))
        {
            mPath[i] = newState;
            return true;
        }
        return false;
    }
}

template <int M>
int Oscilator<M>::manySteps(int steps)
{
    int accepted = 0;
    for (int i = 0; i < steps; ++i)
    {
        if (step())
        {
            ++accepted;
        }
    }
    
    // cout << "Accepted " << accepted << " out of " << steps << endl;
    return accepted;
}

template <int M> 
double Oscilator<M>::energy()
{
    double e = 0;
    for (int j = 0; j < M; ++j)
    {
        int next = (j+1) % M;
        e += energy(mPath[j], mPath[next]);
    }
    return e;
}

template <int M>
double Oscilator<M>::energy(Oscilator<M>::State one, Oscilator<M>::State two)
{
    return 0.5 / beta * (1.0 - M / beta * (one - two) * (one - two)) + 1.0 / M * potential(one);
}

template <int M> 
double Oscilator<M>::logMatrixElement(Oscilator<M>::State one, Oscilator<M>::State two)
{
    return M * 0.5 / beta * (one - two) * (one - two) + beta / M * potential(one);
}

template <int M> 
double Oscilator<M>::logMatrixElementTotal(int start, int end)
{
    double e = 0;
    if (end <= start)
    {
        end += M;
    }
    
    if (start == 0 && end != M)
    {
        e += logMatrixElement(mPath[M-1], mPath[0]);
    }
    else if (start != 0 && (end - start) != M)
    {
        e += logMatrixElement(mPath[start-1], mPath[0]);
    }
    
    for (int i = start; i < end; ++i)
    {
        e += logMatrixElement(mPath[i % M], mPath[(i+1) % M]);
    }
    return e;
}

template <int M> 
typename Oscilator<M>::State Oscilator<M>::randomState()
{
    return gsl_ran_gaussian(rng, epsilon);
}

template <int M>
double Oscilator<M>::potential(State state)
{
    return state * state * (0.5 + lambda * state * state);
}

template <int M>
void Oscilator<M>::measure()
{
    
    E.insertValue(0, energy());
    
    double x = 0;
    for (int j = 0; j < M; ++j)
    {
        x += pow(mPath[j], 2);
    }
    X.insertValue(0, x / M);
    PotE.insertValue(0, potentialEnergy());
}

template <int M> 
double Oscilator<M>::potentialEnergy()
{
    double p = 0;
    for (int i = 0; i < M; ++i)
    {
        p += potential(mPath[i]);
    }
    return p / M;
}



#endif // OSCILATOR_H
