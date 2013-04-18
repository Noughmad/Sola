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

class Oscilator
{
    static const int M = 50;
    
    typedef ::Observable<int, double, double> Observable;
    typedef double State;
    typedef ::Vector<M, State> Path;

public:
    Oscilator();
    virtual ~Oscilator();
    
    double energy(State one, State two);
    double energy();
    double logMatrixElement(State one, State two);
    
    double lambda;
    double beta;
    double epsilon;

    Observable E;
    Observable X;
    Path mPath;
    
    void setRandomPath();
    bool step();
    int manySteps(int steps);
    void measure();
    
    State randomState();
    double potential(State state);
    
    gsl_rng* rng;
};

#endif // OSCILATOR_H
