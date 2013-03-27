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


#ifndef VERIGA_H
#define VERIGA_H

#include <gsl/gsl_odeiv2.h>

class Veriga
{
public:
    Veriga(int N);
    virtual ~Veriga();
    
    virtual void setup() = 0;
    virtual void step() = 0;
    virtual size_t size() const;
    
    double Vprime(double current, double previous, double next);
    double Vprime(double x);
    
    const int N;
    double lambda;

    double K;
    double Q;
    
    double* y;
};

class Hoover : public Veriga
{
public:
    Hoover(int N, double T_L, double T_R);
    virtual ~Hoover();
    
    double T_L;
    double T_R;
    double invTau;
    
    virtual void setup();
    virtual void step();
    virtual size_t size() const;
    
    double t;
    
private:
    gsl_odeiv2_system* sys;
    gsl_odeiv2_driver* driver;
};

#endif // VERIGA_H
