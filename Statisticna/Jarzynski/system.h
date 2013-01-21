/*
     Copyright (C) 2013  Miha ?an?ula <miha@noughmad.eu>
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


#ifndef SYSTEM_H
#define SYSTEM_H

#include <math.h>

typedef char spin_t;

class System
{
public:
    System(int L);
    ~System();

    inline spin_t& value(int i, int j) const;

    void randomState();
    double magnetization();

    double energyChange(int i, int j) const;
    inline void metropolis();
    void metropolisSteps(int N);

    void snapshot(const char* filename);

public:
    inline double h() const {return mH;}
    inline void setH(double h)
    {
        mH = h;
        expBH = exp(2 * beta * mH);
        expmBH = exp(-2 * beta * mH);
    }

    inline double setBeta(double b)
    {
        beta = b;
        for (spin_t i = -4; i <= 4; ++i)
        {
            expBeta[i+4] = exp(2 * i * b);
        }
    }

    inline double expBetaE(spin_t delta) const
    {
        return expBeta[delta + 4];
    }

private:
    spin_t* data;
    const short int L;
    const short int mask;
    const int RandMax;
    double expBH;
    double expmBH;
    double mH;
    double beta;
    double expBeta[9];
};


spin_t& System::value(int i, int j) const
{
    return data[(i & mask) * L + (j & mask)];
}

#endif // SYSTEM_H
