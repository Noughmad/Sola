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


#include "system.h"

#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using namespace std;

System::System(int L)
: L(L)
, data(new spin_t[L*L])
, mask(L-1)
, RandMax(RAND_MAX / L / L)
{
    assert((mask & L) == 0);
}

System::~System()
{
    delete[] data;
}

spin_t& System::value(int i, int j) const
{
    return data[(i & mask) * L + (j & mask)];
}

void System::randomState()
{
    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            value(i, j) = (rand() < (RAND_MAX/2)) ? -1 : 1;
        }
    }
}

double System::energyChange(int i, int j) const
{
    const spin_t factor = value(i+1, j) + value(i-1, j) + value(i, j+1) + value(i, j-1) + h;
    return 2 * value(i, j) * factor;
}

void System::metropolis()
{
    int r = rand();
    int i = r & mask;
    r /= L;
    int j = r & mask;
    r /= L;

    double ee = exp(-beta * energyChange(i, j));

    if (ee > 1 || ee > ((double)rand() / RandMax))
    {
        value(i, j) *= -1;
    }
}

double System::magnetization()
{
    double M = 0;
    for (int i = 0; i < L*L; ++i)
    {
        M += data[i];
    }
    return M;
}

void System::metropolisSteps(int N)
{
    #pragma omp parallel for
    for (int i = 0; i < N; ++i)
    {
        metropolis();
    }
}

void System::snapshot(const char* filename)
{
    ofstream out(filename);
    for (int i = 0; i < L; ++i)
    {
        for (int j = 0; j < L; ++j)
        {
            out << (int)value(i, j) << " ";
        }
        out << endl;
    }
    out.close();
}
