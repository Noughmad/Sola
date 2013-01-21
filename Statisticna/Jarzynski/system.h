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

typedef char spin_t;

class System
{
public:
    System(int L);
    ~System();

    inline spin_t& value(int i, int j) const;

    void randomState();
    double magnetization();

    inline double energyChange(int i, int j) const;
    inline void metropolis();
    void metropolisSteps(int N);

    void snapshot(const char* filename);

public:
    double beta;
    double h;

private:
    spin_t* data;
    const short int L;
    const short int mask;
    const int RandMax;
};

#endif // SYSTEM_H
