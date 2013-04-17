#include <iostream>

#include <NUtils/observable.h>

#include <cstdlib>
#include <boost/graph/graph_concepts.hpp>
#include "oscilator.h"

using namespace std;

inline int mod(int i, int n)
{
    int m = i % n;
    return (m >= 0) ? m : m+n;
}

typedef Observable<int, double> IsingObservable;

class Ising
{
public:
    Ising(int N);
    ~Ising();
    
    IsingObservable E;
    IsingObservable M;
    
    void step();
    void measure();
    
    double beta;
    
    inline int index(int i, int j)
    {
        int l = mod(i, N) * N + mod(j, N);
        assert(l >= 0);
        assert(l < N*N);
        return l;
    }
    
    inline char value(int i, int j)
    {
        return grid[index(i, j)];
    }
    
    inline void setValue(int i, int j, char value)
    {
        grid[index(i, j)] = value;
    }
    
private:
    char* grid;
    int N;
};

Ising::Ising(int N)
: N(N)
{
    grid = new char[N*N];
    
    for (int i = 0; i < N*N; ++i)
    {
        grid[i] = ((double)rand() / RAND_MAX > 0.5) ? -1 : 1;
    }
}

Ising::~Ising()
{
    delete[] grid;
}

void Ising::measure()
{
    double e = 0, m = 0;
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            m += value(i, j);
            e -= value(i, j) * (value(i+1, j) + value(i-1, j) + value(i, j+1) + value(i, j-1));
        }
    }
    
    E.insertValue(0, e / 2); // Because we counted every bond twice
    M.insertValue(0, m);
}

void Ising::step()
{
    int i = rand() % N;
    int j = rand() % N;
    
    double dE = 2 * value(i, j) * (value(i+1, j) + value(i-1, j) + value(i, j+1) + value(i, j-1));
    
    if (dE < 0 || ((double)rand() / RAND_MAX) < exp(-beta * dE))
    {
        setValue(i, j, -value(i, j));
    }
}

void racun(int N, int NMeasure, int NSteps)
{    
    Ising grid(N);
    grid.beta = 0;

    
    while (grid.beta < 10)
    {
        grid.beta += 0.1;
        
        for (int i = 0; i < NMeasure; ++i)
        {
            for (int i = 0; i < NSteps; ++i)
            {
                grid.step();
            }
            grid.measure();
        }
        
        cout << grid.beta << ", " << grid.E.average() << ", " << grid.E.variance() << ", " << grid.M.average() << ", " << grid.M.variance() << endl;
        
    }
}

void oscilator()
{
    Oscilator O;
    O.lambda = 1;
    O.beta = 1000;
    O.epsilon = 1;
    
    for (int i = 0; i < 3000; ++i)
    {
        int a = O.manySteps(100);
        if (a > 70)
        {
            O.epsilon *= pow(2, (a-70) * 0.01);
        }
        else if (a < 30)
        {
            O.epsilon *= pow(2, (a-30) * 0.01);
        }
        cout << O.mPath << endl;
        cout << "Adjusting epsilon to " << O.epsilon << endl;
    }
    
    for (int i = 0; i < 3000; ++i)
    {
        O.manySteps(1000);
        O.measure();
    }
    
    cout << O.beta << ", " << O.E.average() << ", " << O.E.variance() << endl;
    
    cout << O.mPath << endl;
};

int main(int argc, char **argv) {
    
    // racun(128, 10000, 10000);
    oscilator();
    
    
    return 0;
}
