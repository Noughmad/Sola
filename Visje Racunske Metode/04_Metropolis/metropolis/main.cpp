#include <iostream>
#include <fstream>

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

void oscilator(double beta, double lambda, int InitialSteps, int AverageSteps, std::ostream& stream)
{
    Oscilator O;
    O.lambda = lambda;
    O.beta = beta;
    O.epsilon = min(0.1, beta);
    
    cout << endl;
    cout << "Starting with beta = " << beta << endl;
        
    for (int i = 0; i < InitialSteps; ++i)
    {
        int a = O.manySteps(1000);
    }
    
    for (int i = 0; i < AverageSteps; ++i)
    {        
        O.manySteps(1000);
        O.measure();
    }
    
    stream << O.beta << ", " << O.lambda << ", " << O.epsilon << ", " << O.E.average() << ", " << sqrt(O.E.variance()) << ", " << O.X.average() << ", " << sqrt(O.X.variance()) << endl;
};

int main(int argc, char **argv) {
    
    ofstream out_0("g_energija_0.dat");
    ofstream out_03("g_energija_0.3.dat");
    ofstream out_1("g_energija_1.dat");
    ofstream out_3("g_energija_3.dat");
    
    const int I = 10000;
    const int A = 10000;
    
    const double F = sqrt(10);
    
#pragma omp parallel for
    for (int blog = -10; blog < 8; ++blog)
    {
        double beta = pow(F, blog);
        
        oscilator(beta, 0, I, A, out_0);
        oscilator(beta, 0.3, I, A, out_03);
        oscilator(beta, 1, I, A, out_1);
        oscilator(beta, 3, I, A, out_3);
    }
    
    out_0.close();
    out_03.close();
    out_1.close();
    out_3.close();
    
    return 0;
}
