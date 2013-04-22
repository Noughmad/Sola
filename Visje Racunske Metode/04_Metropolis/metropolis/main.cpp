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
            e -= 0.5 * value(i, j) * (value(i+1, j) + value(i-1, j) + value(i, j+1) + value(i, j-1));
        }
    }
    
    E.insertValue(0, e / 2 / N / N); // Because we counted every bond twice
    M.insertValue(0, m / N / N);
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
    double b = 0;

    
    while (b < 3)
    {
        b += 0.1;
        
        Ising grid(N);
        grid.beta = b;
        
        for (int i = 0; i < NMeasure * NSteps; ++i)
        {
            grid.step();
        }
        
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

template <int M>
void oscilatorImpl(double beta, double lambda, int InitialSteps, int AverageSteps, std::ostream& stream)
{
    Oscilator<M> O;
    O.lambda = lambda;
    O.beta = beta;
    O.epsilon = min(1.0, 0.1 * sqrt(beta));
    O.setRandomPath();
    
    cout << endl;
    cout << "Starting with beta = " << beta << endl;
        
    for (int i = 0; i < InitialSteps; ++i)
    {
        int a = O.manySteps(10000);
    }
  
 
    for (int i = 0; i < AverageSteps; ++i)
    {        
        O.manySteps(1000);
        O.measure();
    }
    
    stream << O.beta << ", " << O.lambda << ", " << O.epsilon << ", " 
        << O.E.average() << ", " << sqrt(O.E.variance()) << ", " 
        << O.PotE.average() << ", " << sqrt(O.PotE.variance()) << ", " 
        << O.X.average() << ", " << sqrt(O.X.variance()) << endl;
}


void oscilator(double beta, double lambda, int InitialSteps, int AverageSteps, std::ostream& stream)
{
    oscilatorImpl<5>(beta, lambda, InitialSteps, AverageSteps, stream);
}

void testConvergenceEpsilon()
{
    const char* format = "g_conv_%g.dat";
    for (double eps = 0.1; eps <= 2; eps += 0.1)
    {
        Oscilator<25> O;
        O.beta = 10.0;
        O.lambda = 0;
        O.epsilon = eps;
        O.setRandomPath();
        
        char buf[32];
        sprintf(buf, format, eps);
        
        
        ofstream stream(buf);
        
        for (int j = 0; j < 100; ++j)
        {
            double Eavg = 0;
            const int N = 100;
            for (int i = 0; i < N; ++i)
            {
                O.manySteps(1000);
                const double t = (double)i / N;
                Eavg = Eavg * t + O.energy() * (1-t);
            }
            stream << j * 1000 * 100 << ", " << Eavg << endl;
        }
        
        stream.close();
    }
}

void racun_osc()
{
    ofstream out_0("g_energija_0.dat");
    ofstream out_03("g_energija_0.3.dat");
    ofstream out_1("g_energija_1.dat");
    ofstream out_3("g_energija_3.dat");

    const int I = 1000;
    const int A = 1000;

    const double F = sqrt(10);

    #pragma omp parallel for
    for (int blog = -8; blog < 8; ++blog)
    {
        double beta = pow(F, blog);
        
        cout << "Beta = " << beta << endl;
        
        oscilator(beta, 0, I, A, out_0);
        oscilator(beta, 0.3, I, A, out_03);
        oscilator(beta, 1, I, A, out_1);
        oscilator(beta, 3, I, A, out_3);
    }

    out_0.close();
    out_03.close();
    out_1.close();
    out_3.close();
}

int main(int argc, char **argv)
{
    
    // testConvergenceEpsilon();
    racun_osc();
    // racun(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    
    return 0;
}
