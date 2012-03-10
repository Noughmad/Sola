

#include <QMap>
#include <QPair>
#include <QVector>


typedef QVector<double> State;
typedef QMap<double, State> Solution;
typedef int (*YDot)(double t, const double y[], double f[], void* params);
typedef QPair<double, double> Interval;

// Convenience for four parameters
State createState(double x, double y, double u, double v);

class Integrator
{
public:
  virtual ~Integrator();
  
  virtual Solution integrate(YDot dot, const State& initialState, const Interval& interval, double eps, void* param = 0) = 0;
};

class GslIntegrator : public Integrator
{
public:
  virtual Solution integrate(YDot dot, const State& initialState, const Interval& interval, double eps, void* param = 0);
};

double radius(const State& state);
double energy(const State& state);
double momentum(const State& state);

double odstopanje(const Solution& solution, double (*fun)(const State&)); 

double obhodni_cas(const Solution& solution);
double povratek(const Solution& solution);

double z(double t);
void saveToFile(const QString& filename, const Solution& solution, bool zvezda = false);
void saveFrame(const QString& filename, double t, const State& state, const Interval& xrange, const Interval& yrange, bool zvezda = false);