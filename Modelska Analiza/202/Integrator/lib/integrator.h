

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
  
  virtual Solution integrate(YDot dot, const State& initialState, const Interval& interval) = 0;
};

class GslIntegrator : public Integrator
{
  virtual Solution integrate(YDot dot, const State& initialState, const Interval& interval);
};

double radius(const State& state);
double energy(const State& state);
double momentum(const State& state);