#include "integratortest.h"
#include "../lib/integrator.h"

#include <gsl/gsl_errno.h>

#define ERROR(actual, expected)			\
(expected ? fabs(expected/actual-1) : fabs(expected) )

#define COMPARE_DOUBLE(actual, expected) 	\
if (ERROR(actual,expected) > 1e-8) { QCOMPARE(actual, expected); } else { QVERIFY(true); }

int func(double t, const double y[], double f[], void* params)
{
    const double r = sqrt(y[0] * y[0] + y[1] * y[1]);
    const double rr = 1.0 / (r*r*r);
    f[0] = y[2];
    f[1] = y[3];
    f[2] = -y[0] * rr;
    f[3] = -y[1] * rr;
    
    // std::cout << "Evaluating " << y[0] << " " << y[1] << " " << y[2] << " " << y[3]  << std::endl;

    return GSL_SUCCESS;
}

void IntegratorTest::benchmark(const State& state)
{
  QBENCHMARK(integrator->integrate(func, state, interval, 1e-12));
}


void IntegratorTest::initTestCase()
{
  integrator = new GslIntegrator;
  interval = qMakePair(0.0, 20.0);
}

void IntegratorTest::cleanupTestCase()
{
  delete integrator;
}

void IntegratorTest::benchmarkCircle()
{
  benchmark(createState(1.0, 0.0, 0.0, 1.0));
}

void IntegratorTest::benchmarkLargeEllipse()
{
  benchmark(createState(1.0, 0.0, 0.0, 1.2));
}

void IntegratorTest::benchmarkSmallEllipse()
{
  benchmark(createState(1.0, 0.0, 0.0, 0.2));
}

void IntegratorTest::benchmarkParabola()
{
  benchmark(createState(1.0, 0.0, 0.0, sqrt(2.0)));
}

void IntegratorTest::benchmarkHyperbola()
{
  benchmark(createState(1.0, 0.0, 0.0, 2));
}

void IntegratorTest::testCircle()
{
  State state = createState(1.0, 0.0, 0.0, 1.0);
  Solution solution = integrator->integrate(func, state, interval, 1e-12);
  int i = 0;
  int n = solution.size() / 100;
  if (n == 0)
  {
    n = 1;
  }
  double e = energy(state);
  double l = momentum(state);
  
  qDebug() << solution.size() << n;
  foreach (const State& s, solution)
  {
    if (i % n == 0)
    {
      COMPARE_DOUBLE(radius(s), 1.0);
      COMPARE_DOUBLE(energy(s), e);
      COMPARE_DOUBLE(momentum(s), l);
    }
    ++i;
  }
}

void IntegratorTest::testConservation()
{
  State state = createState(1.0, 0.0, 0.0, 0.1);
  Solution solution = integrator->integrate(func, state, interval, 1e-12);
  int i = 0;
  int n = solution.size() / 100;
  if (n == 0)
  {
    n = 1;
  }
  double e = energy(state);
  double l = momentum(state);
  
  qDebug() << solution.size() << n;
  foreach (const State& s, solution)
  {
    if (i % n == 0)
    {
      COMPARE_DOUBLE(energy(s), e);
      COMPARE_DOUBLE(momentum(s), l);
    }
    ++i;
  }
}

void IntegratorTest::testExtremes()
{
  // TODO: Preveri, pri kaksne kotu dobimo najvecjo in najmanjso oddaljenost
}


#include "integratortest.moc"

QTEST_MAIN(IntegratorTest)


