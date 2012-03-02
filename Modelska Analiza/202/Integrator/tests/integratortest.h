#ifndef INTEGRATORTEST_H
#define INTEGRATORTEST_H

#include <QtTest>

class Integrator;
class IntegratorTest : public QObject
{
  Q_OBJECT
    
private slots:
  void initTestCase();
  void cleanupTestCase();
  void benchmarkCircle();
  void benchmarkSmallEllipse();
  void benchmarkLargeEllipse();
  void benchmarkHyperbola();
  void benchmarkParabola();
  
  void testCircle();
  void testConservation();
  
  void testExtremes();
  
private:
  Integrator* integrator;
  QPair<double,double> interval;
  
  void benchmark(const QVector<double>& state);
};

#endif // INTEGRATORTEST_H
