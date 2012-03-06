#include "integrator.h"

#include <iostream>
#include <iomanip>
#include <math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>

#include <QtCore/QFile>
#include <QtCore/QTextStream>
#include <QtCore/QDebug>
#include <qdir.h>
#include <QtGui/QImage>
#include <QtGui/QPainter>

const double FrameInterval = 0.1;

State createState(double x, double y, double u, double v)
{
  State s(4);
  s[0] = x;
  s[1] = y;
  s[2] = u;
  s[3] = v;
  return s;
}

double z(double t)
{
  return 2*t - 20; 
}

Integrator::~Integrator()
{

}

Solution GslIntegrator::integrate(YDot dot, const State& initialState, const Interval& interval)
{
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;
  const size_t dim = initialState.size();

  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, dim);
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (0.0, 1e-14);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (dim);

  gsl_odeiv2_system sys = {dot, 0, dim, 0};

  double t = 0.0, t1 = 40.0;
  double h = 1e-6;
  double* y = (double*)malloc(sizeof(double) * initialState.size());
  for (int i = 0; i < initialState.size(); ++i)
  {
    y[i] = initialState[i];
  }
  
  Solution solution;
  solution[t] = initialState;

  while (t < t1)
  {
      int status = gsl_odeiv2_evolve_apply (e, c, s, &sys, &t, t1, &h, y);

      if (status != GSL_SUCCESS)
	  break;
      
      State s(dim);
      for (int i = 0; i < dim; ++i)
      {
	s[i] = y[i];
      }
      
      solution[t] = s;
  }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);
  delete[] y;
  
  return solution;
}

double radius(const State& state)
{
  return sqrt(state[0] * state[0] + state[1] * state[1]);
}

double energy(const State& state)
{
  // TODO: Normalization
  // both parts should come with different constant factors
  return 0.5*(state[2] * state[2] + state[3] * state[3]) - 1.0/radius(state);
}

double momentum(const State& state)
{
  // TODO: Calculate
  return 0;
}

void saveToFile(const QString& filename, const Solution& solution, bool zvezda)
{
  const QString dataName = "g_" + filename + ".dat";
  const QString folderName = filename;
  const QString frameName = folderName + "/g_frame_%1.png";
  QDir::current().mkdir(folderName);
  
  int frame = 0;
  double t = 0;
  double nextT = FrameInterval;
  QFile dataFile(dataName);
  dataFile.open(QIODevice::WriteOnly);
  QTextStream dataStream(&dataFile);
  
  Solution::const_iterator it = solution.constBegin();
  const Solution::const_iterator end = solution.constEnd();
  
  Interval xrange = qMakePair(0.0, 0.0);
  Interval yrange = qMakePair(0.0, 0.0);
  
  if (zvezda)
  {
    xrange.first = -3;
    xrange.second = 3;
    yrange.first = -1.6;
  }
  
  for (; it != end; ++it)
  {
    dataStream << it.key() << it.value()[0] << it.value()[1] << it.value()[2] << it.value()[3];
    if (zvezda)
    {
      dataStream << z(it.key());
    }
    dataStream << endl;
    xrange.first = qMin(xrange.first, it.value()[0]);
    xrange.second = qMax(xrange.second, it.value()[0]);
    yrange.first = qMin(yrange.first, it.value()[1]);
    yrange.second = qMax(yrange.second, it.value()[1]);
  }
  
  xrange.first = qMax(xrange.first, -5.0);
  xrange.second = qMin(xrange.second, 5.0);
  yrange.first = qMax(yrange.first, -3.0);
  yrange.second = qMin(yrange.second, 3.0);
  
  dataFile.close();
  
//  qDebug() << xrange << yrange;
  it = solution.constBegin();
  while (it != end)
  {
    saveFrame(frameName.arg(frame, 3, 10, QChar('0')), it.key(), it.value(), xrange, yrange, zvezda);
    
    ++frame;
    t += FrameInterval;
    it = solution.upperBound(t);
  }
}

void saveFrame(const QString& filename, double t, const State& state, const Interval& xrange, const Interval& yrange, bool zvezda)
{
  const int w = zvezda ? 400 : 200;
  const int h = 200;
  QImage image(w, h, QImage::Format_ARGB32);
  image.fill(Qt::black);
  QPainter painter(&image);
  painter.setPen(Qt::NoPen);
  double x,y;
  const double dx = w / (xrange.second - xrange.first);
  const double dy = h / (yrange.second - yrange.first);
  
  // Draw the Sun
  painter.setBrush(QBrush(Qt::yellow));
  x = (- xrange.first) * dx;
  y = (- yrange.first) * dy;
  painter.drawEllipse(QPointF(x, y), 5.0, 5.0);
  
  // Draw the planet
  painter.setBrush(QBrush(Qt::green));
  x = (state[0] - xrange.first) * dx;
  y = (state[1] - yrange.first) * dy;
  painter.drawEllipse(QPointF(x, y), 3.0, 3.0);
  
  // Draw the passing star
  if (zvezda)
  {
    painter.setBrush(QBrush(Qt::red));
    x = (z(t) - xrange.first) * dx;
    y = (-1.5 - yrange.first) * dy;
    painter.drawEllipse(QPointF(x, y), 5.0, 5.0);
  }
  
  image.save(filename);
}

