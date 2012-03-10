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

Solution GslIntegrator::integrate(YDot dot, const State& initialState, const Interval& interval, double eps, void* param)
{
  const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk4;
  const size_t dim = initialState.size();

  gsl_odeiv2_step * s = gsl_odeiv2_step_alloc (T, dim);
  gsl_odeiv2_control * c = gsl_odeiv2_control_y_new (0.0, eps);
  gsl_odeiv2_evolve * e = gsl_odeiv2_evolve_alloc (dim);

  gsl_odeiv2_system sys = {dot, 0, dim, param};

  double t = interval.first, t1 = interval.second;
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
  return state[0] * state[3] - state[1] * state[2];
}

double odstopanje(const Solution& solution, double (*fun)(const State&) )
{
  double o = 0;
  double v = fun(solution.value(0.0));
  for (Solution::const_iterator it = solution.constBegin(); it != solution.constEnd(); ++it)
  {
    o = qMax(o, fabs(fun(it.value()) - v));
  }
  return o;
}

double povratek(const Solution& solution)
{
  double lastY = 0;
  double T = M_PI;
  double err = 0;
  
  int n = solution.size();
  
  
  Solution::const_iterator it = solution.constBegin() + 1;
  Solution::const_iterator end = solution.constEnd() - 1;
  
  for (; it != end; ++it)
  {
    double rp = radius((it-1).value());
    double r = radius(it.value());
    double rn = radius((it+1).value());
    if ( (r - rp) * (r - rn) > 0 )
    {
      err = qMax(err, it.value()[1]);
    }
  }
  return err;
}

double obhodni_cas(const Solution& solution)
{
  double lastY = 0;
  double T = M_PI;
  double err = 0;
  Solution::const_iterator it = solution.constBegin();
  Solution::const_iterator end = solution.constEnd();
  
  for (; it != end; ++it)
  {
    const double y = it.value()[1];
    if (y != 0 && lastY * y < 0)
    {
      err = qMax(err, fabs(it.key() - T));
      T += M_PI;
    }
    lastY = y;
  }
  return err;
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
  dataStream.setFieldAlignment(QTextStream::AlignRight);
  dataStream.setFieldWidth(12);
  dataStream.setPadChar(' ');
  
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
  
  dataFile.close();
  
//  qDebug() << xrange << yrange;
  bool animations = false;
  if (animations)
  {
    it = solution.constBegin();
    while (it != end)
    {
      saveFrame(frameName.arg(frame, 3, 10, QChar('0')), it.key(), it.value(), xrange, yrange, zvezda);
      
      ++frame;
      t += FrameInterval;
      it = solution.upperBound(t);
    }
  }
}

void saveFrame(const QString& filename, double t, const State& state, const Interval& xrange, const Interval& yrange, bool zvezda)
{
  const int w = zvezda ? 800 : 400;
  const int h = 400;
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

