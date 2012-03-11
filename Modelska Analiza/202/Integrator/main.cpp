#include <iostream>
#include <iomanip>
#include <math.h>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include "lib/integrator.h"

#include <QtCore/QDebug>
#include <QFile>
#include <QImage>

int func(double t, const double y[], double f[], void* params)
{
    const double r = sqrt(y[0] * y[0] + y[1] * y[1]);
    const double rr = 1.0 / (r*r*r);
    f[0] = y[2];
    f[1] = y[3];
    f[2] = -y[0] * rr;
    f[3] = -y[1] * rr;

    return GSL_SUCCESS;
}

int fun_zoga(double t, const double y[], double f[], void* params)
{
    f[0] = y[1];
    bool odd = (int)(floor(y[0]) + floor(t)) % 2;
    const double k = *(double*)params;
    f[1] = odd ? k : -k;
    return GSL_SUCCESS;
}

double energija2(double t, const State& state)
{
  Q_ASSERT(state.size() == 4);
  State s(4);
  s[0] = state[0] - z(t);
  s[1] = state[1] + 1.5;
  s[2] = state[2] - 2;
  s[3] = state[3];
  return energy(s);
}

int zvezda(double t, const double y[], double f[], void* params)
{
  const double r = sqrt(y[0] * y[0] + y[1] * y[1]);
  const double rr = 1.0 / (r*r*r);
  f[0] = y[2];
  f[1] = y[3];
  f[2] = -y[0] * rr;
  f[3] = -y[1] * rr;
  
  const double dx = y[0] - z(t);
  const double dy = y[1] + 1.5; 
  
  const double R = sqrt(dx * dx + dy * dy);
  const double RR = 1.0 / (R*R*R);
  
  f[2] -= dx * RR;
  f[3] -= dy * RR;
  
  return GSL_SUCCESS;
}

State zacetni(double phi, double v)
{
  State y(4);
  phi = phi - 10 - M_PI_2;
  y[0] = cos(phi);
  y[1] = sin(phi);
  y[2] = -v * sin(phi);
  y[3] = v * cos(phi);
  
  return y;
}

void planetek()
{
    Integrator* integrator = new GslIntegrator();
    Interval i = qMakePair(0.0, 50.0);
    
    /*
    double eps[] = {1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 1e-15, 1e-16 };
    State s(4);
    s[0] = 1.0;
    s[1] = 0.0;
    s[2] = 0.0;
    s[3] = 0.1;
    for (int j = 0; j < 11; ++j)
    {
      Solution sol = integrator->integrate(func, s, i, eps[j]);
      qDebug() << eps[j] << " " << odstopanje(sol, energy) << " " << odstopanje(sol, momentum)
	<< " " << obhodni_cas(sol) << " " << povratek(sol) << " " << sol.size();
     // saveToFile(QString("test-%1").arg(j+6), sol, false);
    }
    
    return 0;
    
    qDebug() << "Naredil primer brez zvezde";
    
    */
    
    for (double phi = 0; phi < 2*M_PI; phi += 0.02)
    {
      const State state = (integrator->integrate(zvezda, zacetni(phi, 1.0), i, 1e-14).constEnd()-1).value();
      qDebug() << phi << " " << energy(state) << " "<< energija2(i.second, state);
      
      // Testing only
      // break;
    }  
}

void omeji(double* q, double *p)
{
  double r = fmod(fabs(*q), 2);
  if (r > 1)
  {
    r = 2 - r;
    *p = - *p;
  }
  *q = r;
}

void portret(double k, const QString& fileName)
{
  
  // TODO: Direktro uporabi GSL, da ne porablja toliko RAMa. 
  Interval i = qMakePair(0.0, 10.0);
  Integrator* integrator = new GslIntegrator;
  
  QImage image(QSize(800,800), QImage::Format_Indexed8);
  for (int j = 0; j < 256; ++j)
  {
    image.setColor(j, qRgb(j, j, j));
  }
  image.fill(255);
  
  gsl_odeiv2_system sys = {fun_zoga, 0, 2, &k};
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-3, 1, 1);
  double y[2];
  double t;
    
  for (double q = 0; q < 1; q += 0.02)
  {
    qDebug() << "q = " << q;
    for (double p = -5; p < 5; p += 0.05)
    {
      y[0] = q;
      y[1] = p;
      t = 0.0;
      double y[] = {q, p};
      while (t < 10.0)
      {
        int status = gsl_odeiv2_driver_apply_fixed_step(driver, &t, 1e-3, 1, y);
        if (status != GSL_SUCCESS)
        {
          qDebug() << t << gsl_strerror(status);
          break;
        }
        double qq = y[0];
        double pp = y[1];
        omeji(&qq, &pp);
        uchar& pix = image.scanLine( qBound<int>(0, (5+pp) * 80.0, 799) )[ qBound<int>(0, qq * 800.0, 799) ];
        if (pix > 0)
        {
          --pix;
        }
      }
    }
  }
  image.save(fileName);
}

void sled(double k, double q, double p, const QString& fileName)
{
  QFile f(fileName);
  f.open(QIODevice::WriteOnly);
  QTextStream stream(&f);
  auto i = new GslIntegrator();
  double y[] = {q, p};
  
  gsl_odeiv2_system sys = {fun_zoga, 0, 2, &k};
  gsl_odeiv2_driver* driver = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk4, 1e-4, 1e-3, 1e-3);
  double t = 0.0;
  while (t < 20.0)
  {
    int c = fmod(t, 2) > 1 ? 1 : 0;
    double q = y[0], p = y[1];
    omeji(&q, &p);
    stream << t << " " << q << " " << p << " " << c << endl;

    int status = gsl_odeiv2_driver_apply_fixed_step(driver, &t, 1e-4, 1, y);
    if (status != GSL_SUCCESS)
    {
      qDebug() << t << gsl_strerror(status);
      break;
    }
  }
  
  f.close();
}

int main(int argc, char **argv) {
    portret(0.1, "g_potrtret_1.png");
    portret(0.5, "g_potrtret_5.png");
    portret(1.0, "g_potrtret_10.png");
    portret(5.0, "g_potrtret_50.png");
    portret(10.0, "g_potrtret_100.png");
    portret(30.0, "g_potrtret_300.png");
  
    sled(1.0, 0.5, 1.0, "g_sled_10_5_10.dat");
    sled(1.0, 0.5, 0.5, "g_sled_10_5_5.dat");
    sled(1.0, 0.5, 0.0, "g_sled_10_5_0.dat");
    
    sled(0.5, 0.5, 1.0, "g_sled_5_5_10.dat");
    sled(0.5, 0.5, 0.5, "g_sled_5_5_5.dat");
    sled(0.5, 0.5, 0.0, "g_sled_5_5_0.dat");
  
    sled(0.5, 0.5, -1.0, "g_sled_5_5_-10.dat");
    sled(0.5, 0.5, -0.5, "g_sled_5_5_-5.dat");
    sled(1.0, 0.5, -1.0, "g_sled_10_5_-10.dat");
    
    sled(10.0, 0.5, 0.0, "g_sled_100_5_0.dat");
    sled(10.0, 0.5, 2.0, "g_sled_100_5_20.dat");

    return 0;
}
