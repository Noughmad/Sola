#include <iostream>
#include <math.h>
#include <QtGui/QImage>
#include <QtCore/QDebug>

double trk(double a, double b, double c, double* rD)
{
  const double D = b*b - 4*a*c;
  if (D < 0)
  {
    return 2;
  }
  *rD = sqrt(D);
  
  const double t1 = (-b+*rD)/2/a;
  const double t2 = (-b-*rD)/2/a;
    
  if (t1 > 0 && t2 > 0)
  {
    // Oba pozitivna, vrnemo manjsega
    return qMin(t1,t2);
  }
  else if (t1 > 0)
  {
    return t1;
  }
  else if (t2 > 0)
  {
    return t2;
  }
  else
  {
    return 2;
  }
}

int naslednja(double t, double y[], double polje)
{
  Q_ASSERT(t >= 0);
  Q_ASSERT(t <= 1);
  double a = polje * 0.5;
  double b = y[1];
  double c = y[0];
  
  if (qFuzzyCompare(t,1))
  {
    return 0;
  }
    
  if ((qFuzzyCompare(y[0],1) && b > 0) || (qFuzzyIsNull(y[0]) && b < 0))
  {
    y[1] = -y[1];
    return naslednja(t,y,polje);
  }
  
  // parabola: q(t) = at^2 + bt + c
  
  // Kdaj se bomo zaleteli v spodnji rob (q = 0)
  double D1,D2;
  const double t1 = trk(a,b,c,&D1);
  
  // Pa v zgornjega (q = 1)
  const double t2 = trk(a,b,c-1,&D2);
  
  // Koliko casa se imamo do obrata polja
  const double tl = 1-t;
  
  Q_ASSERT(t1 > 0);
  Q_ASSERT(t2 > 0);
  Q_ASSERT(tl > 0);
  
 // qDebug() << "Casi:" << t1 << t2 << tl;
  if (t1 <= t2 && t1 <= tl)
  {
    // Najprej se bomo zaleteli v spodnji rob
    y[0] = 0;
    y[1] = D1;
    return naslednja(t + t1, y, polje);
  }
  else if (t2 <= t1 && t2 <= tl)
  {
    // Najprej se bomo zaleteli v zgornji rob
    y[0] = 1;
    y[1] = -D2;
    return naslednja(t + t2, y, polje);
  }
  else
  {
    // Ne bomo se zaleteli, premikamo se po paraboli do obrata polja
    y[0] = a*tl*tl + b*tl + c;
    y[1] = 2*a*tl + b;
    Q_ASSERT(y[0] >= 0);
    Q_ASSERT(y[0] <= 1);
    return 0;
  }
}

void portret(double q, double p, double polje, QImage* image)
{
  double y[2] = {q, p};
  for (int i = 0; i < 1000; ++i)
  {
    naslednja(0, y, polje);
    naslednja(0, y, -polje);
    
    int qq = qBound<int>(0, y[0] * 800.0, 799);
    int pp = qBound<int>(0, (y[1] + 5) * 80, 799);
    uchar* line = image->scanLine(pp);
    if (line[qq] > 3)
    {
      line[qq] -= 3;
    }
  }
}

void vse(double polje, const QString& fileName)
{
  QImage image(QSize(800,800), QImage::Format_Indexed8);
  for (int j = 0; j < 256; ++j)
  {
    image.setColor(j, qRgb(j, j, j));
  }
  image.fill(255);
  for (double q = 0.005; q < 1; q += 0.01)
  {
    for (double p = -0; p < 5; p += 0.05)
    {
      portret(q,p,polje,&image);
      qDebug() << "Narisal za" << q << p;
    }
  }
  image.save(fileName);
}


int main(int argc, char **argv) {
    vse(0.1, "g_portret_0.1.png");
    vse(1, "g_portret_1.png");
    vse(3, "g_portret_3.png");
    vse(10, "g_portret_10.png");
    return 0;
}
