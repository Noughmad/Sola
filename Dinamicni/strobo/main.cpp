#include <iostream>
#include <math.h>
#include <QtGui/QImage>
#include <QtCore/QDebug>

int naslednja(double t, double y[], double polje)
{
  Q_ASSERT(t >= 0);
  double a = polje * 0.5;
  double b = y[1];
  double c = y[0];
  
  // parabola: q(t) = at^2 + bt + c
  
  // Kdaj se bomo zaleteli v spodnji rob (q = 0)
  const double D1 = b*b-4*a*c;
  double t1 = (D1 > 0) ? qMax( (-b + D1)/2/a, (-b-D1)/2/a) : 2;
  t1 = (t1 > 0) ? t1 : 2;
  
  // Pa v zgornjega (q = 1)
  const double D2 = b*b-4*a*(c-1);
  double t2 = (D2 > 0) ? qMax( (-b + D1)/2/a, (-b-D1)/2/a) : 2;
  t2 = (t2 > 0) ? t2 : 2;
  
  // Koliko casa se imamo do obrata polja
  const double tl = 1-t;
  
  if (t1 < t2 && t1 < tl)
  {
    // Najprej se bomo zaleteli v spodnji rob
    y[0] = 0;
    y[1] = sqrt(D1);
    return naslednja(t + t1, y, polje);
  }
  else if (t2 < t1 && t2 < tl)
  {
    // Najprej se bomo zaleteli v zgornji rob
    y[0] = 1;
    y[1] = -sqrt(D2);
    return naslednja(t + t2, y, polje);
  }
  else
  {
    // Ne bomo se zaleteli, premikamo se po paraboli do obrata polja
    y[0] = a*tl*tl + b*tl + c;
    y[1] = 2*a*tl + b;
    return 0;
  }
}

void portret(double q, double p, double polje, QImage* image)
{
  double y[2] = {q, p};
  double t = 0;
  while (t < 2000)
  {
    naslednja(t, y, polje);
    polje *= -1;
    t += 1;
    
    int qq = qBound<int>(0, y[0] * 800.0, 799);
    int pp = qBound<int>(0, (y[1] + 5) * 80, 799);
    uchar* line = image->scanLine(pp);
    if (line[qq] > 5)
    {
      qDebug() << qq << pp << line[qq];
      line[qq] -= 5;
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
  for (double q = 0.0; q < 1; q += 0.01)
  {
    for (double p = -5; p < 5; p += 0.1)
    {
      portret(q,p,polje,&image);
    }
  }
  image.save(fileName);
}


int main(int argc, char **argv) {
    vse(1, "g_portret_1.png");
    return 0;
}
