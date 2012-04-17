
#include "delitev.h"
#include <mgl/mgl.h>
#include <mgl/mgl_eps.h>
#include <mgl/mgl_zb.h>

#include <QtCore/QDebug>

void Delitev::plot(cholmod_dense* a, int d, double k)
{
    int m = tocke.size();
    
    mreal* x = new mreal[m];
    mreal* y = new mreal[m];
    mreal* z = new mreal[m];
    
    double* ax = (double*)a->x;
    
    for (int i = 0; i < m; ++i)
    {
        x[i] = tocke[i].x;
        y[i] = tocke[i].y;
        z[i] = (tocke[i].noter ? ax[not_indeksi[i]] : 0);
    }
    
    double zRange = 0;
    for (int i = 0; i < a->nzmax; ++i)
    {
        zRange = qMax(zRange, fabs(ax[i]));
    }
    
    qDebug() << zRange;
    
    int j = -1;    
    mreal* t = new mreal[3*trikotniki.size()];

    foreach (const Trikotnik& tr, trikotniki)
    {
        t[++j] = tr.i;
        t[++j] = tr.j;
        t[++j] = tr.k;
    }
    
    mglData tData(trikotniki.size(), 3, t);
    mglData xData(m, x);
    mglData yData(m, y);
    mglData zData(m, z);
    
    mglGraphZB g;
  //  g.Light(true);
  //  g.Rotate(50, 60);
    g.Clf();    
    char filename[63];
    g.SetSize(rect.width() * 800, rect.height() * 800);

    
    g.Axis(mglPoint(rect.left(), rect.top(), -zRange), mglPoint(rect.right(), rect.bottom(), zRange));
    
    g.TriPlot(tData, xData, yData, zData, "d");
    g.TriCont(tData, xData, yData, zData, "k#");
    
    
    sprintf(filename, "g_contour_%s_%d.png", name, d);
    g.WritePNG(filename);
    
    
    mglGraphZB h;
    h.Clf();
 //   h.Rotate(30, -45, 0);
    zRange *= 2.5;
    h.Axis(mglPoint(-1, -1, -zRange), mglPoint(1, 1, zRange));
   // g.Light(true);

    h.TriPlot(tData, xData, yData, zData, "d");
    h.TriPlot(tData, xData, yData, zData, "k#");
    
    sprintf(filename, "g_plot_%s_%d.png", name, d);
    qDebug() << "Plotting" << name;
    h.WritePNG(filename);
    
    
    
    delete[] t;
    delete[] x, y, z;
}