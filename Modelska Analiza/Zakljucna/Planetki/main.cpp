#include <iostream>
#include <QImage>

#include "planets.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc < 5)
    {
        qDebug() << "Usage: planetki mu delta total_time steps [ trajectory_type ]" << endl;
        return 0;
    }

    double mu = atof(argv[1]);
    double delta = atof(argv[2]);
    Planets p(mu, delta);

    double T = atof(argv[3]);
    int N = atoi(argv[4]);
    p.dt = T/N;

    QString type = (argc > 5) ? argv[5] : "direct";
    Trajectory t;
    if (type == "direct")
    {
        t = p.direct_route(N);
    }
    else if (type == "spiral")
    {
        t = p.ellipse(N);
    }
    else if (type == "spiral_imp")
    {
        t = p.ellipse(N, atoi(argv[3])/20);
    }
    else if (type == "spline")
    {
        t = p.spline(N);
    }
    else
    {
        qDebug() << "Unknown orbit type: " << type << endl;
        qDebug() << "Recognized orbit types are 'direct', 'spiral' and 'spline'" << endl;
        return 0;
    }
    
    Q_ASSERT(t.positions[0] == p.planet_one(0));
    Q_ASSERT(t.positions[N] == p.planet_two(T));

    double change;
    int iter = 0;
    p.chebishev = 1;
    do
    {
        ++iter;
        change = p.relax_step(t);
    }
    while (change > 1e-8 && iter < 1e6);
    p.plot(t).save("g_plot_after.png");
    
    QPair<double, double> b = p.burst(t);
    qDebug() << T << iter << b.first << b.second << b.first + b.second;

    QString plot = QString("g_plot_%1_%2_%3_%4.png")
                   .arg(mu).arg(delta).arg(type).arg(T);
    p.plot(t).save(plot);
    
    return 0;
}
