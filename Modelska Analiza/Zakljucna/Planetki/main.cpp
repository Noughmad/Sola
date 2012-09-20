#include <iostream>
#include <QImage>
#include <QPainter>
#include <QSvgGenerator>

#include "planets.h"

using namespace std;

int main(int argc, char **argv)
{
    if (argc < 5)
    {
        qDebug() << "Usage: planetki mu delta total_time steps [ trajectory_type ]" << endl;
        return 0;
    }

    if (strcmp(argv[1], "plot") == 0)
    {
        double delta = atof(argv[2]);
        double T = atof(argv[3]);
        int N = atoi(argv[4]);
        int s = atoi(argv[5]);

        QSvgGenerator image;
        // image.setViewBox(QRectF(0,0  sys, 800, 800));
        image.setSize(QSize(800,800));
        QString file("g_plot_%1_%2.svg");
        image.setFileName(file.arg(round(delta)).arg(T));

        QPainter painter(&image);
        painter.setRenderHints(QPainter::Antialiasing, QPainter::HighQualityAntialiasing);
        Vector star;
        star.x = star.y = 0;

        painter.setPen(Qt::NoPen);
        painter.setBrush(Qt::darkYellow);
        painter.drawEllipse(star, 15, 15);
        
        painter.setBrush(Qt::blue);
        painter.drawEllipse(from_polar(1, 0), 5, 5);
        
        painter.setBrush(Qt::red);
        painter.drawEllipse(from_polar(2, invsqrt8*T + delta), 5, 5);
        
        painter.setBrush(Qt::NoBrush);
        painter.setPen(QPen(Qt::DashLine));
        painter.drawEllipse(star, 100, 100);
        painter.drawEllipse(star, 200, 200);

        double MU[] = {0.0, 0.01, 0.03, 0.1};
        for (int i = 0; i < 4; ++i)
        {
            qDebug() << i << MU[i];
            qDebug() << N << s;
            Planets p(MU[i], delta);
            p.dt = T/N;
            Trajectory t = p.toolbox(N, s);
            
            int iter = 0;
            if (!qFuzzyIsNull(MU[i]))
            {
                double change;
                p.chebishev = 1.0;
                for (iter = 0; iter < 3; ++iter)
                {
                    p.relax_step(t);
                }
                p.chebishev = 2.0/(1 + M_PI/N);
                do
                {
                    ++iter;
                    change = p.relax_step(t);
                }
                while (change > 1e-14 && iter < 1e5);

                
                int f = 250 + 80 * log10(MU[i]);
                painter.setPen(QColor(qRgb(f, 0, 0)));
                
                for (int i = 0; i < N; ++i)
                {
                    painter.drawLine(t.positions[i], t.positions[i+1]);
                }
            }
        }

        painter.end();
    }

    double mu = atof(argv[1]);
    double delta = atof(argv[2]);
    Planets p(mu, delta);

    double T = atof(argv[3]);
    int N = atoi(argv[4]);
    int s = atoi(argv[5]);
    p.dt = T/N;

    Trajectory t = p.toolbox(N, s);
    int iter = 0;

    if (!qFuzzyIsNull(mu))
    {
        double change;
        p.chebishev = 1.0;
        for (iter = 0; iter < 3; ++iter)
        {
            p.relax_step(t);
        }
        p.chebishev = 2.0/(1 + M_PI/N);
        do
        {
            ++iter;
            change = p.relax_step(t);
        }
        while (change > 1e-8 && iter < 1e5);
    }
    // p.plot(t).save("g_plot_after.png");
    
    QPair<double, double> b = p.burst(t);
    qDebug() << T << iter << b.first << b.second << p.adjust_speed(b.first) + p.adjust_speed(b.second);

    /*
    QString plot = QString("g_plot_%1_%2_%3_%4.png")
                   .arg(mu).arg(delta).arg(s).arg(T);
    p.plot(t).save(plot);
    */
    
    return 0;
}
