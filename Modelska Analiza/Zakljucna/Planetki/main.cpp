#include <iostream>
#include <QImage>

#include "planets.h"

using namespace std;

void test()
{
    Planets p(1.0, 0.02);
    p.dt = 0.025;
    p.phase = 0.5;
    Trajectory t = p.ellipse(500, 2);
    double change;
    int iter = 0;
    p.chebishev = 1.9;
    do 
    {
        ++iter;
        change = p.relax_step(t);
    }
    while (change > 1e-8 && iter < 1e5);
    qDebug() << "Final iteration:" << change << iter;
    
    QImage image = p.plot(t);
    image.save("g_plot_circles.png");
}

int main(int argc, char **argv) {
    test();
    return 0;
}
