#include <iostream>
#include <fstream>
#include <math.h>
#include <sstream>

#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_histogram.h>

#include <QtOpenCL/QtOpenCL>
#include <QtOpenCL/qclcontext.h>
#include <CL/cl_platform.h>

#define l(x) y[x-1]
#define n(x) y[x+2]

using namespace std;

const int PoincareIndex = 4;
const double dE = 0.001;

const int ParallelRuns = 100;

QCLContext context;
QCLProgram program;

struct top_params
{
    top_params() {}
    top_params(double l)
    {
        D = 5.0;
        a = 1.0;
        mg = 1.0;
        L = l;
        R = 1.0 - 1.0/D;
    }
    top_params(double D, double a, double L, double mg)
    : D(D), a(a), L(L), mg(mg)
    {
        R = 1.0 - 1.0/D;
    }
    
    float D;
    float a;
    float L;
    float mg;
    float R;
};

inline double energy(const top_params& top, double y[])
{
    return 0.5 * (y[0]*y[0] + y[1]*y[1] + y[2]*y[2]/top.D) + top.mg * (top.L * n(1) + top.a * n(3));
}

void random_state(double y[], const top_params& top, double emax);

void lyapunov(const top_params& top, const double initial[], double* exponent, double* sigma);

inline double interpolate(double t1, double t2, double x1, double x2, double t)
{
    return ((t-t1)*x2 + (t2-t)*x1) / (t2-t1);
}

int odvod(double t, const double y[], double dy[], void* params)
{
    QCLKernel& kernel = *(QCLKernel*)params;

    QCLBuffer yBuffer = context.createBufferCopy(y, 6*ParallelRuns*sizeof(double), QCLMemoryObject::ReadOnly);
    kernel.setArg(6, yBuffer);

    QCLBuffer dyBuffer = context.createBufferDevice(6*ParallelRuns*sizeof(double), QCLMemoryObject::WriteOnly);
    kernel.setArg(7, dyBuffer);

    QCLEvent event = kernel.run();
    event.waitForFinished();

    dyBuffer.read(dy, 6*ParallelRuns*sizeof(double));

    return GSL_SUCCESS;
}

class TopWorkspace
{
public:
    TopWorkspace(const top_params& top, const string& file);
    TopWorkspace(double D, double a, double L);
    ~TopWorkspace();

    void setInitial(const double initial[]);

    double operator-(const TopWorkspace& other);

    void poincare();
    void save();

    inline void apply()
    {
        for (int i = 0; i < 6*ParallelRuns; ++i)
        {
            last[i] = y[i];
        }
        runBuffer = context.createBufferCopy(run, ParallelRuns*sizeof(cl_char), QCLMemoryObject::ReadOnly);
        kernel.setArg(0, runBuffer);
        gsl_odeiv2_driver_apply(driver, &t, t+t1, y);
    }

    inline double energy()
    {
        return ::energy(params, y);
    }

public:
    double y[6*ParallelRuns];
    double t;

private:
    void init(const string& name);
    
    gsl_odeiv2_system system;
    gsl_odeiv2_driver* driver;
    top_params params;
    double t1;
    double last[6*ParallelRuns];
    cl_char run[ParallelRuns];
    ofstream output;

public:
    QCLKernel kernel;
    QCLBuffer runBuffer, yBuffer, dyBuffer;
};

TopWorkspace::TopWorkspace(const top_params& top, const string& file)
{
    params = top;
    init("g_" + file + ".dat");
}


TopWorkspace::TopWorkspace(double D, double a, double L)
{
    params = {D, a, L, 1.0};
    init("g_vrtavka.dat");
}

void TopWorkspace::init(const string& name)
{
    kernel = program.createKernel("odvod");
    QCLWorkSize size(ParallelRuns);
    kernel.setGlobalWorkSize(size);

    cout << kernel.localWorkSize().toString().toStdString() << endl;

    kernel.setArg(1, params.D);
    kernel.setArg(2, params.L);
    kernel.setArg(3, params.mg);
    kernel.setArg(4, params.a);
    kernel.setArg(5, params.R);

    system = {odvod, 0, 6*ParallelRuns, &kernel};
    driver = gsl_odeiv2_driver_alloc_y_new(&system, gsl_odeiv2_step_rk4, 1e-3, 1e-12, 0);

    for (int i = 0; i < ParallelRuns; ++i)
    {
        run[i] = 1;
    }

    t1 = 1e-3;
    output.open(name);
}


TopWorkspace::~TopWorkspace()
{
    gsl_odeiv2_driver_free(driver);
    output.close();
}

void TopWorkspace::setInitial(const double initial[])
{
    for (int i = 0; i < 6*ParallelRuns; ++i)
    {
        y[i] = initial[i];
    }
}

double TopWorkspace::operator-(const TopWorkspace& other)
{
    double d = 0;
    for (int i = 0; i < 6; ++i)
    {
        d += (y[i] - other.y[i]) * (y[i] - other.y[i]);
    }
    return d;
}

void TopWorkspace::save()
{
    output << t;
    for (int i = 0; i < 6; ++i)
    {
        output << " " << y[i];
    }
    output << endl;
}

void TopWorkspace::poincare()
{
    const int p = PoincareIndex;

    for (int i = 0; i < ParallelRuns; ++i)
    {
        run[i] == 0;
    }

    int running = ParallelRuns;
    for (int i = 0; i < ParallelRuns; ++i)
    {
        run[i] = 1;
    }

    apply();
    apply();

    while (running > 0)
    {
        for (int i = 0; i < ParallelRuns; ++i)
        {
            if (run[i] && last[6*i+p] * y[6*i+p] <= 0)
            {
                double tS = t - t1 * fabs(y[6*i+p]) / (fabs(y[6*i+p]) + fabs(last[6*i+p]));
                for (int k = 0; k < 6; ++k)
                {
                    y[6*i+k] = interpolate(t-t1, t, last[6*i+k], y[6*i+k], tS);
                }
                run[i] = 0;
                --running;
                cout << running << endl;
            }
        }
        apply();
    }
}

void lyapunov(const top_params& top, const double initial[], double* exponent, double* sigma)
{
    double factor = 1;
    const double step = 1e-6;

    const int N = 200;
    const int P = 2;

    double z[6*ParallelRuns];
    for (int i = 0; i < 6*ParallelRuns; ++i)
    {
        z[i] = initial[i];
    }

    context.create();
    program = context.buildProgramFromSourceFile("../odvod.cl");

    TopWorkspace top1(top.D, top.a, top.L);
    top1.setInitial(z);

    for (int i = 0; i < 6; ++i)
    {
        z[i] += (double)rand() / RAND_MAX * step;
    }
    

    TopWorkspace top2(top.D, top.a, top.L);
    top2.setInitial(z);

    for (int p = 0; p < P; ++p)
    {
        top1.poincare();
        top2.poincare();
    }

    QCLKernel diff = program.createKernel("razdalja");
    diff.setGlobalWorkSize(QCLWorkSize(ParallelRuns));

    double d0[ParallelRuns];
    QCLBuffer buf1 = context.createBufferCopy(top1.y, 6*ParallelRuns*sizeof(double), QCLMemoryObject::ReadOnly);
    QCLBuffer buf2 = context.createBufferCopy(top2.y, 6*ParallelRuns*sizeof(double), QCLMemoryObject::ReadOnly);
    QCLBuffer buf3 = context.createBufferHost(d0, ParallelRuns*sizeof(double), QCLMemoryObject::WriteOnly);
    diff.setArg(0, buf1);
    diff.setArg(1, buf2);
    diff.setArg(2, buf3);
    diff.run().waitForFinished();

    double** a = new double*[ParallelRuns];
    for (int i = 0; i < ParallelRuns; ++i)
    {
        a[i] = new double[N];
    }


    for (int k = 0; k < N; ++k)
    {
        for (int p = 0; p < P; ++p)
        {
            top1.poincare();
            top2.poincare();
        }

        cout << "One Poincare mapping done: " << k << endl;

        double d[ParallelRuns];
        buf1 = context.createBufferCopy(top1.y, 6*ParallelRuns*sizeof(double), QCLMemoryObject::ReadOnly);
        buf2 = context.createBufferCopy(top2.y, 6*ParallelRuns*sizeof(double), QCLMemoryObject::ReadOnly);
        buf3 = context.createBufferHost(d, ParallelRuns*sizeof(double), QCLMemoryObject::ReadWrite);
        diff.setArg(0, buf1);
        diff.setArg(1, buf2);
        diff.setArg(2, buf3);
        diff.run().waitForFinished();

        for (int i = 0; i < ParallelRuns; ++i)
        {
            a[i][k] = sqrt(d0[i]/d[i]);
            for (int j = 0; j < 6; ++j)
            {
                top2.y[6*i+j] = (1.0-a[i][k]) * top1.y[6*i+j] + a[i][k] * top2.y[6*i+j];
            }

            a[i][k] = -log(a[i][k]);
        }
    }

    double* x = new double[N];
    for (int i = 0; i < N; ++i)
    {
        x[i] = 1.0/(i+1);
    }

    for (int p = 0; p < ParallelRuns; ++p)
    {

        double c0, c1, cov00, cov01, cov11, sumsq;
        gsl_fit_linear(x, 1, a[p], 1, N, &c0, &c1, &cov00, &cov01, &cov11, &sumsq);

        exponent[p] = c0;
        sigma[p] = sqrt(cov00);
    }

    delete[] x;
    for (int i = 0; i < N; ++i)
    {
        delete[] a[i];
    }
    delete[] a;
}

void random_state(double y[], const top_params& top, double emax)
{
    do
    {
        for (int i = 1; i < 4; ++i)
        {
            l(i) = 20 * (2 * (double)rand() / RAND_MAX - 1);
            n(i) = (double)rand() / RAND_MAX;
        }
        double nn = 1.0/sqrt(n(1)*n(1) + n(2)*n(2) + n(3)*n(3));
        for (int i = 1; i < 4; ++i)
        {
            n(i) *= nn;
        }
    }
    while (energy(top, y) > emax);
}

void chaos_part(const top_params& top, double emax)
{
    double y[6*ParallelRuns];

    gsl_histogram* h = gsl_histogram_alloc(14);
    gsl_histogram_set_ranges_uniform(h, -1, 6);

    double exponent[ParallelRuns];
    double sigma[ParallelRuns];
    for (int i = 0; i < ParallelRuns; ++i)
    {
        random_state(&(y[6*i]), top, emax);
    }
    
    lyapunov(top, y, exponent, sigma);

    for (int i = 0; i < ParallelRuns; ++i)
    {
        gsl_histogram_increment(h, exponent[i]/sigma[i]);
    }

    char buf[32];
    sprintf(buf, "g_histogram_%g_%g.dat", emax, top.L);
    FILE* f = fopen(buf, "wt");
    gsl_histogram_fprintf(f, h, "%g", "%g");
    fclose(f);
}

void fazni_prostor(double lambda)
{
    top_params top(lambda);
    stringstream name;
    name << "vrtavka_" << lambda;
    TopWorkspace w(top, name.str());
    for (int i = 0; i < 20; ++i)
    {
        double y[6];
        random_state(y, top, 10.0);
        w.setInitial(y);

        for (int k = 0; k < 2000; ++k)
        {
            w.poincare();
            w.save();
        }

        cout << "Naredil zacetni pogoj " << i << endl;
    }
}

int main(int argc, char **argv) {
    srand(time(0));
    double erg = atof(argv[1]);
    for (double L = 0; L < 2.05; L += 0.1)
    {
        cout << "Starting with lambda=" << L << endl;
        chaos_part(top_params(L), erg);
    }
    return 0;
}
