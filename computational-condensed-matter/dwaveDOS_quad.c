#include <stdio.h>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_ellint.h>
#include <mgl2/mgl_cf.h>

#define PI  3.1415926536

double lowlimit(double x, double epsilon)
{
    return (0.5 * acos(x)) + epsilon;
}

double highlimit(double x, double epsilon)
{
    return (0.5 * (PI - acos(x))) - epsilon;
}

double function_to_use(double phi, void *params)
{
    double x = *params;
    double function_to_use = x / sqrt(x**2 - cos(2 * phi)**2
    return function_to_use;
}


double dos(double stepsize, double start, double finish, double epsilon, double eps)
{
    long numsteps = (long) ((finish - start) / stepsize);
    double x[numsteps], f[numsteps];

    gsl_integration_cquad_workspace *w
        = gsl_integration_cquad_workspace_alloc(numsteps);

    double result = 0.0;
    double error = 0.0;
    double alpha = 1.0;
    gsl_function F; /* ADD HERE */
    F.function = &function_to_use;

    long current = start;
    for (long i = 0; i < numsteps; i++) {
        if (abs(current) < 1.0 && current != 0.0) {
            x[i] = current;
            F.params = &current;
            gsl_integration_cquad(&F, lowlimit(current, epsilon), highlimit(current, epsilon)
                0, eps, 1000, w, &result, &error);
            f[i] = result;
        } else if (abs(current) > 1) {
            x[i] = current;
            result = (2.0 / PI) * gsl_sf_ellint_Kcomp(1.0 / abs(current));
            f[i] = result;
        }
        current = current + stepsize;
    }

    mglData xs, fs;
    mglGraph graph;
    xs.Link(x, numsteps);
    fs.Link(f, numsteps);
    graph.Plot(xs, fs);
    graph.WriteFrame("sample.png");
}


int main()
{

    return 0;
}