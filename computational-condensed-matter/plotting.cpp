#include <matplotlibcpp.h>
#include <cmath>

namespace plt = matplotlibcpp;

int main() 
{
    // Prepare data.
    int n = 5000;
    std::vector<double> x(n,0), y(n,0), z(n,0), w(n,2);
    for(int i=0; i<n; i++) {
        x[i] = i*i;
        y[i] = sin(2*M_PI*i/360.0);
        z[i] = log(i);
    }

    // Set the size of output image to 1200x780 pixels
    // plt::figure_size(1200, 780);
    // // Plot line from given x and y data. Color is selected automatically.
    plt::plot(x, y);
    // // Plot a red dashed line from given x and y data.
    plt::plot(x, w,"r--");
    // // Plot a line whose name will show up as "log(x)" in the legend.
    plt::named_plot("log(x)", x, z);
    // // Set x-axis to interval [0,1000000]
    plt::xlim(0, 1000*1000);
    // // Add graph title
    plt::title("Sample figure");
    // // Enable legend.
    plt::legend();
    // Show the image
    plt::show();
}

/*
 * Compile with:
 *  g++ -I /usr/include/python3.6m -o prog_name plotting.cpp -lm -lpython3.6m
 */