/*
 * C++ operates almost exactly like C. However, OOP was added as well as a bunch of
 * type-safety and other very powerful methods. I recommend referring to the c_example
 * for the basics of C before continuing. C++ filename extensions are .h/.hpp for headers
 * and .cc/.cpp for normal C++ files.
 */

    
#include <iostream>
#include <math.h>


// Let's create a class for our integration
// Usually, class declarations are in a header file and definition happens in a cc/cpp file.
// They can then be included as '#include "Name of header"'
// However, we will do both at the same time here.
class Integration
{
public:
    // Public variables/methods are accessible by anything

    // This is a class constructor
    Integration(double low, double high, long steps);

    /*
     * Custom destructors are declared as:
     * 
     * ~Integration(...)
     */

    // Just in case we want a way to change limits or steps along the way
    void change_limits(double low, double high);
    void change_steps(long steps);
    double integrate(double (*func)(double x));


private:
    // Private variables/methods are only accessible by the class itself

    // Limits: Initialize the limits and steps
    double limits[2] = {0, 0};
    long steps = 0;
};

/*
 * Here, we must DEFINE the funcitons in the class Integration.
 * We use the scope-definition operator '::' as:
 * 
 *      Class_Name::Function_From_Class(...) {do stuff};
 * 
 * This means 'access the function that is defined in the scope of class.'
 */

Integration::Integration(double low, double high, long steps)
{
    // 'this' acts like Python's 'self,' meaning that we want to access only THIS class's
    // member/method.
    this->limits[0] = low;
    this->limits[1] = high;
    this->steps = steps;

    if (this->steps % 2 != 0) {
        // ++ increments by one
        this->steps++;
    }
}


void Integration::change_limits(double low, double high)
{
    this->limits[0] = low;
    this->limits[1] = high;
}

void Integration::change_steps(long steps)
{
    this->steps = steps;
}

double Integration::integrate(double (*func)(double x))
{
    // This follows the same procedure as in c_example.c

    double dx = (this->limits[1] - this->limits[0]) / this->steps;
    double x[this->steps + 1];
    double f[this->steps + 1];

    for (long i = 0; i <= this->steps; i++) {
        x[i] = this->limits[0] + (i * dx);
        f[i] = (*func)(x[i]);
    }
    
    double result = 0.0;
    for (long i = 0; i <= this->steps; i++) {
        if (i == 0 || i == this->steps) {
            result += f[i];
        } else if (i % 2 != 0) {
            result += 4 * f[i];
        } else {
            result += 2 * f[i];
        }
    }
    result = result * (dx / 3);
    return result;
}


// Create test functions

double func1(double x)
{
    return pow(x, 2);
}

double func2(double x)
{
    return cosh(x);
}


// Like in C, the function main() is what is actually run.

int main()
{
    Integration integrate1(0, 5, 100);
    Integration integrate2(0, 10, 500);

    double result1 = integrate1.integrate(&func1);
    double result2 = integrate2.integrate(&func2);

    // To print, the preferred way is using stream operators. These come from
    // the standard library, but we had to include <iostream> to use them.
    // Again, we use the scope-definition operator because "cout" and "endl" are
    // defined in the scope of "std".
    // 'cout' means "C out," as in output.

    std::cout << "Integral 1: x**2 from 0 to 5: " << result1 << std::endl;
    std::cout << "Integral 2: cosh(x) from 0 to 10: " << result2 << std::endl;

    return 0;
}


/*
 * Compile this with g++:
 *      
 *      g++ -o program_name cpp_example.cpp -lm
 * 
 * Then run with "./program_name"
 */