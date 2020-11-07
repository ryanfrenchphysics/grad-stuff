// This is a single-line comment

/*
 * This is a multiline comment. The asterisk between the ends are not necessary, however,
 * this is the preferred style in all C which isn't ANSI (the C used for Linux development).
 * The way this code is written is the preferred standard.
 */

// Libraries are used with the #include directive. Anything with a pound sign is called a
// preprocessor directive. All of these things are run _before_ compilation.
#include<stdio.h>
#include<math.h>


// Create a struct with all of our information. C is a statically-typed language,
// which means we must supply all type information when defining variables and functions
struct IntegralInfo {
    double low, high;       // Double is an 8-byte float
    long steps;             // Long is an 8-byte integer
};


/*
 * Unlike Julia, we can't directly use a function inside a struct. So we will develop a function
 * now that will check if steps is even. C functions are written as:
 *      return-type function-name(type1 var1, ...)
 *      {
 *          do stuff here;
 *      }
 */


// Void return type doesn't return a type. Here, our parameter will be pass-by-pointer,
// which will allow us to modify the value inside the original struct.
// A pass-by-value here would only modify the copied struct, which only exists inside the function.

void check_if_even(struct IntegralInfo *info)
{
    /*
     * Pointers point to a memory location. So if we were to do something like:
     *      pointer_var = 5;
     * then we would be changing the value of the memory location. If we instead want to change
     * the _value_ of the variable pointed to by the pointer, we need to dereference it. We do so
     * by using the * operator:
     *      *pointer_var = 5;
     * Now the variable at memory location pointer_var contains the value 5.
     */

    if ((*info).steps % 2 != 0) {
        (*info).steps += 1;
    }
    return;         // This isn't necessary for void functions, but is suggested by most standards.
}


/*
 * Let's create our Simpson's Rule Integration function.
 * Since we are passing a function to a function, we need to pass a pointer to said
 * function. This is done by
 *      function-return-type (*function-name)(type1 arg1, ...)
 */

// Our simpsons_integrate function MUST be passed a function that has one argument
// of type double and returns a double. We no longer require a pointer to the struct,
// since we only want the values from info and don't wish to modify it.

double simpsons_integrate(double (*func)(double x), struct IntegralInfo info)
{
    double dx = (info.high - info.low) / info.steps;

    // Create empty arrays of type double. We must specify how many elements will be in each array
    long array_len = info.steps;
    double x[array_len];
    double f[array_len];

    // Fill arrays with x and f(x) vals
    for (long i = 0; i <= array_len; i++) {
        x[i] = info.low + (i * dx);
        f[i] = (*func)(x[i]);
    }

    double result = 0.0;
    for (long i = 0; i <= array_len; i++) {
        if (i == 0 || i == array_len) {
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


// Create functions

double func1(double x)
{
    return pow(x, 2);
}

double func2(double x)
{
    return cosh(x);
}


// The thing that is actually run is the function main(). We must ALWAYS have a main() function,
// or else we will get an error. In C, this must return an integer (this goes
// back to the days where return codes were heavily used).

int main()
{
    // Create our info structs
    struct IntegralInfo info1;
    struct IntegralInfo info2;

    // Fill in information
    info1.low = 0;
    info1.high = 5;
    info1.steps = 100;

    info2.low = 0;
    info2.high = 10;
    info2.steps = 500;

    // Since check_if_even requires a pointer, we must pass it the location of our
    // structs. To do this, we use the 'Address of' operator, &:
    check_if_even(&info1);
    check_if_even(&info2);

    // Like before, simpsons_integrate requires the location of the functions,
    // but only the values of the info structs
    double result1 = simpsons_integrate(&func1, info1);
    double result2 = simpsons_integrate(&func2, info2);

    // This is how C's formatted print works. The % operator states what datatype will be
    // passed to it. In our case, we are passing a double, so we use %f for float.
    printf("Integral 1: x**2 from 0 to 5: %f\n", result1);
    printf("Integral 2: cosh(x) from 0 to 10: %f\n", result2);

    // Returning 0 is a formality, it means that main reached the end with no errors.
    // In reality, you can return whatever you want here.
    return 0;
}


/*
 * To run this program, you must use a C compiler. I recommend gcc.
 * To compile, cd to the location of the code and do:
 * 
 *      gcc -o program-name c_example.c -lm
 * 
 * The -o program_name will create your program with the given name. Otherwise,
 * a generic program caled 'a.out' is created.
 * The -lm flag ensures that the math.h header is used (this is necessary for
 * pow() and cosh() functions).
 * 
 * To run the compiled program, do:
 * 
 *      ./program-name
 */