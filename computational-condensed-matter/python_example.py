# Import a module with a pseudonym
import numpy as np
# Import specific submodules
from math import pi
# Import all submodules from a module. This is not recommended, as you can have namespace collisions.
from pprint import *


# Class construction in Python:
class Integration:
    def __init__(self, type='simpsons'):
        # If you'd like to add other methods of integration in the future, one could include other types
        # Use triple quotes (single or double) for multiline comments
        '''
            __init__ is the constructor for the Integration class
            self is a keyword that refers to the object itself. Thus we use it to only use values/methods within
            the object we are concerned with.
        '''

        # Double underscores before a name are the most typical way to set "private" members and methods.
        #   Note: Python doesn't actually have private methods and members like C++. 

        self.__limits = [0, 0]       # Set initial integration limits in a list: [low, high]

    def setlimits(self, low, high):
        # Note: All methods in a class must contain 'self' as the first parameter.
        '''
            Set low and high limits for the integral.
        '''

        # We will check whether or not low and high are numbers:
        try:
            # We will see if both low and high are instances of a float or integer
            if (not(isinstance(low, float) or isinstance(low, int)) or 
                    not(isinstance(high, float) or isinstance(high, int))):
                # If not, raise a value error
                raise ValueError

        except ValueError:
            # If we have a value error, print an error message and return out of the function
            # We will use f-strings, which is the preferred (and quickest) method of formatted printing
            # Notice the 'f' before the opening quotes
            print(f"Error: '{low}' and/or '{high}' is not a float.")
            return

        self.__limits = [low, high]      # If no error, set the limits


    def __simpsons(self, func, low, high, steps=100):
        '''
            Define Simpson's Rule Integration of a function from low to high.
            Default number of steps is 100.
        '''
        n = steps
        if n % 2 != 0:
            # If steps is not even, add one to steps to make it even
            n += 1
        
        deltax = (high - low) / n

        # Create empty list to store values of x & f(x) at each point
        x = []
        f = []

        for i in range(n + 1):
            # Append values of x and f(x) at each point
            x.append(low + (i * deltax))
            f.append(func(x[i]))

        result = 0
        for i in range(n + 1):
            if i == 0 or i == n:
                result += f[i]
            elif i % 2 != 0:
                result += 4 * f[i]
            else:
                result += 2 * f[i]

        result = result * (deltax / 3)
        return result

    def integrate(self, func, steps=100):
        return self.__simpsons(func, self.__limits[0], self.__limits[1], steps)



# To use the Integration class, we need a function to integrate:
def func1(x):
    return x**2

def func2(x):
    return np.cosh(x)


# Create a main function
def main():
    # Create two separate objects

    integral1 = Integration()
    integral2 = Integration(type='simpsons')    # Both of these objects use type of 'simpsons'

    # Set limits of the two integrals
    integral1.setlimits(0, 5)
    integral2.setlimits(0, 10)

    # Run the integrals and store the values.
    result1 = integral1.integrate(func1)
    result2 = integral2.integrate(func2, steps=500)

    # Print out the results of the integrals:
    print(f"Integral 1: x**2 from 0 to 5: {result1}\n")     # \n means 'newline'
    print(f"Integral 2: cosh(x) from 0 to 10: {result2}\n")

    # Now, let's demonstrate what happens if we try to set nonsensical limits on an integral. This will print the error.
    integral1.setlimits('a', "This isn't a float")


# Now, we only want to run the 'main' function if called as a stand-alone script:
if __name__ == "__main__":
    main()

'''
    To run this script, call 'python3 python_example.py' from command line.
    Note: Some operating systems may call Python 3 as just 'python' instead of 'python3'
'''

