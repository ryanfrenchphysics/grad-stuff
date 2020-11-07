#include <matplotlibcpp.h>
#include <cmath>
#include <algorithm>
#include <iostream>


// GLOBAL VARS
int DIMENSIONS = 1;
long NUMBER_K = 1000;
long NUMBER_E = 1000;
double ERROR = 0.01;
double HOP = 1.0;
double E_0 = 2 * DIMENSIONS;
double SPACING = 0.01;

typedef std::vector<std::vector<double>> matdbl;
typedef std::vector<double> vecdbl;


matdbl gen_permutations(matdbl &k_mat)
{
    long num_vals = DIMENSIONS * NUMBER_K);
    vecdbl k_vals_list(num_vals);
    matdbl perm_list;

    for (long i = 0; i < num_vals; ++i) {

    }



}

double dispersion()
{
    return 0.0;
}



vecdbl k_perms(const matdbl &k_mat, long dimensions, long num_k)
{
    vecdbl k_vals(NUMBER_K);
    for (long i = 0; i < NUMBER_K; i++) {
        k_vals[i] = k_mat[i][0];
    }

}
    
int main()
{
    return 0;
}