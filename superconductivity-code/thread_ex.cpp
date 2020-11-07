#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <./eigen/Eigen/Dense>
#include <./eigen/Eigen/Eigenvalues>
#include <./eigen/Eigen/Core>
#include <complex>
#include <vector>
#include <iomanip>
#include <thread>
#include <future>
#include <ctime>
#include <chrono> // Look into this instead of clock from ctime


// Geometric sum of r, N terms
// Sum of r^k from k=0 to k=N
auto g_series(float r, int N)
{
  auto sum = 0.0;
  for(int k = 0; k < N+1; k++){
    sum += pow(r, k);
  }
  return sum;
}


auto partial_g_series(float r, int start, int end)
{
  auto sum = 0.0;
  for(int k = start; k < end; k++){
    sum += pow(r, k);
  }
  return sum;
}


std::vector<int> vec_indices(int vec_len, int chunks)
{
  std::vector<int> ret_vec(chunks);
  int chunk_size = floor(vec_len / chunks);
  int remainder = vec_len % chunks;

  for(int i=0; i<chunks; i++) {
    if (i < remainder) {
      ret_vec[i] = chunk_size + 1;
      continue;
    }
    ret_vec[i] = chunk_size;
  }
  return ret_vec;
}


auto g_series_threaded(float r, int N, int threads)
{
  auto sum = 0.0;
  std::vector<std::future<double>> threadvec;

  auto indices_vec = vec_indices(N, threads);

  int start = 0;
  int end = indices_vec[0] + 1;
  for (int i=0; i<threads; i++) {
    threadvec.emplace_back(std::async(partial_g_series, r, start, end));
    start = end;
    if (i != threads - 1) {
      end += indices_vec[i + 1] + 1;
    }
  }

  for (int i=0; i<threads; i++) {
    sum += threadvec[i].get();
  }

  return sum;
}



/*
 * See board problem
 *
 */



int main(int argc, char *argv[])
{
  std::clock_t start;
  auto r = std::strtod(argv[1], nullptr);
  auto N = std::strtol(argv[2], nullptr, 0);
  auto num_threads = std::strtol(argv[3], nullptr, 0);

  start = std::clock();
  auto s1 = g_series(r, N);
  double duration1 = (std::clock() - start) / (double) CLOCKS_PER_SEC;

  start = std::clock();
  auto s1_threaded = g_series_threaded(r, N, num_threads);
  double duration2 = (std::clock() - start) / (double) CLOCKS_PER_SEC;

  std::cout << "Geometric Series: " << s1 << "\n";
  std::cout << "Time: " << duration1 << "s \n\n";

  std::cout << "Geometric Series, threads-" << num_threads << ": " << s1_threaded << "\n";
  std::cout << "Time: " << duration2 << "s \n\n";

  return 0;
}

/*
 *  Compile with:
 *
 *  gcc -lstdc++ -lm -std=c++14 -Wall -fconcepts -fext-numeric-literals -pthread -g thread_ex.cpp -o thread_test
 */
