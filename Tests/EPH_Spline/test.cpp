
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>
#include <chrono>

#include "eph_spline.h"

using namespace std;

constexpr size_t n {1001};
constexpr double x0 {0.0};
constexpr double dx {0.01};
vector<double> x(n);
vector<double> y(n);

int main(int args, char **argv) {
  // populate function
  // sin
  for(int i = 0; i < n; ++i) {
    x[i] = x0 + i*dx;
    y[i] = sin(x[i]);
  }

  Spline spline(dx, y);

  { // test 0 and end
    spline(x[0]);
    spline(x[x.size() - 1]);
  }

  //~ spline.reverse(0.5);
  { // test reverse lookup
    double l_x = spline.reverse(0.5);
    cout << "sin(" << l_x << ") = " << sin(l_x) << " (0.5)\n";
  }

  std::cout << "Testing interpolation values" << std::endl;
  std::ofstream fn("out.data");
  if(fn.is_open()) {
    double xx = x0;
    while(xx < x0 + (n-1)*dx) {
      fn << xx << " " << sin(xx) << " " << spline(xx) << '\n';
      xx += 0.001;
    }

    fn.close();
  }
  return 0;
}


  /**
   *
   * Do a billion or so calls to GetValue to test timing
   *
   **/
  //~ // test performance
  //~ std::cout << "Testing performance" << std::endl;
  //~ double accumulator = 0.0;
  //~ // random numbers
  //~ std::default_random_engine gen(111);
  //~ std::uniform_real_distribution<double> distr(x0, x0 + dx*(n-1));
  //~ std::vector<double> numbers;
  //~ for(unsigned long int i = 0; i < 1000000000; ++i) {
    //~ numbers.push_back(distr(gen));
  //~ }

  //~ /** timing **/
  //~ auto t1 = std::chrono::system_clock::now();
  //~ for(auto i: numbers) {
    //~ accumulator += spline(i);
  //~ }
  //~ auto t2 = std::chrono::system_clock::now();
  //~ std::cout << "Elapsed time: "
    //~ << ((std::chrono::duration_cast<std::chrono::milliseconds> (t2-t1)).count())/1000.0
    //~ << std::endl;

  //~ std::cout << "Total sum: " << accumulator << " should be: " << 1.83892e+08 << std::endl;

  //~ std::cout << "Testing interpolation values" << std::endl;
  //~ std::ofstream fn("out.data");
  //~ if(fn.is_open()) {
    //~ double xx = x0;
    //~ while(xx < x0 + (n-1)*dx) {
      //~ fn << xx << " " << sin(xx) << " " << spline(xx) << '\n';
      //~ xx += 0.001;
    //~ }

    //~ fn.close();
  //~ }

  //~ return 0;
//~ }

