
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <chrono>

#include "eph_spline.h"

const unsigned int n = 1001;
const double x0 = 0.0;
const double dx = 0.01;
double x[n];
double y[n];

int main(int args, char **argv) {
  // populate function
  // sin
  for(int i = 0; i < n; ++i) {
    x[i] = x0 + i*dx;
    y[i] = sin(x[i]);
  }
  
  EPH_Spline spline(x0, dx);
  std::cout << "Testing call to GetValue on an empty spline" << std::endl;
  try {
    spline.GetValue(x0);
    std::cout << "FAIL" << std::endl;
  }
  catch (...) {
    std::cout << "PASS" << std::endl;
  }
  
  for(int i = 0; i < n; ++i) {
    spline << y[i];
  }
  
  std::cout << "Testing call to GetValue before initing" << std::endl;
  try {
    spline.GetValue(x0);
    std::cout << "FAIL" << std::endl;
  }
  catch (...) {
    std::cout << "PASS" << std::endl;
  }
  
  spline << true;
  
  // test borders
  std::cout << "Testing borders" << std::endl;
  std::cout << "Left : " << x0 << " " 
    << spline.GetValue(x0) << std::endl;
  std::cout << "Right: " << x0 + (n-1)*dx << " " 
    << spline.GetValue(x0 + (n-1)*dx) << std::endl;
    
  // test outside the range
  std::cout << "Testing outside range" << std::endl;
  try {
    spline.GetValue(x0-0.1);
    std::cout << "FAIL" << std::endl;
  }
  catch (...) {
    std::cout << "PASS" << std::endl;
  }
  try {
    spline.GetValue(x0 + (n-1)*dx + 0.1);
    std::cout << "FAIL" << std::endl;
  }
  catch (...) {
    std::cout << "PASS" << std::endl;
  }
  
  /**
   *
   * Do a billion or so calls to GetValue to test timing 
   * 
   **/
  // test performance
  std::cout << "Testing performance" << std::endl;
  double accumulator = 0.0;
  double daccumulator = 0.0;
  double ddaccumulator = 0.0;
  
  // random numbers
  std::default_random_engine gen(111);
  std::uniform_real_distribution<double> distr(x0, x0 + dx*(n-1));
  std::vector<double> numbers;
  for(unsigned long int i = 0; i < 1000000000; ++i)
    numbers.push_back(distr(gen));
  
  
  /** timing **/
  auto t1 = std::chrono::system_clock::now();
  for(auto i: numbers) {
        accumulator += spline.GetValue(i);
        daccumulator += spline.GetDValue(i);
        ddaccumulator += spline.GetDDValue(i);
  }
  auto t2 = std::chrono::system_clock::now();
  std::cout << "Elapsed time: " 
    << ((std::chrono::duration_cast<std::chrono::milliseconds> (t2-t1)).count())/1000.0 
    << std::endl; 
  
  std::cout << "Total sum: " << accumulator << " should be: " << 1.83892e+08 << std::endl;
  std::cout << "Total dsum: " << daccumulator << " should be: " << -5.43821e+07 << std::endl;
  std::cout << "Total ddsum: " << ddaccumulator << " should be: " << -1.84457e+08 << std::endl;
  
  std::cout << "Testing interpolation values" << std::endl;
  std::ofstream fn("out.data");
  if(fn.is_open()) {
    double xx = x0;
    while(xx < x0 + (n-1)*dx) {
      fn << xx << " " << sin(xx) << " " << spline.GetValue(xx) << " " 
        << cos(xx) << " " << spline.GetDValue(xx) << " "
        << -sin(xx) << " " << spline.GetDDValue(xx) <<'\n';
      xx += 0.001;
    }
    
    fn.close();
  }
  
  return 0;
}

