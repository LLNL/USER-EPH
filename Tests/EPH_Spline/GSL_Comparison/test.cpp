
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <chrono>

#include <gsl/gsl_interp.h>

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
  
  for(int i = 0; i < n; ++i) {
    spline << y[i];
  }
  
  spline << true;
  
  /**
   *
   * Do a billion or so calls to GetValue to test timing 
   * 
   **/
  // test performance
  std::cout << "Testing performance" << "\n";
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
    << "\n";
  
  std::cout << "Total sum: " << accumulator << " should be: " << 1.83892e+08 << "\n";
  std::cout << "Total dsum: " << daccumulator << " should be: " << -5.43821e+07 << std::endl;
  std::cout << "Total ddsum: " << ddaccumulator << " should be: " << -1.84457e+08 << std::endl;
  
  /* GSL interpolation */
  gsl_interp *interp;
  gsl_interp_accel *accel;
  
  /* possible options
   * gsl_interp_linear
   * gsl_interp_polynomial
   * gsl_interp_cspline
   * gsl_interp_akima
   * gsl_interp_steffen
   */
  
  interp = gsl_interp_alloc(gsl_interp_akima, n);
  accel = gsl_interp_accel_alloc();
  
  gsl_interp_init(interp, x, y, n);
  
  accumulator = 0;
  /** timing **/
  t1 = std::chrono::system_clock::now();
  for(auto i: numbers) {
        accumulator += gsl_interp_eval(interp, x, y, i, accel);
        daccumulator += gsl_interp_eval_deriv(interp, x, y, i, accel);
        ddaccumulator += gsl_interp_eval_deriv2(interp, x, y, i, accel);;
  }
  t2 = std::chrono::system_clock::now();
  std::cout << "Elapsed time (GSL): "
    << ((std::chrono::duration_cast<std::chrono::milliseconds> (t2-t1)).count())/1000.0 
    << "\n";
  
  std::cout << "Total sum: " << accumulator << " should be: " << 1.83892e+08 << "\n";
  std::cout << "Total dsum: " << daccumulator << " should be: " << -5.43821e+07 << std::endl;
  std::cout << "Total ddsum: " << ddaccumulator << " should be: " << -1.84457e+08 << std::endl;
  
  gsl_interp_free(interp);
  gsl_interp_accel_free(accel);
  
  return 0;
}

