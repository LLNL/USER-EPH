/*
 * Authors of the extension Artur Tamm, Alfredo Correa
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_LINEAR
#define EPH_LINEAR

#include <stddef.h> // bring size_t into scope

#include <vector>
#include <iterator> // std::distance
#include <algorithm> // std::lower_bound

#include <iostream>

/*
 * This is a linear interpolation that supports reverse lookup.
 * Therefore the Ce(T) used has to be monotonic (growing or constant).
 */

struct EPH_Linear {
  double dx;
  std::vector<double> y;
  std::vector<double> dy;
  
  EPH_Linear() = default;
  
  template<typename y_it>
  EPH_Linear(double _dx, y_it y_begin, y_it y_end)
      : dx {_dx} {
    while(y_begin != y_end) {
      y.push_back(*(y_begin++));
    }
    
    for(size_t i = 0; i < y.size() - 1; ++i) {
      dy.push_back((y[i+1] - y[i]) / dx);
    }
  }
  
  // given an x find the y value
  double operator() (double _x) { 
    size_t idx = static_cast<size_t>(_x / dx); // we ignore the end point
    if(idx < y.size()) {
      double delta = _x - idx * dx;
      return y[idx] + dy[idx] * delta;
    }
    return 0.;
  }
  
  // given a y find the appropriate x
  double reverse_lookup(double _y) {
    auto it {std::upper_bound(y.begin(), y.end(), _y)};
    if(it != y.end()) {
      size_t idx {static_cast<size_t>(std::distance(y.begin(), it))};
      idx--;
      return idx * dx + 1. / dy[idx] * (_y - y[idx]);
    }
    
    return 0.;
  }
  
  double derivative(double _x) {
    size_t idx = static_cast<size_t>(_x / dx); // we ignore the end point
    if(idx < y.size()) { return dy[idx]; }
    
    return dy[dy.size() - 1];
  }
  
};

#endif
