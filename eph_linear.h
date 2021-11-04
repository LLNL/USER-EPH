/*
 * Authors of the extension Artur Tamm, Alfredo Correa
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_LINEAR
#define EPH_LINEAR

#include <vector>

/*
 * This is a linear interpolation that supports reverse lookup.
 * Therefore the Ce(T) used has to be monotonic (growing or constant).
 */

struct EPH_Linear {
  std::vector<double> x;
  std::vector<double> y;
  
  
  EPH_Linear() = default;
  
  template<typename x_it, typename y_it>
  EPH_Linear(x_it x_begin, x_it x_end, y_it y_begin) {
    
  }
  
};

