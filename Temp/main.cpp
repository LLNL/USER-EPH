
#include <iostream>

#include "../eph_spline.h"

int main() {
  std::vector<double> y = {10, 9, 8, 7, 6};
  double dx = 1; 
  
  EPH_Spline<> spl(dx, y);
  //EPH_Spline<> spl;
  
  //std::cout << spl(0) << " " << spl(1) << " " << spl(2) << " " << spl(3) 
  //  << " " << spl(4) << " " << spl(1.5) << '\n';
  
  //std::cout << spl(4.1) << '\n';
  
  for(size_t i = 0; i < 100; ++i) {
    double x = 0.05 * i;
    std::cout << x << " " << spl(x) << '\n';
  }
  return 0;
}

