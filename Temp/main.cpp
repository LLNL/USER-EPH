
#include <iostream>

#include "../eph_spline.h"

int main() {
  std::vector<double> x = { 1, 2, 3, 4, 5};
  std::vector<double> y = {10, 9, 8, 7, 6};
  auto dx = x[1] - x[0]; 
  
  EPH_Spline<> spl(dx, y);
  //EPH_Spline<> spl;
  
  std::cout << spl(1) << " " << spl(2) << " " << spl(1.5) << " " << spl(3) 
    << " " << spl(4) << " " << spl(4.999) << '\n';
  
  //std::cout << spl(5) << '\n';
  return 0;
}

