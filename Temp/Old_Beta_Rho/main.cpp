
#include <iostream>

#include "eph_spline.h"
#include "eph_beta.h"

#if 1
int main() {
  EPH_Beta beta("Beta_Rho.beta");
  
  double r_cutoff = beta.getCutoff();
  double dr = 0.001;
  size_t N = r_cutoff / dr;
  
  std::cout << "# " << beta.getElementsNumber() << 
    ' ' << beta.getCutoff() << ' ' << beta.getRhoCutoff() << '\n';
  
  for(size_t i = 0; i < beta.getElementsNumber(); ++i) {
    std::cout << "# " << beta.getName(i) << ' ' << beta.getNumber(i) << '\n';
  }
  
  for(size_t i = 0; i < N; ++i) {
    double r = i * dr;
    
    std::cout << r << " " << beta.getRho(0, r) << '\n';
  }
  
  return 0;
}

#endif
