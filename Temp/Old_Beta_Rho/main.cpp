
#include <iostream>

#include "eph_spline.h"
#include "eph_beta.h"

#if 0
int main() {
  EPH_Beta beta("Beta_Rho.beta");
  
  double r_cutoff = beta.get_r_cutoff();
  double dr = 0.001;
  size_t N = r_cutoff / dr;
  
  std::cout << "# " << beta.get_n_elements() << 
    ' ' << beta.get_r_cutoff() << ' ' << beta.get_rho_cutoff() << '\n';
  
  for(size_t i = 0; i < beta.get_n_elements(); ++i) {
    std::cout << "# " << beta.get_element_name(i) << ' ' << beta.get_element_number(i) << '\n';
  }
  
  for(size_t i = 0; i < N; ++i) {
    double r = i * dr;
    
    std::cout << r << " " << beta.get_rho(0, r) << '\n';
  }
  
  return 0;
}

#endif
