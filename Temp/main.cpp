
#include <iostream>

#include "../eph_spline.h"
#include "../eph_beta.h"

#if 0
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
#endif

#if 1
int main() {
  EPH_Beta<> beta;
  beta = EPH_Beta<>("Beta_Rho.beta");
  
  double r_cutoff = beta.get_r_cutoff();
  double dr = 0.001;
  size_t N = r_cutoff / dr;
  
  std::cout << "# " << beta.get_n_elements() << 
    ' ' << beta.get_r_cutoff() << ' ' << beta.get_rho_cutoff() << '\n';
  
  for(size_t i = 0; i < beta.get_n_elements(); ++i) {
    std::cout << "# " << beta.get_element_name(i) << ' ' << beta.get_element_number(i) << '\n';
  }
  
  /*
  for(size_t i = 0; i < N; ++i) {
    double r = i * dr;
    
    std::cout << r << " " << beta.get_rho(0, r) << ' ' << beta.get_rho_r_sq(0, r * r) << '\n';
  }
  */
  
  double rho_cutoff = beta.get_rho_cutoff();
  double drho = 0.01;
  size_t M = rho_cutoff / drho;
  
  for(size_t i = 0; i < M; ++i) {
    double rho = i * drho;
    
    std::cout << rho << ' ' << beta.get_beta(0, rho) << ' ' << beta.get_alpha(0, rho) << '\n';
  }
  
  return 0;
}


#endif

#if 0

using Float = double;

template <typename _F = Float>
using Allocator = std::allocator<_F>;

template<typename _F = Float, typename _A = Allocator<_F>>
using Container = std::vector<_F, _A>;

//using Beta = EPH_Beta<Container<Float, Allocator<Float>>>;
using Beta = EPH_Beta<Float, Allocator, Container>;


int main() {
  Container<uint8_t> a = {'1', '2', '3'};
  Container<> b = {1.0, 2.0, 3.0};
  
  std::cout << a[0] << ' ' << a[1] << ' ' << a[2] << '\n';
  std::cout << b[0] << ' ' << b[1] << ' ' << b[2] << '\n';
  
  Beta beta();
  
  return 0;
}

#endif
