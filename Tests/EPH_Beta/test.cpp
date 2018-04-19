
#include <iostream>

#include "eph_beta.h"

int main(int args, char **argv) {
  std::cout << "EPH_Beta class tests" << std::endl;
  
  std::cout << "Load file with beta(rho)" << std::endl;
  EPH_Beta beta("NiFe.beta");
  std::cout << "done" << std::endl;
  
  std::cout << "Number of elements " << beta.getElementsNumber() << 
    " should be 2" << std::endl;
  
  std::cout << "Element names: " << std::endl;
  for(int i = 0; i < beta.getElementsNumber(); ++i) {
    std::cout << "  " << beta.getName(i) << std::endl;
  }
  std::cout << "Cutoff: " << beta.getCutoff() << " should be 5.0" << std::endl;
  std::cout << "Rho Cutoff: " << beta.getRhoCutoff() << " should be 10.0" << std::endl;
  
  std::cout << "Element numbers: " << std::endl;
  for(int i = 0; i < beta.getElementsNumber(); ++i) {
    std::cout << "  " << beta.getNumber(i) << std::endl;
  }
  
  // test some density values
  std::cout << "Some rhos" << std::endl;
  for(int i = 0; i < beta.getElementsNumber(); ++i) {
    std::cout << "  " << beta.getRho(i, 0.0) << std::endl;
  }
  
  std::cout << "Some betas" << std::endl;
  for(int i = 0; i < beta.getElementsNumber(); ++i) {
    std::cout << "  " << beta.getBeta(i, 0.0) << std::endl;
  }
  
  // test values outside the range
  try {
    beta.getName(beta.getElementsNumber());
  } catch(...) {
    std::cout << "Caught out of range element" << std::endl;
  }
  
  return 0;
}

