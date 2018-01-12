/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_BETA
#define EPH_BETA

// external headers 
#include <limits>

// internal headers
#include "eph_spline.h"

class EPH_Beta {
  public:
    EPH_Beta(const char* file); // initialise EPH_Beta based on file
    ~EPH_Beta();
    
    unsigned int getElementsNumber() {
      return elements;
    }
    
    char* getName(const unsigned int element) {
      if(element > elements) return nullptr;
      
      return elementName[element];
    }
    
    double getCutoff() {
      return cutoff;
    }
    
    unsigned int getNumber(const unsigned int element) {
      if(element > elements) return std::numeric_limits<int>::max();
      
      return elementNumber[element];
    }
    
    double getRho(const unsigned int element, const double r) {
      if(element > elements) return std::numeric_limits<double>::quiet_NaN();
      
      if(r > cutoff) return 0.0;
      return (rho[element])->GetValue(r);
    }
    
    double getDRho(const unsigned int element, const double r) {
      if(element > elements) return std::numeric_limits<double>::quiet_NaN();
      
      if(r > cutoff) return 0.0;
      return (rho[element])->GetDValue(r);
    }
    
    double getBeta(const unsigned int element, const double rho) {
      if(element > elements) return std::numeric_limits<double>::quiet_NaN();
      
      return (beta[element])->GetValue(rho);
    }
    
    double getDBeta(const unsigned int element, const double rho) {
      if(element > elements) return std::numeric_limits<double>::quiet_NaN();
      
      return (beta[element])->GetDValue(rho);
    }
  
  private:
    static constexpr unsigned int lineLength = 1024;
    static constexpr unsigned int nameLength = 8;
    
    unsigned int elements;
    
    unsigned int* elementNumber;
    char** elementName; // list of element names
    double cutoff;
    EPH_Spline** rho; // rho(r)
    EPH_Spline** beta; // beta(rho)
};

#endif
