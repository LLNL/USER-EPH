/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_BETA
#define EPH_BETA

// external headers 
#include <vector>
#include <string>

// internal headers
#include "eph_spline.h"

class EPH_Beta {
  public:
    EPH_Beta(const char* file); // initialise EPH_Beta based on file
    EPH_Beta() = delete;
    
    int getElementsNumber() const {
      return elementName.size();
    }
    
    std::string getName(const unsigned int element) const {
      if(element >= elementName.size()) throw;
      
      return elementName[element];
    }
    
    double getCutoff() const {
      return rcut;
    }

    double getRhoCutoff() const {
      return rhocut;
    }
    
    unsigned int getNumber(const unsigned int element) const {
      if(element >= elementNumber.size()) throw;
      
      return elementNumber[element];
    }
    
    double getRho(const unsigned int element, const double r) const {
      if(element >= rho.size()) throw;
      
      //if(r > rcut) return 0.0; // TODO: this might be unnecessary
      return rho[element].GetValue(r);
    }
    
    double getDRho(const unsigned int element, const double r) const {
      if(element >= rho.size()) throw;
      
      //if(r > rcut) return 0.0;
      return rho[element].GetDValue(r);
    }
    
    double getBeta(const unsigned int element, const double rho) const {
      if(element > beta.size()) throw;
      return beta[element].GetValue(rho);
    }
    
    double getDBeta(const unsigned int element, const double rho) const {
      if(element > beta.size()) throw;
      return beta[element].GetDValue(rho);
    }
  
  private:
    static constexpr unsigned int lineLength = 1024; // this is for parsing
    
    double rcut;
    double rhocut;
    std::vector<int> elementNumber;
    std::vector<std::string> elementName;
    std::vector<EPH_Spline> rho;
    std::vector<EPH_Spline> beta;
};

#endif
