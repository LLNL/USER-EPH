/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

//#define DEBUG_EPH

// external headers
#ifdef DEBUG_EPH
#include <iostream>
#endif

#include <fstream>
#include <string>
#include <sstream>
#include <cstring>

// internal headers
#include "eph_beta.h"
#include "eph_spline.h"

EPH_Beta::EPH_Beta(const char* file) {
  std::ifstream fd(file);
  
  if(!fd.is_open()) {
    return;
  }
  
  char line[lineLength];
  
  // read first three lines
  fd.getline(line, lineLength);
  fd.getline(line, lineLength);
  fd.getline(line, lineLength);
  
  unsigned int nElements;
  unsigned int nPointsRho;
  unsigned int nPointsBeta;
  
  double dr;
  double drho;
  
  fd >> nElements;
  
  if(nElements < 1) {
    throw;
  }
  
  elementName.resize(nElements);
  elementNumber.resize(nElements);
  rho.resize(nElements);
  beta.resize(nElements);
  
  fd.getline(line, lineLength);
  std::string str(line);
  std::istringstream strstream(str);
  
  for(unsigned int i = 0; i < nElements; ++i) {
    std::string elem;
    strstream >> elem;
    elementName[i] = elem;
  }
  
  fd >> nPointsRho;
  fd >> dr;
  fd >> nPointsBeta;
  fd >> drho;
  fd >> rcut;
  rhocut = drho * (nPointsBeta-1);
  
  for(unsigned int i = 0; i < nElements; ++i) {
    fd >> elementNumber[i];
    
    rho[i].SetDiscretisation(0.0, dr);
    for(unsigned int j = 0; j < nPointsRho; ++j) {
      double rho_r;
      fd >> rho_r;
      rho[i] << rho_r;
    }
    rho[i] << true;
    
    beta[i].SetDiscretisation(0.0, drho);
    for(unsigned int j = 0; j < nPointsBeta; ++j) {
      double beta_rho;
      fd >> beta_rho;
      beta[i] << beta_rho;
    }
    beta[i] << true;
  }

  fd.close();
}
