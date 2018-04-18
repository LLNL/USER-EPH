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
  elements = 0;
  
  std::ifstream fd(file);
  
  if(!fd.is_open()) {
    return;
  }
  
  unsigned int inputType; // 0 old; 1 new json (TODO)
  char symbol;
  char line[lineLength];
  
  // read first free lines
  inputType = 1;
  
  symbol = fd.peek();
  if(symbol == '#') {
    fd.getline(line, lineLength);
  
    symbol = fd.peek();
    if(symbol == '#') {
      fd.getline(line, lineLength);
      
      symbol = fd.peek();
      fd.getline(line, lineLength);
  
      if(symbol == '#') {
        inputType = 0;
      }
    }
  }
  
  if(inputType == 0) {
    unsigned int nElements;
    unsigned int nPointsRho;
    unsigned int nPointsBeta;
    
    double dr;
    double drho;
    
    fd >> nElements;
    
    if(nElements < 1) {
      return;
    }
    
    elements = nElements;
    elementName = new char*[nElements];
    elementNumber = new unsigned int[nElements];
    this->rho = new EPH_Spline*[nElements];
    this->beta = new EPH_Spline*[nElements];
    
    fd.getline(line, lineLength);
    std::string str(line);
    std::istringstream strstream(str);
    
    for(unsigned int i = 0; i < nElements; ++i) {
      elementName[i] = new char[nameLength];
      this->rho[i] = new EPH_Spline();
      this->beta[i] = new EPH_Spline();
      
      std::string elem;
      strstream >> elem;
      
      std::strcpy(elementName[i], elem.c_str());
    }
    
    fd >> nPointsRho;
    fd >> dr;
    fd >> nPointsBeta;
    fd >> drho;
    fd >> cutoff;
    
    for(unsigned int i = 0; i < elements; ++i) {
      fd >> elementNumber[i];
      
      (*(this->rho[i])).SetDiscretisation(0.0, dr);
      for(unsigned int j = 0; j < nPointsRho; ++j) {
        double rho_r;
        fd >> rho_r;
        *(this->rho[i]) << rho_r;
      }
      *(this->rho[i]) << true;
      
      (*(this->beta[i])).SetDiscretisation(0.0, drho);
      for(unsigned int j = 0; j < nPointsBeta; ++j) {
        double beta_rho;
        fd >> beta_rho;
        *(this->beta[i]) << beta_rho;
      }
      *(this->beta[i]) << true;
      
      
    }
  }
  else {
    // TODO: json type input
  }
  
  
  fd.close();
}

EPH_Beta::~EPH_Beta() {
  for(unsigned int i = 0; i < elements; ++i) {
    delete[] elementName[i];
    delete rho[i];
    delete beta[i];
  }
  
  delete[] elementName;
  delete[] elementNumber;
  delete[] rho;
  delete[] beta;
}

