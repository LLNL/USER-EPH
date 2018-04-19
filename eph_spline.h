/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_SPLINE
#define EPH_SPLINE

#include <vector>

class EPH_Spline {
  public:
  private:
    double x_First; // coordinate of the first element
    double x_Last; // coordinate of the last element
    double dx; // discretisation step
    
    //unsigned int points; // number of points in the spline
    std::vector<double> y; // y values
    
    // coeffiecients for function values
    // splines are always in a + b*x + c*x**2 + d*x**3 form
    std::vector<double> a;
    std::vector<double> b;
    std::vector<double> c;
    std::vector<double> d;
    
    // coefficients for the derivative
    // da + db * x + dc * x**2
    std::vector<double> da;
    std::vector<double> db;
    std::vector<double> dc;
    
    // coefficients for the second derivative
    std::vector<double> dda;
    std::vector<double> ddb;
    
  public:
    EPH_Spline() : EPH_Spline(0.0, 0.0) {} // default constructor
    EPH_Spline(const double x0, const double dx) { 
        x_First = x0;
        this->dx = dx;
        x_Last = x0 - 0.1; // this ensures that GetValue() will return NaN until FindCoefficients has been called
    }
    EPH_Spline(const double x0, const double dx, const double *y, const unsigned int points); // initialise spline
    
    // special type of initialisation (show off)
    void SetDiscretisation(const double x0, const double dx) {
      *this << false;
      this->x_First = this->x_Last = x0;
      this->dx = dx;
    }
    
    EPH_Spline& operator<< (const double y);
    EPH_Spline& operator<< (const bool init);
    
    double GetValue(const double x) const { // get function value at x
      double result = 0.0;
      unsigned int index = 0;
      
      if(x < this->x_First)
        throw;
      else if(x > this->x_Last)
        throw;
      
      index = FindIndex(x);
      
      double x2 = x*x;
      double x3 = x*x2;
  
      result = a[index] + b[index] * x + c[index] * x2 + d[index] * x3;
  
      return result;
    }
    
    double GetDValue(const double x) const { // get function derivative at x
      double result = 0.0;
      unsigned int index = 0;
      
      if(x < this->x_First)
        throw;
      else if(x > this->x_Last)
        throw;
      
      index = FindIndex(x);
      
      double x2 = x*x;
      double x3 = x*x2;
      
      result = da[index] + db[index] * x + dc[index] * x2;
      
      return result;
    }
      
    double GetDDValue(const double x) const { // get second derivative at x (discontinuous!)
      double result = 0.0;
      unsigned int index = 0;
      
      if(x < this->x_First)
        throw;
      else if(x > this->x_Last)
        throw;
      
      index = FindIndex(x);
      
      double x2 = x*x;
      double x3 = x*x2;
      
      result = dda[index] + ddb[index] * x;
      
      return result;
    }
  
  private: // some private functions for spline initialisation and value calculation
    unsigned int FindIndex(const double x) const {
      return ((x-this->x_First)/dx);
    }
    
    void FindCoefficients();
};

#endif
