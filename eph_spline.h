/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_SPLINE
#define EPH_SPLINE

#include <vector>
#include <stdexcept>

class EPH_Spline {
  private:
    double x_First; // coordinate of the first element
    double x_Last; // coordinate of the last element
    double dx; // discretisation step
    
    std::vector<double> y; // y values
    
    // coeffiecients for splines
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
    
    // current version will brake if less there are less than 4 points
    constexpr static int min_size = 4;
    
  public:
    EPH_Spline() : EPH_Spline(0.0, 0.0) {} // default constructor
    
    EPH_Spline(const double x0, const double dx) : // constructor for use with << operator
      x_First {x0}, 
      dx {dx}, 
      x_Last {x0 - 0.1} 
      { }
      
    EPH_Spline(const double x0, const double dx, const double *y, const unsigned int points); // constructor to initialise the spline fully
    
    EPH_Spline(const double x0, const double dx, std::vector<double>&& y) : // constructor to initialise the spline fully
      x_First {x0}, 
      dx {dx}, 
      x_Last {x0 + y.size()*dx},
      y {std::move(y)}
    {
      FindCoefficients();
    }
    
    // special type of initialisation; useful when initialising from file
    void SetDiscretisation(const double x0, const double dx) {
      *this << false;
      this->x_First = this->x_Last = x0;
      this->dx = dx;
    }
    
    // add a point to spline
    EPH_Spline& operator<< (const double y);
    // initialise the spline based or points or reset
    EPH_Spline& operator<< (const bool init);
    
    // get the value of the function at x
    inline double GetValue(const double x) const {
      double result = 0.0;
      unsigned int index = 0;
      
      // maybe this should be a debug option and asserted instead?
      #ifndef EPH_UNSAFE
      if(x < this->x_First)
        throw std::runtime_error("eph_spline: argument smaller than the lower bound");
      else if(x > this->x_Last)
        throw std::runtime_error("eph_spline: argument larger than the upper bound");
      #endif
      
      index = FindIndex(x);
      
      double x2 = x*x;
      double x3 = x*x2;
  
      result = a[index] + b[index] * x + c[index] * x2 + d[index] * x3;
  
      return result;
    }
    
    // get a derivative of the function at x
    inline double GetDValue(const double x) const { 
      double result = 0.0;
      unsigned int index = 0;
      
      // maybe this should be a debug option and asserted instead?
      #ifndef EPH_UNSAFE
      if(x < this->x_First)
        throw std::runtime_error("eph_spline: argument smaller than the lower bound");
      else if(x > this->x_Last)
        throw std::runtime_error("eph_spline: argument larger than the upper bound");
      #endif
      
      index = FindIndex(x);
      
      double x2 = x*x;
      
      result = da[index] + db[index] * x + dc[index] * x2;
      
      return result;
    }
    
    // get the second derivative of the function at x (discontinuous!)
    inline double GetDDValue(const double x) const { 
      double result = 0.0;
      unsigned int index = 0;
      
      // maybe this should be a debug option and asserted instead?
      #ifndef EPH_UNSAFE
      if(x < this->x_First)
        throw std::runtime_error("eph_spline: argument smaller than the lower bound");
      else if(x > this->x_Last)
        throw std::runtime_error("eph_spline: argument larger than the upper bound");
      #endif
      
      index = FindIndex(x);
      
      result = dda[index] + ddb[index] * x;
      
      return result;
    }
  
  private: // some private functions for spline initialisation and value calculation
    inline unsigned int FindIndex(const double x) const {
      return ((x-this->x_First)/dx);
    }
    
    void FindCoefficients();
};

#endif
