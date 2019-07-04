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
    
    EPH_Spline(double x0, double dx) : // constructor for use with << operator
      x_First {x0}, 
      dx {dx}, 
      x_Last {x0 - 0.1} 
      { }
      
    EPH_Spline(double x0, double dx, const double *y, size_t points); // constructor to initialise the spline fully
    
    EPH_Spline(double x0, double dx, const std::vector<double>& a_y) : // constructor to initialise the spline fully
      x_First {x0}, 
      dx {dx}, 
      x_Last {x0 + a_y.size()*dx},
      y {a_y}
    {
      a.resize(y.size());
      b.resize(y.size());
      c.resize(y.size());
      d.resize(y.size());

      da.resize(y.size());
      db.resize(y.size());
      dc.resize(y.size());

      dda.resize(y.size());
      ddb.resize(y.size());

      FindCoefficients();
    }
    
    EPH_Spline(double x0, double dx, std::vector<double>&& y) : // constructor to initialise the spline fully
      x_First {x0}, 
      dx {dx}, 
      x_Last {x0 + y.size()*dx},
      y {std::move(y)}
    {
      a.resize(y.size());
      b.resize(y.size());
      c.resize(y.size());
      d.resize(y.size());

      da.resize(y.size());
      db.resize(y.size());
      dc.resize(y.size());

      dda.resize(y.size());
      ddb.resize(y.size());

      FindCoefficients();
    }
    
    // special type of initialisation; useful when initialising from file
    void SetDiscretisation(double x0, double dx) {
      *this << false;
      this->x_First = this->x_Last = x0;
      this->dx = dx;
    }
    
    // add a point to spline
    EPH_Spline& operator<< (const double y);
    // initialise the spline based or points or reset
    EPH_Spline& operator<< (const bool init);
    
    // get the value of the function at x
    double GetValue(double x) const {
      // maybe this should be a debug option and asserted instead?
      #ifndef DNDEBUG
      if(x < this->x_First)
        throw std::runtime_error("eph_spline: argument smaller than the lower bound");
      else if(x > this->x_Last)
        throw std::runtime_error("eph_spline: argument larger than the upper bound");
      #endif
      
      auto const index = FindIndex(x);
      
      return a[index] + x * (b[index] + x * (c[index] + x * d[index]));
    }
    
    // get a derivative of the function at x
    inline double GetDValue(double x) const { 
      // maybe this should be a debug option and asserted instead?
      #ifndef DNDEBUG
      if(x < this->x_First)
        throw std::runtime_error("eph_spline: argument smaller than the lower bound");
      else if(x > this->x_Last)
        throw std::runtime_error("eph_spline: argument larger than the upper bound");
      #endif
      
      auto index = FindIndex(x);
      
      return  da[index] + x * (db[index] + x * dc[index]);
    }
    
    // get the second derivative of the function at x (discontinuous!)
    double GetDDValue(double x) const { 
      // maybe this should be a debug option and asserted instead?
      #ifndef DNDEBUG
      if(x < this->x_First)
        throw std::runtime_error("eph_spline: argument smaller than the lower bound");
      else if(x > this->x_Last)
        throw std::runtime_error("eph_spline: argument larger than the upper bound");
      #endif
      
      auto index = FindIndex(x);
      
      return dda[index] + ddb[index] * x;
    }
  
  private: // some private functions for spline initialisation and value calculation
    size_t FindIndex(double x) const {
      return (x-x_First)/dx;
    }
    
    void FindCoefficients();
};

#endif
