/*
 * Authors of the extension Artur Tamm, Alfredo Correa
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_SPLINE
#define EPH_SPLINE

#include <vector>
#include <cmath>
#include <cassert>

/* 
 * Stripped down version of the EPH_Spline class
 * 
 * 
 */
#if 0
template<typename Container = std::vector<double>>
class EPH_Spline {
  using Float = typename Container::value_type;
  
  public:
    EPH_Spline() {}
    EPH_Spline(Float i_dx, Container &i_y) :
      dx {i_dx},
      y {i_y},
      a(i_y.size()),
      b(i_y.size()),
      c(i_y.size()),
      d(i_y.size())
    {
      size_t points = y.size();
      
      assert(dx > 0); // dx has to be positive
      assert(points > min_size); // EPH_Spline needs at least 4 points
  
      // we use b, c, and d as temporary buffers
      Float z0; // z_-2
      Float z1; // z_-1
      
      Float z2; // z_k-1
      Float z3; // z_k
      
      for(size_t i = 0; i < points-1; ++i) {
        //b -> z
        b[i] = (y[i+1]-y[i])/dx;
      }
    
      z1 = 2.0*b[0] - b[1];
      z0 = 2.0*z1 - b[0];
      
      z2 = 2.0*b[points-2] - b[points-1];
      z3 = 2.0*z2 - b[points-1];
      
      b[points-1] = z2;
    
      for(size_t i = 2; i < points-2; ++i) {
        //c -> w_i-1 ; d -> w_i

        c[i] = fabs(b[i+1] - b[i]);
        d[i] = fabs(b[i-1] - b[i-2]);
      }
    
      // special cases
      c[0] = fabs(b[1]-b[0]);
      d[0] = fabs(z1-z0);
      
      c[1] = fabs(b[2]-b[1]);
      d[1] = fabs(b[0]-z1);
      
      c[points-2] = fabs(z2-b[points-2]);
      d[points-2] = fabs(b[points-3]-b[points-4]);
      
      c[points-1] = fabs(z3-z2);
      d[points-1] = fabs(b[points-2]-b[points-3]);
      
      //derivatives
      for(size_t i = 0; i < points; ++i) {
        Float w0, w1;
        Float d_2, d_1, d0, d1;
        
        if(i == 0) {
          d_2 = z0; d_1 = z1; d1 = b[i+1];
        }
        else if(i == 1) {
          d_2 = z1; d_1 = b[i-1]; d1 = b[i+1];
        }
        else {
          d_2 = b[i-2]; d_1 = b[i-1]; d1 = b[i+1];
        }
        
        d0 = b[i]; w1 = c[i]; w0 = d[i];
        
        // special cases
        if(d_2 == d_1 && d0 != d1) {
          a[i] = d_1;
        }
        else if(d0 == d1 && d_2 == d_1) {
          a[i] = d0;
        }
        else if(d_1 == d0) {
          a[i] = d0;
        }
        else if(d_2 == d_1 && d0 == d1 && d0 != d_1) {
          a[i] = 0.5 * (d_1 + d0);
        }
        else {
          a[i] = (d_1*w1 + d0*w0) / (w1+w0);
        }
      }
  
      // solve the equations
      for(size_t i = 0; i < points-1; ++i) {
        Float dx3 = dx*dx*dx;
        
        Float x0_1 = i * dx;
        Float x0_2 = i * dx * x0_1;
        Float x0_3 = i * dx * x0_2;
        
        Float x1_1 = (i+1) * dx;
        Float x1_2 = (i+1) * dx * x1_1;
        Float x1_3 = (i+1) * dx * x1_2;
        
        d[i] = (-a[i]*x0_1 - a[i+1]*x0_1 + a[i]*x1_1 + a[i+1]*x1_1 + 2.0*y[i] - 2.0*y[i+1]) / dx3;
        c[i] = (-a[i] + a[i+1] + 3.0*d[i]*x0_2 - 3.0*d[i]*x1_2) / 2.0 / dx;
        b[i] = (c[i]*x0_2 + d[i]*x0_3 - c[i]*x1_2 - d[i]*x1_3 - y[i] + y[i+1]) / dx;
        a[i] = y[i] - b[i]*x0_1 - c[i]*x0_2 - d[i]*x0_3;
      }
  
      a[points-1] = y[points-1];
      b[points-1] = 0.0;
      c[points-1] = 0.0;
      d[points-1] = 0.0;
    }
    
    Float operator() (Float i_x) const {
      assert(i_x >= 0);
      
      size_t index = i_x/dx;
      assert(index < y.size());
      
      return a[index] + i_x * (b[index] + i_x * (c[index] + i_x * d[index]));
    }
    
  private:
    constexpr static size_t min_size {3};
    
    Float dx;
    Container y;
    Container a, b, c, d;
};
#endif

//using EPH_Spline = EPH_Spline_Generic<double, std::allocator, std::vector>;
//using EPH_Spline = EPH_Spline_Generic<>;

template<typename Float = double, template<typename> class Allocator = std::allocator, template <typename _F = Float, typename _A = Allocator<Float>> class Container = std::vector>
class EPH_Spline {
  public:
    EPH_Spline() {}
    EPH_Spline(Float i_dx, Container<> &i_y) :
      dx {i_dx},
      y {i_y},
      c(i_y.size())
    {
      size_t points = y.size();
      
      assert(dx > 0); // dx has to be positive
      assert(points > min_size); // EPH_Spline needs at least 4 points
  
      // we use b, c, and d as temporary buffers
      Float z0; // z_-2
      Float z1; // z_-1
      
      Float z2; // z_k-1
      Float z3; // z_k
      
      for(size_t i = 0; i < points-1; ++i) {
        //b -> z
        c[i].b = (y[i+1]-y[i])/dx;
      }
    
      z1 = 2.0*c[0].b - c[1].b;
      z0 = 2.0*z1 - c[0].b;
      
      z2 = 2.0*c[points-2].b - c[points-1].b;
      z3 = 2.0*z2 - c[points-1].b;
      
      c[points-1].b = z2;
    
      for(size_t i = 2; i < points-2; ++i) {
        //c -> w_i-1 ; d -> w_i

        c[i].c = fabs(c[i+1].b - c[i].b);
        c[i].d = fabs(c[i-1].b - c[i-2].b);
      }
    
      // special cases
      c[0].c = fabs(c[1].b-c[0].b);
      c[0].d = fabs(z1-z0);
      
      c[1].c = fabs(c[2].b-c[1].b);
      c[1].d = fabs(c[0].b-z1);
      
      c[points-2].c = fabs(z2-c[points-2].b);
      c[points-2].d = fabs(c[points-3].b-c[points-4].b);
      
      c[points-1].c = fabs(z3-z2);
      c[points-1].d = fabs(c[points-2].b-c[points-3].b);
      
      //derivatives
      for(size_t i = 0; i < points; ++i) {
        Float w0, w1;
        Float d_2, d_1, d0, d1;
        
        if(i == 0) {
          d_2 = z0; d_1 = z1; d1 = c[i+1].b;
        }
        else if(i == 1) {
          d_2 = z1; d_1 = c[i-1].b; d1 = c[i+1].b;
        }
        else {
          d_2 = c[i-2].b; d_1 = c[i-1].b; d1 = c[i+1].b;
        }
        
        d0 = c[i].b; w1 = c[i].c; w0 = c[i].d;
        
        // special cases
        if(d_2 == d_1 && d0 != d1) {
          c[i].a = d_1;
        }
        else if(d0 == d1 && d_2 == d_1) {
          c[i].a = d0;
        }
        else if(d_1 == d0) {
          c[i].a = d0;
        }
        else if(d_2 == d_1 && d0 == d1 && d0 != d_1) {
          c[i].a = 0.5 * (d_1 + d0);
        }
        else {
          c[i].a = (d_1*w1 + d0*w0) / (w1+w0);
        }
      }
  
      // solve the equations
      for(size_t i = 0; i < points-1; ++i) {
        Float dx3 = dx*dx*dx;
        
        Float x0_1 = i * dx;
        Float x0_2 = i * dx * x0_1;
        Float x0_3 = i * dx * x0_2;
        
        Float x1_1 = (i+1) * dx;
        Float x1_2 = (i+1) * dx * x1_1;
        Float x1_3 = (i+1) * dx * x1_2;
        
        c[i].d = (-c[i].a*x0_1 - c[i+1].a*x0_1 + c[i].a*x1_1 + c[i+1].a*x1_1 + 2.0*y[i] - 2.0*y[i+1]) / dx3;
        c[i].c = (-c[i].a + c[i+1].a + 3.0*c[i].d*x0_2 - 3.0*c[i].d*x1_2) / 2.0 / dx;
        c[i].b = (c[i].c*x0_2 + c[i].d*x0_3 - c[i].c*x1_2 - c[i].d*x1_3 - y[i] + y[i+1]) / dx;
        c[i].a = y[i] - c[i].b*x0_1 - c[i].c*x0_2 - c[i].d*x0_3;
      }
  
      c[points-1].a = y[points-1];
      c[points-1].b = 0.0;
      c[points-1].c = 0.0;
      c[points-1].d = 0.0;
    }
    
    Float operator() (Float i_x) const {
      assert(i_x >= 0);
      
      size_t index = i_x/dx;
      assert(index < y.size());
      
      return c[index].a + i_x * (c[index].b + i_x * (c[index].c + i_x * c[index].d));
    }
    
  private:
    constexpr static size_t min_size {3};
    
    struct Coefficients {
      Float a, b, c, d;
    };
    
    Float dx;
    Container<> y;
    Container<Coefficients, Allocator<Coefficients>> c;
};

#if 0
class EPH_Spline_OLD {
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
    EPH_Spline& operator<< (double y);
    // initialise the spline based or points or reset
    EPH_Spline& operator<< (bool init);
    
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
    double GetDValue(double x) const { 
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
#endif
