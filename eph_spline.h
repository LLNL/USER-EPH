/*
 * Authors of the extension Artur Tamm, Alfredo Correa
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_SPLINE
#define EPH_SPLINE

#include <vector>
#include <cmath>
#include <cassert>
#include <cstddef>

/* 
 * Stripped down version of the EPH_Spline class
 * 
 * 
 */

template<typename Float = double, template<typename> class Allocator = std::allocator, template <typename _F = Float, typename _A = Allocator<Float>> class Container = std::vector>
class EPH_Spline {
  public:
    EPH_Spline() {}
    EPH_Spline(Float dx, Container<> &y) :
      inv_dx {1./dx},
      c(y.size())
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
        c[i].b = (y[i+1]-y[i]) / dx;
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
      c[0].c = fabs(c[1].b - c[0].b);
      c[0].d = fabs(z1-z0);
      
      c[1].c = fabs(c[2].b - c[1].b);
      c[1].d = fabs(c[0].b-z1);
      
      c[points-2].c = fabs(z2 - c[points-2].b);
      c[points-2].d = fabs(c[points-3].b - c[points-4].b);
      
      c[points-1].c = fabs(z3 - z2);
      c[points-1].d = fabs(c[points-2].b - c[points-3].b);
      
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
        if(d_2 == d_1 && d0 != d1)
          c[i].a = d_1;
        else if(d0 == d1 && d_2 == d_1)
          c[i].a = d0;
        else if(d_1 == d0)
          c[i].a = d0;
        else if(d_2 == d_1 && d0 == d1 && d0 != d_1)
          c[i].a = 0.5 * (d_1 + d0);
        else
          c[i].a = (d_1*w1 + d0*w0) / (w1+w0);
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
    
    Float operator() (Float x) const {
      assert(x >= 0);
      
      size_t index = x * inv_dx; // use inv_dx = 1./dx instead
      assert(index < c.size());
      
      return c[index].a + x * (c[index].b + x * (c[index].c + x * c[index].d));
    }
    
  protected:
    constexpr static size_t min_size {3};
    
    struct Coefficients {
      Float a, b, c, d;
    };
    
    Float inv_dx;
    Container<Coefficients, Allocator<Coefficients>> c;
};

#endif
