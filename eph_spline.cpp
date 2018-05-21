/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

// external headers
#include <algorithm>
#include <cmath>
#include <stdexcept>

// internal headers
#include "eph_spline.h"

EPH_Spline::EPH_Spline(const double x0, const double dx, const double* y, const unsigned int points) : EPH_Spline(x0, dx) {
  if(!y) throw std::runtime_error("eph_spline: null pointer supplied");
  if(points <= min_size) throw std::runtime_error("eph_spline: not enough points supplied");
  
  // resize vectors so they would fit enough data
  this->y.resize(points);
  
  a.resize(points);
  b.resize(points);
  c.resize(points);
  d.resize(points);
  
  da.resize(points);
  db.resize(points);
  dc.resize(points);
  
  dda.resize(points);
  ddb.resize(points);
  
  for(int i = 0; i < points; ++i) {
    this->y[i] = y[i];
  }
  
  FindCoefficients();
}

EPH_Spline& EPH_Spline::operator<< (const double y) {
  this->y.push_back(y);
  a.push_back(0.0);
  b.push_back(0.0);
  c.push_back(0.0);
  d.push_back(0.0);
  
  da.push_back(0.0);
  db.push_back(0.0);
  dc.push_back(0.0);
  
  dda.push_back(0.0);
  ddb.push_back(0.0);
  
  return *this;
}

EPH_Spline& EPH_Spline::operator<< (const bool init) {
  if(!init) {
    x_Last = x_First - 0.1;
    y.clear();
    
    a.clear();
    b.clear();
    c.clear();
    d.clear();
    
    da.clear();
    db.clear();
    dc.clear();
    
    dda.clear();
    ddb.clear();
    
    return *this;
  }
  
  if((this->y).size() > min_size) FindCoefficients();
  else throw std::runtime_error("eph_spline: not enough points supplied");
  
  return *this; 
}

// calculate coeffiecients based on x0, x1, ,y0, y1, y0' and y1'
void EPH_Spline::FindCoefficients() {
  unsigned int points = y.size();
  
  if(!(dx > 0.0)) throw std::runtime_error("eph_spline: negative incerement step provided");
  
  x_Last = x_First + points * dx;
  
  // we use da, db, and dc as temporary buffers
  double z0; // z_-2
  double z1; // z_-1
  
  double z2; // z_k-1
  double z3; // z_k
  
  for(int i = 0; i < points-1; ++i) {
    //da -> z
    da[i] = (y[i+1]-y[i])/dx;
  }
  
  z1 = 2.0*da[0] - da[1];
  z0 = 2.0*z1 - da[0];
  
  z2 = 2.0*da[points-2] - da[points-1];
  z3 = 2.0*z2 - da[points-1];
  
  da[points-1] = z2;
  
  for(int i = 2; i < points-2; ++i) {
    //db -> w_i-1 ; dc -> w_i

    db[i] = fabs(da[i+1] - da[i]);
    dc[i] = fabs(da[i-1] - da[i-2]);
  }
  
  // special cases
  db[0] = fabs(da[1]-da[0]);
  dc[0] = fabs(z1-z0);
  
  db[1] = fabs(da[2]-da[1]);
  db[1] = fabs(da[0]-z1);
  
  db[points-2] = fabs(z2-da[points-2]);
  dc[points-2] = fabs(da[points-3]-da[points-4]);
  
  db[points-1] = fabs(z3-z2);
  dc[points-1] = fabs(da[points-2]-da[points-3]);
  
  //derivatives
  for(unsigned int i = 0; i < points; ++i) {
    double w0,w1;
    double d_2, d_1, d0, d1;
    
    if(i == 0) {
      d_2 = z0;
      d_1 = z1;
      d1 = da[i+1];
    }
    else if(i == 1) {
      d_2 = z1;
      d_1 = da[i-1];
      d1 = da[i+1];
    }
    else {
      d_2 = da[i-2];
      d_1 = da[i-1];
      d1 = da[i+1];
    }
    
    d0 = da[i];
    w1 = db[i];
    w0 = dc[i];
    
    // special cases
    if(d_2 == d_1 && d0 != d1) {
      ddb[i] = d_1;
    }
    else if(d0 == d1 && d_2 == d_1) {
      ddb[i] = d0;
    }
    else if(d_1 == d0) {
      ddb[i] = d0;
    }
    else if(d_2 == d_1 && d0 == d1 && d0 != d_1) {
      ddb[i] = 0.5*(d_1 + d0);
    }
    else {
      ddb[i] = (d_1*w1 + d0*w0)/(w1+w0);
    }
  }
  
  // solve the equations
  for(unsigned int i = 0; i < points-1; ++i) {
    double dx3 = dx*dx*dx;
    
    double x0_1 = (x_First + i*dx);
    double x0_2 = (x_First + i*dx)*x0_1;
    double x0_3 = (x_First + i*dx)*x0_2;
    
    double x1_1 = (x_First + (i+1)*dx);
    double x1_2 = (x_First + (i+1)*dx)*x1_1;
    double x1_3 = (x_First + (i+1)*dx)*x1_2;
    
    d[i] = (-ddb[i]*x0_1-ddb[i+1]*x0_1+ddb[i]*x1_1+ddb[i+1]*x1_1+2.0*y[i]-2.0*y[i+1])/dx3;
    c[i] = (-ddb[i]+ddb[i+1]+3.0*d[i]*x0_2-3.0*d[i]*x1_2)/2.0/dx;
    b[i] = (c[i]*x0_2+d[i]*x0_3-c[i]*x1_2-d[i]*x1_3-y[i]+y[i+1])/dx;
    a[i] = y[i] - b[i]*x0_1 - c[i]*x0_2 - d[i]*x0_3;
  }
  
  a[points-1] = y[points-1];
  b[points-1] = 0.0;
  c[points-1] = 0.0;
  d[points-1] = 0.0;
  
  // initialise other parts
  for(unsigned int i = 0; i < points; ++i) {
    da[i] = b[i];
    db[i] = 2.0*c[i];
    dc[i] = 3.0*d[i];
    
    dda[i] = db[i];
    ddb[i] = 2.0*dc[i];
  }
}
