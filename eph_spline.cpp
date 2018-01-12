/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: artur.tamm.work@gmail.com
 */

//#define DEBUG_EPH

// external headers
#include <algorithm>
#include <limits>
#include <cmath>

#ifdef DEBUG_EPH
#include <iostream>
#endif

// internal headers
#include "eph_spline.h"

EPH_Spline::EPH_Spline() {
  this->inited = false;
  this->lastIndex = 0;
  
  this->points = 0;
  
  this->x = nullptr;
  this->y = nullptr;
  
  this->a = nullptr;
  this->b = nullptr;
  this->c = nullptr;
  this->d = nullptr;
  
  this->da = nullptr;
  this->db = nullptr;
  this->dc = nullptr;
  
  this->dda = nullptr;
  this->ddb = nullptr;
}

EPH_Spline::~EPH_Spline() {
  if(x) delete[] x;
  if(y) delete[] y;
  
  if(a) delete[] a;
  if(b) delete[] b;
  if(c) delete[] c;
  if(d) delete[] d;
  
  if(da) delete[] da;
  if(db) delete[] db;
  if(dc) delete[] dc;
  
  if(dda) delete[] dda;
  if(ddb) delete[] ddb;
}

EPH_Spline::EPH_Spline(const EPH_Spline& spline) {
  inited = false;
  dx = 0;
  points = 0;
  lastIndex = 0;
  
  x = nullptr;
  y = nullptr;
  
  a = nullptr;
  b = nullptr;
  c = nullptr;
  d = nullptr;
  
  da = nullptr;
  db = nullptr;
  dc = nullptr;
  
  dda = nullptr;
  ddb = nullptr;
  
  *this = spline;
}

EPH_Spline& EPH_Spline::operator= (const EPH_Spline& spline) {
  points = spline.points;
  inited = spline.inited;
  dx = spline.dx;
  lastIndex = spline.lastIndex;
  
  if(spline.points > 0) {
    if(x) delete[] x;
    x = new double[points];
    std::copy_n(spline.x, points, x);
    
    if(y) delete[] y;
    y = new double[points];
    std::copy_n(spline.y, points, y);
    
    if(a) delete[] a;
    a = new double[points];
    std::copy_n(spline.a, points, a);
    
    if(b) delete[] b;
    b = new double[points];
    std::copy_n(spline.b, points, b);
    
    if(c) delete[] c;
    c = new double[points];
    std::copy_n(spline.c, points, c);
    
    if(d) delete[] d;
    d = new double[points];
    std::copy_n(spline.d, points, d);
    
    if(da) delete[] da;
    da = new double[points];
    std::copy_n(spline.da, points, da);
    
    if(db) delete[] db;
    db = new double[points];
    std::copy_n(spline.db, points, db);
    
    if(dc) delete[] dc;
    dc = new double[points];
    std::copy_n(spline.dc, points, dc);
        
    if(dda) delete[] dda;
    dda = new double[points];
    std::copy_n(spline.dda, points, dda);
    
    if(ddb) delete[] ddb;
    ddb = new double[points];
    std::copy_n(spline.ddb, points, ddb);
  }
  else {
    x = nullptr;
    y = nullptr;
    
    a = nullptr;
    b = nullptr;
    c = nullptr;
    d = nullptr;
    
    da = nullptr;
    db = nullptr;
    dc = nullptr;
    
    dda = nullptr;
    ddb = nullptr;
  }
}

bool EPH_Spline::InitSpline(const double* x, const double* y, const unsigned int points) {
  if(!x || !y) return false;
  if(points < 2) return false;
  
  this->points = points;
  
  if(this->x) delete[] this->x;
  if(this->y) delete[] this->y;
  
  this->x = new double[points];
  this->y = new double[points];
  
  std::copy_n(x, points, this->x);
  dx = this->x[1] - this->x[0];
  
  std::copy_n(y, points, this->y);
  
  if(a) delete[] a;
  if(b) delete[] b;
  if(c) delete[] c;
  if(d) delete[] d;
  
  if(da) delete[] da;
  if(db) delete[] db;
  if(dc) delete[] dc;
  
  if(dda) delete[] dda;
  if(ddb) delete[] ddb;
  
  a = new double[points];
  b = new double[points];
  c = new double[points];
  d = new double[points];
  
  da = new double[points];
  db = new double[points];
  dc = new double[points];
  
  dda = new double[points];
  ddb = new double[points];
  
  // do the magic
  // we use da, db, and dc as temporary buffers
  double z0; // z_-2
  double z1; // z_-1
  
  double z2; // z_k-1
  double z3; // z_k
  
  for(int i = 0; i < points-1; ++i) {
    //da -> z
    da[i] = (y[i+1]-y[i])/(x[i+1]-x[i]);
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
  
  FindCoefficients();
  return true;
}

// cubic spline
double EPH_Spline::GetValue(const double x) {
  double result = 0.0;
  unsigned int index = 0;
  
  if(x < this->x[0])
    return std::numeric_limits<double>::quiet_NaN();
  else if(x > this->x[points-1])
    return std::numeric_limits<double>::quiet_NaN();
  
  index = FindIndex(x);
  
  double x2 = x*x;
  double x3 = x*x2;
  
  result = a[index] + b[index] * x + c[index] * x2 + d[index] * x3;
  
  return result;
}

double EPH_Spline::GetDValue(const double x) {
  double result = 0.0;
  unsigned int index = 0;
  
  if(x < this->x[0])
    return std::numeric_limits<double>::quiet_NaN();
  else if(x > this->x[points-1])
    return std::numeric_limits<double>::quiet_NaN();
  
  index = FindIndex(x);
  
  double x2 = x*x;
  double x3 = x*x2;
  
  result = da[index] + db[index] * x + dc[index] * x2;
  
  return result;
}

double EPH_Spline::GetDDValue(const double x) {
  double result = 0.0;
  unsigned int index = 0;
  
  if(x < this->x[0])
    return std::numeric_limits<double>::quiet_NaN();
  else if(x > this->x[points-1])
    return std::numeric_limits<double>::quiet_NaN();
  
  index = FindIndex(x);
  
  double x2 = x*x;
  double x3 = x*x2;
  
  result = dda[index] + ddb[index] * x;
  
  return result;
}

// calculate coeffiecients based on x0, x1, ,y0, y1, y0' and y1'
void EPH_Spline::FindCoefficients() {
  // solve the equations
  for(unsigned int i = 0; i < points-1; ++i) {
    double dx = x[i+1]-x[i];
    double dx3 = dx*dx*dx;
    
    double x0_2 = x[i]*x[i];
    double x0_3 = x[i]*x0_2;
    
    double x1_2 = x[i+1]*x[i+1];
    double x1_3 = x[i+1]*x1_2;
    
    d[i] = (-ddb[i]*x[i]-ddb[i+1]*x[i]+ddb[i]*x[i+1]+ddb[i+1]*x[i+1]+2.0*y[i]-2.0*y[i+1])/dx3;
    c[i] = (-ddb[i]+ddb[i+1]+3.0*d[i]*x0_2-3.0*d[i]*x1_2)/2.0/dx;
    b[i] = (c[i]*x0_2+d[i]*x0_3-c[i]*x1_2-d[i]*x1_3-y[i]+y[i+1])/dx;
    a[i] = y[i] - b[i]*x[i] - c[i]*x0_2 - d[i]*x0_3;
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

