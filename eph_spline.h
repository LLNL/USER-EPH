/*
 * Authors of the extension Artur Tamm, Alfredo Caro, Alfredo Correa, Mattias Klintenberg
 * e-mail: arturt@ut.ee
 */

#ifndef EPH_SPLINE
#define EPH_SPLINE

class EPH_Spline {
  public:
  private:
    bool inited;
    
    unsigned int points; // number of points in the spline
    double* x; // x values
    double* y; // y values
    
    double* a; // splines are always in a + b*x + c*x^2 + d*x^3 form
    double* b;
    double* c;
    double* d;
    
    double* da; // this is for the derivative
    double* db;
    double* dc;
    
    double* dda; // this is for the second derivative
    double* ddb;
    
    double dx; 
    
    unsigned int lastIndex;
    
  public:
    EPH_Spline(); // default constructor
    ~EPH_Spline(); // destructor
    
    EPH_Spline& operator= (const EPH_Spline& spline);
    
    //TODO: copy constructor; maybe not needed
    EPH_Spline(const EPH_Spline& spl);
    
    bool InitSpline(const double *x, const double *y, const unsigned int points); // initialise spline
    double GetValue(const double x);
    double GetDValue(const double x);
    double GetDDValue(const double x);
  
  private: // some private functions for spline initialisation and value calculation
    unsigned int FindIndex(const double x) {
      unsigned int index = 0;
      index = (unsigned int) ((x-this->x[0])/dx);
      return index;
    }
    
    void FindCoefficients();
};

#endif
