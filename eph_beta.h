/*
 * Authors of the extension Artur Tamm, Alfredo Correa
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_BETA
#define EPH_BETA

// external headers 
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <fstream>

// internal headers
#include "eph_spline.h"

/* 
 * Stripped down version of beta(rho) class
 * 
 */

// TODO: consider storing data as alpha instead of beta to reduce the number of sqrt

#if 0
template<typename Container = std::vector<double>>
class EPH_Beta {
  using Float = typename Container::value_type;
  using Spline = EPH_Spline<Container>;
  
  public:
    EPH_Beta() : 
      n_elements {0},
      r_cutoff {0},
      rho_cutoff {0}
      {}
    
    EPH_Beta(const char* file) {
      std::ifstream fd(file);
  
      assert(fd.is_open());
      
      char line[max_line_length];
      
      // read first three lines
      // these are comments so we ignore them
      fd.getline(line, max_line_length);
      fd.getline(line, max_line_length);
      fd.getline(line, max_line_length);
      
      // read the header
      fd >> n_elements;
      
      assert(n_elements > 0);
      
      element_name.resize(n_elements);
      element_number.resize(n_elements);
      rho.resize(n_elements);
      beta.resize(n_elements);
      
      // read the number of elements and their names
      fd.getline(line, max_line_length);
      std::string str(line);
      std::istringstream strstream(str);
      
      for(size_t i = 0; i < n_elements; ++i) {
        std::string elem;
        strstream >> elem;
        element_name[i] = elem;
      }
      
      // read spline parameters
      size_t n_points_rho;
      size_t n_points_beta;
      
      Float dr;
      Float drho;
      
      fd >> n_points_rho;
      fd >> dr;
      fd >> n_points_beta;
      fd >> drho;
      fd >> r_cutoff;
      rho_cutoff = drho * (n_points_beta-1);
      
      // read spline knots for rho and beta for each element
      for(size_t i = 0; i < n_elements; ++i) {
        fd >> element_number[i];
        
        Container l_rho(n_points_rho);
        for(size_t j = 0; j < n_points_rho; ++j)
          fd >> l_rho[j];
        
        rho[i] = Spline(dr, l_rho);
        
        Container l_beta(n_points_beta);
        for(size_t j = 0; j < n_points_beta; ++j)
          fd >> l_beta[i];
        
        beta[i] = Spline(drho, l_beta);
      }

      fd.close();
    }
    
    size_t get_n_elements() const {
      return n_elements;
    }
    
    Float get_r_cutoff() const {
      return r_cutoff;
    }
    
    Float get_rho_cutoff() const {
      return rho_cutoff;
    }
    
    size_t get_element_number(size_t i_index) const {
      assert(i_index < n_elements);
      return element_number[i_index];
    }
    
    std::string get_element_name(size_t i_index) const {
      assert(i_index < n_elements);
      return element_name[i_index];
    }
    
    Float get_rho(size_t i_index, Float i_r) const {
      assert(i_index < n_elements);
      assert(i_r < r_cutoff);
      
      return rho[i_index](i_r);
    }
    
    Float get_beta(size_t i_index, Float i_rho) const {
      assert(i_index < n_elements);
      assert(i_rho < rho_cutoff);
      
      return beta[i_index](i_rho);
    }
  
  private:
    static constexpr unsigned int max_line_length = 1024; // this is for parsing
    
    Float r_cutoff; // cutoff for locality
    Float rho_cutoff; // cutoff for largest site density
    size_t n_elements; // number of elements
    
    std::vector<uint8_t> element_number;
    std::vector<std::string> element_name;
    std::vector<Spline> rho;
    std::vector<Spline> beta;
};
#endif

#if 1
template<typename Float = double, template<typename> class Allocator = std::allocator, template <typename _F = Float, typename _A = Allocator<Float>> class Container = std::vector>
class EPH_Beta {
  using Spline = EPH_Spline<Float, Allocator, Container>;
  using Container_Float = Container<Float, Allocator<Float>>;
  public:
    EPH_Beta() : 
      n_elements {0},
      r_cutoff {0},
      rho_cutoff {0}
      {}
    
    EPH_Beta(const char* file) {
      std::ifstream fd(file);
  
      assert(fd.is_open());
      
      char line[max_line_length];
      
      // read first three lines
      // these are comments so we ignore them
      fd.getline(line, max_line_length);
      fd.getline(line, max_line_length);
      fd.getline(line, max_line_length);
      
      // read the header
      fd >> n_elements;
      
      assert(n_elements > 0);
      
      element_name.resize(n_elements);
      element_number.resize(n_elements);
      rho.resize(n_elements);
      beta.resize(n_elements);
      
      // read the number of elements and their names
      fd.getline(line, max_line_length);
      std::string str(line);
      std::istringstream strstream(str);
      
      for(size_t i = 0; i < n_elements; ++i) {
        std::string elem;
        strstream >> elem;
        element_name[i] = elem;
      }
      
      // read spline parameters
      size_t n_points_rho;
      size_t n_points_beta;
      
      Float dr;
      Float drho;
      
      fd >> n_points_rho;
      fd >> dr;
      fd >> n_points_beta;
      fd >> drho;
      fd >> r_cutoff;
      rho_cutoff = drho * (n_points_beta-1);
      
      // read spline knots for rho and beta for each element
      for(size_t i = 0; i < n_elements; ++i) {
        fd >> element_number[i];
        
        Container_Float l_rho(n_points_rho);
        for(size_t j = 0; j < n_points_rho; ++j)
          fd >> l_rho[j];
        
        rho[i] = Spline(dr, l_rho);
        
        Container_Float l_beta(n_points_beta);
        for(size_t j = 0; j < n_points_beta; ++j)
          fd >> l_beta[i];
        
        beta[i] = Spline(drho, l_beta);
      }

      fd.close();
    }
    
    size_t get_n_elements() const {
      return n_elements;
    }
    
    Float get_r_cutoff() const {
      return r_cutoff;
    }
    
    Float get_rho_cutoff() const {
      return rho_cutoff;
    }
    
    size_t get_element_number(size_t i_index) const {
      assert(i_index < n_elements);
      return element_number[i_index];
    }
    
    std::string get_element_name(size_t i_index) const {
      assert(i_index < n_elements);
      return element_name[i_index];
    }
    
    Float get_rho(size_t i_index, Float i_r) const {
      assert(i_index < n_elements);
      assert(i_r < r_cutoff);
      
      return rho[i_index](i_r);
    }
    
    Float get_beta(size_t i_index, Float i_rho) const {
      assert(i_index < n_elements);
      assert(i_rho < rho_cutoff);
      
      return beta[i_index](i_rho);
    }
  
  private:
    static constexpr unsigned int max_line_length = 1024; // this is for parsing
    
    Float r_cutoff; // cutoff for locality
    Float rho_cutoff; // cutoff for largest site density
    size_t n_elements; // number of elements
    
    Container<uint8_t, Allocator<uint8_t>> element_number;
    Container<std::string, Allocator<std::string>> element_name;
    Container<Spline, Allocator<Spline>> rho;
    Container<Spline, Allocator<Spline>> beta;
};
#endif

#if 0
#include <stdexcept>

class EPH_Beta {
  public:
    EPH_Beta(const char* file); // initialise EPH_Beta based on file
    EPH_Beta() = delete;
    
    // how many types are described by this class
    inline int getElementsNumber() const {
      return elementName.size();
    }
    
    // get element name described by the class
    inline std::string getName(const unsigned int element) const {
      #ifndef DNDEBUG
      if(element >= elementName.size()) throw std::runtime_error("eph_beta: out of scope index provided");
      #endif
      return elementName[element];
    }
    
    // get a cutoff used in rho(r)
    inline double getCutoff() const {
      return rcut;
    }
    
    // get a cutoff used in beta(rho)
    inline double getRhoCutoff() const {
      return rhocut;
    }
    
    // get the periodic table number for an element (this is not used)
    inline unsigned int getNumber(const unsigned int element) const {
      #ifndef DNDEBUG
      if(element >= elementNumber.size()) throw std::runtime_error("eph_beta: out of scope index provided");
      #endif
      return elementNumber[element];
    }
    
    // return the density value
    inline double getRho(const unsigned int element, const double r) const {
      #ifndef DNDEBUG
      if(element >= rho.size()) throw std::runtime_error("eph_beta: out of scope index provided");
      #endif
      return rho[element].GetValue(r);
    }
    
    // return the derivative of the density
    inline double getDRho(const unsigned int element, const double r) const {
      #ifndef DNDEBUG
      if(element >= rho.size()) throw std::runtime_error("eph_beta: out of scope index provided");
      #endif
      return rho[element].GetDValue(r);
    }
    
    // get coupling parameter value
    inline double getBeta(const unsigned int element, const double rho) const {
      #ifndef EPH_UNSAFE
      if(element > beta.size()) throw std::runtime_error("eph_beta: out of scope index provided");
      #endif
      return beta[element].GetValue(rho);
    }
    
    // get the derivative of the coupling parameter
    inline double getDBeta(const unsigned int element, const double rho) const {
      #ifndef DNDEBUG
      if(element > beta.size()) throw std::runtime_error("eph_beta: out of scope index provided");
      #endif
      return beta[element].GetDValue(rho);
    }
  
  private:
    static constexpr unsigned int lineLength = 1024; // this is for parsing
    double rcut;
    double rhocut;
    
    std::vector<int> elementNumber;
    std::vector<std::string> elementName;
    std::vector<EPH_Spline> rho;
    std::vector<EPH_Spline> beta;
};
#endif

#endif
