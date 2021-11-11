/*
 * Authors of the extension Artur Tamm, Alfredo Correa
 * e-mail: artur.tamm.work@gmail.com
 */

#ifndef EPH_KAPPA
#define EPH_KAPPA

// external headers
#include <memory>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <fstream>
#include <cstddef>

// internal headers
#include "eph_spline.h"
#include "eph_linear.h"

/*
 * this class reads the per atom electronic properties
 *
 */

template<typename Float = double, template<typename> class Allocator = std::allocator, template <typename _F = Float, typename _A = Allocator<Float>> class Container = std::vector>
struct EPH_kappa {
  using Spline = EPH_Spline<Float, Allocator, Container>;
  using Container_Float = Container<Float, Allocator<Float>>;
  using Linear = EPH_Linear;

  static constexpr unsigned int max_line_length = 1024; // this is for parsing

  Float r_cutoff; // cutoff for locality
  Float r_cutoff_sq; // sq version for convienience
  Float T_max;
  size_t n_elements; // number of elements
  size_t n_pairs;

  Container<int, Allocator<int>> element_number;
  Container<std::string, Allocator<std::string>> element_name;
  Container<Spline, Allocator<Spline>> rho_r; // spatial correlation rho(r) [n_elements]
  Container<Spline, Allocator<Spline>> rho_r_sq; // rho(r**2) [n_elements]
  Container<Linear, Allocator<Linear>> E_T_atomic; // T(E) temperature [n_elements]
  Container<Linear, Allocator<Linear>> K_T_atomic; // K(E) conductivity [n_pairs]

  EPH_kappa() :
    n_elements {0},
    n_pairs {0},
    r_cutoff {0},
    r_cutoff_sq {0},
    T_max {0}
    {}

  EPH_kappa(char const* file) {
    std::ifstream fd(file);
    assert(fd.is_open() && "Unable to open input file");
    char line[max_line_length];

    // read first three lines
    // these are comments so we ignore them
    fd.getline(line, max_line_length);
    fd.getline(line, max_line_length);
    fd.getline(line, max_line_length);

    // read the header
    fd >> n_elements;

    assert(n_elements > 0 && "File contains zero elements");

    n_pairs = (n_elements > 1) ? (n_elements + 1) * (n_elements - 1) / 2 : 1;

    element_name.resize(n_elements);
    element_number.resize(n_elements);

    rho_r.resize(n_elements);
    rho_r_sq.resize(n_elements);
    E_T_atomic.resize(n_elements);
    K_T_atomic.resize(n_pairs);

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
    size_t n_points_r, n_points_T;

    Float dr, dr_sq, dT;

    // read general tabulation properties
    fd >> n_points_r;
    fd >> dr;
    fd >> r_cutoff;

    fd >> n_points_T;
    fd >> dT;
    fd >> T_max;

    r_cutoff_sq = r_cutoff * r_cutoff;
    dr_sq = r_cutoff_sq / (static_cast<double>(n_points_r - 1));

    Container_Float _rho_r(n_points_r);
    Container_Float _C_T(n_points_T);
    Container_Float _K_T(n_points_T);

    // read spline knots for rho and beta for each element
    for(size_t i = 0; i < n_elements; ++i) {
      // workaround to read an uint8_t
      int val; // read element number
      fd >> val;
      element_number[i] = val;

      // read locality rho(r)
      for(size_t j = 0; j != n_points_r; ++j) {
        fd >> _rho_r[j];
      }
      rho_r[i] = Spline(dr, _rho_r);

      // create square version
      for(size_t j = 0; j < n_points_r; ++j) {
        _rho_r[j] = rho_r[i](sqrt(j * dr_sq));
      }
      rho_r_sq[i] = Spline(dr_sq, _rho_r);

      // read thermal properties
      for(size_t j = 0; j != n_points_T; ++j) {
        fd >> _C_T[j];
      }
      
      // create Ee(Te) mapping
      _C_T[0] = 0.;
      for(size_t i = 1; i < n_points_T; ++i) {
        _C_T[i] = _C_T[i - 1] + _C_T[i] * dT;
      }
      E_T_atomic[i] = EPH_Linear(dT, _C_T.begin(), _C_T.end());
    }

    // read kappa(E)
    for(size_t i = 0; i < n_pairs; ++i) {
      for(size_t j = 0; j < n_points_T; ++j) {
        fd >> _K_T[j];
      }
      
      K_T_atomic[i] = EPH_Linear(dT, _K_T.begin(), _K_T.end());
    }

    fd.close();
  }

  Spline& get_K_T(int i_type, int j_type) { // this is a convenience function
    int k = i_j_to_k(i_type, j_type, n_elements);
    return K_T_atomic[k];
  }

  static int i_j_to_k(int i_type, int j_type, int n) { // temporary solution
    int k = 0;

    for(int i = 0; i < n; ++i) {
      for(int j = i; j < n; ++j) {
        if(i == i_type && j == j_type) { break; }
        ++k;
      }

      if(i == i_type) { break; }
    }

    return k;
  }
};

using Kappa = EPH_kappa<Float, Allocator, Container>;

#endif
