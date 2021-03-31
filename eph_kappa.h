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

/*
 * this class reads the per atom electronic properties
 *
 */

template<typename Float = double, template<typename> class Allocator = std::allocator, template <typename _F = Float, typename _A = Allocator<Float>> class Container = std::vector>
struct EPH_kappa {
  using Spline = EPH_Spline<Float, Allocator, Container>;
  using Container_Float = Container<Float, Allocator<Float>>;

  static constexpr unsigned int max_line_length = 1024; // this is for parsing

  Float r_cutoff; // cutoff for locality
  Float r_cutoff_sq; // sq version for convienience
  Float T_max;
  size_t n_elements; // number of elements

  Container<int, Allocator<int>> element_number;
  Container<std::string, Allocator<std::string>> element_name;
  Container<Spline, Allocator<Spline>> rho;
  Container<Spline, Allocator<Spline>> rho_sq;
  Container<Spline, Allocator<Spline>> E_e_atomic;
  Container<Spline, Allocator<Spline>> C_e_atomic;
  Container<Spline, Allocator<Spline>> kappa_e_atomic;

  EPH_kappa() :
    n_elements {0},
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

    element_name.resize(n_elements);
    element_number.resize(n_elements);

    rho.resize(n_elements);
    E_e_atomic.resize(n_elements);
    C_e_atomic.resize(n_elements);
    kappa_e_atomic.resize(n_elements);

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
    size_t n_points_rho, n_points_T;

    Float dr, dr_sq, dT;

    // read general tabulation properties
    fd >> n_points_rho;
    fd >> dr;
    fd >> r_cutoff;
    fd >> n_points_T;
    fd >> dT;
    fd >> T_max;

    r_cutoff_sq = r_cutoff * r_cutoff;
    dr_sq = r_cutoff_sq / (n_points_rho - 1);

    Container_Float _rho_r(n_points_rho);
    Container_Float _E_e_T(n_points_T);
    Container_Float _C_e_T(n_points_T);
    Container_Float _kappa_e_T(n_points_T);

    // read spline knots for rho and beta for each element
    for(size_t i = 0; i < n_elements; ++i) {
      // workaround to read an uint8_t
      int val; // read element number
      fd >> val;
      element_number[i] = val;

      // read locality rho(r)
      for(size_t j = 0; j != n_points_rho; ++j) {
        fd >> _rho_r[j];
      }
      rho[i] = Spline(dr, _rho_r);

      // create square version
      for(size_t j = 0; j != n_points_rho; ++j) {
        _rho_r[j] = rho[i](sqrt(j * dr_sq));
      }
      rho_sq[i] = Spline(dr_sq, _rho_r);

      // read thermal properties
      for(size_t j = 0; j != n_points_T; ++j) {
        fd >> _E_e_T[j] >> _C_e_T[j] >> _kappa_e_T[j];
      }

      E_e_atomic[i] = Spline(dT, _E_e_T);
      C_e_atomic[i] = Spline(dT, _C_e_T);
      kappa_e_atomic[i] = Spline(dT, _kappa_e_T);
    }
    fd.close();
  }

};

using Kappa = EPH_kappa<Float, Allocator, Container>;

#endif
