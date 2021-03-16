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
 * Stripped down version of beta(rho) class
 *
 */

// TODO: consider storing data as alpha instead of beta to reduce the number of sqrt
template<typename Float = double, template<typename> class Allocator = std::allocator, template <typename _F = Float, typename _A = Allocator<Float>> class Container = std::vector>
class EPHKappa {
  public:
    using Spline = EPH_Spline<Float, Allocator, Container>;
    using Container_Float = Container<Float, Allocator<Float>>;

    EPHKappa() :
      n_elements {0},
      r_cutoff {0},
      rho_cutoff {0}
      {}

    EPHKappa(const char* file) {
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
      alpha.resize(n_elements);
      beta.resize(n_elements);
      rho_r_sq.resize(n_elements);

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
      Float dr_sq;
      Float drho;

      fd >> n_points_rho;
      fd >> dr;
      fd >> n_points_beta;
      fd >> drho;
      fd >> r_cutoff;

      r_cutoff_sq = r_cutoff * r_cutoff;
      rho_cutoff = drho * (n_points_beta-1);

      dr_sq = r_cutoff_sq / (n_points_rho - 1);

      // read spline knots for rho and beta for each element
      for(size_t i = 0; i < n_elements; ++i) {
        // workaround to read an uint8_t
        unsigned short val; // TODO: change this;
        fd >> val;
        element_number[i] = val;

        Container_Float l_rho(n_points_rho);
        for(size_t j = 0; j != n_points_rho; ++j)
          fd >> l_rho[j];

        rho[i] = Spline(dr, l_rho);

        // create square version
        for(size_t j = 0; j != n_points_rho; ++j)
          l_rho[j] = rho[i](sqrt(j * dr_sq));

        rho_r_sq[i] = Spline(dr_sq, l_rho);

        Container_Float l_beta(n_points_beta);
        for(size_t j = 0; j != n_points_beta; ++j)
          fd >> l_beta[j];

        beta[i] = Spline(drho, l_beta);

        // create alpha from beta
        for(size_t j = 0; j != n_points_beta; ++j)
          l_beta[j] = sqrt(l_beta[j]);

        alpha[i] = Spline(drho, l_beta);
      }

      fd.close();
    }

    size_t get_n_elements() const {
      return n_elements;
    }

    Float get_r_cutoff() const {
      return r_cutoff;
    }

    Float get_r_cutoff_sq() const {
      return r_cutoff_sq;
    }

    Float get_rho_cutoff() const {
      return rho_cutoff;
    }

    uint8_t get_element_number(size_t index) const {
      assert(index < n_elements);
      return element_number[index];
    }

    std::string get_element_name(size_t index) const {
      assert(index < n_elements);

      return element_name[index];
    }

    Float get_rho(size_t index, Float r) const {
      assert(index < n_elements);
      assert(r < r_cutoff);

      return rho[index](r);
    }

    Float get_rho_r_sq(size_t index, Float r_sq) const {
      assert(index < n_elements);
      assert(r_sq < r_cutoff_sq);

      return rho_r_sq[index](r_sq);
    }

    Float get_beta(size_t index, Float rho_i) const {
      assert(index < n_elements);
      assert(rho_i < rho_cutoff);

      return beta[index](rho_i);
    }

    Float get_alpha(size_t index, Float rho_i) const {
      assert(index < n_elements);
      assert(rho_i < rho_cutoff);

      return alpha[index](rho_i);
    }

  protected:
    static constexpr unsigned int max_line_length = 1024; // this is for parsing

    Float r_cutoff; // cutoff for locality
    Float r_cutoff_sq; // cutoff sq for locality mostly unused
    Float rho_cutoff; // cutoff for largest site density
    size_t n_elements; // number of elements

    Container<uint8_t, Allocator<uint8_t>> element_number;
    Container<std::string, Allocator<std::string>> element_name;
    Container<Spline, Allocator<Spline>> rho;
    Container<Spline, Allocator<Spline>> rho_r_sq;
    Container<Spline, Allocator<Spline>> alpha;
    Container<Spline, Allocator<Spline>> beta;
};

using Kappa = EPHKappa<Float, Allocator, Container>;

#endif
