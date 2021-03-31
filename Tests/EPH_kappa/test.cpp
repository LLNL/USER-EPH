
#include <iostream>
#include <fstream>
#include <cstdio>

#include "eph_kappa.h"

int main(int args, char **argv) {
  double kB = 8.6173303e-5;

  std::cout << "EPH_kappa class tests" << std::endl;
  { // create an empty eph_kappa object
    Kappa kappa;
  }

  { // load a file and test contents
    EPH_kappa<> kappa("Cu.kappa");
    assert(kappa.n_elements == 1 && "number of elements not 1");

    std::cout << "Elements in the file: " << std::endl;
    for(int i = 0; i < kappa.n_elements; ++i) {
      std::cout << "  " << kappa.element_name[i] << std::endl;
    }
    std::cout << '\n';

    // save interpolated functions into a file
    {
      for(size_t i = 0; i < kappa.n_elements; ++i) { // rho(r)
        char fn[128];
        sprintf(fn, "out/rho_r_%ld.data", i);
        std::ofstream fout(fn);

        double dr = 0.01;
        double r = 0.0;

        while(r < kappa.r_cutoff) {
          fout << r << ' ' << kappa.rho_r[i](r) << '\n';
          r += dr;
        }
      }

      for(size_t i = 0; i < kappa.n_elements; ++i) { // T(E)
        char fn[128];
        sprintf(fn, "out/T_E_%ld.data", i);
        std::ofstream fout(fn);

        double dE = 0.1;
        double E = 0.0;

        while(E < kappa.E_max) {
          fout << E << ' ' << kappa.T_E_atomic[i](E) << '\n';

          E += dE;
        }
      }

      for(size_t i = 0; i < kappa.n_elements; ++i) { // kappa(E)
        for(size_t j = 0; j < kappa.n_elements; ++j) {
          char fn[128];
          sprintf(fn, "out/K_E_%ld_%ld.data", i, j);
          std::ofstream fout(fn);

          double dE = 0.1;
          double E = 0.0;

          while(E < kappa.E_max) {
            fout << E << ' ' << E / kB << ' ' << kappa.get_K_E(i, j)(E) << '\n';

            E += dE;
          }
        }
      }
    }
  }
  return 0;
}

