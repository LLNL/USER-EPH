
#include <iostream>
#include <fstream>

#include "eph_kappa.h"

int main(int args, char **argv) {
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
      { // rho(r)
        std::ofstream fout("out/rho_r.data");

        double dr = 0.01;
        double r = 0.0;

        while(r < kappa.r_cutoff) {
          fout << r << ' ' << kappa.rho[0](r) << '\n';

          r += dr;
        }
      }

      { // E_e(T)
        std::ofstream fout("out/E_e_T.data");

        double dT = 0.1;
        double T = 0.0;

        while(T < kappa.T_max) {
          fout << T << ' ' << kappa.E_e_atomic[0](T) << '\n';

          T += dT;
        }
      }

      { // C_e(T)
        std::ofstream fout("out/C_e_T.data");

        double dT = 0.1;
        double T = 0.0;

        while(T < kappa.T_max) {
          fout << T << ' ' << kappa.C_e_atomic[0](T) << '\n';

          T += dT;
        }
      }

      { // kappa_e(T)
        std::ofstream fout("out/kappa_e_T.data");

        double dT = 0.1;
        double T = 0.0;

        while(T < kappa.T_max) {
          fout << T << ' ' << kappa.kappa_e_atomic[0](T) << '\n';

          T += dT;
        }
      }
    }
  }
  return 0;
}

