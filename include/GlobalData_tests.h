#ifndef GLOBALDATA_TESTS_H_
#define GLOBALDATA_TESTS_H_

#include <array>
#include <iostream>
#include <string>

#include "Operators.h"
#include "typedefs.h"

namespace GlobalData_tests {
  
  // simplify and clean read_parameters function
  void lattice_input_data_handling(const std::string path_output,
                                   const std::string name_lattice, 
                                   const std::string path_config, int Lt, 
                                   int Lx, int Ly, int Lz);
  void eigenvec_perambulator_input_data_handling(const int number_of_eigen_vec, 
      const std::string path_eigenvectors, const std::string name_eigenvectors, 
      const std::string path_perambulators, 
      const std::string name_perambulators);
  void config_input_data_handling(const int start_config, const int end_config, 
                                  const int delta_config);
  void quark_check(quark quarks);

}

#endif // GLOBALDATA_TESTS_H_
