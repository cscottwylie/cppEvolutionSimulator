#include <sstream>
#include <string>
#include <vector>

#include "paths.hpp"
#include TEMP_TEMPLATES
#include PARAMETERS

using namespace std;

namespace evolve{
Parameters::Parameters(std::string filename) {
  std::ifstream param_file(filename.c_str());
  if (!param_file) {
    std::cout << "Couldn't open parameter file.";
    abort();
  };
  while (param_file) {
    char next_char;
    param_file.get(next_char);
    param_string += next_char;
  };
};

double Parameters::get_double(std::string param_name) const {
  double param;
  get_param(param_name, param);
  return param;
};

int Parameters::get_int(std::string param_name) const {
  int param;
  get_param(param_name, param);
  return param;
};

bool Parameters::get_bool(std::string param_name) const {
  bool param;
  get_param(param_name, param);
  return param;
};

std::string Parameters::get_string(std::string param_name) const {
  std::string param;
  get_param(param_name, param);
  return param;
};

std::ostream& operator<<(std::ostream& out, const Parameters& prm) {
  out << prm.param_string << std::endl;
  return out;
};

} // end namespace block
