// This is a fairly generic, nifty class for obtaining parameters from a file.  The file's contents
// are stored in a string, from which the actual values can be extracted with member functions.  


#ifndef _PARAMETERS_
#define _PARAMETERS_

#include <assert.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "paths.hpp"
#include TEMP_TEMPLATES
#include RV_GENERATORS

using namespace std;

namespace evolve{

class Parameters;
std::ostream& operator<<(std::ostream&, const Parameters&);

class Parameters {
public:
  explicit Parameters(std::string filename);

  double get_double(std::string param_name) const;
  int get_int(std::string param_name) const;
  bool get_bool(std::string param_name) const;
  std::string get_string(std::string param_name) const;

  friend std::ostream& operator<<(std::ostream& out, const Parameters&);
private:
  std::string param_string;

  template<typename T>
  void get_param(std::string param_name, T& param) const; 
};

template<typename T>
void Parameters::get_param(std::string param_name, T& param) const {
  std::string::size_type loc = param_string.find(param_name, 0);
  if( loc != std::string::npos ) {
    std::stringstream oss(param_string);
    oss.seekg(loc);
    std::string chunk;
    oss >> chunk >> chunk >> param;
  } 
  else {
    std::cout << "Couldn't find parameter: " << param_name 
	      << " in parameter file." << std::endl << std::endl;
    abort();
  };
};


} // end namespace block

#endif
