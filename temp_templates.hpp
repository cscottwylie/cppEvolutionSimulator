// Declarations *and definitions* of template functions used throughout evolution code.  Mutiple
// defintions (via includes in separate compilation units) are ok with g++ b/c these are templates.  

#ifndef _TEMP_TEMPLATES_
#define _TEMP_TEMPLATES_

#include <assert.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

namespace evolve{

template<typename T>
void swap_pop(std::vector<T>& vec, int i) {
  assert(i >= 0);
  assert(i < (int) vec.size());
  vec[i] = vec.back();
  vec.pop_back();
};

template<typename T>
inline void make_write_safe(boost::shared_ptr<T>& ptr) {
  if (not ptr.unique()) ptr = boost::shared_ptr<T>(new T(*ptr));
};


template<class T>                                // Get const reference
T const& as_const( T const& r ) { return r; };

// Helper function for getting parameters
template<class p_type>
void get_param(p_type& param, std::string p_name, std::ifstream& p_file) {
  std::string s;
  p_file >> s;  //takes the first chunk from p_file
  if (s != p_name) {
    std::cout << "In parameter file, expected to see parameter " 
	 << "'" << p_name << "'" << " in parameter file." << std::endl;
    abort();
  };
  p_file >> s;
  if (s != "=") {
    std::cout << "In parameter file, expected to see an '=' after "
	 << "'" << p_name << "'" << " in parameter file." << std::endl;
    abort();
  };
  p_file >> param;
};


template <typename T>
void string_format(std::string & val, const T & t) {
  std::ostringstream oss; 
  oss << t; 
  val=oss.str(); // extract string and assign it to val
}

template <typename T>
void append_param(std::string& filename,
		  const std::string& label,
		  const T & val) {
  // Omit parameters with log-value less than the following constant
  const double omit_below = -998.0;
  if (val > omit_below) {
    std::string val_string;
    string_format(val_string, val); 
    filename += label + val_string;
  };
};


}

#endif
