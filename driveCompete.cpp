//#define NDEBUG
#include <ctime>
#include <cmath>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "paths.hpp"     
#include EXPERIMENT
#include POPULATION
#include ORGANISM
#include RV_GENERATORS
#include TEMP_TEMPLATES
#include PARAMETERS

using namespace evolve;
using namespace std;
 
namespace {          
  using namespace evolve;
  evolve::Parameters prm  ("parameters_compete.txt");      // Create parameter object from file
}
 
int main() {

  //clock_t start_time= clock();                        // Start program timing clock
  long seed = time( NULL)+ getpid();                  // Get random number generator seed 
  srand48( seed);                                     // Seed random number generator
  //std::cout<< "seed= "<< seed<< std::endl;
  //cout<<prm;
  
  Organism::add_states(3);   
  
  Organism::set_state_params(0, prm);                           // connect Parameters to Organism

  Organism::set_state_params(1, prm);
  Organism::set_state_params(2, prm);
  int numFix= 0;
  
  for( int itrial= 0; itrial< prm.get_int("trials"); ++itrial){
	  Organism org_w;                                                // Empty genome, pnat_product= 1 
	  Organism org_t;
	  org_t.set_tracked(1);
	  
	  Population pop;                                               // create Population      
	  pop.set_pop_capacity( prm.get_int( "pop_capacity") );
   
		for( int i= 0; i< prm.get_int( "pop_capacity")- prm.get_int( "cells_init_tracked"); ++i)                
		  pop.add_org( org_w, 0);          
		for( int i= 0; i< prm.get_int( "cells_init_tracked"); ++i)
		  pop.add_org( org_t, 1);   
		// ------------------------------------------------------------------------
		Experiment exp;                                                    // Create experiment
		exp.set_population( pop).set_stop_cond( fixed_or_lost);
		exp.start();
		if (exp.population().num_wld_orgs() == 0) ++numFix;
    };
   
    //cout << "Pfix = "<< (double)numFix/prm.get_int("trials")<< endl;
    cout<< (double)numFix/prm.get_int("trials")<< endl;
return 0;

}