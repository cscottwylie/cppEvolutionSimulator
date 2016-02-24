// function definitions for Org_state, Organism_data, and Organism classes

#include <iostream>
#include <cmath>
#include <vector>
#include <assert.h>
#include "paths.hpp"
#include ORGANISM
#include RV_GENERATORS

namespace evolve{

ptrdiff_t myrandom (ptrdiff_t i) { return rnd_int(i);} //namespace scope function used below
std::vector<Org_state> Organism::state_list;           // namespace scope static object

Org_state& Organism::state(int st) {
  assert(state_list.size() > 0);
  assert(st >= 0);
  assert(st < (int) state_list.size()); 
  return state_list[st];
};


// *********************** Organism writing functions ***********************
// ** each calls make_write_safe() which reserves memory for new Org_data  **
// ** this means you should only call these when a new lineage in created  **
// **************************************************************************


Organism& Organism::set_tracked(bool trk) {
  if (trk != tracked()) {
    make_write_safe(data_ptr);
    data_ptr->set_tracked(trk);
  };           
  return *this;   
};

// ******************************** Mutation ******************************** 
/*void Organism::mutate(int st) {
  int up_muts   = rnd_binomial(state(st).up_mut_prob()  ,num_zeros());
  int down_muts = rnd_binomial(state(st).down_mut_prob(),num_ones());
  if ((up_muts > 0) or (down_muts > 0)) {
    make_write_safe(data_ptr);                    // Reserve memory for new lineage if needed
    data_ptr->flip_genes(up_muts, down_muts);     // Now perform genomic changes
  };
};*/

bool Organism::mutate( int st) {
  double mut_prob_ben= 1- exp( -state( st).mut_rate_ben() );
  double mut_prob_del= 1- exp( -state( st).mut_rate_del() );
  int allele_change= 0;
  
  if( allele_state() == 1)
    allele_change-= rnd_uniform() < mut_prob_del;
  
  else if ( allele_state() == 0)  // deleterious genotypes can't mutate in Ben Allen's model
    allele_change+= ( rnd_uniform() < mut_prob_ben) - ( rnd_uniform() < mut_prob_del );
    
  if( allele_change != 0) {
    make_write_safe( data_ptr);
    data_ptr->set_allele_state( allele_state() + allele_change); 
  };
  return 0;   // this is a relic from when lethal mutations were implemented
};


// ************************** Organism constructor  **************************
// ** new memory reserved each time called.  copy constructor will often be called
// ** implicitly by vector::push_back(), which will NOT reserve new memory

Organism::Organism() : data_ptr(new Organism_data()) {};  


// ***********************  Organism data constructor ***********************
Organism_data::Organism_data() 
  : is_tracked(false),
    n_in_lineage(0),
    line_index(-1),                   // -1's below mean "not yet assigned by pop."
    n_in_state(Organism::num_states()),
    allele(0) {}; // zero vector of length num_states() 


// ********************** Organism-property constructor ***********************
Org_state::Org_state() 
  : mt_rate_b(0.0),
    mt_rate_d(0.0),
    b_pre(0.0),
    c_rate (0.0),
    f_adjust(0.0),
    d_rate(0.0){};

// ******************* Organism-property setting functions  *******************
// all of these called together in Organism::set_state_params()

Org_state& Org_state::set_mut_rate_ben(double ben_mut_rate) {
  assert(ben_mut_rate >= 0.0);
  mt_rate_b = ben_mut_rate;
  return *this;
};

Org_state& Org_state::set_mut_rate_del(double del_mut_rate) {
  assert(del_mut_rate >= 0.0);
  mt_rate_d = del_mut_rate;
  return *this;
};

Org_state& Org_state::set_sel_coeff_ben(double sel_coeff_ben){
  assert(sel_coeff_ben >= 0.0);
  s_b= sel_coeff_ben;
  return *this;
};

Org_state& Org_state::set_sel_coeff_del(double sel_coeff_del){
  assert(sel_coeff_del >= 0.0);
  assert(sel_coeff_del <= 1.0);
  s_d= sel_coeff_del;
  return *this;
};

Org_state& Org_state::set_birth_prefactor(double birth_prefactor) {
  assert(birth_prefactor >= 0.0);
  b_pre = birth_prefactor;
  return *this;
};

Org_state& Org_state::set_chg_rate(double change_rate) {
  assert(change_rate >= 0.0);
  c_rate = change_rate;
  return *this;
};


Org_state& Org_state::set_death_rate(double dth_rate){
  d_rate = dth_rate;
  return *this;
};

void Organism::add_states(int num_states) {
  for(int i=0; i<num_states; ++i) state_list.push_back(Org_state());
};


// ********************** Organism state screen output  ***********************
std::ostream& operator<<(std::ostream& out, const Org_state& props) {
  out << "mut_rate_ben    = " << props.mut_rate_ben()    << std::endl
      << "mut_rate_del    = " << props.mut_rate_del()    << std::endl
      << "sel_coeff_ben   = " << props.sel_coeff_ben()   << std::endl
      << "sel_coeff_del   = " << props.sel_coeff_del()   << std::endl
      << "birth_prefactor = " << props.birth_prefactor() << std::endl
      << "change_rate     = " << props.chg_rate()        << std::endl;
  return out;
};


// ********************* Equality operators for organisms *********************
bool operator==(const Organism& org_a, const Organism& org_b) {
  return org_a.get_ptr() == org_b.get_ptr();
};

bool operator!=(const Organism& org_a, const Organism& org_b) {
  return org_a.get_ptr() != org_b.get_ptr();
};

// ************************ Print organism to output ************************
std::ostream& operator<<(std::ostream& out, const Organism& org) {
  out << "org data_ptr   = " << org.get_ptr()   << std::endl
      << "allele         = " << org.allele_state() <<std::endl;
  out << "tracked        = " << org.tracked()   << std::endl
      << "num_in_lineage = " << org.num_in_lineage() << std::endl
      << "num_in_state   = [";
  for (int st = 0; st < Organism::num_states(); ++st) 
    out << " " << org.num_in_state(st);
  out << " ]" << std::endl
      << "lineage_index  = " << org.lineage_index() << std::endl;
  return out;
};

void Organism::set_state_params(int st, const Parameters& prm) {
  std::string state_number;
  string_format(state_number, st);
  std::string suffix = "_s";
  suffix.append(state_number);

  std::string pname_mut_ben     = "mut_ben";
  std::string pname_mut_del     = "mut_del";
  std::string pname_s_ben           = "s_ben";
  std::string pname_s_del           = "s_del";
  std::string pname_birth_prefactor = "birth_prefactor";
  std::string pname_log_chg_rate    = "log_chg_rate";
  std::string pname_death_rate      = "death_rate";
  
  pname_mut_ben.append(suffix);
  pname_mut_del.append(suffix);
  pname_s_ben.append(suffix);
  pname_s_del.append(suffix);
  pname_birth_prefactor.append(suffix);
  pname_log_chg_rate.append(suffix);
  pname_death_rate.append(suffix);

  double mut_ben     = prm.get_double(pname_mut_ben); 
  double mut_del     = prm.get_double(pname_mut_del);      
  double s_ben           = prm.get_double(pname_s_ben);
  double s_del           = prm.get_double(pname_s_del);
  double birth_prefactor = prm.get_double(pname_birth_prefactor);
  double log_chg_rate    = prm.get_double(pname_log_chg_rate);
  double death_rate      = prm.get_double(pname_death_rate);
  
  Organism::state(st)      
    .set_mut_rate_ben   (mut_ben  )
    .set_mut_rate_del   (mut_del  )
    .set_sel_coeff_ben  (s_ben)
    .set_sel_coeff_del  (s_del)
    .set_birth_prefactor( birth_prefactor)
    .set_chg_rate       (pow(10.0, log_chg_rate))
    .set_death_rate     (death_rate);
};

void Organism::write_states() {
  for(int st=0; st < Organism::num_states(); ++st) {
    std::cout << "Organism State " << st << " :" << std::endl
	      << Organism::state(st) << std::endl;
  };
};



}
