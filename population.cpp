#include <iostream>
#include <vector>
#include <cmath>
#include <assert.h>

#include "paths.hpp"
#include POPULATION
#include ORGANISM
#include RV_GENERATORS
#include TEMP_TEMPLATES

using namespace std;
namespace evolve{

Population::Population()
  : b_rate_tots  (Organism::num_states()),
    sum_sq_b_rates(Organism::num_states()),
    b_rate_ubnds (Organism::num_states()),
    tot_rates    (Organism::num_states()),
    orgs               (Organism::num_states()),
    tot_event_rate     (0.0),
    n_orgs             (0),
    n_births           (0),
    gens               (0.0),
    n_deaths           (0),
    n_state_chg        (0),
    leth_muts          (0),
    n_trk_orgs         (0),
    n_trk_births       (0),
    n_trk_deaths       (0),
    n_trk_state_chg    (0) {};
    
double Population::birth_rate() const {
  double tot = 0;
  for (int i=0; i<Organism::num_states(); ++i) {
    tot += b_rate_tots[i];
  };
  return tot;
};

double Population::sum_squared_birth_rate() const {
  double tot= 0;
  for( int st= 0; st< Organism::num_states(); ++st)
    tot += sum_sq_b_rates[ st];
  return tot;
};



double Population::death_rate() const {
  double tot = 0;
  for (int st=0; st < Organism::num_states(); ++st)
    tot += (orgs[st].size() * Organism::state(st).death_rate());
  return tot;
};

void Population::reset_counts() {
  n_births    = 0;
  n_deaths    = 0;
  n_state_chg = 0;
  n_trk_births    = 0;
  n_trk_deaths    = 0;
  n_trk_state_chg = 0;
};

void Population::add_rates(const Organism& org, int st) {
  assert(st >= 0);
  assert(st < Organism::num_states());

  const Org_state& os = Organism::state(st);
  double fit = org_birth_rate(org, st);
  double tot = ( os.chg_rate() + os.death_rate() + fit);
  tot_event_rate += tot;
  tot_rates[st] += tot;
  b_rate_tots[st] += fit;
  sum_sq_b_rates[ st]+= fit* fit;
  //n_ones        += org.num_ones();
  if (fit > b_rate_ubnds[st]) 
    b_rate_ubnds[st] = fit;
};

void Population::remove_rates(const Organism& org, int st) {
  assert(st >= 0);
  assert(st < Organism::num_states());

  const Org_state& os = Organism::state(st);
  double fit = org_birth_rate(org, st);
  double tot = (os.chg_rate() + os.death_rate() + fit);
  tot_event_rate -= tot;
  b_rate_tots[st] -= fit;
  sum_sq_b_rates[ st]-= fit* fit;
  tot_rates[st] -= tot;

};

void Population::add_org(Organism& org, int st) {
  assert(st >= 0);
  assert(st < (int) orgs.size());
  assert(st < (int) Organism::num_states());

  orgs[st].push_back(org);
  
  add_rates(org, st);         // two important helper functions called here
  add_to_lineage_data(org, st);
  
  ++n_orgs;
  if (org.tracked()) ++n_trk_orgs;
};

void Population::death(int st){
assert(st >=0);
assert(st < Organism::num_states());
assert(orgs[st].size() > 0);  // someone here to kill

uint ch = rnd_int(orgs[st].size());
// Remove dead organisms rates/lineage info
  remove_rates(orgs[st][ch], st);    
  remove_from_lineage_data(orgs[st][ch], st);
  
  --n_orgs;
  ++n_deaths;
  if (orgs[st][ch].tracked()) {
    --n_trk_orgs;
    ++n_trk_deaths;
  };
  
  swap_pop(orgs[st], ch); 
};


void Population::state_changer(int st) {
  // This changer designed for 3 state system
  assert(Organism::num_states() == 3);     
  assert(orgs.size() == 3);                   // really, again checking there're 3 states 
  assert(st != 0);                                // State zero shouldn't switch (by fiat)
  assert(num_in_state( st) > 0);
  
  int ch = rnd_int(orgs[st].size());
  // std::cout << "organism number " << ch << std::endl;
  Organism org = orgs[st][ch];  // same pointer

  // Essentially kill orgs[st][ch], but w/out possibility of removing lineage from progenitor list
  org.dec_num_in_state(st);     
  org.dec_num_in_lineage();     
  remove_rates(org, st);

  assert(ch >= 0);
  assert(ch < (int) orgs[st].size());
  swap_pop(orgs[st], ch);
  --n_orgs;
  
  // this necessary b/c add_org increments n_trk_orgs if tracked
  if (org.tracked()) --n_trk_orgs;

  // Add in new organism in new state
  if (st == 1) add_org(org, 2);     
  else {
    assert(st == 2);
    add_org(org, 1);  
  };

  ++n_state_chg;
  if (org.tracked()) ++n_trk_state_chg;
};

void Population::hack_st_change(int num_to_switch){  
  for (int i = 0; i < num_to_switch; ++i){
    int ind_ch = rnd_int(orgs[0].size());
    Organism newguy = orgs[0][ind_ch];
    
    remove_rates(orgs[0][ind_ch], 0);
    remove_from_lineage_data(orgs[0][ind_ch],0);  
    swap_pop(orgs[0], ind_ch);
    --n_orgs;
    
    newguy.reset_lineage_counts();
    newguy.set_tracked(1);
    add_org(newguy, 1);  
   };
};

void Population::update_birth_ub() {           // explicitly update birth rate upper bounds by state
  for( int st= 0; st< Organism::num_states(); ++st) {
    b_rate_ubnds[ st]= 0.0;
    for( int who= 0; who< orgs[ st].size(); ++who)
      if( org_birth_rate( orgs[ st][ who], st) > b_rate_ubnds[ st] )
        b_rate_ubnds[ st]= org_birth_rate( orgs[ st][ who], st );
  };
};
        

  
void Population::do_event() {
    
  double ch = rnd_uniform() * event_rate();
  int st = 0;
  while ((ch -= tot_rates[st]) > 0) { ++st; };
  assert(st >= 0);  
  assert(st < Organism::num_states());  

  if ((ch += b_rate_tots[st]) > 0) {
    birth(st);
    
  int death_st= 0;
  int death_ch= rnd_int( num_orgs() );
  while( ( death_ch -= orgs[ death_st].size() ) >= 0 ) ++death_st;  
  death( death_st);                       // Moran process: call death after every birth
    
  }
  else if ((ch +=(orgs[st].size() * Organism::state(st).death_rate() )) > 0)
    death(st);
  
  else {
    //std::cout << "state-change in state " << st << " ... ";
    //assert((ch += (orgs[st].size() * Organism::state(st).chg_rate())) > 0);
    state_changer(st);
  };
};

std::ostream& operator<<(std::ostream& out, const Population& pop) {
  out << "|---  Population  -------------------------------------------------|"
      << std::endl
      << "tot_event_rate = " << pop.event_rate()        << std::endl
      << "tot_death_rate   = " << pop.death_rate()          << std::endl
      << "num_orgs       = " << pop.num_orgs()          << std::endl
      << "wld_lineages   = " << pop.num_wld_lineages()  << std::endl
      << "trk_lineages   = " << pop.num_trk_lineages()  << std::endl
      << "num_trk_orgs   = " << pop.num_trk_orgs()      << std::endl
      << "state_b_rates  = [ ";
  for (int i=0; i<Organism::num_states(); ++i)
    out << pop.b_rate_tots[i] << " ";
  out << "]" << std::endl;
  out << "state_b_ubnds  = [ ";
  for (int i=0; i<Organism::num_states(); ++i)
    out << pop.b_rate_ubnds[i] << " ";
  out << "]" << std::endl;
  out << "tot_rates= [ ";
  for (int i=0; i<Organism::num_states(); ++i)
    out << pop.tot_rates[i] << " ";
  out << "]" << std::endl;
  for(uint i=0; i<pop.orgs.size(); ++i) {
    out << "--------  "
	<< endl<<"Organisms in state " << i  
	<< "  -----------------------------" << std::endl;
    for(uint j=0; j<pop.orgs[i].size(); ++j) { 
      out << pop.orgs[i][j];
      out << "----------------------------------------------------------"
	  << std::endl;
    };
  };
  out << "---- Tracked Lineages ---------" << std::endl;
  for(int i=0; i < pop.num_trk_lineages(); ++i) {
    out << pop.trk_lines[i] << std::endl;
  };
  out << "---- Wild Lineages ---------" << std::endl;
  for(int i=0; i < pop.num_wld_lineages(); ++i) {
    out << pop.wld_lines[i] << std::endl;
  };
  
  return out;
};

}
