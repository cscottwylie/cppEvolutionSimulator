//  This header defines the Population class, which holds Organism objects.  The number of member
//  Organisms equals the census size.  However, the number of underlying Organism_data objects, i.e.
//  the number of distinct lineages, is generally smaller than this.  Representatives of each lineage
//  are stored trk_lines, wld_lines, depending on their tracking flag status.  
//  
//  The total rates associated with each type of event are stored as data members.  
//
//  Member functions include do_event(), imlementing Gillespie's algorithm for stochastically 
//  choosing which Poisson process occurs.  Also, there are functions for birth, death, mutation,
//  and phenotypic switching.   


#ifndef _POPULATION_
#define _POPULATION_

#include <vector>
#include <iostream>
#include <assert.h>
#include <boost/serialization/vector.hpp>

#include "paths.hpp"
#include ORGANISM
#include RV_GENERATORS
#include TEMP_TEMPLATES

using namespace std;
namespace evolve {

typedef unsigned int uint;
class Population;                                       // defined below
ostream& operator<<(std::ostream&, const Population&);  //namespace scope function defined in .cpp

// ****************************************************************************
// ***********************          Population          ***********************
// ****************************************************************************
class Population {
public:
  //void state_changer(int state);
  Population();                                  // Construct empty population
  
  void set_pop_capacity(int);                   // could be fixed N or logistic carrying capacity
  
  void do_event();                         // Chooses which Poisson process occurs (birth/death,etc)  
  void hack_st_change(int num_to_switch);       // quick fix for adding tracked orgs to burned pop
  void update_birth_ub();                  // explicitly update upper bound of birth rate

  void add_org(Organism& org, int state);

  double event_rate() const;                 // birth + death  + change state
 
  double org_birth_rate(const Organism&, int state) const;
  double birth_rate()             const;      // birth rate depends on genome, so must be calculated
  double sum_squared_birth_rate() const;
  double death_rate()             const;     
  
  int num_orgs()          const;                // Get info on population
  int num_lineages()      const;
  int num_in_state(int)   const;
  int num_births()        const;
  int pop_capacity()      const;                 // carrying capacity of population
  int num_deaths()        const;
  int num_state_chg()     const;
  double generations()    const;

  int num_trk_orgs()      const;
  int num_trk_lineages()  const;
  int num_trk_births()    const;
  int num_trk_deaths()    const;
  int num_trk_state_chg() const;
  int num_lethal_muts()   const;

  int num_wld_orgs()      const;
  int num_wld_lineages()  const;
  int num_wld_births()    const;
  int num_wld_deaths()    const;
  int num_wld_state_chg() const;

  void reset_counts();            // Resets birth/death/state-chg counts.  called if
                                  // population loaded from file

  const Organism& org(int st, int i) const; // if you MUST deal directly with Organism's interface
  const Organism& rnd_org() const;
  const Organism& trk_prog(int) const;
  const Organism& wld_prog(int) const;

  friend std::ostream& operator<<(std::ostream& out, const Population& pop);
private:  
  std::vector<double> b_rate_tots;          // Birth-rate totals by state.  other rate totals
                                            // easily calculated, thus not stored
  std::vector<double> sum_sq_b_rates;       // sum of squared birth rates for computing variance[br]
  std::vector<double> b_rate_ubnds;         // Birth-rate upper bounds by state
  std::vector<double> tot_rates;            // Each states tot event rate 
  std::vector<std::vector<Organism> > orgs; // Orgs in pop. organized by state
  std::vector<Organism> trk_lines;          // Tracked lineages progenitors
  std::vector<Organism> wld_lines;          // Tracked lineages progenitors
  
  double tot_event_rate;     // Total rate an internally handled event happens 

  int n_orgs;                // Counters
  int n_births;              // number of births mod n_orgs... gets too high otherwise
  double gens;
  int pop_cap;
  int n_deaths;
  int n_state_chg;
  int leth_muts;

  int n_trk_orgs;         
  int n_trk_births;
  int n_trk_deaths;
  int n_trk_state_chg;
  
  void death(int state);
  void birth(int state);     // Basic functions by state
  void state_changer(int state);

  void add_rates(const Organism&, int state);         // Helper functions      
  void remove_rates(const Organism&, int state);
  void add_to_lineage_data(Organism&, int state);
  void remove_from_lineage_data(Organism&, int state);

  // Enable reading/writing of object to archive file
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
  ar & b_rate_tots;  
  ar & sum_sq_b_rates;
  ar & b_rate_ubnds; 
  ar & tot_rates;
  ar & gens;
  ar & orgs;  
  ar & trk_lines;       
  ar & wld_lines;       
  ar & tot_event_rate;
  ar & n_orgs;
  ar & n_births;
  ar & pop_cap;
  ar & n_deaths;
  ar & n_state_chg;
  ar & n_trk_orgs;
  ar & n_trk_births;
  ar & n_trk_deaths;
  ar & n_trk_state_chg;
  };
};

inline double Population::event_rate()    const {return tot_event_rate;  };
inline int    Population::num_orgs()      const {return n_orgs;          };
inline int    Population::num_births()    const {return n_births;        };
inline double Population::generations()   const {return gens;            };
inline int    Population::pop_capacity()  const {return pop_cap;         };
inline int    Population::num_deaths()    const {return n_deaths;        };
inline int    Population::num_state_chg() const {return n_state_chg;     };
inline int    Population::num_trk_orgs()  const {return n_trk_orgs;      };

inline int Population::num_in_state(int st) const {
  assert(st >= 0);
  assert(st < Organism::num_states());
  return orgs[st].size();
};

inline int Population::num_trk_lineages() const {return trk_lines.size(); };

inline int Population::num_lineages()     const {
  return wld_lines.size() + trk_lines.size();
};

inline int Population::num_wld_lineages()  const {return wld_lines.size();       };
inline int Population::num_trk_births()    const {return n_trk_births;           };
inline int Population::num_trk_deaths()    const {return n_trk_deaths;           };
inline int Population::num_trk_state_chg() const {return n_trk_state_chg;        };
inline int Population::num_lethal_muts()   const {return leth_muts;              };
inline int Population::num_wld_orgs()      const {return n_orgs - n_trk_orgs;    };
inline int Population::num_wld_births()    const {return n_births - n_trk_births;};
inline int Population::num_wld_deaths()    const {return n_deaths - n_trk_deaths;};



inline int Population::num_wld_state_chg() const {
  return n_state_chg - n_trk_state_chg;
};

inline const Organism& Population::wld_prog(int index) const {
  assert(index < (int) wld_lines.size());
  return wld_lines[index];
};

inline const Organism& Population::trk_prog(int index) const {
  assert(index < (int) trk_lines.size());
  return trk_lines[index];
};

inline void Population::set_pop_capacity(int p_cap){
  pop_cap = p_cap;
};


inline double Population::org_birth_rate( const Organism& org, int st) const {
  double fit= Organism::state( st).birth_prefactor();
  if( org.allele_state() == 1) 
    fit*= (1+ Organism::state( st).sel_coeff_ben() );
  else if( org.allele_state() == -1)
    fit*= (1- Organism::state( st).sel_coeff_del() );
  return fit;
};
  

inline void Population::add_to_lineage_data(Organism& org, int st) {
  assert(st >= 0);
  assert(st < Organism::num_states());

  org.inc_num_in_state(st);  
  org.inc_num_in_lineage();
  if (org.lineage_index() == -1) {
    if (org.tracked()) {
      trk_lines.push_back(org);
      org.set_lineage_index(trk_lines.size() - 1);
    }
    else {
      wld_lines.push_back(org);
      org.set_lineage_index(wld_lines.size() - 1);
    };
  };
};

inline void Population::remove_from_lineage_data(Organism& org, int st) {
  assert(st >= 0);
  assert(st < Organism::num_states());
  assert(org.lineage_index() != -1);               // lineage index assigned positive # when added
  assert(org.num_in_state(st) > 0);                // can't remove an org that isn't there
  assert(org.num_in_lineage() > 0);

  org.dec_num_in_state(st);
  org.dec_num_in_lineage();                        // merely reduces a counter
  if (org.num_in_lineage() == 0) {
    if (org.tracked()) {
      assert(org.lineage_index() >= 0);
      assert(org.lineage_index() < (int) trk_lines.size());
      trk_lines.back().set_lineage_index(org.lineage_index());
      swap_pop(trk_lines, org.lineage_index());
    }
    else {
      assert(org.lineage_index() >= 0);
      if( org.lineage_index() >= (int) wld_lines.size() ){ 
        cout<<"index= "<<org.lineage_index()<<", wld_lines.size= "<<wld_lines.size()<<endl; };
        assert(org.lineage_index() < (int) wld_lines.size());
      wld_lines.back().set_lineage_index(org.lineage_index());
      swap_pop(wld_lines, org.lineage_index());
    };
  };
};

inline const Organism& Population::org(int st, int i) const {
  assert(st >= 0);
  assert(st <= Organism::num_states());
  assert(i < (int)orgs[st].size());
  
  return orgs[st][i];
};

inline const Organism& Population::rnd_org() const {
  assert(num_orgs() > 0);
  assert(orgs.size() > 0);
  int st  = 0;
  uint ch = rnd_int(n_orgs);
  while (ch >= orgs[st].size()) {
    ch -= orgs[st].size();
    ++st;
  };
  return orgs[st][ch];
};

inline void Population::birth(int st) {
  assert(b_rate_tots[st] > 0);  
  assert(b_rate_ubnds[st] > 0);
  assert(tot_rates[st] > 0);
  
  const double ubound = b_rate_ubnds[st];
  int ch;
  do ch = rnd_int(orgs[st].size());
  while ((rnd_uniform() * ubound) > org_birth_rate(orgs[st][ch], st));
  //  std::cout << "organism number " << ch << std::endl;
  
  Organism parent= orgs[st][ch];
  Organism child=  parent;                 // points to same data as parent
  
  bool is_lethal=  child.mutate(st);       // Organism::mutate(int) made new pointer if necessary
  if (is_lethal){
    ++leth_muts;
    return;               // If lethal mutation occurred, don't add child (was killed)   
  };
  
  if (child != parent) child.reset_lineage_counts();        // i.e. was make_write_safe() called
  
  add_org(child, st);
  ++n_births;
  gens += (double)1/n_orgs;
  
  if (parent.tracked()) ++n_trk_births;
};

/*inline void Population::birth( int st) {
  assert(b_rate_tots[ st]> 0);  
  assert(b_rate_ubnds[ st]> 0);
  assert(tot_rates[ st]> 0);
  
  const double ubound= b_rate_ubnds[ st];
  int ch;
  
  do ch = rnd_int( orgs[ st].size() );
  while ( ( rnd_uniform()* ubound)> org_birth_rate( orgs[ st][ ch], st) );
  
  Organism original= orgs[ st][ ch];            // all 3 have identical pointers
  Organism parent  = original;
  Organism child   = original;  
  if ( parent.tracked() ) ++n_trk_births;
  
  // *********    ALWAYS call mutate function, which determines # muts     ********* ///
  bool is_lethal_child=  child.mutate(st);      // Organism::mutate() made new pointer if necessary
  bool is_lethal_parent= parent.mutate(st);     // Organism::mutate() made new pointer if necessary
  
  if( not is_lethal_child) {
    if( child!= original) child.reset_lineage_counts();  // first reset counts
    add_org( child, st);                               // THEN add to orgs                      
  }
  else ++leth_muts;
  
  if( parent!= original) {                 // otherwise, if parent didn't change, leave parent alone  
    parent.reset_lineage_counts();
    if( not is_lethal_parent)  add_org( parent, st);
    else ++leth_muts;                                   // increment counter, but DON'T add parent 

    /////////////// essentially kill original in this chunk   ///////////////////////////
    remove_rates( orgs[st][ch], st);              // next 3 lines effectively kill the DNA recipient
    remove_from_lineage_data( orgs[st][ch], st);
    swap_pop( orgs[st], ch);
    --n_orgs;             //n_orgs was incremented in add_org(), so undo that here and next 2 lines
    if ( orgs[st][ch].tracked() ) --n_trk_orgs;
    /////////////////////////////////////////////////////////////////////////////////////  
  };
  
  
  ++n_births;
  gens+= (double)1/n_orgs;
  
};*/





} //end of evolve namespace


#endif

