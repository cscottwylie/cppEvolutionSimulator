//  version for Ben Allen's simple 3-level fitness landscape.

//  This header contains definitions for classes Org_state, Organism_data, and Organism, as well as
//  definitions of member inline functions.  Other function definitions are in organism.cpp
//
//  Org_state is essentially a set of parameters governing a phenotype. 
// 
//  Organism is the copy-on-write facade for Organism_data.  These objects represent *lineages* 
//  of possibly several biological cells identical by descent.  The objects contain a genome, a tracking 
//  flag, and the phenotypic distribution of the clones they represent.  They also have member 
//  functions for manipulating the genome (mutate etc.).

// The genome is an integer from {+1, 0, -1}


#ifndef _ORGANISM_
#define _ORGANISM_

#include <vector>
#include <iostream>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "paths.hpp"
#include RV_GENERATORS
#include PARAMETERS
#include TEMP_TEMPLATES

namespace evolve{

class Parameters;  // defined elsewhere, needed by Org_state::set_state_params(...)

//these classes defined here
class Org_state;      
class Organism_data;
class Organism;

//namespace scope function prototypes, defined in .cpp
bool          operator!= (const Organism&, const Organism& );
bool          operator== (const Organism&, const Organism& );
std::ostream& operator<< (std::ostream&,   const Organism& );
std::ostream& operator<< (std::ostream&,   const Org_state&);


// ****************************************************************************
// *********************          Organism State          *********************
// ****************************************************************************
class Org_state {
public:
  Org_state();                        // Constructor sets properties to 0.0
  
  double mut_rate_del()    const;       // Get properties 
  double mut_rate_ben()    const;
  double sel_coeff_ben()   const;
  double sel_coeff_del()   const;
  double birth_prefactor() const;
  double chg_rate()        const; 
  double death_rate()      const;
  
  // Functions for setting the state's properties, all return a reference
  // to the current state to facilitate chaining of fuction calls.
  Org_state& set_mut_rate_del    (double); 
  Org_state& set_mut_rate_ben    (double);  
  Org_state& set_sel_coeff_ben   (double);
  Org_state& set_sel_coeff_del   (double);  
  Org_state& set_birth_prefactor (double); 
  Org_state& set_chg_rate        (double);   
  Org_state& set_death_rate      (double);
private:
  double mt_rate_b;   // probability of beneficial mutation per replication
  double mt_rate_d;   // probability of deleterious mutation per replication
  double s_b;         // selection coefficient of beneficial mutations
  double s_d;         // selection coefficient of deleterious mutations
  double b_pre;       // fitness of "neutral" genotype  
  double c_rate;      // Rate organsim can switch it's state
  double f_adjust;    // Additive adjustement to birth rate
  double d_rate;      // Death rate per unit time

  // Enable reading/writing of object to archive file
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & mt_rate_b;
    ar & mt_rate_d;
    ar & s_b;
    ar & s_d;
    ar & b_pre;
    ar & c_rate;
    ar & f_adjust;
    ar & d_rate;
  };
};  

// ******************* Organism-property reading functions  *******************
inline double Org_state::mut_rate_ben()    const {return mt_rate_b;};
inline double Org_state::mut_rate_del()    const {return mt_rate_d;};
inline double Org_state::sel_coeff_ben()   const {return s_b;      };
inline double Org_state::sel_coeff_del()   const {return s_d;      };
inline double Org_state::birth_prefactor() const {return b_pre;    };
inline double Org_state::chg_rate()        const {return c_rate;   };
inline double Org_state::death_rate()      const {return d_rate;   };


// ****************************************************************************
// *********************          Organism Data           *********************
// ****************************************************************************  
class Organism_data {
public:
  friend class Organism;          // Organism needs access to its data
private:
  Organism_data();                // Starts with empty genome, not tracked

  bool is_tracked;                // Marks org. as tracked. No physical effect
  int  n_in_lineage;              // Num. orgs with identical data via descent
  int  line_index;                // Position of lineage's progenitor in list 
  std::vector<int> n_in_state;    // Num. orgs in lineage with particular state
  int allele;                     // +1, 0, or -1: entire "genome"
  
  // -----------------------   Data reading functions   -----------------------
  int     allele_state() const;    
  int     num_in_lineage()  const;     // Number of orgs identical by descent
  int     num_in_state(int) const;     // Number of orgs in lineage in each state
  int     lineage_index()   const;     
  bool    tracked()         const; 
  
  // ------------------   Genome modifying functions   ------------------

  void set_tracked(bool); 
  void set_allele_state(int);

  // ----------   Helper functions for changing lineage info   ----------
  void set_lineage_index(int);
  void inc_num_in_lineage();
  void dec_num_in_lineage();
  void inc_num_in_state(int);   // Increment num_in_state counter
  void dec_num_in_state(int);   // Decrement num_in_state counter
  
  // -------   Enable reading/writing of object to archive file   -------
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & is_tracked;
    ar & n_in_lineage;
    ar & line_index;
    ar & n_in_state;
    ar & allele;
  };
};


inline bool Organism_data::tracked()        const {return is_tracked;     };
inline int  Organism_data::allele_state()   const {return allele;         };
inline int  Organism_data::num_in_lineage() const {return n_in_lineage;   };
inline int  Organism_data::lineage_index()  const {return line_index;     };

inline void Organism_data::set_allele_state(int g)     {allele = g;       };
inline void Organism_data::set_tracked(bool trk)       {is_tracked = trk; }; 
inline void Organism_data::set_lineage_index(int idx)  {line_index = idx; }; 
inline void Organism_data::inc_num_in_lineage()        {++n_in_lineage;   };
inline void Organism_data::inc_num_in_state(int st)    {++n_in_state[st]; };
inline void Organism_data::dec_num_in_lineage()        {--n_in_lineage;   };
inline void Organism_data::dec_num_in_state(int st)    {--n_in_state[st]; };


// ****************************************************************************
// ***************                 Organism                 *******************
// ****************          (copy-on-write facade)         *******************
// ****************************************************************************

class Organism {
public:
  Organism();                             // Org starts with "neutral" genotype (0)

  bool mutate(int state);                  // returns 0 if "lethal" mutation occurred
  
  // Data reading functions
  int  num_in_lineage()  const;
  int  num_in_state(int) const;
  int  lineage_index()   const;
  bool tracked()         const; 
  int allele_state()     const;
 
  // Data setting functions

  Organism& set_tracked(bool);

  void inc_num_in_lineage();
  void dec_num_in_lineage();
  void set_lineage_index(int);
  void inc_num_in_state(int);     // Increment num_in_state counter
  void dec_num_in_state(int);     // Decrement num_in_state counter
  void reset_lineage_counts();
                                   
  // Change/access possible organism states, all orgs share set of states.
  static void add_states(int);          // Adds "all 0.0" states
  static void set_state_params(int state, const Parameters&);       
  static Org_state& state(int i);       // Allows access to i-th state
  static void write_states(); 
  static int num_states(); 

  
  // Gets underlying raw pointer to data. DON'T DO ANYTHING WITH THIS!
  const Organism_data* get_ptr() const {return boost::get_pointer(data_ptr);};

  // Convenient to give ouput operator access
  friend std::ostream& operator<<(std::ostream&, const Organism&);
  
  // Enable reading/writing of object to archive file
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    ar & data_ptr;
    ar & state_list;
  }
private:
  boost::shared_ptr<Organism_data> data_ptr; // Pointer to actual organism data
  static std::vector<Org_state> state_list;  // Contains org's possible states
};

// *********************** Organism reading functions  ***********************
inline int    Organism::allele_state() const {return data_ptr->allele;       };
inline bool   Organism::tracked()      const {return data_ptr->tracked();    };
inline int    Organism::num_states()         {return state_list.size();      };


inline int  Organism::lineage_index() const {
  return data_ptr->lineage_index();
};

inline int  Organism::num_in_state(int st) const {
  return data_ptr->num_in_state(st);
};

inline int  Organism::num_in_lineage() const {
  return data_ptr->num_in_lineage();
};

inline int  Organism_data::num_in_state(int st) const {
  assert(st >= 0);
  assert(st < (int) n_in_state.size());
  assert(st < Organism::num_states());
  return n_in_state[st];
};

inline void Organism::reset_lineage_counts() {
  // reserve memory for new org_data w/ 0 in lineage, initialize allele to 0, flag line_index w/ -1
  make_write_safe(data_ptr);
  data_ptr->n_in_lineage = 0;
  fill_n(data_ptr->n_in_state.begin(), Organism::num_states(), 0);
  data_ptr->line_index = -1;
};

inline void Organism::inc_num_in_lineage() {data_ptr->inc_num_in_lineage(); };
inline void Organism::dec_num_in_lineage() {data_ptr->dec_num_in_lineage(); };

inline void Organism::inc_num_in_state(int st) {
  assert(state_list.size() > 0);
  assert(st >= 0);
  assert(st < (int) state_list.size()); 
  data_ptr->inc_num_in_state(st);
};

inline void Organism::dec_num_in_state(int st) {
  assert(state_list.size() > 0);
  assert(st >= 0);
  assert(st < (int) state_list.size());
  assert(data_ptr->num_in_state(st) > 0);
  data_ptr->dec_num_in_state(st);
};

inline void Organism::set_lineage_index(int i) {
  assert(i >= 0);
  data_ptr->set_lineage_index(i);
};


} // closing namespace block


#endif

