#include <iostream>

#include "paths.hpp"
#include EXPERIMENT
#include RV_GENERATORS

using namespace std;
namespace evolve{

void Experiment::start() {
  pre_snapshot(*this);                         // (function) value of pre_snapshot is set in driver
  t_last_snapshot = t_elapsed;
  g_last_snapshot= pop.generations();
  /// ********************** main loop here  **************************//
  while(not stop_cond(*this)) {
  
    pop.do_event();
    if( pop.num_orgs()== 0) break;             // extinction occurred.  handle this case in driver
    
    t_elapsed += rnd_expo(pop.event_rate() );
    
    if (snapshot_cond(*this)) {
      snapshot(*this);
      t_last_snapshot = t_elapsed;
      g_last_snapshot= pop.generations();
      pop.update_birth_ub();                  // there's probably a better place to put this, but...
    }; 
  };
  /// ******************************************************************//
  post_snapshot(*this);
  t_last_snapshot = t_elapsed;
  g_last_snapshot= pop.generations();
};

void Write_snapshot::operator()(const Experiment& exp){
  const Population& p = exp.population();
  
  //std::cout<<"meanfit = "<< p.tot_ones()/p.num_orgs()<<std::endl;
  o_file<< p.generations()    << "\t" << 
           exp.time_elapsed() << "\t" <<
           p.birth_rate()/p.num_orgs() << "\t";
           
  for (int i = 0; i < Organism::num_states(); ++i)
    o_file << p.num_in_state(i) << "\t";
            
  o_file<< p.num_trk_orgs()   << "\t" <<
           p.num_lineages()   << "\t" <<
         //  p.event_rate()     << "\t" << 
           p.sum_squared_birth_rate()/ p.num_orgs()<< "\t"<< std::endl;
};

void load_pop(Population& pop, std::string filename) {
  std::ifstream file(filename.c_str());  
  boost::archive::text_iarchive arch(file);
  arch >> pop;
  pop.reset_counts();
};

void save_pop(const Population& pop, std::string filename) {
  // Open output file and archive the population
  std::ofstream out_file(filename.c_str());    
  boost::archive::text_oarchive out_archive(out_file);    
  out_archive << as_const(pop);  
};


const  Population& Experiment::population()    const {return pop;                       };
double Experiment::time_last_snapshot()        const {return t_last_snapshot;            };
double Experiment::generations_last_snapshot() const {return g_last_snapshot;            };
double Experiment::time_elapsed()              const {return t_elapsed;                  };
double Experiment::generations_elapsed()       const {return population().generations(); };


Experiment& Experiment::set_population(const Population& p) {
  pop = p;
  return *this;
};


Experiment& Experiment::set_stop_cond(Exp_cond cond){
  stop_cond = cond;
  return *this;
};
  
Experiment& Experiment::set_snapshot_cond(Exp_cond cond)  {
  snapshot_cond = cond;
  return *this;
};

Experiment& Experiment::set_pre_snapshot(Exp_snap snap)  {
  pre_snapshot = snap;
  return *this;
};

Experiment& Experiment::set_snapshot(Exp_snap snap) {
  snapshot = snap;
  return *this;
};

Experiment& Experiment::set_post_snapshot(Exp_snap snap) {
  post_snapshot = snap;
  return *this;
};

Experiment::Experiment() 
  : t_elapsed(0.0),
    t_last_snapshot(-1.0),     // Indicates no snapshot taken yet
    g_last_snapshot(-1.0),
    pre_snapshot(nothing),
    snapshot(nothing),
    post_snapshot(nothing),
    snapshot_cond(never),
    stop_cond(always) {};

// *********   Boolean tests for experimental conditions.  Namespace scope   *********
bool fixed(const Experiment& exp) {
  return (exp.population().num_trk_orgs() == exp.population().num_orgs());
};	  
bool lost(const Experiment& exp) {
  return (exp.population().num_trk_orgs() == 0);
};	  
bool fixed_or_lost(const Experiment& exp) {return (lost(exp) or fixed(exp)); };
bool never(const Experiment& exp)         {return false; };
bool always(const Experiment& exp)        {return true;  };
void nothing(const Experiment& exp)       {};
void test(const Experiment& exp) {
  std::cout << "time_elapsed = " << exp.time_elapsed() 
       << "  num_lineages = " << exp.population().num_lineages() << std::endl;
};





}
