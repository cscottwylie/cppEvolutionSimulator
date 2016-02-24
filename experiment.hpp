// The Experiment class contains a population that evolves, once Experiment.start() is called.
//
// Other aspects of the experiment include conditions for terminating evolution and how/when to 
// record data.  These conditions are represented as boost::functions of type Exp_snap, Exp_cond, 
// and data members in Experiment.  (double) t_elapsed is also a data member.
//
// The condition functions, e.g. Time_since_start(double) are actually classes, whose data members
// can be interpreted as the function's argument.  


#ifndef _EXPERIMENT_
#define _EXPERIMENT_

#include <iostream>
#include <boost/function.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include "paths.hpp"
#include POPULATION

using namespace std;

namespace evolve {

class Experiment;                                             // class defined below

// namespace scope functions for (de)archiving populaitons with boost:: library
void load_pop(Population&, std::string);
void save_pop(const Population&, std::string);

// function place holders: will be associated with, e.g. fixed_or_lost
typedef boost::function<void (const Experiment&)> Exp_snap;   // how to record data
typedef boost::function<bool (const Experiment&)> Exp_cond;   // when to record data, start, quit

// namespace scope functions.  Will be assigned to Exp_cond
void test         (const Experiment&);
bool fixed        (const Experiment&);
bool lost         (const Experiment&);
bool fixed_or_lost(const Experiment&); 
bool never        (const Experiment&);       
bool always       (const Experiment&);      
void nothing      (const Experiment&);    

// Will be assigned to Exp_cond, just as functions above, e.g. fixed_or_lost.  The only
// purpose of these "classes" is to define a function that holds a value.
class Mean_fit_at_least;
class Time_since_start;
class Time_since_last_snapshot;

class Write_snapshot; 					// Important: this specifies which data is written down

// ****************************************************************************
// ***********************          Experiment          ***********************
// ****************************************************************************

class Experiment {
public:
  Experiment();
  void start();

  Experiment& set_population( const Population&);
  Experiment& set_stop_cond    ( Exp_cond);
  Experiment& set_snapshot_cond( Exp_cond);
  Experiment& set_pre_snapshot ( Exp_snap);
  Experiment& set_snapshot     ( Exp_snap);
  Experiment& set_post_snapshot( Exp_snap);
  
  const  Population& population()    const;
  double time_elapsed()              const; 
  double generations_elapsed()       const;
  double time_last_snapshot()        const; 
  double generations_last_snapshot() const;
private:
  Population  pop;
  double t_elapsed;
  double t_last_snapshot;
  double g_last_snapshot;
  
  Exp_snap pre_snapshot;
  Exp_snap snapshot;
  Exp_snap post_snapshot;
  Exp_cond snapshot_cond;
  Exp_cond stop_cond;
};

// really this is a function that holds a file, more than a "class"
class Write_snapshot {
public:
  Write_snapshot(std::ofstream& out_file) : o_file(out_file) {};
  void operator()(const Experiment&);
private:
  std::ofstream& o_file;
};


class Mean_fit_at_least {
public:
  explicit Mean_fit_at_least(double fit) : mean_fit(fit) {assert(fit >= 0.0);};
  
  bool operator()(const Experiment& exp) const {
    return ((exp.population().birth_rate() / exp.population().num_orgs())
	    > mean_fit);
  };
private:
  double mean_fit;
};


class Time_since_start {
public:
  explicit Time_since_start(double t) : time(t) {assert(time >= 0.0); };
  bool operator()(const Experiment& exp) const {
    return exp.time_elapsed() > time; 
  };
private:
  double time;
};

class Generations_since_start {
public:
  explicit Generations_since_start( double t) : gens( t) {assert( gens>= 0.0);};
  bool operator()( const Experiment& exp) const{
    return exp.population().generations()> gens;
  };
private:
  double gens;
};
  
class Time_since_last_snapshot {
public:
  explicit Time_since_last_snapshot(double time_interval) 
    : t_interval(time_interval) { assert(time_interval > 0.0); };
  
  bool operator()(const Experiment& exp) {
    return exp.time_elapsed() > (exp.time_last_snapshot() + t_interval);
  };
private:
  double t_interval;
};

class Generations_since_last_snapshot {
public:
  explicit Generations_since_last_snapshot( double gen_interval)
    : g_interval( gen_interval) { assert( gen_interval> 0.0 ); };
    
  bool operator()( const Experiment& exp){
    return exp.population().generations()> (exp.generations_last_snapshot()+ g_interval);
  };
private:
  double g_interval;
};


}  // end namespace block

#endif
