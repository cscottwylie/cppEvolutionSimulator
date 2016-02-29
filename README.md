# This repo contains c++ code for simulating biological evolution.  
### The simulation contains several classes, compartmentalized into corresponding files.  The .hpp files contain declarations, and are a good place to start browsing the functionality offered.  The .cpp files contain the implementation of classes and their member functions. 

In order from lowest to highest level, the files/classes are:
- **Organism.hpp**.  Each Organism object contains a pointer to an Organism_data object, which contains a representation of the genome.  This scheme saves memory, since several Organisms within a population may be genetically identical.  New data/pointers are only created upon genetic changes; thus Organism is a copy-on-write facade for Organism_data
- **Population.hpp**.  A Population contains a std::vector of Organism, along with member functions for modifying the Population in Kosher ways.   The main dynamical function is Population::do_event(), which selects among birth, death, etc. according to Gillespie's algorithm (which exactly simulates multi-type Poisson processes).  
- **Experiment.hpp** allows user to set up common evolutionary scenarios, such as competition experiments (which terminate when 1 of 2 competitors go extinct), or running for a fixed number of generations.
- **rv_generators.hpp and temp_templates.hpp** contain  random number generators and miscallaneous helper functions.
- **parameters*.txt** contain the parameters needed to run various experiments

#### This code was written by Aaron Trout and Scott Wylie, and was used in the following publications:
- "Optimal Strategy for Competence Differentiation in Bacteria", PLoS Genetics, 2010, C. Scott Wylie, Aaron Trout, et al.
- "A Biophysical Protein  Folding Model Accounts for Most Mutational Fitness Effects in Viruses", PNAS, 2011, C. Scott Wylie et al.
- "Mutation-Induced Extinction in Finite Populations: Lethal Mutagenesis and Lethal Isolation", PLoS Computational Biology, 2012, C. Scott Wylie et al.
