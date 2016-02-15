#ifndef ODNEATCONTROLLER_H
#define ODNEATCONTROLLER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Utilities/Graphics.h"
#include "Controllers/Controller.h"
#include "WorldModels/RobotWorldModel.h"
#include "odNeat/include/odNeatAgentObserver.h"
#include "odNeat/include/odNeatSpecies.h"
#include <odneatgc/network.h>
#include <odneatgc/genome.h>
#include <set>
#include <tuple>
#include <map>
#include <iomanip>

using namespace ODNEATGC;


//Population (map int to species) including genome and species
//Species are local, received genomes have to be assigned a new species
typedef std::vector<odNeatSpecies> local_population;

class odNeatController : public Controller
{
private:
    int _iteration; int _birthdate; // evaluation when this controller was initialized
    double _currentFitness; double _energy;
    //OdNeat-------------------------------------------
    Network *nn;

    //Fitness is measured every gFitnessFreq steps
    int _fitnessUpdateCounter;
    double update_energy_level();
    void updateFitness();
    bool in_maturation_period();
    Genome* selectBestGenome();

    //TabuList: set of dropped genomes, with corresponding timeout
    std::map<Genome*, int> _tabu;
    void add_to_tabu_list(Genome* g);
    bool tabu_list_approves(Genome* g);
    int tabu_contains(Genome* g);

    //Population
    bool findInPopulation(Genome* g);
    void add_to_population(message msg);
    bool population_accepts(double f);
    int computeSpeciesOfGenome(Genome* g);

    void add_to_species(message msg);
    double speciesFitness(int sp);
    void adjustSpeciesFitness();
    /*
    int findInSpecies(Genome* g);
    void cleanPopAndSpecies();
    int computeSpeciesId(Genome* g);
    void recompute_all_species();
    void adjust_population_size();
    void adjust_species_fitness();
    void adjust_active_species_fitness(int species);
    int selectSpecies();*/

    /*Genome* selectParent(int spId);
    void update_population(Genome* offspring);*/
    Genome* generate_offspring();
    bool doBroadcast();
    //EndOdNeat--------------------------------------------------------------------

    void createNN();
    void stepBehaviour(); void stepEvolution();
    
    void broadcastGenome(); void storeGenome(message m);

    GC _genomeId; local_population _pop; int _newSpId;

    int _lifetime;    
    // ANN
    double _minValue; double _maxValue;
    unsigned int _nbInputs; unsigned int _nbOutputs;
    //gene (link) and neuron local counters
    int _g_count; int _n_count;
        
    void resetRobot();
    
public:

    odNeatController(RobotWorldModel *wm);
    ~odNeatController();
    Genome *_genome; // current genome in evaluation

    void reset();
    void step();

    int getBirthdate() { return _birthdate; }
    double getFitness(){ return _currentFitness;}
    void logGenome(std::string s);
    void pickItem();
};


#endif

