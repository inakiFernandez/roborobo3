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
typedef std::vector<odNeatSpecies*> local_population;

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
    //Genome* selectBestGenome();

    //TabuList: set of dropped genomes, with corresponding timeout
    std::map<Genome*, int> _tabu;
    void add_to_tabu_list(Genome* g);
    bool tabu_list_approves(Genome* g);
    int tabu_contains(Genome* g);
    bool update_tabu_list(Genome* g);

    //Population
    bool findInPopulation(GC gId);
    void add_to_population(Genome* g, double f);
    bool population_accepts(double f);
    odNeatSpecies* computeSpeciesOfGenome(Genome* g);

    void add_to_species(message msg);
    double speciesFitness(int sp);
    void adjustSpeciesFitness();
    odNeatSpecies* selectSpecies();
    genomeFitness selectParent(odNeatSpecies* sp);
    void recompute_all_species();
    /*void cleanPopAndSpecies();*/

    Genome* generate_offspring();
    bool doBroadcast();
    //EndOdNeat--------------------------------------------------------------------

    void createNN();
    void stepBehaviour(); void stepEvolution();
    
    void broadcastGenome(); void storeGenome(message m);

    int _newSpId;

    int _lifetime;    
    // ANN
    double _minValue; double _maxValue;
    unsigned int _nbInputs; unsigned int _nbOutputs;
    //gene (link) and neuron local counters
    int _g_count; int _n_count;
    //Number of gathered items
    int _items;
        
    void resetRobot();
    
public:

    odNeatController(RobotWorldModel *wm);
    ~odNeatController();
    Genome *_genome; // current genome in evaluation
    GC _genomeId; local_population _pop;
    void reset();
    void step();

    int getBirthdate() { return _birthdate; }
    double getFitness(){ return _currentFitness;}
    void logGenome(std::string s);
    void pickItem();
    int getItems();
};


#endif

