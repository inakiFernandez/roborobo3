/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#include <string>
#ifndef ODNEATSHAREDDATA_H
#define ODNEATSHAREDDATA_H
#include "odneatgc/genome.h"

using namespace ODNEATGC;
//Genome, energy, node counter, link gene counter
typedef std::tuple<Genome*, double, int,int> message;

class odNeatSharedData {
	
	public: 
	
	// -----
	
    static double gSigmaRef; //! reference value of sigma
    static int gPopulationSize; //!size of population per robot
    static bool gClearPopulation; //! empty population after replacement?
    static bool gStoreOwn; //!store own genome in local population?
    static int gEvaluationTime; //! theoretical duration of a generation (ie. maximum time a controller will be evaluated on a robot)
	static int gIteration; //! used by every class to know what is the current iteration step of roborobo
    static int gSelectionMethod; //! Random by default

	static bool gPropertiesLoaded;

    static bool gCommunicationBySensors; // comm. with sensors or by distance (slower)
    static int gFitness; //indicates which fitness to use

    static std::string gOutGenomeFile; //filename for last genome (text format)

    static bool gSaveGenome;

    //OdNeat parameters
    static double gDefaultInitialEnergy;
    static double gEnergyThreshold;
    static double gMaxEnergy;

    static int gMaturationPeriod;

    static double gCompatThreshold;

    static int gFitnessFreq;
    static int gTabuTimeout;
    static double gTabuThreshold;

    static double gEnergyItemValue;
    // -----
    static double gEnergyConsumption;

    static bool gUpdateGC;
};


#endif
