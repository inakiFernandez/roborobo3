/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#include <string>
#include <vector>
#ifndef ORIGINALSHAREDDATA_H
#define ORIGINALSHAREDDATA_H

class OriginalSharedData {
	
	public: 
	
    static double gSigmaRef; //! reference value of sigma
    static int gPopulationSize; //!size of population per robot
    static bool gClearPopulation; //! empty population after replacement?
    static bool gStoreOwn; //!store own genome in local population?
    static int gEvaluationTime; //! theoretical duration of a generation (ie. maximum time a controller will be evaluated on a robot)
	static int gIteration; //! used by every class to know what is the current iteration step of roborobo
    static int gSelectionMethod; //! Random by default

	static bool gPropertiesLoaded;

    static int gControllerType; // MLP, Perceptron, Elman
    static int gNbHiddenLayers; // default: 1
    static int gNbNeuronsPerHiddenLayer; // default: 5
    static int gNeuronWeightRange; // default: 800.0 (ie. weights are in [-400,+400[
    static bool gWithBias;

    static bool gCommunicationBySensors; // comm. with sensors or by distance (slower)
    static int gFitness; //indicates which fitness to use
    static std::vector<int> gTaskSeq;

    static int gTaskIdx;
    static std::vector<int> gTimeSeq;

    static std::string gOutGenomeFile; //filename for last genome (text format)
    static std::string gEvolutionLogFile; // base filename for evolution data (fitness, etc...)
    static bool gSaveGenome;
    static bool gIsLoadGenome;

    static bool gWithCollectColorEffector;
    static int gBrait;
    static int gForgetMethod;
};


#endif
