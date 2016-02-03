/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#ifndef INCREMENTSHAREDDATA_H
#define INCREMENTSHAREDDATA_H

class IncrementSharedData {
	
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

    static int gNbHiddenLayers; // default: 1
    static int gNbNeuronsPerHiddenLayer; // default: 5
    static int gNeuronWeightRange; // default: 800.0 (ie. weights are in [-400,+400[
    static bool gWithBias;

    static int gFitness; //indicates fitness to use
};


#endif
