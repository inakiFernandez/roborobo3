/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#include <string>
#ifndef TESTNNSHAREDDATA_H
#define TESTNNSHAREDDATA_H

class TestNNSharedData {
	
	public: 
	
	// -----
    static int gEvaluationTime; //! theoretical duration of a generation (ie. maximum time a controller will be evaluated on a robot)
	static int gIteration; //! used by every class to know what is the current iteration step of roborobo

	static bool gPropertiesLoaded;

    static int gNbHiddenLayers; // default: 1
    static int gNbNeuronsPerHiddenLayer; // default: 5
    static int gNeuronWeightRange; // default: 800.0 (ie. weights are in [-400,+400[
    static bool gWithBias;

    static int gFitness; //indicates fitness to use

    static std::string gOutGenomeFile; //filename for last genome (binary format)
};


#endif
