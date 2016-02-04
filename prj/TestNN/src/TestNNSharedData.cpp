/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "TestNN/include/TestNNSharedData.h"

int TestNNSharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

// global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver
bool TestNNSharedData::gPropertiesLoaded = false;

int TestNNSharedData::gNbHiddenLayers = 1;
int TestNNSharedData::gNbNeuronsPerHiddenLayer = 5;
int TestNNSharedData::gNeuronWeightRange = 800;
bool TestNNSharedData::gWithBias = false;

int TestNNSharedData::gFitness = -1;

std::string TestNNSharedData::gOutGenomeFile = "";
