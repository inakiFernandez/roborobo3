/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "Original/include/OriginalSharedData.h"

double OriginalSharedData::gSigmaRef = 0.0; // reference value of sigma
int OriginalSharedData::gPopulationSize = 1; // size of population per robot
bool OriginalSharedData::gClearPopulation; // empty population after replacement?
bool OriginalSharedData::gStoreOwn; // store own genome in local population?
int OriginalSharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

int OriginalSharedData::gSelectionMethod = 0; // Random by default

// global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver
bool OriginalSharedData::gPropertiesLoaded = false;

int OriginalSharedData::gControllerType = 0; //MLP, Perceptron, Elman
int OriginalSharedData::gNbHiddenLayers = 1;
int OriginalSharedData::gNbNeuronsPerHiddenLayer = 5;
int OriginalSharedData::gNeuronWeightRange = 800;
bool OriginalSharedData::gWithBias = false;

bool OriginalSharedData::gCommunicationBySensors = true;
int OriginalSharedData::gFitness = -1;
std::vector<int> OriginalSharedData::gTaskSeq;
std::vector<int> OriginalSharedData::gTimeSeq;
int OriginalSharedData::gTaskIdx = 0;

std::string OriginalSharedData::gOutGenomeFile = "";
bool OriginalSharedData::gSaveGenome = true;
bool OriginalSharedData::gIsLoadGenome = true;
std::string OriginalSharedData::gEvolutionLogFile = "";

bool OriginalSharedData::gWithCollectColorEffector = false;
