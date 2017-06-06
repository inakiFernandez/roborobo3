/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "OriginalEA2017/include/OriginalEA2017SharedData.h"

double OriginalEA2017SharedData::gSigmaRef = 0.0; // reference value of sigma
int OriginalEA2017SharedData::gPopulationSize = 1; // size of population per robot
bool OriginalEA2017SharedData::gClearPopulation; // empty population after replacement?
bool OriginalEA2017SharedData::gStoreOwn; // store own genome in local population?
int OriginalEA2017SharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

int OriginalEA2017SharedData::gSelectionMethod = 0; // Random by default

// global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver
bool OriginalEA2017SharedData::gPropertiesLoaded = false;

int OriginalEA2017SharedData::gControllerType = 0; //MLP, Perceptron, Elman
int OriginalEA2017SharedData::gNbHiddenLayers = 1;
int OriginalEA2017SharedData::gNbNeuronsPerHiddenLayer = 5;
int OriginalEA2017SharedData::gNeuronWeightRange = 800;
bool OriginalEA2017SharedData::gWithBias = false;

bool OriginalEA2017SharedData::gCommunicationBySensors = true;
int OriginalEA2017SharedData::gFitness = -1;
std::vector<int> OriginalEA2017SharedData::gTaskSeq;
std::vector<int> OriginalEA2017SharedData::gTimeSeq;
int OriginalEA2017SharedData::gTaskIdx = 0;

std::string OriginalEA2017SharedData::gOutGenomeFile = "";
bool OriginalEA2017SharedData::gSaveGenome = true;
bool OriginalEA2017SharedData::gIsLoadGenome = true;
std::string OriginalEA2017SharedData::gEvolutionLogFile = "";

bool OriginalEA2017SharedData::gWithCollectColorEffector = false;

int OriginalEA2017SharedData::gBrait = -1;

int OriginalEA2017SharedData::gForgetMethod = -1;


int OriginalEA2017SharedData::gUpdateGC = true;

double OriginalEA2017SharedData::gSelPressure= 1.0;
