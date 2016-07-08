/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "Collect2/include/Collect2SharedData.h"

double Collect2SharedData::gSigmaRef = 0.0; // reference value of sigma
int Collect2SharedData::gPopulationSize = 1; // size of population per robot
bool Collect2SharedData::gClearPopulation; // empty population after replacement?
bool Collect2SharedData::gStoreOwn; // store own genome in local population?
int Collect2SharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

int Collect2SharedData::gSelectionMethod = 0; // Random by default

// global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver
bool Collect2SharedData::gPropertiesLoaded = false;

int Collect2SharedData::gControllerType = 0; //MLP, Perceptron, Elman
int Collect2SharedData::gNbHiddenLayers = 1;
int Collect2SharedData::gNbNeuronsPerHiddenLayer = 5;
int Collect2SharedData::gNeuronWeightRange = 800;
bool Collect2SharedData::gWithBias = false;

bool Collect2SharedData::gCommunicationBySensors = true;
int Collect2SharedData::gFitness = -1;
std::vector<int> Collect2SharedData::gTaskSeq;
std::vector<int> Collect2SharedData::gTimeSeq;
int Collect2SharedData::gTaskIdx = 0;

std::string Collect2SharedData::gOutGenomeFile = "";
bool Collect2SharedData::gSaveGenome = true;
