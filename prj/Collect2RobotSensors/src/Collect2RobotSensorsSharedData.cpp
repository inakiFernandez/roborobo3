/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "Collect2RobotSensors/include/Collect2RobotSensorsSharedData.h"

double Collect2RobotSensorsSharedData::gSigmaRef = 0.0; // reference value of sigma
int Collect2RobotSensorsSharedData::gPopulationSize = 1; // size of population per robot
bool Collect2RobotSensorsSharedData::gClearPopulation; // empty population after replacement?
bool Collect2RobotSensorsSharedData::gStoreOwn; // store own genome in local population?
int Collect2RobotSensorsSharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

int Collect2RobotSensorsSharedData::gSelectionMethod = 0; // Random by default

// global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver
bool Collect2RobotSensorsSharedData::gPropertiesLoaded = false;

int Collect2RobotSensorsSharedData::gControllerType = 0; //MLP, Perceptron, Elman
int Collect2RobotSensorsSharedData::gNbHiddenLayers = 1;
int Collect2RobotSensorsSharedData::gNbNeuronsPerHiddenLayer = 5;
int Collect2RobotSensorsSharedData::gNeuronWeightRange = 800;
bool Collect2RobotSensorsSharedData::gWithBias = false;

bool Collect2RobotSensorsSharedData::gCommunicationBySensors = true;
int Collect2RobotSensorsSharedData::gFitness = -1;
std::vector<int> Collect2RobotSensorsSharedData::gTaskSeq;
std::vector<int> Collect2RobotSensorsSharedData::gTimeSeq;
int Collect2RobotSensorsSharedData::gTaskIdx = 0;

std::string Collect2RobotSensorsSharedData::gOutGenomeFile = "";
bool Collect2RobotSensorsSharedData::gSaveGenome = true;
