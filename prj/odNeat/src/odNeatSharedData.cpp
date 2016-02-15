/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#include "odNeat/include/odNeatSharedData.h"

double odNeatSharedData::gSigmaRef = 0.0; // reference value of sigma
int odNeatSharedData::gPopulationSize = 1; // size of population per robot
bool odNeatSharedData::gClearPopulation; // empty population after replacement?
bool odNeatSharedData::gStoreOwn; // store own genome in local population?
int odNeatSharedData::gEvaluationTime = 0; // how long a controller will be evaluated on a robot

int odNeatSharedData::gSelectionMethod = 0; // Random by default

// global variable local to file -- TODO: move specific properties loader in dedicated WorldObserver
bool odNeatSharedData::gPropertiesLoaded = false;

bool odNeatSharedData::gCommunicationBySensors = true;
int odNeatSharedData::gFitness = -1;

std::string odNeatSharedData::gOutGenomeFile = "";
bool odNeatSharedData::gSaveGenome = true;

double odNeatSharedData::gDefaultInitialEnergy = 100;
double odNeatSharedData::gEnergyThreshold = 0;
double odNeatSharedData::gMaxEnergy = 2*odNeatSharedData::gDefaultInitialEnergy;

int odNeatSharedData::gMaturationPeriod = 50;

double odNeatSharedData::gCompatThreshold = 10.0;

int odNeatSharedData::gFitnessFreq = 10;
int odNeatSharedData::gTabuTimeout = 15;
double odNeatSharedData::gTabuThreshold = 1.0;


double odNeatSharedData::gEnergyItemValue = 10.0;
double odNeatSharedData::gEnergyConsumption = 0.1;

bool odNeatSharedData::gUpdateGC = true;
