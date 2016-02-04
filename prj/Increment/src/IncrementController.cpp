/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Increment/include/IncrementController.h"
#include "Increment/include/IncrementWorldObserver.h"

#include "World/World.h"
#include "Utilities/Misc.h"
#include <math.h>
#include <string>

#include <neuralnetworks/MLP.h>
#include <neuralnetworks/Perceptron.h>
#include <neuralnetworks/Elman.h>

using namespace Neural;

IncrementController::IncrementController( RobotWorldModel *wm )
{
    _wm = wm; nn = NULL;

    // evolutionary engine
    _minValue = -3.0;
    _maxValue = 3.0;
	_currentSigma = IncrementSharedData::gSigmaRef;
    _currentFitness = 0.0;
    _populationSize = IncrementSharedData::gPopulationSize;
    _genomeId.robot_id = _wm->_id;
    _genomeId.gene_id = 0;
	resetRobot();
    _lifetime = -1;
    // behaviour
    _iteration = 0; _birthdate = 0;
    _wm->updateLandmarkSensor();
    _wm->setRobotLED_colorValues(255, 0, 0);
    if(_wm->_id == 0)
    {
        std::cout << "W:" << computeRequiredNumberOfWeights() << std::endl;
        std::cout << "(I:" << nn->getNbInputs();
        if (IncrementSharedData::gWithBias)
            std::cout << ", B)";
        else
            std::cout << ")";
        std::cout << ", O: " << nn->getNbOutputs() << std::endl;
    }
}

IncrementController::~IncrementController()
{
    _parameters.clear();
    delete nn;
    nn = NULL;
}

void IncrementController::reset()
{   
    _parameters.clear();
    _parameters = _genome;
}

void IncrementController::resetRobot()
{
    //Number of sensors
    _nbInputs = 0;

    //proximity sensors
    _nbInputs += _wm->_cameraSensorsNb;

    //If task=collect, add object sensors
    if (IncrementSharedData::gFitness == 1)
    {
        // gathering object distance
        _nbInputs +=  _wm->_cameraSensorsNb;
    }

    //Number of effectors
    _nbOutputs = 2;

    //NN structure
    _nbHiddenLayers = IncrementSharedData::gNbHiddenLayers;
    _nbNeuronsPerHiddenLayer = new std::vector<unsigned int>(_nbHiddenLayers);
    for(unsigned int i = 0; i < _nbHiddenLayers; i++)
        (*_nbNeuronsPerHiddenLayer)[i] = IncrementSharedData::gNbNeuronsPerHiddenLayer;

    createNN();

    unsigned int const nbGene = computeRequiredNumberOfWeights();

    _genome.clear();
    double w;
    //Random genome (weights) initialization
    for ( unsigned int i = 0 ; i != nbGene ; i++ )
    {
        // weights: random init between -1 and +1
        w = (double)(rand() % 10000)/5000.0 - 1.0;
        _genome.push_back(w);
    }
    _currentGenome = _genome;
    setNewGenomeStatus(true);
    _genomesList.clear();
    _fitnessList.clear();

}

void IncrementController::createNN()
{
    delete nn;
    // MLP
    nn = new MLP(_parameters, _nbInputs, _nbOutputs, *(_nbNeuronsPerHiddenLayer),
                 IncrementSharedData::gWithBias);
}

unsigned int IncrementController::computeRequiredNumberOfWeights()
{
    unsigned int res = nn->getRequiredNumberOfWeights();
    return res;
}

void IncrementController::step()
{
	_iteration++;
    stepEvolution();
    stepBehaviour();
    double distance_sensor;
    double coef_obstacle = 1.0;

    double trans = _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
    double rot = _wm->_desiredRotationalVelocity  / gMaxRotationalSpeed;
    double deltaFitness;
    //Fitness measurement and update
    switch (IncrementSharedData::gFitness) {
    case 0:
        //Floreano's increment
        for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
        {
            distance_sensor = _wm->getDistanceValueFromCameraSensor(i) /
                    _wm->getCameraSensorMaximumDistanceValue(i);
            if (distance_sensor < coef_obstacle)
                coef_obstacle = distance_sensor;
        }

        //deltaFitness: fitness contribution at current time-step
        //abs(trans) in [0,1], (1 - abs(rot)) in [0,1], coeffObstacle is in [0,1]
        deltaFitness = fabs(trans) * (1 - fabs(rot)) * coef_obstacle;
        //Incrementally averaging deltas
        _currentFitness = (_currentFitness * _lifetime + deltaFitness) / (_lifetime + 1);
        break;
    case 1:
        //Counting items: already done in agent observer
        break;
    default:
        break;
    }
}


//################BEHAVIOUR METHODS################
void IncrementController::stepBehaviour()
{

    // ---- Build inputs ----

    std::vector<double>* inputs = new std::vector<double>(_nbInputs);
    int inputToUse = 0;
    
    // distance sensors
    for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
    {
        (*inputs)[inputToUse] = _wm->getDistanceValueFromCameraSensor(i) /
                _wm->getCameraSensorMaximumDistanceValue(i);
        inputToUse++;
        
        //If task=collect, add object sensors
        if (IncrementSharedData::gFitness == 1)
        {
            int objectId = _wm->getObjectIdFromCameraSensor(i);
            // input: physical object?
            //sensing distance to energy item, 1.0 if not energy item
            int type = 1;
            if ( PhysicalObject::isInstanceOf(objectId) )
            {
                if ( type == gPhysicalObjects[objectId - gPhysicalObjectIndexStartOffset]->getType() )
                  (*inputs)[inputToUse] = _wm->getDistanceValueFromCameraSensor(i) /
                        _wm->getCameraSensorMaximumDistanceValue(i);
                else
                    (*inputs)[inputToUse] = 1.0;
                inputToUse++;

            }
            else
            {
                //Not a physical object.
                //But: should still fill in the inputs (max distance, 1.0)
                (*inputs)[inputToUse] = 1.0;
                inputToUse++;
            }
       }
    }
    

    // ---- compute and read out ----
    nn->setWeigths(_parameters); // set genome
    nn->setInputs(*inputs);
    nn->step();
    std::vector<double> outputs = nn->readOut();

    //Direct Kinematic model
    //_wm->_desiredTranslationalValue = outputs[0]; _wm->_desiredRotationalVelocity = outputs[1];

    //Differential model
    double lW = outputs[0];
    double rW = outputs[1];
    _wm->_desiredTranslationalValue = (rW + lW) / 2;
    _wm->_desiredRotationalVelocity = (lW - rW) / 2;
    
    // normalize to motor interval values
    _wm->_desiredTranslationalValue =  _wm->_desiredTranslationalValue * gMaxTranslationalSpeed;
    _wm->_desiredRotationalVelocity = _wm->_desiredRotationalVelocity * gMaxRotationalSpeed;
    
    delete (inputs);
}



//################ EVOLUTION ENGINE METHODS ################
void IncrementController::stepEvolution()
{
    //broadcasting genome : robot broadcasts its genome
    //to all neighbors (contact-based wrt proximity sensors)
    broadcastGenome();
    _lifetime++;
    //agent's lifetime ended: replace genome (if possible)
    if(_lifetime >= IncrementSharedData::gEvaluationTime)
	{
        if (IncrementSharedData::gStoreOwn)
            storeOwnGenome();

        loadNewGenome();
        _currentFitness = 0.0;
        _lifetime = 0;

       if (IncrementSharedData::gClearPopulation)
       {
           _genomesList.clear();
           _fitnessList.clear();
       }
    }
    // check for new NN parameters
    if ( getNewGenomeStatus() )
	{
		reset();
        _genomeId.gene_id += 1;
		setNewGenomeStatus(false);        
	}
}

void IncrementController::loadNewGenome()
{
    switch(IncrementSharedData::gSelectionMethod)
    {
        case 0:
            // case: 1+ genome(s) imported, random pick.
            if (_genomesList.size() > 0)
            {
                selectRandomGenome();
            }
            else
            {
                // case: no imported genome
                resetRobot(); // destroy then create a new NN
            }
            break;
       case 1:
            // case: 1+ genome(s) imported, select best.
            if (_genomesList.size() > 0)
            {
                selectBestGenome();
            }
            else
            {
                // case: no imported genome do nothing
            }
            break;
        case 2:
             // case: 1+ genome(s) imported, select rank-based.
             if (_genomesList.size() > 0)
             {
                 selectRankBasedGenome();
             }
             else
             {
                 // case: no imported genome do nothing
             }
             break;
       default:
            std::cout << "Error: unknown selection method: "
                      << IncrementSharedData::gSelectionMethod << std::endl;
            exit(-1);
    }
}

void IncrementController::selectRandomGenome()
{
    if(_genomesList.size() != 0)
    {
        int randomIndex = rand()%_genomesList.size();
        std::map<GC, std::vector<double> >::iterator it = _genomesList.begin();
        while (randomIndex !=0 )
        {
            it ++;
            randomIndex --;
        }
        
        _currentGenome = (*it).second;
        
        mutate(_currentSigma);
        
        setNewGenomeStatus(true);
        
        _birthdate = gWorld->getIterations();
    }
}

void IncrementController::selectBestGenome()
{
    if(_genomesList.size() != 0)
    {
        double maxFitness = -1.0;
        struct GC indexBest;
        indexBest.robot_id = -1; indexBest.gene_id = -1;
        std::map<GC, double>::iterator it = _fitnessList.begin();
        for(; it != _fitnessList.end(); it++)
        {
            if ((*it).second > maxFitness)
            {
                maxFitness = (*it).second;
                indexBest = (*it).first;
            }
        }

        //if best individual found => replace
        if (indexBest.robot_id != -1)
        {
            _currentGenome = _genomesList[indexBest];
        }
        else
        {
            //This should never happen
            std::cout << "ERROR in best selection" << std::endl;
            exit(-1);
            //Do nothing (current genome will be mutated)
        }

        mutate(_currentSigma);

        setNewGenomeStatus(true);
        _birthdate = gWorld->getIterations();
    }
}

void IncrementController::selectRankBasedGenome()
{
    unsigned int n = _genomesList.size();
    if(n != 0)
    {
        double total_prob = n * (n + 1) / 2;

        int i = 0;
        double p = (rand() / static_cast<double>(RAND_MAX)) * total_prob;
        //choose i-th index with i/total_prob
        while((p -= (i+1)) > 0)
            i++;
        std::vector<std::pair<GC, double>> pairs;
        for (auto itr = _fitnessList.begin(); itr != _fitnessList.end(); ++itr)
            pairs.push_back(*itr);
        
        sort(pairs.begin(), pairs.end(),
             [=](const std::pair<GC, double> &a,
                 const std::pair<GC, double> &b){ return a.second < b.second;});

        _currentGenome = _genomesList[pairs[i].first];

        mutate(_currentSigma);

        setNewGenomeStatus(true);
        _birthdate = gWorld->getIterations();
    }
}

void IncrementController::mutate( float sigma) // mutate within bounds.
{
	_genome.clear();
    
	_currentSigma = sigma;
	for (unsigned int i = 0 ; i != _currentGenome.size() ; i++ )
	{
		double value = _currentGenome[i] + getGaussianRand(0,_currentSigma);
		// bouncing upper/lower bounds
		if ( value < _minValue )
		{
			double range = _maxValue - _minValue;
			double overflow = - ( (double)value - _minValue );
			overflow = overflow - 2*range * (int)( overflow / (2*range) );
			if ( overflow < range )
				value = _minValue + overflow;
			else // overflow btw range and range*2
				value = _minValue + range - (overflow-range);
		}
		else if ( value > _maxValue )
		{
			double range = _maxValue - _minValue;
			double overflow = (double)value - _maxValue;
			overflow = overflow - 2*range * (int)( overflow / (2*range) );
			if ( overflow < range )
				value = _maxValue - overflow;
			else // overflow btw range and range*2
				value = _maxValue - range + (overflow-range);
		}
        
		_genome.push_back(value);
	}
    
	_currentGenome = _genome;
    
    //gLogManager->write(s); gLogManager->flush();
}

void IncrementController::updateFitness(double delta)
{
    _currentFitness += delta;
}

// ################ COMMUNICATION METHODS ################
void IncrementController::broadcastGenome()
{
    for( int i = 0 ; i < _wm->_cameraSensorsNb; i++)
    {
        int targetIndex = _wm->getObjectIdFromCameraSensor(i);

        // sensor ray bumped into a robot : communication is possible
        if ( targetIndex >= gRobotIndexStartOffset )
        {
            // convert image registering index into robot id.
            targetIndex = targetIndex - gRobotIndexStartOffset;

            IncrementController* targetRobotController =
                    dynamic_cast<IncrementController*>
                    (gWorld->getRobot(targetIndex)->getController());

            if ( ! targetRobotController )
            {
                std::cerr << "Error: observer not compatible" << std::endl;
                exit(-1);
            }

            // other agent stores my genome.
            //Genome id as gene clock (different for each agent)
            targetRobotController->storeGenome(_currentGenome, _genomeId, _currentFitness);
        }
    }
}

void IncrementController::storeGenome(std::vector<double> genome, GC senderId, double fitness)
{
    if (_genomesList.size() < _populationSize)
    {
        _genomesList[senderId] = genome;
        _fitnessList[senderId] = fitness;
    }
    else
    {
        //TOTEST replace worse individual in population
        double minFitness = std::numeric_limits<double>::max();
        struct GC indexWorse;
        indexWorse.robot_id = -1; indexWorse.gene_id = -1;
        std::map<GC, double>::iterator it = _fitnessList.begin();
        for(; it != _fitnessList.end(); it++)
        {
            if ((*it).second < minFitness)
            {
                minFitness = (*it).second;
                indexWorse = (*it).first;
            }
        }
        if ((indexWorse.robot_id != -1) && (minFitness < fitness))
        {
            _fitnessList.erase(indexWorse);
            _genomesList.erase(indexWorse);

            _fitnessList[senderId] = fitness;
            _genomesList[senderId] = genome;
        }
    }

}

void IncrementController::storeOwnGenome()
{
    //std::cout << _currentFitness << std::endl;
    storeGenome(_currentGenome, _genomeId, _currentFitness);
    //std::cout << "Pop: " << _fitnessList.size() << " # " << _genomesList.size() << std::endl;
}

void IncrementController::logGenome(std::string s)
{
    std::ofstream genomeF;
    genomeF.open(s);
    genomeF << "[W\n";
    genomeF << std::fixed << std::setprecision(10);

    for(auto it = _currentGenome.begin(); it != _currentGenome.end(); ++it)
    {
        genomeF << (*it) << "\n";
    }

    genomeF << "[F\n" << _currentFitness << "\n";
    //Structure (NN)
    genomeF << "[I\n" << _nbInputs << "\n";//Inputs
    if(IncrementSharedData::gWithBias)
        genomeF << "[B\n" << 1 << "\n";//Bias
    else
        genomeF << "[B\n" << 0 << "\n";
    genomeF << "[H\n" << _nbHiddenLayers << "\n"
            << "[N\n" << IncrementSharedData::gNbNeuronsPerHiddenLayer << "\n";//Neur per hidden layer
    genomeF << "[O\n" << _nbOutputs << "\n";

    genomeF.close();
}

//#################AUXILIARY#############################
std::ostream& operator<<(std::ostream& os, const GC& gene_clock)
{
    os << "(R" << gene_clock.robot_id << "; G" << gene_clock.gene_id << ")";
    return os;
}

