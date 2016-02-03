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
    _lifetime = 0;
    // behaviour
    _iteration = 0; _birthdate = 0;
    _wm->updateLandmarkSensor();
    _wm->setRobotLED_colorValues(255, 0, 0);
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

    // gathering object distance
    _nbInputs =  _wm->_cameraSensorsNb;
    //proximity sensors
    _nbInputs += _wm->_cameraSensorsNb;

    //Number of effectors
    _nbOutputs = 2;

    //NN structure
    _nbHiddenLayers = IncrementSharedData::gNbHiddenLayers;
    _nbNeuronsPerHiddenLayer = new std::vector<unsigned int>(_nbHiddenLayers);
    for(unsigned int i = 0; i < _nbHiddenLayers; i++)
        (*_nbNeuronsPerHiddenLayer)[i] = IncrementSharedData::gNbNeuronsPerHiddenLayer;

    createNN();

    unsigned int const nbGene = computeRequiredNumberOfWeights();
    if(_wm->_id == 0)
    {
        std::cout << "W:" << nbGene << std::endl;
        std::cout << "I:" << nn->getNbInputs()<< ", O: " << nn->getNbOutputs() << std::endl;
    }
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
    nn = new MLP(_parameters, _nbInputs, _nbOutputs, *(_nbNeuronsPerHiddenLayer), IncrementSharedData::gWithBias);
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
    _lifetime++;
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
        //(1+trans)/2 in [0,1], (1 - abs(rot)) in [0,1], coeffObstacle is in [0,1]
        //deltaFitness: fitness contribution at current time-step
        deltaFitness = (1+trans)/2 * (1 - abs(rot)) * coef_obstacle;
        //std::cout << "delta:" << deltaFitness << ", life: " << _lifetime << std::endl;

        _currentFitness = (_currentFitness * (_lifetime - 1) + deltaFitness) / _lifetime;
        std::cout << _currentFitness << std::endl;
        break;
    case 1:
        //Counting items: already done in agent observer
        break;
    default:
        break;
    }
}


// ################ BEHAVIOUR METHOD(S)################
void IncrementController::stepBehaviour()
{

    // ---- Build inputs ----

    std::vector<double>* inputs = new std::vector<double>(_nbInputs);
    int inputToUse = 0;
    
    // distance sensors
    for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
    {
        (*inputs)[inputToUse] = _wm->getDistanceValueFromCameraSensor(i) / _wm->getCameraSensorMaximumDistanceValue(i);
        inputToUse++;
        
        int objectId = _wm->getObjectIdFromCameraSensor(i);

        // input: physical object?
        //sensing distance to energy item, 1.0 if not energy item
        int type = 1;
        if ( PhysicalObject::isInstanceOf(objectId) )
        {
            if ( type == gPhysicalObjects[objectId - gPhysicalObjectIndexStartOffset]->getType() )
              (*inputs)[inputToUse] = _wm->getDistanceValueFromCameraSensor(i) / _wm->getCameraSensorMaximumDistanceValue(i);
            else
                (*inputs)[inputToUse] = 1.0;
            inputToUse++;

        }
        else
        {
            // not a physical object. But: should still fill in the inputs (max distance, 1.0)
            (*inputs)[inputToUse] = 1.0;
            inputToUse++;
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
    _wm->_desiredTranslationalValue = (rW + lW) / 2; _wm->_desiredRotationalVelocity = (lW - rW) / 2;
    
    // normalize to motor interval values
    _wm->_desiredTranslationalValue =  _wm->_desiredTranslationalValue * gMaxTranslationalSpeed;
    _wm->_desiredRotationalVelocity = _wm->_desiredRotationalVelocity * gMaxRotationalSpeed;
    
    delete (inputs);
}



//################ EVOLUTION ENGINE METHODS ################
void IncrementController::stepEvolution()
{
    // * broadcasting genome : robot broadcasts its genome to all neighbors (contact-based wrt proximity sensors)

    broadcastGenome();
    
	// * lifetime ended: replace genome (if possible)
    if( dynamic_cast<IncrementWorldObserver*>(gWorld->getWorldObserver())->getLifeIterationCount()
            >= IncrementSharedData::gEvaluationTime-1 )
	{
        if (IncrementSharedData::gStoreOwn)
            storeOwnGenome();

        loadNewGenome();
        //std::cout << "R" << _wm -> _id << " Fit: " << _currentFitness << std::endl;
        //std::cout << _currentFitness << std::endl;
        _currentFitness = 0.0;
       //std::cout << "Local pop size: " <<_fitnessList.size() << std::endl;
       if (IncrementSharedData::gClearPopulation)
       {
           _genomesList.clear();
           _fitnessList.clear();
       }
    }
    
    if ( getNewGenomeStatus() ) // check for new NN parameters
	{
		reset();
        _genomeId.gene_id += 1;
		setNewGenomeStatus(false);
        _lifetime = 0;
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
       default:
            std::cout << "Error: unknown selection method: " << IncrementSharedData::gSelectionMethod << std::endl;
            exit(-1);
    }



}

void IncrementController::selectRandomGenome()
{
    if(_genomesList.size() != 0)
    {
        int randomIndex = rand()%_genomesList.size();
        //std::map<int, std::vector<double> >::iterator it = _genomesList.begin();
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
        //std::map<int, double>::iterator it = _fitnessList.begin();
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
    if ( _wm->isAlive() == true )  	// only if agent is active (ie. not just revived) and deltaE>0.
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
                        dynamic_cast<IncrementController*>(gWorld->getRobot(targetIndex)->getController());
                
                if ( ! targetRobotController )
                {
                    std::cerr << "Error: observer not compatible" << std::endl;
                    exit(-1);
                }

                // other agent stores my genome.
                //TODO sender ID should be genomeId not robotId. TODO change also in storeOwnGenome
                //targetRobotController->storeGenome(_currentGenome, _wm->getId(), _currentFitness);
                targetRobotController->storeGenome(_currentGenome, _genomeId, _currentFitness);
            }
        }
    }
}

void IncrementController::storeGenome(std::vector<double> genome, GC senderId, double fitness)
{
    /*std::map<GC, double>::iterator it = _fitnessList.begin();
    for(; it != _fitnessList.end(); it++)
    {
        std::cout << (*it).first << ", fit: " << (*it).second << std::endl;
    }*/
    if (_genomesList.size() < _populationSize)
    {
        _genomesList[senderId] = genome;
        _fitnessList[senderId] = fitness;
    }
    else
    {
        //TODO replace worse individual in population
        double minFitness = std::numeric_limits<double>::max();
        struct GC indexWorse;
        indexWorse.robot_id = -1; indexWorse.gene_id = -1;
        //std::map<int, double>::iterator it = _fitnessList.begin();
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
            std::cout << "Replacing" << std::endl;
        }
    }

}

void IncrementController::storeOwnGenome()
{
    //std::cout << "Storing own: " << _genomeId << std::endl;
    storeGenome(_currentGenome, _genomeId, _currentFitness);
}


//#################AUXILIARY#############################
std::ostream& operator<<(std::ostream& os, const GC& gene_clock)
{
    os << "(R" << gene_clock.robot_id << "; G" << gene_clock.gene_id << ")";
    return os;
}

