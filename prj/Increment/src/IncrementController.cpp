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
    _minValue = -1.0;
	_maxValue = 1.0;
	_currentSigma = IncrementSharedData::gSigmaRef;
    _currentFitness = 0.0;
    _populationSize = IncrementSharedData::gPopulationSize;
	resetRobot();
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
    //std::cout << "W:" << nbGene << std::endl;
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
    stepBehaviour();

    //_wm->_desiredTranslationalValue = 0.0;
    //_wm->_desiredRotationalVelocity = 0.0;
}


// ################ BEHAVIOUR METHOD(S)      ################
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
        int type = 1; //sensing distance to energy item, 1.0 if not energy item
        if ( PhysicalObject::isInstanceOf(objectId) )
        {
            if ( 1 == gPhysicalObjects[objectId - gPhysicalObjectIndexStartOffset]->getType() )
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
    
    nn->setWeigths(_parameters); // create NN

    nn->setInputs(*inputs);

    nn->step();
    
    std::vector<double> outputs = nn->readOut();
    /*
    if ( gWorld->getIterations()%10 == 0 )
    {
        for(unsigned int i=0; i < inputs->size(); i++)
        {
            std::cout << (*inputs)[i] << ", ";
        }
        std::cout << " # ";

        for(unsigned int i=0; i< outputs.size(); i++)
        {
            std::cout << outputs[i] << ", ";
        }
        std::cout << std::endl << std::flush;
    }*/

    //Direct Kinematic model
    //_wm->_desiredTranslationalValue = outputs[0]; _wm->_desiredRotationalVelocity = outputs[1];

    //Differential model
    double lW = outputs[0];
    double rW = outputs[1];
    _wm->_desiredTranslationalValue = (rW + lW) / 2; _wm->_desiredRotationalVelocity = (lW - rW) / 2;
    
    // normalize to motor interval values
    _wm->_desiredTranslationalValue =  _wm->_desiredTranslationalValue * gMaxTranslationalSpeed;//gMaxTranslationalSpeed;
    _wm->_desiredRotationalVelocity = _wm->_desiredRotationalVelocity * gMaxRotationalSpeed; //0.0;
    
    delete (inputs);
}



// ################ EVOLUTION ENGINE METHODS ################
void IncrementController::stepEvolution()
{
    // * broadcasting genome : robot broadcasts its genome to all neighbors (contact-based wrt proximity sensors)

    broadcastGenome();
    
	// * lifetime ended: replace genome (if possible)
    if( dynamic_cast<IncrementWorldObserver*>(gWorld->getWorldObserver())->getLifeIterationCount()
            >= IncrementSharedData::gEvaluationTime-1 )
	{
        loadNewGenome();
        std::cout << "R" << _wm -> _id << " Fit: " << _currentFitness << std::endl;
        _currentFitness = 0.0;
    }
    
    if ( getNewGenomeStatus() ) // check for new NN parameters
	{
		reset();
		setNewGenomeStatus(false);
	}
}

void IncrementController::loadNewGenome()
{
    // case: 1+ genome(s) imported, random pick.
    if (_genomesList.size() > 0)
        selectRandomGenome();
    else
        // case: no imported genome
        resetRobot(); // destroy then create a new NN
}

void IncrementController::selectRandomGenome()
{
    if(_genomesList.size() != 0)
    {
        int randomIndex = rand()%_genomesList.size();
        std::map<int, std::vector<double> >::iterator it = _genomesList.begin();
        while (randomIndex !=0 )
        {
            it ++;
            randomIndex --;
        }
        
        _currentGenome = (*it).second;
        
        mutate(_currentSigma);
        
        setNewGenomeStatus(true);
        
        _birthdate = gWorld->getIterations();

        _genomesList.clear();
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
                targetRobotController->storeGenome(_currentGenome, _wm->getId(), _birthdate);
            }
        }
    }
}

void IncrementController::storeGenome(std::vector<double> genome, int senderId, int senderBirthdate)
{
    _genomesList[senderId] = genome;
    _birthdateList[senderId] = senderBirthdate;
}
