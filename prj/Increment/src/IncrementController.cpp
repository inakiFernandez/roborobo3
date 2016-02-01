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

	resetRobot();

    // behaviour
    _iteration = 0; _birthdate = 0;
    _wm->updateLandmarkSensor();
    
    _wm->setAlive(true); _wm->setRobotLED_colorValues(255, 0, 0);

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
    _nbInputs = 0;

    if ( gExtendedSensoryInputs )     {
        //  isItAnAgent? + agentAngleDifference? + isItAWall? + gathering object distance
        _nbInputs = (1 + 1 + 1 + PhysicalObjectFactory::getNbOfTypes()) * _wm->_cameraSensorsNb;
    }

    _nbInputs += _wm->_cameraSensorsNb; // proximity sensors

    _nbOutputs = 2;

    _nbHiddenLayers = IncrementSharedData::gNbHiddenLayers;
    _nbNeuronsPerHiddenLayer = new std::vector<unsigned int>(_nbHiddenLayers);
    for(unsigned int i = 0; i < _nbHiddenLayers; i++)
        (*_nbNeuronsPerHiddenLayer)[i] = IncrementSharedData::gNbNeuronsPerHiddenLayer;

    createNN();

    unsigned int const nbGene = computeRequiredNumberOfWeights();

    if ( gVerbose )
        std::cout << std::flush ;

    _genome.clear();
    double w;
    for ( unsigned int i = 0 ; i != nbGene ; i++ )
    {
        // weights: random init between -1 and +1
        w = (double)(rand() % 1000)/500.0 - 1.0;
        _genome.push_back(w);
    }
    _currentGenome = _genome;
    setNewGenomeStatus(true);
    _genomesList.clear();

}

void IncrementController::createNN()
{
    delete nn;
    // MLP
    nn = new MLP(_parameters, _nbInputs, _nbOutputs, *(_nbNeuronsPerHiddenLayer));
}


void IncrementController::step()
{
	_iteration++;
    stepEvolution();
    if ( _wm->isAlive() )
	{
        stepBehaviour();
    }
    else
    {
        _wm->_desiredTranslationalValue = 0.0;
        _wm->_desiredRotationalVelocity = 0.0;
    }

}


// ################ ######################## ################
// ################ BEHAVIOUR METHOD(S)      ################
// ################ ######################## ################

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
        
        if ( gExtendedSensoryInputs )
        {
            //      - inputs[N_physicalobjecttypes]: 0 or 1 is active, other are set to zero.
            //      - inputs[2]: (a) is it a robot? (b) is it from the same group? (c) what is its relative orientation wrt. current robot
            //      - inputs[1]: is it a wall (=1), or nothing (=0)
            //      Comment: from 0 (nothing) to 2 (robot, same group) active inputs.
            
            int objectId = _wm->getObjectIdFromCameraSensor(i);

            // input: physical object? which type?
            if ( PhysicalObject::isInstanceOf(objectId) )
            {
                int nbOfTypes = 1 + PhysicalObjectFactory::getNbOfTypes();
                for ( int i = 0 ; i != nbOfTypes ; i++ )                 
                {
                    // if ( i == ( objectId - gPhysicalObjectIndexStartOffset ) ) // [bug]: discovered by Inaki F. -- solved 2014-09-21
                    if ( i == gPhysicalObjects[objectId - gPhysicalObjectIndexStartOffset]->getType() )
                        (*inputs)[inputToUse] = 1; // match
                    else
                        (*inputs)[inputToUse] = 0;
                    inputToUse++;
                }
            }
            else
            {
                // not a physical object. But: should still fill in the inputs (with zeroes)
                int nbOfTypes = PhysicalObjectFactory::getNbOfTypes();
                for ( int i = 0 ; i != nbOfTypes ; i++ )
                {
                    (*inputs)[inputToUse] = 0;
                    inputToUse++;
                }
            }

            // input: another agent? If yes: same group?
            if ( Agent::isInstanceOf(objectId) )
            {
                // this is an agent
                (*inputs)[inputToUse] = 1;
                inputToUse++;

                // relative orientation? (ie. angle difference wrt. current agent)
                double srcOrientation = _wm->_agentAbsoluteOrientation;
                double tgtOrientation = gWorld->getRobot(objectId-gRobotIndexStartOffset)->getWorldModel()->_agentAbsoluteOrientation;

                double delta_orientation = - ( srcOrientation - tgtOrientation );
                if ( delta_orientation >= 180.0 )
                    delta_orientation = - ( 360.0 - delta_orientation );
                else
                    if ( delta_orientation <= -180.0 )
                        delta_orientation = - ( - 360.0 - delta_orientation );
                (*inputs)[inputToUse] = delta_orientation/180.0;
                inputToUse++;
                
            }
            else
            {
                (*inputs)[inputToUse] = 0; // not an agent...
                inputToUse++;
                (*inputs)[inputToUse] = 0; // ...and no orientation.
                inputToUse++;

            }
            
            // input: wall or empty?
            // not empty, but cannot be identified: this is a wall.
            if ( objectId >= 0 && objectId < gPhysicalObjectIndexStartOffset )
                (*inputs)[inputToUse] = 1;
            else
                // nothing. (objectId=-1)
                (*inputs)[inputToUse] = 0;
            inputToUse++;
        }
    }
    

    // ---- compute and read out ----
    
    nn->setWeigths(_parameters); // create NN
    
    nn->setInputs(*inputs);
    /*for(unsigned int i=0; i < inputs->size(); i++)
    {
        std::cout << (*inputs)[i] << ", ";
    }
    std::cout << " # ";*/

    nn->step();
    
    std::vector<double> outputs = nn->readOut();

    /*for(unsigned int i=0; i< outputs.size(); i++)
    {
        std::cout << outputs[i] << ", ";
    }
    std::cout << std::endl << std::flush;*/

    _wm->_desiredTranslationalValue = outputs[0]; _wm->_desiredRotationalVelocity = outputs[1];

    
    // normalize to motor interval values
    _wm->_desiredTranslationalValue = _wm->_desiredTranslationalValue * gMaxTranslationalSpeed;
    _wm->_desiredRotationalVelocity = _wm->_desiredRotationalVelocity * gMaxRotationalSpeed;
    
    delete (inputs);
}





unsigned int IncrementController::computeRequiredNumberOfWeights()
{
    unsigned int res = nn->getRequiredNumberOfWeights();
    return res;
}

// ################ ######################## ################
// ################ EVOLUTION ENGINE METHODS ################
// ################ ######################## ################

void IncrementController::stepEvolution()
{
    // * broadcasting genome : robot broadcasts its genome to all neighbors (contact-based wrt proximity sensors)

    broadcastGenome();
    
	// * lifetime ended: replace genome (if possible)
    if( dynamic_cast<IncrementWorldObserver*>(gWorld->getWorldObserver())->getLifeIterationCount()
            >= IncrementSharedData::gEvaluationTime-1 )
	{
        loadNewGenome();
    }
    
    if ( getNewGenomeStatus() ) // check for new NN parameters
	{
		reset();
		setNewGenomeStatus(false);
	}
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

        // Logging
        std::string s = std::string("");
        s += "{" + std::to_string(gWorld->getIterations()) + "} [" +
                std::to_string(_wm->getId()) + "::" + std::to_string(_birthdate) +
                "] descends from [" + std::to_string((*it).first) +
                "::" + std::to_string(_birthdateList[(*it).first]) + "]\n";
        gLogManager->write(s);
        gLogManager->flush();
        
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
    
    // Logging
    std::string s = std::string("");
    s += "{" + std::to_string(gWorld->getIterations()) + "} [" + std::to_string(_wm->getId()) + "::" +
            std::to_string(_birthdate) + "] [sigma=" + std::to_string(_currentSigma) + "]\n";
    gLogManager->write(s);
    gLogManager->flush();
}





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
                    std::cerr << "Error from robot " << _wm->getId() <<
                                 " : the observer of robot " << targetIndex <<
                                 " is not compatible" << std::endl;
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


void IncrementController::loadNewGenome()
{
    if ( _wm->isAlive() || gEnergyRefill )
    {
        // Logging
        std::string s = "{" + std::to_string(gWorld->getIterations()) + "} [" + std::to_string(_wm->getId()) +
                "::" + std::to_string(_birthdate) + "] [genomeList:" + std::to_string(_genomesList.size()) + "]\n";
        gLogManager->write(s);
        gLogManager->flush();
        
        if (_genomesList.size() > 0)
        {
            // case: 1+ genome(s) imported, random pick.
            
            selectRandomGenome();
            
            _wm->setAlive(true);
            if ( _wm->getEnergyLevel() == 0 )
                _wm->setEnergyLevel(gEnergyInit);
            _wm->setRobotLED_colorValues(255, 0, 0);
            
        }
        else
        {
            // case: no imported genome - wait for new genome.
            
            // Logging
            std::string s = std::string("");
            s += "{" + std::to_string(gWorld->getIterations()) + "} [" + std::to_string(_wm->getId())
                    + "::" + std::to_string(_birthdate) + "] no_genome.\n";
            gLogManager->write(s);
            gLogManager->flush();
            
            resetRobot(); // destroy then create a new NN

            _wm->setAlive(true);
        }
        
        // log the genome
        
        if ( _wm->isAlive() )
        {
            // Logging
            std::string s = std::string("");
            s += "{" + std::to_string(gWorld->getIterations()) + "} [" + std::to_string(_wm->getId()) + "::" +
                    std::to_string(_birthdate) + "] new_genome: ";
            for(unsigned int i=0; i<_genome.size(); i++)
            {
                s += std::to_string(_genome[i]) + " ";
            }
            s += "\n";
            gLogManager->write(s);
            gLogManager->flush();
        }
    }
}

