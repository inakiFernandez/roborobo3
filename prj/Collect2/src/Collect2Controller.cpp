/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Collect2/include/Collect2Controller.h"
#include "Collect2/include/Collect2WorldObserver.h"

#include "World/World.h"
#include "Utilities/Misc.h"
#include <math.h>
#include <string>

#include <neuralnetworks/MLP.h>
#include <neuralnetworks/Perceptron.h>
#include <neuralnetworks/Elman.h>

using namespace Neural;

Collect2Controller::Collect2Controller( RobotWorldModel *wm )
{
    _wm = wm; nn = NULL;

    // neural weights limits
    _minValue = -Collect2SharedData::gNeuronWeightRange/2;
    _maxValue = Collect2SharedData::gNeuronWeightRange/2;
	_currentSigma = Collect2SharedData::gSigmaRef;

    //Evo variables
    _currentFitness = 0.0;
    _genomeId.robot_id = _wm->_id;
    _genomeId.gene_id = 0;

    resetRobot();

    _lifetime = -1; _iteration = 0; _birthdate = 0;

    _wm->updateLandmarkSensor();
    _wm->setRobotLED_colorValues(255, 0, 0);

    if(_wm->_id == 0)
    {
        std::cout << "W:" << computeRequiredNumberOfWeights() << std::endl;
        std::cout << "(I:" << nn->getNbInputs();
        if (Collect2SharedData::gWithBias)
            std::cout << ", B)";
        else
            std::cout << ")";
        std::cout << ", O: " << nn->getNbOutputs() << std::endl;       
    }
}

Collect2Controller::~Collect2Controller()
{
    _parameters.clear();
    delete nn;
    nn = NULL;
}

void Collect2Controller::reset()
{   
    _parameters.clear();
    _parameters = _currentGenome;
}

void Collect2Controller::resetRobot()
{
    //Number of sensors
    _nbInputs = 0;

    //proximity sensors
    _nbInputs += _wm->_cameraSensorsNb;

    //If task=collect, add object sensors
    if ((Collect2SharedData::gFitness == 1) || (Collect2SharedData::gFitness == 2))
    {
        // gathering object distance
        _nbInputs +=  _wm->_cameraSensorsNb;
    }

    if(gInitialNumberOfRobots > 1)
        _nbInputs += 2; //closest agent relative orientation and is it collecting?

    //Current translational and rotational speeds
    //_nbInputs += 2;

    //Number of effectors
    _nbOutputs = 2;

    //NN structure
    _nbHiddenLayers = Collect2SharedData::gNbHiddenLayers;
    _nbNeuronsPerHiddenLayer = new std::vector<unsigned int>(_nbHiddenLayers);
    for(unsigned int i = 0; i < _nbHiddenLayers; i++)
        (*_nbNeuronsPerHiddenLayer)[i] = Collect2SharedData::gNbNeuronsPerHiddenLayer;

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

void Collect2Controller::createNN()
{
    delete nn;
    switch ( Collect2SharedData::gControllerType )
    {
        case 0:
            nn = new MLP(_parameters, _nbInputs, _nbOutputs,
                         *(_nbNeuronsPerHiddenLayer), Collect2SharedData::gWithBias);
            break;
        case 1:
            nn = new Perceptron(_parameters, _nbInputs, _nbOutputs);
            break;
        case 2:
            nn = new Elman(_parameters, _nbInputs, _nbOutputs, *(_nbNeuronsPerHiddenLayer));
            break;
        default: // default: no controller
            std::cerr << "[ERROR] gController type unknown (value: "
                      << Collect2SharedData::gControllerType << ").\n";
            exit(-1);
    };

}

unsigned int Collect2Controller::computeRequiredNumberOfWeights()
{
    unsigned int res = nn->getRequiredNumberOfWeights();
    return res;
}

void Collect2Controller::step()
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
    switch (Collect2SharedData::gFitness) {
    case 0:
        //Navigation fitness instant increment
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
    case 2:
        //Counting items cooperatively:  already done in agent observer
        break;
    default:
        break;
    }
}


//################BEHAVIOUR METHODS################
void Collect2Controller::stepBehaviour()
{
    // ---- Build inputs ----
    std::vector<double>* inputs = new std::vector<double>(_nbInputs);
    int inputToUse = 0;
    
    int type = 1; //object type= energy
    for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
    {
        // distance sensors
        //If there is object but not "Physical object"
            //=> wall (get distance sensor) else 0.0
        int objId = _wm->getObjectIdFromCameraSensor(i);
        if ( PhysicalObject::isInstanceOf(objId) )
        {
            //(*inputs)[inputToUse] = 1.0;
            (*inputs)[inputToUse] = 0.0;
            inputToUse++;
        }
        else
        {
            //(*inputs)[inputToUse] = _wm->getDistanceValueFromCameraSensor(i) /
             //_wm->getCameraSensorMaximumDistanceValue(i);
            (*inputs)[inputToUse] = 1.0 - _wm->getDistanceValueFromCameraSensor(i) /
                          _wm->getCameraSensorMaximumDistanceValue(i);
            inputToUse++;
        }
        
        //If task=collect, add object sensors
        if ((Collect2SharedData::gFitness == 1) ||(Collect2SharedData::gFitness == 2))
        {
            int objectId = _wm->getObjectIdFromCameraSensor(i);
            // input: physical object?
            //sensing distance to energy item, 0.0 [not 1.0] if not energy item
            if ( PhysicalObject::isInstanceOf(objectId) )
            {
                if ( type == gPhysicalObjects[objectId - gPhysicalObjectIndexStartOffset]->getType() )
                  //(*inputs)[inputToUse] = _wm->getDistanceValueFromCameraSensor(i) /
                  //      _wm->getCameraSensorMaximumDistanceValue(i);
                  (*inputs)[inputToUse] = 1.0 - _wm->getDistanceValueFromCameraSensor(i) /
                      _wm->getCameraSensorMaximumDistanceValue(i);
                else
                    //(*inputs)[inputToUse] = 1.0;
                    (*inputs)[inputToUse] = 0.0;
                inputToUse++;
            }
            else
            {
                //Not a physical object. But: should still fill in the inputs 0.0
                //(*inputs)[inputToUse] = 1.0;
                (*inputs)[inputToUse] = 0.0;
                inputToUse++;
            }
       }
    }

    //Angle wrt closest agent, relative orientation
    if(gInitialNumberOfRobots > 1)
    {
        int closestAgent = -1;
        double minDistance = 10000.0;
        for(int a = 0; a < gNumberOfRobots; a++)
        {
            if(a != _wm->getId())
            {
                double distance = sqrt(
                        pow(_wm->getXReal() - gWorld->getRobot(a)->getWorldModel()->getXReal(),2)
                        + pow(_wm->getYReal() - gWorld->getRobot(a)->getWorldModel()->getYReal() ,2));
                if(distance < minDistance)
                {
                   minDistance = distance;
                   closestAgent = a;
                }
            }
        }
        double srcOrientation = _wm->_agentAbsoluteOrientation + 90.0;

        if(srcOrientation > 180.0)
            srcOrientation = srcOrientation - 360.0;
        srcOrientation = srcOrientation / 180.0;

        double x = _wm->getXReal();
        double y = _wm->getYReal();
        double aX = gWorld->getRobot(closestAgent)->getWorldModel()->getXReal();
        double aY = gWorld->getRobot(closestAgent)->getWorldModel()->getYReal();

        double dX = (aX - x);
        double dY = (y - aY);
        double angle;
        if(y > aY)
            angle = atan(dX / dY);
        else
        {
            if(x < aX)
                angle = atan(dX / dY) + M_PI;
            else
                angle = atan(dX / dY) - M_PI;
        }

        angle = angle / M_PI;
        angle = angle - (srcOrientation);

        //angle stores the orientation towards the nearest robot
        (*inputs)[inputToUse] = angle;
        inputToUse++;

        //closest agent : is collecting?
        int targetIndexCloseAgent = gWorld->getRobot(closestAgent)->getWorldModel()->getGroundSensorValue();
        bool isCollecting = PhysicalObject::isInstanceOf(targetIndexCloseAgent);
        if (isCollecting)
           (*inputs)[inputToUse] = 1.0;
        else
           (*inputs)[inputToUse] = -1.0;

        inputToUse++;
    }

    //Previous translational and rotational speeds (acts as recurrent connections from last step)
    //(*inputs)[inputToUse] = _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
    //inputToUse++;
    //(*inputs)[inputToUse] = _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
    //inputToUse++;
    /*if(_wm->_id == 0)
    {
        for(auto it = inputs->begin(); it < inputs->end(); it++)
            std::cout << (*it) << " | ";
        std::cout << std::endl;
    }*/

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
void Collect2Controller::stepEvolution()
{
    //broadcasting genome : robot broadcasts its genome
    //to all neighbors (contact-based wrt proximity sensors)
    broadcastGenome();
    _lifetime++;
    //agent's lifetime ended: replace genome (if possible)
    if(_lifetime >= Collect2SharedData::gEvaluationTime)
    {
        if (Collect2SharedData::gStoreOwn)
            storeOwnGenome();

        loadNewGenome();
        //std::cout << _currentFitness << std::endl;
        //for(auto it = _fitnessList.begin(); it != _fitnessList.end(); it++)
        //    std::cout << it->second << std::endl;
        _currentFitness = 0.0;
        _lifetime = 0;

       if (Collect2SharedData::gClearPopulation)
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

void Collect2Controller::loadNewGenome()
{
   // If 1+ genome(s) imported, select best.
   if (_genomesList.size() > 0)
   {
       selectBestGenome();
       //selectRankBasedGenome();
   }
   else
   {
        // case: no imported genome: mutate current
   }
   mutate(_currentSigma);
   setNewGenomeStatus(true);
   _birthdate = gWorld->getIterations();
}

void Collect2Controller::selectBestGenome()
{

    double maxFitness = std::numeric_limits<double>::lowest();

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
    }

}
void Collect2Controller::selectRankBasedGenome()
{
    unsigned int n = _genomesList.size();
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
}
void Collect2Controller::mutate(float sigma) // mutate within bounds.
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

void Collect2Controller::updateFitness(double delta)
{
    _currentFitness += delta;
}

// ################ COMMUNICATION METHODS ################
void Collect2Controller::broadcastGenome()
{
    if (Collect2SharedData::gCommunicationBySensors)
    {
        for( int i = 0 ; i < _wm->_cameraSensorsNb; i++)
        {
            int targetIndex = _wm->getObjectIdFromCameraSensor(i);

            // sensor ray bumped into a robot : communication is possible
            if ( targetIndex >= gRobotIndexStartOffset )
            {
                // convert image registering index into robot id.
                targetIndex = targetIndex - gRobotIndexStartOffset;

                Collect2Controller* targetRobotController =
                        dynamic_cast<Collect2Controller*>
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
    else
    {
        for (int i = 0; i < gNumberOfRobots; i++)
        {
            //Do not send to self
            if (i != _wm->_id)
            {
                Collect2Controller* targetRobotController =
                        dynamic_cast<Collect2Controller*>
                        (gWorld->getRobot(i)->getController());

                if ( ! targetRobotController )
                {
                    std::cerr << "Error: observer not compatible" << std::endl;
                    exit(-1);
                }
                double distance = sqrt((_wm->_xReal - targetRobotController->_wm->_xReal)*
                        (_wm->_xReal - targetRobotController->_wm->_xReal) +
                        (_wm->_yReal - targetRobotController->_wm->_yReal)*
                        (_wm->_yReal - targetRobotController->_wm->_yReal));

                if (distance <= gMaxRadioDistance)
                {
                    //std::cout << "Close range => send from R" << _wm->_id << " to R"
                    //          << targetRobotController->_wm->_id << std::endl;
                    // other agent stores my genome.
                    targetRobotController->storeGenome(_currentGenome, _genomeId, _currentFitness);
                }
            }
        }
    }
}

void Collect2Controller::storeGenome(std::vector<double> genome, GC senderId, double fitness)
{
    if(_genomesList.find(senderId) != _genomesList.end())
    {
        //Update fitness
        _fitnessList[senderId] = fitness;
    }
    else
    {
        if (_genomesList.size() < (unsigned int)Collect2SharedData::gPopulationSize)
        {
            _genomesList[senderId] = genome;
            _fitnessList[senderId] = fitness;
        }
        else
        {
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
            if (minFitness <= fitness)
            {
                _fitnessList.erase(indexWorse);
                _genomesList.erase(indexWorse);

                _fitnessList[senderId] = fitness;
                _genomesList[senderId] = genome;
            }
        }
    }
}

void Collect2Controller::storeOwnGenome()
{
    storeGenome(_currentGenome, _genomeId, _currentFitness);
}

void Collect2Controller::logGenome(std::string s)
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
    if(Collect2SharedData::gWithBias)
        genomeF << "[B\n" << 1 << "\n";//Bias
    else
        genomeF << "[B\n" << 0 << "\n";
    genomeF << "[H\n" << _nbHiddenLayers << "\n"
            << "[N\n" << Collect2SharedData::gNbNeuronsPerHiddenLayer << "\n";//Neur per hidden layer
    genomeF << "[O\n" << _nbOutputs << "\n";

    genomeF.close();
}


