/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Original/include/OriginalController.h"
#include "Original/include/OriginalWorldObserver.h"

#include "World/World.h"
#include "Utilities/Misc.h"
#include <math.h>
#include <string>

#include <neuralnetworks/MLP.h>
#include <neuralnetworks/Perceptron.h>
#include <neuralnetworks/Elman.h>
#include <boost/algorithm/string.hpp>


using namespace Neural;

OriginalController::OriginalController( RobotWorldModel *wm )
{
    _wm = wm; nn = NULL;

    // neural weights limits
    _minValue = -OriginalSharedData::gNeuronWeightRange/2;
    _maxValue = OriginalSharedData::gNeuronWeightRange/2;
	_currentSigma = OriginalSharedData::gSigmaRef;

    //Evo variables
    _currentFitness = 0.0;
    _genomeId.robot_id = _wm->_id;
    _genomeId.gene_id = 0;

    //_wm->setRobotLED_colorValues(255, 0, 0);
    //Set initial color sensor to 0.0
    _wm->setRobotLED_colorValues(128, 128, 0);

    _withCollectColorEffector = OriginalSharedData::gWithCollectColorEffector;


    resetRobot();
    //braitenberg genome (8sensors)
    static const double arr[] = {0.5, -0.5, -1.0, 1.0, 0.75, -0.25, //W
                                 0.5, -0.5, -1.0, 1.0, 0.75, -0.25, //NW
                                 -0.5, -0.5, 1.0, 1.0, -0.5, -0.5, //N
                                 -0.5, 0.5, 1.0,-1.0, -0.25, 0.75, //NE
                                 -0.5, 0.5, 1.0,-1.0, -0.25, 0.75, //E
                                 -0.5, 0.5, 1.0,-1.0, -0.25, 0.75, //SE
                                 0.5, 0.5, -1.0, -1.0, 0.5, 0.5, //S
                                 0.5, -0.5, -1.0, 1.0, 0.75, 0.0, //SW
                                 0.5, 0.5, //bias
                                 0.0, 0.0, //lW recurrent
                                 0.0, 0.0, //rW recurrent
                                };


    std::vector<double> vec (arr, arr + sizeof(arr) / sizeof(arr[0]) );
    _braitWeights = vec;

    _lifetime = -1; _iteration = 0; _birthdate = 0;

    _wm->updateLandmarkSensor();

    if(_wm->_id == 0)
    {
        std::cout << "W:" << computeRequiredNumberOfWeights() << std::endl;
        std::cout << "(I:" << nn->getNbInputs();
        if (OriginalSharedData::gWithBias)
            std::cout << ", B)";
        else
            std::cout << ")";
        std::cout << ", O: " << nn->getNbOutputs() << std::endl;       
    }
}

OriginalController::~OriginalController()
{
    _genome.clear();
    delete nn;
    nn = NULL;
}

void OriginalController::reset()
{   
}

void OriginalController::resetRobot()
{
    //Number of effectors
    _nbOutputs = 2;
    //Number of sensors
    _nbInputs = 0;
    if (OriginalSharedData::gWithBias)
        _nbInputs += 1;
    //proximity sensors
    _nbInputs += _wm->_cameraSensorsNb;

    //If task=collect, add object sensors
    if ((OriginalSharedData::gFitness == 1) || (OriginalSharedData::gFitness == 2))
    {
        // gathering object distance
        _nbInputs +=  _wm->_cameraSensorsNb;
        // agents
        _nbInputs +=  _wm->_cameraSensorsNb;

        if(_withCollectColorEffector)
        {
            //COLOR DISPLAYED BY ROBOT
            _nbInputs += _wm->_cameraSensorsNb;
            _nbOutputs += 1;
        }
    }

    //Current translational and rotational speeds
    _nbInputs += 2;



    //NN structure
    _nbHiddenLayers = OriginalSharedData::gNbHiddenLayers;
    _nbNeuronsPerHiddenLayer = new std::vector<unsigned int>(_nbHiddenLayers);
    for(unsigned int i = 0; i < _nbHiddenLayers; i++)
        (*_nbNeuronsPerHiddenLayer)[i] = OriginalSharedData::gNbNeuronsPerHiddenLayer;

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
    if(OriginalSharedData::gIsLoadGenome)
        readGenome(OriginalSharedData::gOutGenomeFile + std::to_string(_wm->_id) + ".log");

    _genomesList.clear();
    _fitnessList.clear();

}

void OriginalController::createNN()
{
    delete nn;
    switch ( OriginalSharedData::gControllerType )
    {
        case 0:
            nn = new MLP(_genome, _nbInputs, _nbOutputs,
                         *(_nbNeuronsPerHiddenLayer),false);
            break;
        case 1:
            nn = new Perceptron(_genome, _nbInputs, _nbOutputs);
            break;
        case 2:
            nn = new Elman(_genome, _nbInputs, _nbOutputs, *(_nbNeuronsPerHiddenLayer));
            break;
        default: // default: no controller
            std::cerr << "[ERROR] gController type unknown (value: "
                      << OriginalSharedData::gControllerType << ").\n";
            exit(-1);
    };

}

unsigned int OriginalController::computeRequiredNumberOfWeights()
{
    unsigned int res = nn->getRequiredNumberOfWeights();
    return res;
}

void OriginalController::step()
{
	_iteration++;
    if(!OriginalSharedData::gIsLoadGenome)
        stepEvolution();
    stepBehaviour();
    double distance_sensor;
    double coef_obstacle = 1.0;
    double trans = _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
    double rot = _wm->_desiredRotationalVelocity  / gMaxRotationalSpeed;
    double deltaFitness;
    //Fitness measurement and update
    switch (OriginalSharedData::gFitness) {
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
        //abs(trans) in [0,1], (1 - abs(rot)) in [0,1], coeffObstacle in [0,1]
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
        std::cerr << "[ERROR] Wrong fitness ID" << std::endl;
        exit(-1);
        break;
    }
}


//################BEHAVIOUR METHODS################
void OriginalController::stepBehaviour()
{
    // ---- Build inputs ----
    std::vector<double>* inputs = new std::vector<double>(_nbInputs);
    int inputToUse = 0;
    
    int type = 1; //object type= energy
    for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
    {
        // distance sensors
        //If object but not "Physical object"
        //=> wall (get distance sensor) else 0.0
        int objId = _wm->getObjectIdFromCameraSensor(i);
        if ( PhysicalObject::isInstanceOf(objId) || Agent::isInstanceOf(objId))
        {
            (*inputs)[inputToUse] = 0.0;
            inputToUse++;
        }
        else
        {
            (*inputs)[inputToUse] = 1.0 - _wm->getDistanceValueFromCameraSensor(i) /
                          _wm->getCameraSensorMaximumDistanceValue(i);
            inputToUse++;
        }
        
        //If task=collect, add object sensors
        if ((OriginalSharedData::gFitness == 1) ||(OriginalSharedData::gFitness == 2))
        {
            int objectId = _wm->getObjectIdFromCameraSensor(i);
            // input: physical object?
            //sensing distance to energy item, 0.0 [not 1.0] if not energy item
            if ( PhysicalObject::isInstanceOf(objectId) )
            {
                if ( type == gPhysicalObjects[objectId - gPhysicalObjectIndexStartOffset]->getType() )
                  (*inputs)[inputToUse] = 1.0 - _wm->getDistanceValueFromCameraSensor(i) /
                      _wm->getCameraSensorMaximumDistanceValue(i);
                else
                    (*inputs)[inputToUse] = 0.0;
                inputToUse++;
            }
            else
            {
                //Not a physical object. But: should still fill in the inputs 0.0
                (*inputs)[inputToUse] = 0.0;
                inputToUse++;
            }
       }
        //Agent sensors
       if ( Agent::isInstanceOf(objId) )
       {
           (*inputs)[inputToUse] = 1.0 - _wm->getDistanceValueFromCameraSensor(i) /
               _wm->getCameraSensorMaximumDistanceValue(i);
           inputToUse++;

           if(_withCollectColorEffector)
           {
               OriginalController* c = dynamic_cast<OriginalController*>(
                           (gWorld->getRobot(objId - gRobotIndexStartOffset))->getController());
               (*inputs)[inputToUse] = c->getColorEffector();
               inputToUse++;
           }
       }
       else
       {
           (*inputs)[inputToUse] = 0.0;
           inputToUse++;

          if(_withCollectColorEffector)
          {
              (*inputs)[inputToUse] = 0.0;
              inputToUse++;
          }
       }
    }

    if (OriginalSharedData::gWithBias)
        (*inputs)[inputToUse] = 1.0;
        inputToUse++;

    //Previous translational and rotational speeds (acts as recurrent connections from last step)
    //TODO Should take wheel speed instead
    double lW = _wm->_desiredTranslationalValue / gMaxTranslationalSpeed
            + _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
    double rW = _wm->_desiredTranslationalValue / gMaxTranslationalSpeed
            - _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
    (*inputs)[inputToUse] = lW;
    inputToUse++;
    (*inputs)[inputToUse] = rW;
    inputToUse++;
    /*if(_wm->_id == 1)
    {
        for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
        {
            std::cout << (*inputs)[4*i+3] << " | ";
        }
        //for(auto it = inputs->begin(); it < inputs->end(); it++)
        //    std::cout << (*it) << " | ";
        std::cout << std::endl;
    }*/

    // ---- compute and read out ----
    nn->setWeigths(_genome); // set genome
    bool doBraitenberg = false; //true; //
    if (doBraitenberg)
        nn->setWeigths(_braitWeights);
    nn->setInputs(*inputs);
    nn->step();
    std::vector<double> outputs = nn->readOut();

    //Direct Kinematic model
    //_wm->_desiredTranslationalValue = outputs[0]; _wm->_desiredRotationalVelocity = outputs[1];

    //Differential model
    lW = outputs[0];
    rW = outputs[1];

    if(_withCollectColorEffector)
    {
        //Use a discrete set of color values
        int redValue = roundDown((int)((outputs[2] + 1.0)/2.0 * 256.0),32);
        _wm->setRobotLED_colorValues(redValue, 255 - redValue, 0);
    }

    _wm->_desiredTranslationalValue = (rW + lW) / 2;
    _wm->_desiredRotationalVelocity = (lW - rW) / 2;

    // normalize to motor interval values
    _wm->_desiredTranslationalValue =  _wm->_desiredTranslationalValue * gMaxTranslationalSpeed;
    _wm->_desiredRotationalVelocity = _wm->_desiredRotationalVelocity * gMaxRotationalSpeed;
    
    delete (inputs);
}



//################ EVOLUTION ENGINE METHODS ################
void OriginalController::stepEvolution()
{
    //broadcasting genome : robot broadcasts its genome
    //to all neighbors (contact-based wrt proximity sensors)
    broadcastGenome();
    _lifetime++;

    //agent's lifetime ended: replace genome (if possible)
    if(_lifetime >= OriginalSharedData::gEvaluationTime)
    {
        if (OriginalSharedData::gStoreOwn)
            storeOwnGenome();

        loadNewGenome();

        _currentFitness = 0.0;
        _lifetime = 0;

       if (OriginalSharedData::gClearPopulation)
       {
           _genomesList.clear();
           _fitnessList.clear();
       }
    }

}

void OriginalController::loadNewGenome()
{
   // If 1+ genome(s) imported, select best.
   if (_genomesList.size() > 0)
   {       
       if (_genomesList.size() != _fitnessList.size())
       {
           std::cerr << "[ERROR] Wrong pop size" << std::endl;
           exit(-1);
       }
       //selectBestGenome();
       //selectRankBasedGenome();
       selectTournament(1.0); // argument selection pressure (1.0 = elitist)
   }
   else
   {
        // case: no imported genome: mutate current
   }
   mutate(_currentSigma);
   reset();
   _genomeId.gene_id += 1;

   _birthdate = gWorld->getIterations();
}

void OriginalController::selectBestGenome()
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
        _genome = _genomesList[indexBest];
    }
    else
    {
        //This should never happen
        std::cout << "ERROR in best selection" << std::endl;
        exit(-1);
    }

}
void OriginalController::selectRankBasedGenome()
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

    _genome = _genomesList[pairs[i].first];
}
void OriginalController::selectTournament(double sp){
    /* the size of the tournament */
    int inspected = sp * (double) _genomesList.size();

    /* shuffle indexes */
    std::vector<GC> v;
    for(auto i: _genomesList)
        v.push_back(i.first);
    std::random_shuffle(v.begin(), v.end());

    /* get the best from the inspected */
    double max_fit =  _fitnessList[v[0]];
    GC    best_g  =  (*v.begin());
    int    j=1; /* index in v */
    for (int i=1 ; i<inspected; i++)
    {
        double f  = _fitnessList[v[j]];
        if(f > max_fit)
        {
            max_fit = f;
            best_g = v[j] ;
        }
        j++;
    }

    /*for (int i=0 ; i<inspected ; i++)
        std::cout << _fitnessList[v[i]] << ", ";
    std::cout << std::endl;*/

    _genome = _genomesList[best_g];
}
void OriginalController::mutate(float sigma) // mutate within bounds.
{
    std::vector<double> g;
    g.clear();
    
	_currentSigma = sigma;
    for (unsigned int i = 0 ; i != _genome.size() ; i++ )
	{
        double value = _genome[i] + getGaussianRand(0,_currentSigma);
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
        
        g.push_back(value);
	}
    
    _genome = g;
    
}

void OriginalController::updateFitness(double delta)
{
    _currentFitness += delta;
}

// ################ COMMUNICATION METHODS ################
void OriginalController::broadcastGenome()
{

    for( int i = 0 ; i < _wm->_cameraSensorsNb; i++)
    {
        int targetIndex = _wm->getObjectIdFromCameraSensor(i);

        // sensor ray bumped into a robot : communication is possible
        if ( targetIndex >= gRobotIndexStartOffset )
        {
            // convert image registering index into robot id.
            targetIndex = targetIndex - gRobotIndexStartOffset;

            OriginalController* targetRobotController =
                    dynamic_cast<OriginalController*>
                    (gWorld->getRobot(targetIndex)->getController());

            if ( ! targetRobotController )
            {
                std::cerr << "Error: observer not compatible" << std::endl;
                exit(-1);
            }

            // other agent stores my genome.
            //Genome id as gene clock (different for each agent)
            targetRobotController->storeGenome(_genome, _genomeId, _currentFitness);
        }
    }
}

void OriginalController::storeGenome(std::vector<double> genome, GC senderId, double fitness)
{
    if(_genomesList.find(senderId) != _genomesList.end())
    {
        //Update fitness
        _fitnessList[senderId] = fitness;
    }
    else
    {
        //Local population size ignored
        if (_genomesList.size() < (unsigned int)OriginalSharedData::gPopulationSize)
        {
            _genomesList[senderId] = genome;
            _fitnessList[senderId] = fitness;

        }
        else
        {
            //if medea this should not happen
            //std::cerr << "[ERROR] If medea this should not happen" << std::endl;
            //exit(-1);
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

void OriginalController::storeOwnGenome()
{
    storeGenome(_genome, _genomeId, _currentFitness);
}

void OriginalController::logGenome(std::string s)
{
    std::ofstream genomeF;
    genomeF.open(s);

    genomeF << "[F " << _currentFitness << "\n";
    //Structure (NN)
    genomeF << "[I " << _nbInputs << "\n";//Inputs
    if(OriginalSharedData::gWithBias)
        genomeF << "[B " << 1 << "\n";//Bias
    else
        genomeF << "[B " << 0 << "\n";
    genomeF << "[H " << _nbHiddenLayers << "\n"
            << "[N " << OriginalSharedData::gNbNeuronsPerHiddenLayer << "\n";//Neur per hidden layer
    genomeF << "[O " << _nbOutputs << "\n";


    genomeF << "[W ";
    genomeF << std::fixed << std::setprecision(10);

    for(auto it = _genome.begin(); it != _genome.end(); ++it)
    {
        genomeF << (*it) << " ";
    }


    genomeF.close();
}

void OriginalController::readGenome(std::string s)
{
    std::ifstream genomeF;
    std::string line;

    genomeF.open(s);
    std::vector<std::vector<std::string>> allTokens;
    while(getline(genomeF, line))
    {
        //Extract all separate tokens in the line
        std::vector<std::string> tokens;
        boost::split(tokens,line, boost::is_any_of("\t "));
        allTokens.push_back(tokens);
    }
    genomeF.close();
    double fitness = 0.0;
    std::vector<double> weights;
    int inN = 0,hLayers = 0,hN = 0,oN = 0;
    int b = 0;

    for(unsigned int i=0;i< allTokens.size();i++)
    {
        std::vector<std::string> tokens = allTokens[i];
        if(tokens[0] == "[F")
            fitness =  stod(tokens[1]);
        else if(tokens[0] == "[I")
            inN =  stoi(tokens[1]);
        else if(tokens[0] == "[B")
            b =  stoi(tokens[1]);
        else if(tokens[0] == "[H")
            hLayers =  stoi(tokens[1]);
        else if(tokens[0] == "[N")
            hN =  stoi(tokens[1]);
        else if(tokens[0] == "[O")
            oN =  stoi(tokens[1]);
        else if(tokens[0] == "[W")
        {
            for(unsigned int j=1; j < tokens.size();j++)
            {
               boost::trim(tokens[j]);
               if(tokens[j] != "")
                weights.push_back(stod(tokens[j]));
            }
        }
    }
    //NN structure
    _nbHiddenLayers = hLayers;
    _nbNeuronsPerHiddenLayer = new std::vector<unsigned int>(_nbHiddenLayers);
    for(unsigned int i = 0; i < _nbHiddenLayers; i++)
        (*_nbNeuronsPerHiddenLayer)[i] = hN;

    _nbInputs = inN; //bias already included
    _nbOutputs = oN;

    createNN();
    _genome = weights;
    std::cout << fitness << ", " << b<< std::endl;
}
int OriginalController::roundDown(int numToRound, int multiple)
{
    int result;

    if (multiple == 0)
        return numToRound;

    int remainder = abs(numToRound) % multiple;
    if (remainder == 0)
    {
        result = numToRound;
    }
    else
    {
        if (numToRound < 0)
            result = -(abs(numToRound) - remainder + multiple);
        else
            result = numToRound - remainder ;
    }
    //std::cout << numToRound << " | " << multiple << " | "<< result << std::endl;
    return result;
}