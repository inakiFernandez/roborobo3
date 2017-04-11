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
    _wm = wm; nn = NULL; nnF = NULL;

    // neural weights limits
    _minValue = -OriginalSharedData::gNeuronWeightRange/2;
    _maxValue = OriginalSharedData::gNeuronWeightRange/2;
	_currentSigma = OriginalSharedData::gSigmaRef;

    _doEvoTopo = (OriginalSharedData::gControllerType == 3);

    //Evo variables
    _currentFitness = 0.0;
    _genomeId.robot_id = _wm->_id;
    _genomeId.gene_id = 0;

    //Set initial color sensor to 0.0
    _wm->setRobotLED_colorValues(128, 128, 0);

    _withCollectColorEffector = OriginalSharedData::gWithCollectColorEffector;
    _previousColor = 0;
    _nbColorChanges = 0;
    resetRobot();


    _braitWeights = getBraitenberg();

    _lifetime = -1; _iteration = 0; _birthdate = 0;

    _wm->updateLandmarkSensor();

    //Stored previous controller. Initialized to NO PREVIOUS
    _storedF = std::vector<double>();


    if((_wm->_id == 0) && !_doEvoTopo)
    {
        std::cout << "W:" << computeRequiredNumberOfWeights() << std::endl;
        std::cout << "(I:" << nnF->getNbInputs();
        if (OriginalSharedData::gWithBias)
            std::cout << ", B)";
        else
            std::cout << ")";
        std::cout << ", O: " << nnF->getNbOutputs() << std::endl;
    }
}

OriginalController::~OriginalController()
{
    _genomeF.clear();
    if(nnF != NULL)
    {
        delete nnF;
        nnF = NULL;
    }
    if(nn != NULL)
    {
        delete nn;
        nn = NULL;
    }
    delete _genome;    
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
    if ((OriginalSharedData::gFitness == 1)
            || (OriginalSharedData::gFitness == 2))
    {
        // gathering object distance
        _nbInputs +=  _wm->_cameraSensorsNb;
        // agents
        _nbInputs +=  _wm->_cameraSensorsNb;

        if(_withCollectColorEffector)
        {
            //COLOR DISPLAYED BY OTHER ROBOT
            //OR
            //COLOR OF THE DETECTED ITEM
            //cf. stepBehavior()
            _nbInputs += _wm->_cameraSensorsNb;
            _nbOutputs += 1;
        }
    }

    //Current translational and rotational speeds
    _nbInputs += 2;


    if(!_doEvoTopo)
    {
        //NN structure
        _nbHiddenLayers = OriginalSharedData::gNbHiddenLayers;
        _nbNeuronsPerHiddenLayer = new std::vector<unsigned int>(_nbHiddenLayers);
        for(unsigned int i = 0; i < _nbHiddenLayers; i++)
            (*_nbNeuronsPerHiddenLayer)[i] = OriginalSharedData::gNbNeuronsPerHiddenLayer;

        createNN();
        unsigned int const nbGene = computeRequiredNumberOfWeights();

        _genomeF.clear();
        double w;
        //Random genome (weights) initialization
        for ( unsigned int i = 0 ; i != nbGene ; i++ )
        {
            // weights: random init between -1 and +1
            w = (double)(rand() % 10000)/5000.0 - 1.0;
            _genomeF.push_back(w);
        }
        if(OriginalSharedData::gIsLoadGenome)
            readGenomeF(OriginalSharedData::gOutGenomeFile + std::to_string(_wm->_id) + ".log");
        _genomesFList.clear();
    }
    else
    {
        // Inputs, outputs, 0 hidden neurons, fully connected.
        //Initial Genes=>common to all agents, thus identified by a common historical marker
        _genome = new Genome (_genomeId,_nbInputs, _nbOutputs);
        //Initialize weights in [-1,1]
        _genome->initialize_link_weights();
        //Fully connected
        _g_count = _nbInputs * _nbOutputs + 1; //next Gene Clock
        _n_count = _nbInputs + _nbOutputs + 1; //next Neuron Gene Clock
        _genomesList.clear();
        createNN();
    }
    _fitnessList.clear();
}

void OriginalController::createNN()
{
    if(nn != NULL)
        delete nn;
    if(nnF != NULL)
        delete nnF;
    switch ( OriginalSharedData::gControllerType )
    {
        case 0:
            nnF = new MLP(_genomeF, _nbInputs, _nbOutputs,
                         *(_nbNeuronsPerHiddenLayer),false);
            break;
        case 1:
            nnF = new Perceptron(_genomeF, _nbInputs, _nbOutputs);
            break;
        case 2:
            nnF = new Elman(_genomeF, _nbInputs, _nbOutputs, *(_nbNeuronsPerHiddenLayer));
            break;
        case 3:
            nn = _genome->genesis();
            break;
        default: // default: no controller
            std::cerr << "[ERROR] gController type unknown (value: "
                      << OriginalSharedData::gControllerType << ").\n";
            exit(-1);
    };

}

unsigned int OriginalController::computeRequiredNumberOfWeights()
{
    unsigned int res = nnF -> getRequiredNumberOfWeights();
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
    switch(OriginalSharedData::gFitness)
    {
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
    //If number color changes measured on a per timestep basis
    _nbColorChanges = 0;

    // ---- Build inputs ----
    std::vector<double>* inputs = new std::vector<double>(_nbInputs);
    int inputToUse = 0;
    
    int type = 1; //object type= energy
    for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
    {
        // distance sensors
        //If object but not "Physical object"
        //=> wall (get proximity sensor) else 0.0
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
                if(_withCollectColorEffector)
                {                    
                    /* Sensing the color of the item*/
                    double colorObj = gPhysicalObjects[objectId - gPhysicalObjectIndexStartOffset]->getColorValue();
                    (*inputs)[inputToUse] = colorObj;
                    inputToUse++;
                }
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
               /* TORESET ALSO IN THE ELSE PART
                * Sensing the color displayed by other robot*/
               /*OriginalController* c = dynamic_cast<OriginalController*>(
                           (gWorld->getRobot(objId - gRobotIndexStartOffset))->getController());
               (*inputs)[inputToUse] = c->getColorEffector();
               inputToUse++;*/
           }
       }
       else
       {
           (*inputs)[inputToUse] = 0.0;
           inputToUse++;

          if(_withCollectColorEffector)
          {
              /*(*inputs)[inputToUse] = 0.0;
              inputToUse++;*/
          }
       }
    }

    if (OriginalSharedData::gWithBias)
    {
        (*inputs)[inputToUse] = 1.0;
        inputToUse++;
    }

    //Previous left and right wheel speeds (acts as recurrent connections from last step)
    double lW = _wm->_desiredTranslationalValue / gMaxTranslationalSpeed
            + _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
    double rW = _wm->_desiredTranslationalValue / gMaxTranslationalSpeed
            - _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
    (*inputs)[inputToUse] = lW;
    inputToUse++;
    (*inputs)[inputToUse] = rW;
    inputToUse++;
    std::vector<double> outputs;
    // compute and read out
    if(!_doEvoTopo)
    {
        nnF->setWeigths(_genomeF); // set genome
        bool doBraitenberg = (OriginalSharedData::gBrait != 0)
                && (OriginalSharedData::gBrait != 4)
                && (OriginalSharedData::gBrait != 5)
                && (OriginalSharedData::gBrait != 6);
        if (doBraitenberg)
            nnF->setWeigths(_braitWeights);
        nnF->setInputs(*inputs);
        nnF->step();
        outputs = nnF->readOut();
    }
    else
    {
        nn->load_sensors (&((*inputs)[0]));
        if (!(nn->activate ()))
        {
            std::cerr << "[ERROR] Activation of ANN not correct: genome R"
                      << _genome->genome_id.robot_id << ", G" << _genome->genome_id.gene_id << std::endl;
            //save_genome();
            exit (-1);
        }

        for (auto out_iter  = nn->outputs.begin();
             out_iter != nn->outputs.end();
             out_iter++)
            outputs.push_back((*out_iter)->activation);
    }

    //Direct Kinematic model
    //_wm->_desiredTranslationalValue = outputs[0];
    //_wm->_desiredRotationalVelocity = outputs[1];

    //Differential model
    lW = outputs[0];
    rW = outputs[1];

    if(_withCollectColorEffector)
    {
        //Random color effector (neural output ignored). Control experiment
        //outputs[2] = (double)(rand() % 10000)/5000.0 - 1.0;


        //Use a discrete set of color values
        int redValue = roundDown((int)((outputs[2] + 1.0)/2.0 * 256.0),32);
        _wm->setRobotLED_colorValues(redValue, 255 - redValue, 0);
    }

    _wm->_desiredTranslationalValue = (rW + lW) / 2;
    _wm->_desiredRotationalVelocity = (lW - rW) / 2;

    //_wm->_desiredTranslationalValue = 0;
    //_wm->_desiredRotationalVelocity = 0;

    // normalize to motor interval values
    _wm->_desiredTranslationalValue =  _wm->_desiredTranslationalValue * gMaxTranslationalSpeed;
    _wm->_desiredRotationalVelocity = _wm->_desiredRotationalVelocity * gMaxRotationalSpeed;

    //Logging number of color changes
    int currentColor = roundDown((int)((outputs[2] + 1.0)/2.0 * 256.0),32);
    if(_previousColor != currentColor)
        _nbColorChanges++;
    _previousColor = currentColor;

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
        {
            if(!_doEvoTopo)
                storeOwnGenomeF();
            else
                storeOwnGenome();
        }

        loadNewGenome();

        _currentFitness = 0.0;
        _nbColorChanges = 0;
        _lifetime = 0;

       if (OriginalSharedData::gClearPopulation)
       {
           _genomesList.clear();
           _genomesFList.clear();
           _fitnessList.clear();
       }
    }

}

void OriginalController::loadNewGenome()
{
   bool popNotEmpty;
   if(!_doEvoTopo)
       popNotEmpty = (_genomesFList.size() > 0);
   else
       popNotEmpty = (_genomesList.size() > 0);
   // If 1+ genome(s) imported, select best.
   if (popNotEmpty)
   {       
       if ((_genomesList.size() != _fitnessList.size()) && (_genomesFList.size() != _fitnessList.size()))
       {
           std::cerr << "[ERROR] Wrong pop size" << std::endl;
           exit(-1);
       }
       //selectBestGenome();
       //selectRankBasedGenome();
       selectTournament(OriginalSharedData::gSelPressure); // argument selection pressure (1.0 = elitist)
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
        if(!_doEvoTopo)
            _genomeF = _genomesFList[indexBest];
        else
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
    unsigned int n = (_doEvoTopo?_genomesList.size():_genomesFList.size());
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

    if(!_doEvoTopo)
        _genomeF = _genomesFList[pairs[i].first];
    else
        _genome = _genomesList[pairs[i].first];
}
void OriginalController::selectTournament(double sp){
    /* the size of the tournament */
    int size = (_doEvoTopo?_genomesList.size():_genomesFList.size());
    int inspected = sp * (double) size;

    /* shuffle indexes */
    std::vector<GC> v;
    for(auto i: _fitnessList)
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
    if(!_doEvoTopo)
        _genomeF = _genomesFList[best_g];
    else
        _genome = _genomesList[best_g];
}
void OriginalController::mutate(float sigma) // mutate within bounds.
{
    if(!_doEvoTopo)
        mutateFixed(sigma);
    else
        mutateEvoTopo(sigma);
}

void OriginalController::mutateFixed(float s)
{
    std::vector<double> g;
    g.clear();
    
    _currentSigma = s;
    for (unsigned int i = 0 ; i != _genomeF.size() ; i++ )
	{
        double value = _genomeF[i] + getGaussianRand(0,_currentSigma);
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
    
    _genomeF = g;
    
}
void OriginalController::mutateEvoTopo(float s)
{
    /*if((Helper::randFloat() < Helper::mateOnlyProb)
            && !(g1->genome_id == g2->genome_id))
    {
        result = g1 -> mate(g2,_genomeId, f1,f2);
    }
    else
    {
        result = g1;
    }*/
    _genomeId.gene_id++;
    if(Helper::randFloat() < Helper::mutateProb)//Mutate
    {
        _genome = _genome-> mutate(s, _wm->_id,_genomeId,_n_count,_g_count);
    }
    _genome->genome_id = _genomeId;
    if(_doEvoTopo)
        createNN();

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
            if(!_doEvoTopo)
            {
                targetRobotController->storeGenomeF(_genomeF, _genomeId, _currentFitness);
            }
            else
                targetRobotController->storeGenome(_genome, _genomeId, _currentFitness, _g_count, _n_count);
        }
    }
}

void OriginalController::storeGenomeF(std::vector<double> genome, GC senderId, double fitness)
{
    if(_genomesFList.find(senderId) != _genomesFList.end())
    {
        //Update fitness
        _fitnessList[senderId] = fitness;
    }
    else
    {
        //Local population size ignored
        if (_genomesFList.size() < (unsigned int)OriginalSharedData::gPopulationSize)
        {
            _genomesFList[senderId] = genome;
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
                _genomesFList.erase(indexWorse);

                _fitnessList[senderId] = fitness;
                _genomesFList[senderId] = genome;
            }
        }
    }
}
void OriginalController::storeGenome(Genome* genome, GC senderId, double fitness, int g_count, int n_count)
{
    if(OriginalSharedData::gUpdateGC)
    {
        //Update gene clocks for nodes and links:
        //minimize the number of arbitrary sorting orders in genome alignment
        //due to concurrent mutations in different agents
        _n_count = std::max(_n_count,n_count);
        _g_count = std::max(_g_count,g_count);
    }
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
    storeGenome(_genome, _genomeId, _currentFitness, _g_count, _n_count);
}
void OriginalController::storeOwnGenomeF()
{
    storeGenomeF(_genomeF, _genomeId, _currentFitness);
}

void OriginalController::logGenomeF(std::string s)
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

    for(auto it = _genomeF.begin(); it != _genomeF.end(); ++it)
    {
        genomeF << (*it) << " ";
    }


    genomeF.close();
}

void OriginalController::readGenomeF(std::string s)
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
    _genomeF = weights;
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

//Measures a distance between previous stored neural controller and current one
double OriginalController::forget()
{
    //TODO for evolving topo
    double result = 0.0;
    switch(OriginalSharedData::gForgetMethod)
    {
        //Structural difference
        case 1:
            result = forgetStructural();
            break;
        //Behavioral distance (?)
        case 2:
            //TODO
            break;
    }
    return result;
}

double OriginalController::forgetStructural()
{
    //euclidean distance of weight vector (fixed-topology)//TODO evo topo
    double result = 0.0;
    for(unsigned int i = 0; i <_genomeF.size(); i++)
    {
        result += (_genomeF[i] - _storedF[i]) * (_genomeF[i] - _storedF[i]);
    }
    return sqrt(result);
}

std::vector<double> OriginalController::getBraitenberg()
{
    //braitenberg genome (8sensors no color)
    /*static const double arr[] = {0.5, -0.5, -1.0, 1.0, 0.75, -0.25, //W
                                 0.5, -0.5, -1.0, 1.0, 0.75, -0.25, //NW
                                 -0.5, -0.5, 1.0, 1.0, -0.5, -0.5, //N
                                 -0.5, 0.5, 1.0,-1.0, -0.25, 0.75, //NE
                                 -0.5, 0.5, 1.0,-1.0, -0.25, 0.75, //E
                                 -0.5, 0.5, 1.0,-1.0, -0.25, 0.75, //SE
                                 0.5, 0.5, -1.0, -1.0, 0.5, 0.5, //S
                                 0.5, -0.5, -1.0, 1.0, 0.75, -0.25, //SW
                                 0.5, 0.5, //bias
                                 0.0, 0.0, //lW recurrent
                                 0.0, 0.0, //rW recurrent
                                };*/

    //braitenberg genome Colors intermediate (8sensors x (obst,item,robot,colorrobot) color=T1, always red)
    std::vector<double> arr1 = {0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, //0.0, 0.0, 0.0, //W
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, //0.0, 0.0, 0.0, //NW
                                   -0.5, -0.5, 0.0, 1.0, 1.0, 0.0, -0.5, -0.5, 0.0, //0.0, 0.0, 0.0, //N
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, //0.0, 0.0, 0.0, //NE
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, //0.0, 0.0, 0.0, //E
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, //0.0, 0.0, 0.0, //SE
                                   0.5, 0.5, 0.0, -1.0, -1.0, 0.0,  0.5, 0.5, 0.0, //0.0, 0.0, 0.0, //S
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, //0.0, 0.0, 0.0, //SW
                                   0.5, 0.5, 0.5, //bias
                                   0.0, 0.0, 0.0, //lW recurrent
                                   0.0, 0.0, 0.0 //rW recurrent
                                };

    //braitenberg genome Colors intermediate (8sensors x (obst,item,robot,colorrobot) color=T2, always green)
    std::vector<double> arr2 = {0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, //0.0, 0.0, 0.0, //W
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, //0.0, 0.0, 0.0, //NW
                                   -0.5, -0.5, 0.0, 1.0, 1.0, 0.0, -0.5, -0.5, 0.0, //0.0, 0.0, 0.0, //N
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, //0.0, 0.0, 0.0, //NE
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, //0.0, 0.0, 0.0, //E
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, //0.0, 0.0, 0.0, //SE
                                   0.5, 0.5, 0.0, -1.0, -1.0, 0.0,  0.5, 0.5, 0.0, //0.0, 0.0, 0.0, //S
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, //0.0, 0.0, 0.0, //SW
                                   0.5, 0.5, -0.5, //bias
                                   0.0, 0.0, 0.0, //lW recurrent
                                   0.0, 0.0, 0.0 //rW recurrent
                                };

    //braitenberg genome Colors intermediate (8sensors x (obst,item,robot,colorrobot) color=T1, red if alone, green if item + robot)
    std::vector<double>  arr3 = {0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, //0.0, 0.0, 0.0, //W
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, -1.0, //0.0, 0.0, 0.0, //NW
                                   -0.5, -0.5, 0.0, 1.0, 1.0, 0.0, -0.5, -0.5, 0.0, //0.0, 0.0, 0.0, //N
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, -1.0, //0.0, 0.0, 0.0, //NE
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, //0.0, 0.0, 0.0, //E
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, //0.0, 0.0, 0.0, //SE
                                   0.5, 0.5, 0.0, -1.0, -1.0, 0.0,  0.5, 0.5, 0.0, //0.0, 0.0, 0.0, //S
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, //0.0, 0.0, 0.0, //SW
                                   0.5, 0.5, 0.5, //bias
                                   0.0, 0.0, 0.0, //lW recurrent
                                   0.0, 0.0, 0.0 //rW recurrent
                                };

    /*//braitenberg genome (8sensors x (obst,item,robot,colorrobot) color=T1, always red)
    std::vector<double> arr1 = {0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, 0.0, 0.0, 0.0, //W
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, 0.0, 0.0, 0.0, //NW
                                   -0.5, -0.5, 0.0, 1.0, 1.0, 0.0, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, //N
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, 0.0, 0.0, 0.0, //NE
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, 0.0, 0.0, 0.0, //E
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, 0.0, 0.0, 0.0, //SE
                                   0.5, 0.5, 0.0, -1.0, -1.0, 0.0,  0.5, 0.5, 0.0, 0.0, 0.0, 0.0, //S
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, 0.0, 0.0, 0.0, //SW
                                   0.5, 0.5, 5.0, //bias
                                   0.0, 0.0, 0.0, //lW recurrent
                                   0.0, 0.0, 0.0 //rW recurrent
                                };

    //braitenberg genome (8sensors x (obst,item,robot,colorrobot) color=T2, always green)
    std::vector<double> arr2 = {0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, 0.0, 0.0, 0.0, //W
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, 0.0, 0.0, 0.0, //NW
                                   -0.5, -0.5, 0.0, 1.0, 1.0, 0.0, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, //N
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, 0.0, 0.0, 0.0, //NE
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, 0.0, 0.0, 0.0, //E
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, 0.0, 0.0, 0.0, //SE
                                   0.5, 0.5, 0.0, -1.0, -1.0, 0.0,  0.5, 0.5, 0.0, 0.0, 0.0, 0.0, //S
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, 0.0, 0.0, 0.0, //SW
                                   0.5, 0.5, -5.0, //bias
                                   0.0, 0.0, 0.0, //lW recurrent
                                   0.0, 0.0, 0.0 //rW recurrent
                                };

    //braitenberg genome (8sensors x (obst,item,robot,colorrobot) color=T1, red if alone, green if item + robot)
    std::vector<double>  arr3 = {0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, 0.0, 0.0, 0.0, //W
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, -10.0, 0.0, 0.0, 0.0, //NW
                                   -0.5, -0.5, 0.0, 1.0, 1.0, 0.0, -0.5, -0.5, 0.0, 0.0, 0.0, 0.0, //N
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, -10.0, 0.0, 0.0, 0.0, //NE
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, 0.0, 0.0, 0.0, //E
                                   -0.5, 0.5, 0.0, 1.0,-1.0, 0.0, -0.25, 0.75, 0.0, 0.0, 0.0, 0.0, //SE
                                   0.5, 0.5, 0.0, -1.0, -1.0, 0.0,  0.5, 0.5, 0.0, 0.0, 0.0, 0.0, //S
                                   0.5, -0.5, 0.0, -1.0, 1.0, 0.0, 0.75, -0.25, 0.0, 0.0, 0.0, 0.0, //SW
                                   0.5, 0.5, 5.0, //bias
                                   0.0, 0.0, 0.0, //lW recurrent
                                   0.0, 0.0, 0.0 //rW recurrent
                                };*/
    std::vector<double> vec ;

    switch (OriginalSharedData::gBrait)
    {
        case 1:
            vec = arr1;
            break;
        case 2:
            vec = arr2;
            break;
        case 3:
            vec = arr3;
            break;
        case 4:
            //Initialize population to mutation of B1
            _genomeF = arr1;
            mutate(_currentSigma);
        case 5:
            //Initialize population to mutation of B2
            _genomeF = arr2;
            mutate(_currentSigma);
            break;
        case 6:
            //Initialize population to mutation of B3
            _genomeF = arr3;
            mutate(_currentSigma);
            break;
    }

    return vec;
}
