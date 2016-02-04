/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "TestNN/include/TestNNController.h"
#include "TestNN/include/TestNNWorldObserver.h"

#include "World/World.h"
#include "Utilities/Misc.h"
#include <math.h>
#include <string>

#include <neuralnetworks/MLP.h>
#include <neuralnetworks/Perceptron.h>
#include <neuralnetworks/Elman.h>

using namespace Neural;

TestNNController::TestNNController( RobotWorldModel *wm )
{
    _wm = wm; nn = NULL;

    _minValue = -3.0;
    _maxValue = 3.0;
    _currentFitness = 0.0;
    _genomeId.robot_id = _wm->_id;
    _genomeId.gene_id = 0;
    loadTestGenome(TestNNSharedData::gOutGenomeFile + std::to_string(_wm->_id) + ".log");

    resetRobot();
    _lifetime = -1;
    // behaviour
    _iteration = 0; _birthdate = 0;
    _wm->updateLandmarkSensor();
    _wm->setRobotLED_colorValues(255, 0, 0);
}

TestNNController::~TestNNController()
{
    _parameters.clear();
    delete nn;
    nn = NULL;
}

void TestNNController::reset()
{   
    _parameters.clear();
    _parameters = _genome;
}

void TestNNController::resetRobot()
{    
    createNN();
    _currentGenome = _genome;
    reset();
}

void TestNNController::createNN()
{
    delete nn;
    // MLP
    nn = new MLP(_parameters, _nbInputs, _nbOutputs, *(_nbNeuronsPerHiddenLayer),
                 TestNNSharedData::gWithBias);
}

unsigned int TestNNController::computeRequiredNumberOfWeights()
{
    unsigned int res = nn->getRequiredNumberOfWeights();
    return res;
}

void TestNNController::step()
{
	_iteration++;
    _lifetime++;
    if(_lifetime >= TestNNSharedData::gEvaluationTime)
    {
        _currentFitness = 0.0;
        _lifetime = 0;
    }
    stepBehaviour();

    double distance_sensor;
    double coef_obstacle = 1.0;

    double trans = _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
    double rot = _wm->_desiredRotationalVelocity  / gMaxRotationalSpeed;
    double deltaFitness;
    //Fitness measurement and update
    switch (TestNNSharedData::gFitness) {
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
        //TestNNally averaging deltas
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
void TestNNController::stepBehaviour()
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
        if (TestNNSharedData::gFitness == 1)
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
void TestNNController::updateFitness(double delta)
{
    _currentFitness += delta;
}
void TestNNController::loadTestGenome(std::string filename)
{
    std::ifstream infile;
    infile.open(filename, std::ios::in);
    std::string s;
    std::vector<std::string> vString;
    vString.clear();
    std::vector<double> weights;
    weights.clear();
    double d, f = 0.0;
    int in = -1, out = -1, hLayers = -1, nHidLayer = -1, bias = -1;

    while (std::getline(infile, s))
    {
        vString.push_back(s);
    }
    infile.close();

    for (auto it = vString.begin(); it != vString.end(); ++it)
    {
        s = (*it);
        if (s[0] == '[')
        {
            switch(s[1])
            {
                case 'W':
                    ++it; s = (*it);
                    while(s[0] != '[')
                    {
                        d = std::stod(s);
                        weights.push_back(d);
                        ++it;
                        s = (*it);
                    }
                    --it;
                    break;
                case 'F':
                    ++it; s = (*it);
                   f = std::stod(s);
                   break;
                case 'I':
                   ++it; s = (*it);
                   in = std::stoi(s);
                   break;
                case 'O':
                   ++it; s = (*it);
                   out = std::stoi(s);
                   break;
                case 'H':
                   ++it; s = (*it);
                   hLayers = std::stoi(s);
                   break;
                case 'N':
                   ++it; s = (*it);
                   nHidLayer = std::stoi(s);
                   break;
                case 'B':
                    ++it; s = (*it);
                    bias = std::stoi(s);
                    break;
                default:
                    std::cout << "Wrong line in file?" << std::endl;
            }
       }
       else
            std::cout << "Reading parameter went wrong" << std::endl;
    }
    _nbInputs = in;
    _nbOutputs = out;
    _nbHiddenLayers = hLayers;
    if (bias == 1)
        TestNNSharedData::gWithBias = true;
    else
        TestNNSharedData::gWithBias = false;
    _nbNeuronsPerHiddenLayer = new std::vector<unsigned int>(_nbHiddenLayers);

    for (unsigned int i = 0; i < _nbHiddenLayers; i++)
         (*_nbNeuronsPerHiddenLayer)[i] = nHidLayer;
    _genome = weights;
    _recordedFitness = f;
    std::cout << "Recorded fitness: " << _recordedFitness << std::endl;
}
