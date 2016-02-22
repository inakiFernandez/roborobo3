/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "odNeat/include/odNeatController.h"
#include "odNeat/include/odNeatWorldObserver.h"

#include "odneatgc/helper.h"

#include "World/World.h"
#include "Utilities/Misc.h"
#include <math.h>
#include <string>
#include <set>

#include <neuralnetworks/MLP.h>
#include <neuralnetworks/Perceptron.h>
#include <neuralnetworks/Elman.h>

using namespace Neural;

odNeatController::odNeatController( RobotWorldModel *wm )
{
    _wm = wm; nn = NULL;
    // evolutionary engine
    _minValue = -3.0;
    _maxValue = 3.0;

    _genomeId.robot_id = _wm->_id;
    _genomeId.gene_id = 0;

    resetRobot();
    _genome->genome_id = _genomeId;
    _lifetime = -1;
    // behaviour
    _iteration = 0; _birthdate = 0;
    _wm->updateLandmarkSensor();
    _wm->setRobotLED_colorValues(255, 0, 0);
    _newSpId = 0;

    add_to_population(_genome->duplicate(),_energy);
}

odNeatController::~odNeatController()
{
    //TODO invoke delete on each object
    _pop.clear();
    _tabu.clear();
    delete nn;
    nn = NULL;
    delete _genome;
}

void odNeatController::reset()
{

    _fitnessUpdateCounter = 0;
    _genome -> nbFitnessUpdates = 0;
    _energy = odNeatSharedData::gDefaultInitialEnergy;
    _currentFitness = _energy;
    if (odNeatSharedData::gFitness == 1)
    {
        _items = 0;
    }
    _genome->species = -1;
    _genome->nbFitnessUpdates ++;
    _birthdate = gWorld->getIterations();
    _lifetime = -1;

    /*std::cout << _genome->genes.size() << std::endl;
    for(int i = 0; i < 1; i++)
    {
        add_to_population(_genome->duplicate(), _energy + i);
        _genomeId.gene_id++;
        createNN();
        _genome = _genome->mutate (odNeatSharedData::gSigmaRef,
                                   _wm->_id,_genomeId,_n_count,_g_count);
        _genome->genome_id = _genomeId;
        Helper::mutateAddNodeProb = 0.0;
        //_genome->print_to_file(std::cout);

        std::cout << _genome->genes.size() << std::endl;
        //std::cout<< "------------------------------" << std::endl;
    }
    for(int i = 0; i < 50; i++)
    {
        _genome->mutate_link_weights(odNeatSharedData::gSigmaRef);
        _genome -> print_to_file(std::cout);
    }
    exit(-1);*/

    createNN();

    //Add own genome (even if fitness is not yet measured)
    add_to_population(_genome->duplicate(),_energy);
}

void odNeatController::resetRobot()
{
    //Number of sensors
    _nbInputs = 0;
    //proximity sensors
    _nbInputs += _wm->_cameraSensorsNb;
    //If task=collect, add object sensors
    if (odNeatSharedData::gFitness == 1)
    {
        // gathering object distance
        _nbInputs +=  _wm->_cameraSensorsNb;
    }
    //Energy
    _nbInputs += 1;
    //Bias
    _nbInputs += 1;

    //Current translational and rotational speeds
    //_nbInputs += 2;

    //Number of effectors
    _nbOutputs = 2;

    // Inputs, outputs, 0 hidden neurons, fully connected.
    //Initial Genes=>common to all agents, thus identified by a common historical marker    
    _genome = new Genome (_genomeId,_nbInputs, _nbOutputs);
    //Weights at 0.0, so mutate
    _genome->mutate_link_weights(odNeatSharedData::gSigmaRef);
    //Fully connected
    _g_count = _nbInputs * _nbOutputs + 1;
    _n_count = _nbInputs + _nbOutputs + 1;
}

void odNeatController::createNN()
{
    if (nn != NULL)
        delete nn;
    nn = _genome->genesis();
}

void odNeatController::step()
{
    _iteration++;
    _lifetime++;
    //If Inter-robot reproduction event
    //robot broadcasts genome to all neighbors
    if(doBroadcast()) broadcastGenome();

    stepBehaviour ();

    _energy = update_energy_level(); updateFitness();    

    _wm->setRobotLED_colorValues((int) _energy / odNeatSharedData::gMaxEnergy, 0, 0);

    if ((_energy <=odNeatSharedData::gEnergyThreshold)
            && !(in_maturation_period()))
    {
        //TODO LOGS
        stepEvolution();
        reset();
    }
}

void odNeatController::pickItem()
{
    _energy +=odNeatSharedData::gEnergyItemValue;
    _items++;
}
int odNeatController::getItems()
{
    return _items;
}

//################ EVOLUTION ENGINE METHODS ################
void odNeatController::stepEvolution()
{
    _genome->nbFitnessUpdates ++;
    //Update own genome's energy in population with final estimate
    add_to_population(_genome->duplicate(), _currentFitness);
    add_to_tabu_list(_genome->duplicate());

    Genome* offspring =  generate_offspring();
    delete _genome;
    _genome = offspring; //genesis() and add_to_pop are called later, on reset()
}


Genome* odNeatController::generate_offspring()
{
    _genomeId.gene_id++;
    Genome* result;
    //Mutate already calls duplicate
    odNeatSpecies* sp = selectSpecies();
    Genome* g1 = selectParent(sp); Genome* g2 = selectParent(sp);

    //TOTEST Selection on species

    /*//Mate with probability, only if selected genomes are not the same
    if((Helper::randFloat() < Helper::mateOnlyProb) && !(g1->genome_id == g2->genome_id))
        result = g1->mate(g2,_genomeId, fitnessG1,fitnessG2);
    else
    {*/
    result = g1;

    /*}*/
    if(Helper::randFloat() < Helper::mutateProb)//Mutate
    {
        result = result-> mutate(odNeatSharedData::gSigmaRef,
                                  _wm->_id,_genomeId,_n_count,_g_count);
    }
    if((result->mom_id.gene_id != -1)
            && (result->genome_id.gene_id != -1)
            && (result->dad_id.gene_id != -1))
    {
        //offspring comes from mating with the right id's
    }
    else
    {
        /*result->genome_id = newId; result->mom_id = g1->genome_id;
        result->dad_id.gene_id = -1; result->dad_id.robot_id = -1;*/
    }
    result->genome_id = _genomeId;
    return result;
}

bool odNeatController::in_maturation_period()
{
    if(gWorld->getIterations () <=
            _birthdate + odNeatSharedData::gMaturationPeriod)
        return true;
    return false;
}

// ##############POPULATION############################
void odNeatController::add_to_population(Genome* g, double f)
{
    GC receivedId = g->genome_id;
    double receivedF = f;
    bool found = false;
    if(findInPopulation(receivedId))
    {
        //Update fitness estimate
        for(auto itP = _pop.begin(); itP != _pop.end(); itP++)
        {
            for(auto itS = (*itP)->_lGenomes.begin();
                itS != (*itP)->_lGenomes.end(); itS++)
            {
                if(itS->first == receivedId)
                {
                    Genome* existingG = itS->second.g; double existingF = itS->second.f;

                    existingG->nbFitnessUpdates++;
                    itS->second.f = existingF + (receivedF - existingF)
                            /existingG->nbFitnessUpdates;
                    delete g;

                    found = true;
                    break;
                }
            }
            if(found)
                break;
        }
    }
    else
    {        
        int popSize = 0;
        for (auto itP = _pop.begin(); itP != _pop.end(); itP++)
        {
             popSize += (*itP)->_lGenomes.size();
        }
        g->nbFitnessUpdates++;
        if(popSize < odNeatSharedData::gPopulationSize)
        {
            int species = computeSpeciesOfGenome(g);
            if(species != -1)
            {
                //Add to existing species
                for(auto itP = _pop.begin(); itP != _pop.end(); itP++)
                {
                    if((*itP)->_id == species)
                    {
                        (*itP)->add(g, f);
                        break;
                    }
                }
            }
            else
            {
                //Create new species
                odNeatSpecies* newSpecies = new odNeatSpecies(_newSpId);
                _newSpId++;
                newSpecies->add(g, f);
                _pop.push_back(newSpecies);
            }
        }
        else
        {            
            double worseF = 100000.0; GC worseId;             
             std::vector<odNeatSpecies*>::iterator itSpeciesToEraseFrom;
             Genome* gToErase = NULL;
            //Search the worse genome (the one with worse adjusted fitness)
            for(auto itP = _pop.begin();itP != _pop.end();itP++)
            {
                std::map<GC, genomeFitness> lGenomes = (*itP)->_lGenomes;
                for (auto itS = lGenomes.begin(); itS != lGenomes.end(); itS++)
                {
                    //Compare to adjusted fitness
                    if( (itS->second.f/lGenomes.size() < worseF))
                    {
                        worseF = itS->second.f/lGenomes.size();
                        worseId = itS->first;
                        gToErase = itS->second.g;
                        itSpeciesToEraseFrom = itP;
                    }
                }
            }
            add_to_tabu_list(gToErase->duplicate());

            if((*itSpeciesToEraseFrom)->_lGenomes.erase(worseId) != 1)
            {
                std::cerr << "Bad worse id" << std::endl;
                exit(-1);
            }

            if((*itSpeciesToEraseFrom)->_lGenomes.size() == 0)
            {
                _pop.erase(itSpeciesToEraseFrom);
            }
            int species = computeSpeciesOfGenome(g);
            if(species != -1)
            {
                //Add to existing species
                for(auto itP = _pop.begin(); itP != _pop.end(); itP++)
                {
                    if((*itP)->_id == species)
                    {
                        (*itP)->add(g, f);
                        break;
                    }
                }
            }
            else
            {
                //Create new species
                odNeatSpecies* newSpecies = new odNeatSpecies(_newSpId);
                _newSpId++;
                newSpecies->add(g, f);
                _pop.push_back(newSpecies);
            }
        }
    }
    adjustSpeciesFitness();
}

bool odNeatController::findInPopulation(GC gId)
{
    for (auto itP = _pop.begin(); itP != _pop.end(); itP++)
    {
        if((*itP)->has(gId))
            return true;
    }
    return false;
}
bool odNeatController::population_accepts(double f)
{
    int popSize = 0;
    for (auto itP = _pop.begin(); itP != _pop.end(); itP++)
    {
         popSize += (*itP)->_lGenomes.size();
    }
    if(popSize < odNeatSharedData::gPopulationSize)
        return true;
    else
    {
        //Test if there is genome with lower adjusted fitness
        for(auto itP = _pop.begin();itP != _pop.end();itP++)
        {
            std::map<GC, genomeFitness> lGenomes = (*itP)->_lGenomes;
            for (auto itS = lGenomes.begin(); itS != lGenomes.end(); itS++)
            {
                //Compare to adjusted fitness
                if( (itS->second.f/lGenomes.size() < f)
                        && !(itS->first == _genome->genome_id))
                    return true;
            }
        }
    }
    return false;
}

int odNeatController::computeSpeciesOfGenome(Genome* g)
{
    int result = -1;
    for(auto itP = _pop.begin(); itP != _pop.end(); itP++)
    {
        //Representative genome
        Genome* repr = (*itP)->_lGenomes.begin()->second.g;
        //Find first compatible representative genome
        if(g->dissimilarity(repr) < odNeatSharedData::gCompatThreshold)
        {
            result = (*itP)->_id;
            break;
        }
    }
    return result;
}
void odNeatController::adjustSpeciesFitness()
{
    for(auto itP = _pop.begin(); itP != _pop.end(); itP++)    
        (*itP)->computeSpeciesFitness();
}
odNeatSpecies* odNeatController::selectSpecies()
{
    odNeatSpecies* result = NULL;
    double totalAdjFitness = 0.0;

    for (auto it= _pop.begin(); it != _pop.end(); it++)
    {
        totalAdjFitness += (*it)->_speciesFitness;
    }

    double random = Helper::randFloat() * totalAdjFitness;

    auto it = _pop.begin();
    while (random > 0.0)
    {
        random -= (*it)->_speciesFitness;
        it++;
    }
    //This test should not be necessary
    if(random <=0.0)
    {
        it--;
        result = *it;
    }
    if(result == NULL)
    {
        std::cerr << "No species selected..." << std::endl; exit(-1);
    }
    return result;
}
Genome* odNeatController::selectParent(odNeatSpecies* sp)
{
    Genome* result = NULL;
    /*//Intraspecies binary tournament
    auto randomIt1 = sp->_lGenomes.begin(), randomIt2 = sp->_lGenomes.begin();

    if(sp->_lGenomes.size() > 1)
    {
        int ind1 =  rand () % sp->_lGenomes.size();
        int ind2 =  rand () % sp->_lGenomes.size();

        while(ind1 == ind2)
        {
            ind2 =  rand () % sp->_lGenomes.size();
        }

        std::advance(randomIt1,ind1);
        std::advance(randomIt2,ind2);

        if(randomIt1->second.f >= randomIt2->second.f)
        {
            result = randomIt1 -> second.g;
        }
        else
        {
            result = randomIt2-> second.g;
        }
    }
    else
        result = sp->_lGenomes.begin() -> second.g;*/
    //TOERASE BEST
    double bestF = -1.0;
    for(auto it = sp->_lGenomes.begin(); it != sp->_lGenomes.end(); it++)
    {
        if(it -> second.f > bestF)
        {
            result = it -> second.g;
        }
    }

    return result;
}

// ##############TABU LIST############################
bool odNeatController::tabu_list_approves(Genome* g)
{
    bool result = true;

    for(auto it = _tabu.begin();it != _tabu.end();it++)
    {
        if(it->first->dissimilarity(g) < odNeatSharedData::gTabuThreshold)
        {
            result = false; break;
        }
    }

    return result;
}
bool odNeatController::update_tabu_list(Genome* g)
{
    //Returns "is it different than all in tabu?"
    bool result = true;
    auto it = _tabu.begin(), tabuEnd = _tabu.end();

    for(;it != tabuEnd;it++)
    {
        //std::cout << it->first->dissimilarity(g) << std::endl;
        if(it->first->dissimilarity(g) < odNeatSharedData::gTabuThreshold)
        {
            it->second = odNeatSharedData::gTabuTimeout;
            result = false;
        }
        else
        {
            //Decrease timeout of genome on tabu list
            it -> second -= 1;
            //If timeout over, erase genome from tabu list
            if(it->second <= 0)
            {
                delete it->first;
                //End iterator changes on erase
                it = _tabu.erase(it); tabuEnd =  _tabu.end();
                if( it == tabuEnd)
                    break;
            }
        }
    }
    return result;
}

void odNeatController::add_to_tabu_list(Genome* g)
{
    //First, update the timeouts of all non-similar elements
    bool toAdd = update_tabu_list(g);
    //Then add it if it's different than all in tabu
    if(toAdd)
        _tabu[g] = odNeatSharedData::gTabuTimeout;
}

// ################ COMMUNICATION METHODS ################
bool odNeatController::doBroadcast()
{
    bool result = true; //TODO false FALSE!!
    /*
    double adjustedFitness = std::get<1>(species[_genome->species]);
    double totalAdjFitness = 0.0;
    std::map<GC,double>::iterator it = species.begin();
    for(; it != species.end(); it++)
    {
        totalAdjFitness += it->second;
    }

    if(!(totalAdjFitness == 0.0))
        if(randFloat() < (adjustedFitness / totalAdjFitness))
            result = true;*/
    return result;
}
void odNeatController::broadcastGenome()
{
    //Communication based on distance
    for (int i = 0; i < gNumberOfRobots; i++)
    {
        //Do not send to self
        if (i != _wm->_id)
        {
            odNeatController* targetRobotController =
                    dynamic_cast<odNeatController*>
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
                Genome* dupl = _genome->duplicate();
                targetRobotController->
                        storeGenome(message(dupl, _energy, _n_count, _g_count));
            }
        }
    }
}

void odNeatController::storeGenome(message m)
{
    if(odNeatSharedData::gUpdateGC)
    {
        //Update gene clocks for nodes and links:
        //minimize the number of arbitrary sorting orders in genome alignment
        //due to concurrent mutations in different agents
        _n_count = std::max(_n_count,std::get<2>(m));
        _g_count = std::max(_g_count,std::get<3>(m));
    }

    //Starting to measure energy on received genome
    std::get<0>(m) -> nbFitnessUpdates = 1;
    //If genome different from the ones in tabu list
    //and there is space in population (or it is better than the existing)
    if(tabu_list_approves(std::get<0>(m))
            && population_accepts(std::get<1>(m)))
    {
        add_to_population(std::get<0>(m), std::get<1>(m));
    }
    else
    {
        //Delete genome to free memory
        delete std::get<0>(m);
    }
}

void odNeatController::logGenome(std::string s)
{
    //std::ofstream genomeF;
    //genomeF.open(s);
    //TODO
    //genomeF.close();
}

//################BEHAVIOUR METHODS################
void odNeatController::stepBehaviour()
{
    // ---- Build inputs ----
    std::vector<double>* inputs =
            new std::vector<double>(_nbInputs); int inputToUse = 0;
    int type = 1; //object type= energy
    for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
    {
        // distance sensors
        //If there is object but not energy item => wall (get distance sensor) else 0.0 [not1.0]
        int objId = _wm->getObjectIdFromCameraSensor(i);
        if ( PhysicalObject::isInstanceOf(objId) )
        {

            //(*inputs)[inputToUse] = 1.0;
            (*inputs)[inputToUse] = 0.0; inputToUse++;
        }
        else
        {
            //(*inputs)[inputToUse] = _wm->getDistanceValueFromCameraSensor(i) /
            // _wm->getCameraSensorMaximumDistanceValue(i);
            (*inputs)[inputToUse] = 1.0 - _wm->getDistanceValueFromCameraSensor(i) /
                    _wm->getCameraSensorMaximumDistanceValue(i); inputToUse++;
        }

        //If task=collect, add object sensors
        if (odNeatSharedData::gFitness == 1)
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
                //Not a physical object.
                //But: should still fill in the inputs 0.0 //(max distance, 1.0)
                //(*inputs)[inputToUse] = 1.0;
                (*inputs)[inputToUse] = 0.0; inputToUse++;
            }
        }
    }

    (*inputs)[inputToUse++] = _energy / odNeatSharedData::gMaxEnergy;
    //Bias
    (*inputs)[inputToUse++] = 1.0;

    //Previous translational and rotational speeds (acts as recurrent connections from last step)
    /*(*inputs)[inputToUse] = _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
    inputToUse++;
    (*inputs)[inputToUse] = _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
    inputToUse++;*/
    /*for(auto it = inputs->begin(); it < inputs->end(); it++)
        std::cout << (*it) << " ";
    std::cout << std::endl;*/

    nn->load_sensors (&((*inputs)[0]));

    // ---- compute and read out ----
    if (!(nn->activate ()))
    {
        std::cerr << "[ERROR] Activation of ANN not correct: genome R"
                  << _genome->genome_id.robot_id << ", G" << _genome->genome_id.gene_id << std::endl;
        //save_genome();
        exit (-1);
    }

    // Read the output
    std::vector<double> outputs;
    for (auto out_iter  = nn->outputs.begin();
         out_iter != nn->outputs.end();
         out_iter++)
        outputs.push_back((*out_iter)->activation);
    /*for(auto it = outputs.begin(); it < outputs.end(); it++)
        std::cout << (*it) << " ";
    std::cout << std::endl;*/
    //Set the outputs to the right effectors, and rescale the intervals
    //Translational velocity in [-1,+1]. Robots can move backwards. Neat uses sigmoid [0,+1]

    //Direct Kinematic model
    //_wm->_desiredTranslationalValue = outputs[0]; _wm->_desiredRotationalVelocity = outputs[1];

    //Differential Kinematic model
    double lW = outputs[0];
    double rW = outputs[1];
    _wm->_desiredTranslationalValue = (rW + lW) / 2;
    _wm->_desiredRotationalVelocity = (lW - rW) / 2;

    // normalize to motor interval values
    _wm->_desiredTranslationalValue =  _wm->_desiredTranslationalValue * gMaxTranslationalSpeed;
    _wm->_desiredRotationalVelocity = _wm->_desiredRotationalVelocity * gMaxRotationalSpeed;

    delete (inputs);
}
double odNeatController::update_energy_level()
{
    double result = 0.0;
    double transV, rotV, coef_obstacle;
    switch(odNeatSharedData::gFitness)
    {
    case 0:
        transV =  _wm->_desiredTranslationalValue / gMaxTranslationalSpeed;
        rotV = _wm->_desiredRotationalVelocity / gMaxRotationalSpeed;
        coef_obstacle = 1.0;
        for(int i  = 0; i < _wm->_cameraSensorsNb; i++)
        {
            if (_wm->getDistanceValueFromCameraSensor(i)/_wm->getCameraSensorMaximumDistanceValue(i)
                    < coef_obstacle)
                coef_obstacle = _wm->getDistanceValueFromCameraSensor(i)/_wm->getCameraSensorMaximumDistanceValue(i);
        }

        //deltaFitness: energy contribution at current time-step
        //abs(trans) in [0,1], (1 - abs(rot)) in [0,1], coeffObstacle is in [0,1]
        result = fabs(transV) * (1 - fabs(rotV)) * coef_obstacle;
        result =  2.0 * result -1.0;
        break;

    case 1:
        //Fixed rate of energy consumption
        //Energy gathering at energy point done in agent observer
        result = -odNeatSharedData::gEnergyConsumption;
        break;
    default:
        std::cerr << "[ERROR] Unknown fitness function selected." <<
                     "Check gFitness parameter in properties file." << std::endl;
        exit(-1);
    }
    //cap energy (in [0,maxEnergy])
    return std::max(0.0,std::min(result + _energy,odNeatSharedData::gMaxEnergy));
}

//Incrementally update fitness based on energy
void odNeatController::updateFitness ()
{
    _fitnessUpdateCounter++;
    if(_fitnessUpdateCounter >= odNeatSharedData::gFitnessFreq)
    {
        _genome->nbFitnessUpdates++;
        _currentFitness = (_currentFitness) + ((_energy -  _currentFitness)/_genome->nbFitnessUpdates);
        _fitnessUpdateCounter =  0;        
    }
    for(auto itP = _pop.begin(); itP != _pop.end(); itP++)
    {
        if((*itP)->_lGenomes.find(_genomeId) != (*itP)->_lGenomes.end())
        {
            (*itP)->_lGenomes.find(_genomeId)->second.f = _currentFitness;
        }
    }
}
/*
Genome* odNeatController::selectBestGenome()
{
    Genome* result = NULL;
    if(_pop.size() != 0)
    {
        double maxFitness = -1.0, curFit;
        struct GC indexBest;
        indexBest.robot_id = -1; indexBest.gene_id = -1;
        local_population::iterator itP = _pop.begin();
        for(; itP != _pop.end(); itP++)
        {
            std::map<GC, message> lGenomes = std::get<0>(itP->second);
            for(auto itS = lGenomes.begin(); itS != lGenomes.end(); itS++)
            {
                curFit = std::get<1>(itS->second);
                if (curFit > maxFitness)
                {
                    maxFitness = curFit;
                    result = std::get<0>(itS->second);
                }
            }
        }

        //if best individual found, return it
        if (result != NULL)
        {
           //Do nothing (current genome will be mutated)
        }
        else
        {
            //This should never happen
            std::cout << "ERROR in best selection" << std::endl;
            exit(-1);
        }
    }
    else
    {
        //This should never happen
        std::cout << "ERROR in best selection: empty pop" << std::endl;
        exit(-1);
    }
    return result;
}
*/
/*
void odNeatController::cleanPopAndSpecies()
{
    bool ok = true;
    std::map<int,std::pair<std::set<Genome*>,double>>::iterator itTestSp = species.begin();
    std::set<Genome*>::iterator itG;
    std::map<int,message>::iterator itNull = population.begin();
    //Erase Null genomes  (?)
    for(;itNull!= population.end();itNull++)
    {
        if(std::get<0>((itNull->second)) == NULL)
        {
            //Iterator is not invalidated
            population.erase(itNull->first);
        }
    }
    for(;itTestSp != species.end();itTestSp++)
    {
        if(std::get<0>(itTestSp->second).size() == 0 )
            itTestSp = species.erase(itTestSp);
        else
        {
            itG = ((std::get<0>(itTestSp->second)).begin());
            for(; itG != (std::get<0>(itTestSp->second)).end(); itG++)
            {
                if(!findInPopulation(*itG))
                {
                    std::cerr << "[ERROR] Indiv in species not in pop" << std::endl;
                    ok = false;
                }
                if((*itG)->species == -1)
                {
                    std::cerr << "[ERROR] Id species == -1 on cleaning species. Genome == " << (*itG)->genome_id << std::endl;
                    exit(-1);
                }
            }
        }
    }

    std::map<int,message>::iterator itTest = population.begin();

    for(;itTest != population.end();itTest++)
    {
        if(findInSpecies(std::get<0>(itTest->second)) == -1)
        {
            std::cerr << "[ERROR] Cleaning, indiv. in pop not in species" << std::endl;
            ok = false;
        }
        if(std::get<0>(itTest->second)->species == -1)
        {
            std::cerr << "[ERROR] Id species == -1 on cleaning pop. Genome == " << std::get<0>(itTest->second)->genome_id << std::endl;
            exit(-1);
        }

    }
    if(!ok)
    {
        std::cerr << "[ERROR] Something happened on cleaning" << std::endl;
        exit(-1);
    }
}*/
