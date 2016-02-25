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
    _currentFitness = 0.0;
    _energy  = 0.0;
    _items = 0;

    resetRobot();
    /*for (int i=0; i < 20; i++)
    {
        Genome* off = _genome->
                mutate(odNeatSharedData::gSigmaRef, _wm->_id,_genomeId,_n_count,_g_count);
        for (int j = 0; j < 10; j++)
        {
            off = off-> mutate(odNeatSharedData::gSigmaRef,
                                   _wm->_id,_genomeId,_n_count,_g_count);
        }
        std::cout << _genome->dissimilarity(off) << std::endl;
    }
    exit(-1);*/
    _genome->genome_id = _genomeId;
    _lifetime = -1;
    // behaviour
    _iteration = 0; _birthdate = 0;
    _wm->updateLandmarkSensor();
    _wm->setRobotLED_colorValues(255, 0, 0);
    _newSpId = 0;
}

odNeatController::~odNeatController()
{
    //TODO invoke delete on each object?
    _pop.clear();
    _tabu.clear();
    delete nn;
    nn = NULL;
    delete _genome;
}

void odNeatController::reset()
{
    //std::cout << "Changed genome at it. " << gWorld->getIterations()
    //          << " with fitness: " << _currentFitness << std::endl;
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

    createNN();

    recompute_all_species();
    //Add new genome (even if fitness is not yet measured)
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
    //_nbInputs += 1;
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
    if(doBroadcast())
        broadcastGenome();

    stepBehaviour ();

    _energy = update_energy_level(); updateFitness();

    _wm->setRobotLED_colorValues((int) (256 * _energy / odNeatSharedData::gMaxEnergy),
                                 0, (int) (256 * _energy / odNeatSharedData::gMaxEnergy));

    if ((_energy <=odNeatSharedData::gEnergyThreshold)
            && !(in_maturation_period()))
    {
        //std::cout << "Changed" << std::endl;
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

    /*for(auto itP = _pop.begin(); itP != _pop.end(); itP++)
    {
        for(auto itS = (*itP)->_lGenomes.begin();
            itS != (*itP)->_lGenomes.end(); itS++)
        {
            for(auto itS2 = (*itP)->_lGenomes.begin();
                itS2 != (*itP)->_lGenomes.end(); itS2++)
            {
                std::cout << itS->second.g->dissimilarity(itS2->second.g) << std::endl;
            }

        }
        std::cout << "####################" << std::endl;
    }*/

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
    genomeFitness gF1 = selectParent(sp); genomeFitness gF2 = selectParent(sp);
    Genome* g1 = gF1.g; double f1 = gF1.f;
    Genome* g2 = gF1.g; double f2 = gF2.f;
    //Mate with probability, only if selected genomes are not the same
    if((Helper::randFloat() < Helper::mateOnlyProb) && !(g1->genome_id == g2->genome_id))
        result = g1 -> mate(g2,_genomeId, f1,f2);
    else
    {
        result = g1;
    }
    if(Helper::randFloat() < Helper::mutateProb)//Mutate
    {
        result = result-> mutate(odNeatSharedData::gSigmaRef,
                                  _wm->_id,_genomeId,_n_count,_g_count);
    }

    //TODO check mom and dad id's

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
            odNeatSpecies* species = computeSpeciesOfGenome(g);

            if(species != NULL)
            {
                species->add(g, f);
                g->species = species->_id;
            }
            else
            {
                //Create new species
                odNeatSpecies* newSpecies = new odNeatSpecies(_newSpId);
                _newSpId++;
                newSpecies->add(g, f);
                g->species = newSpecies->_id;
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
                    if( (itS->second.f/lGenomes.size() < worseF)
                            && !(itS->first == _genome->genome_id))
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
            odNeatSpecies* species = computeSpeciesOfGenome(g);
            if(species != NULL)
            {
                //Add to existing species
                species->add(g, f);
                g->species = species->_id;
            }
            else
            {
                //Create new species
                odNeatSpecies* newSpecies = new odNeatSpecies(_newSpId);
                _newSpId++;
                newSpecies->add(g, f);
                g->species = newSpecies->_id;
                _pop.push_back(newSpecies);
            }
        }
    }
    adjustSpeciesFitness();
    if(_genome->genome_id == g->genome_id)
        _genome->species = g->species;
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

odNeatSpecies* odNeatController::computeSpeciesOfGenome(Genome* g)
{
    odNeatSpecies* result = NULL;
    if(g->species == -1)
    {
        for(auto itP = _pop.begin(); itP != _pop.end(); itP++)
        {
            //Representative genome
            Genome* repr = (*itP)->_lGenomes.begin()->second.g;
            //Find first compatible representative genome
            if(g->dissimilarity(repr) < odNeatSharedData::gCompatThreshold)
            {
                result = (*itP);
                break;
            }
        }
    }
    else
    {
        for(auto itP = _pop.begin(); itP != _pop.end(); itP++)
        {
            if((*itP)->_id == g->species)
            {
                result = (*itP);
                break;
            }
        }
    }
    return result;
}
void odNeatController::recompute_all_species()
{
    std::map<GC, genomeFitness> allGenomes;
    //Regroup all genomes in population
    for(auto itP = _pop.begin(); itP != _pop.end(); )
    {
        for(auto itS = (*itP)->_lGenomes.begin();
            itS != (*itP)->_lGenomes.end(); itS++)
        {
            allGenomes.insert((*itS));
        }
        delete (*itP);
        itP = _pop.erase(itP);
    }
    _pop.clear();
    //Reset species counter
    _newSpId = 0;
    for(auto itAll = allGenomes.begin(); itAll != allGenomes.end();
        itAll++)
    {
        add_to_population(itAll->second.g, itAll->second.f);
    }
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
genomeFitness odNeatController::selectParent(odNeatSpecies* sp)
{
    Genome* result = NULL; double rF = -1.0;

    //Intraspecies binary tournament
    auto randomIt1 = sp->_lGenomes.begin(), randomIt2 = sp->_lGenomes.begin();

    if(sp->_lGenomes.size() > 1)
    {
        int ind1 =  rand () % sp->_lGenomes.size();
        int ind2 =  rand () % sp->_lGenomes.size();

        //Sampling without replacement
        while(ind1 == ind2)
        {
            ind2 =  rand () % sp->_lGenomes.size();
        }

        std::advance(randomIt1,ind1);
        std::advance(randomIt2,ind2);

        //Pick the best of the two sampled genomes
        if(randomIt1->second.f >= randomIt2->second.f)
        {
            result = randomIt1 -> second.g;
            rF = randomIt1 -> second.f;
        }
        else
        {
            result = randomIt2-> second.g;
            rF = randomIt2 -> second.f;
        }
    }
    else
    {
        result = sp->_lGenomes.begin()-> second.g;
        rF = sp->_lGenomes.begin()-> second.f;
    }
    //TOERASE BEST in species
    /*double bestF = -1.0;
    for(auto it = sp->_lGenomes.begin(); it != sp->_lGenomes.end(); it++)
    {
        if(it -> second.f > bestF)
        {
            result = it -> second.g;
        }
    }*/
    genomeFitness r; r.g = result; r.f = rF;
    return r;
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
    //Updates timeout of different genomes
    bool result = true;
    auto it = _tabu.begin(), tabuEnd = _tabu.end();    

    for(;it != tabuEnd;it++)
    {
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
    bool result = false;

    double adjustedFitness = computeSpeciesOfGenome(_genome)->_speciesFitness;

    double totalAdjFitness = 0.0;
    for(auto it = _pop.begin(); it != _pop.end(); it++)
    {
        totalAdjFitness += (*it)->_speciesFitness;
    }

    if(!(totalAdjFitness == 0.0))
        if(Helper::randFloat() < (adjustedFitness / totalAdjFitness))
            result = true;
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
        //Delete genome because not used
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

    //(*inputs)[inputToUse++] = _energy / odNeatSharedData::gMaxEnergy;
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
                coef_obstacle = _wm->getDistanceValueFromCameraSensor(i)
                        /_wm->getCameraSensorMaximumDistanceValue(i);
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
    //Update genome's fitness in population
    for(auto itP = _pop.begin(); itP != _pop.end(); itP++)
    {
        //Search in current species
        if((*itP)->_lGenomes.find(_genomeId) != (*itP)->_lGenomes.end())
        {
            (*itP)->_lGenomes.find(_genomeId)->second.f = _currentFitness;
            break;
        }
    }
    adjustSpeciesFitness();
}
