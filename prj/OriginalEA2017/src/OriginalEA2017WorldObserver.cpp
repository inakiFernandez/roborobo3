/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Observers/AgentObserver.h"
#include "Observers/WorldObserver.h"
#include "OriginalEA2017/include/OriginalEA2017WorldObserver.h"
#include "OriginalEA2017/include/OriginalEA2017Controller.h"
#include "World/World.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <vector>
#include <fstream>

OriginalEA2017WorldObserver::OriginalEA2017WorldObserver( World* world ) : WorldObserver( world )
{
    _world = world;

	// ==== loading project-specific properties

	gProperties.checkAndGetPropertyValue("gSigmaRef",&OriginalEA2017SharedData::gSigmaRef,true);
    gProperties.checkAndGetPropertyValue("gPopSize",&OriginalEA2017SharedData::gPopulationSize,true);
    gProperties.checkAndGetPropertyValue("gClearPopulation",&OriginalEA2017SharedData::gClearPopulation,true);
    gProperties.checkAndGetPropertyValue("gStoreOwn",&OriginalEA2017SharedData::gStoreOwn,true);
    //gProperties.checkAndGetPropertyValue("gSelectionMethod",&OriginalEA2017SharedData::gSelectionMethod,true);

    gProperties.checkAndGetPropertyValue("gCommunicationBySensors",&OriginalEA2017SharedData::gCommunicationBySensors,true);
    gProperties.checkAndGetPropertyValue("gFitness",&OriginalEA2017SharedData::gFitness,true);

    std::string tasks;
    gProperties.checkAndGetPropertyValue("gTaskSeq",&tasks,true);
    std::vector<std::string> taskV;
    boost::algorithm::split(taskV, tasks, boost::algorithm::is_any_of(","));
    for(auto it = taskV.begin(); it != taskV.end(); it++)
    {
        int t = stoi((*it));
        OriginalEA2017SharedData::gTaskSeq.push_back(t);
    }
    std::cout << std::endl;
    if(OriginalEA2017SharedData::gTaskSeq.size() > 0)
    {
        OriginalEA2017SharedData::gFitness = OriginalEA2017SharedData::gTaskSeq[0];
        OriginalEA2017SharedData::gTaskIdx = 1;
    }
    listCollected.resize(gNbOfPhysicalObjects);

    std::string transitionTime;
    gProperties.checkAndGetPropertyValue("gTimeChange",&transitionTime,true);
    std::vector<std::string> timeV;
    boost::algorithm::split(timeV, transitionTime, boost::algorithm::is_any_of(","));
    for(auto it = timeV.begin(); it != timeV.end(); it++)
    {
        int t = stoi((*it));
        OriginalEA2017SharedData::gTimeSeq.push_back(t);
    }
    if(OriginalEA2017SharedData::gTimeSeq.size() != OriginalEA2017SharedData::gTaskSeq.size())
    {
        std::cerr << "Task sequence and time of transition sequence not same size." << std::endl;
        exit(-1);
    }

    gProperties.checkAndGetPropertyValue("gEvaluationTime",&OriginalEA2017SharedData::gEvaluationTime,true);

    gProperties.checkAndGetPropertyValue("gControllerType",&OriginalEA2017SharedData::gControllerType,true);

    gProperties.checkAndGetPropertyValue("gNbHiddenLayers",&OriginalEA2017SharedData::gNbHiddenLayers,true);
	gProperties.checkAndGetPropertyValue("gNbNeuronsPerHiddenLayer",&OriginalEA2017SharedData::gNbNeuronsPerHiddenLayer,true);
	gProperties.checkAndGetPropertyValue("gNeuronWeightRange",&OriginalEA2017SharedData::gNeuronWeightRange,true);
    gProperties.checkAndGetPropertyValue("gWithBias",&OriginalEA2017SharedData::gWithBias,true);

    gProperties.checkAndGetPropertyValue("gOutGenomeFile",&OriginalEA2017SharedData::gOutGenomeFile,true);

    gProperties.checkAndGetPropertyValue("gEvolutionLogFile",&OriginalEA2017SharedData::gEvolutionLogFile,true);

    gProperties.checkAndGetPropertyValue("gIsLoadGenome",&OriginalEA2017SharedData::gIsLoadGenome,true);
    gProperties.checkAndGetPropertyValue("gLogGenome",&OriginalEA2017SharedData::gSaveGenome,true);

    gProperties.checkAndGetPropertyValue("gWithCollectColorEffector",&OriginalEA2017SharedData::gWithCollectColorEffector,true);

    gProperties.checkAndGetPropertyValue("gBrait",&OriginalEA2017SharedData::gBrait,true);
    gProperties.checkAndGetPropertyValue("gSelPressure",&OriginalEA2017SharedData::gSelPressure, true);
    gProperties.checkAndGetPropertyValue("gForgetMethod",&OriginalEA2017SharedData::gForgetMethod,true);

    //OdNeat----------------------------------
    //Mutations
    gProperties.checkAndGetPropertyValue("mutate_only_prob",&Helper::mutateProb,true);
    gProperties.checkAndGetPropertyValue("mutate_link_weights_prob",&Helper::mutateLinkWeightsProb,true);
    gProperties.checkAndGetPropertyValue("mutate_individual_weight_prob",&Helper::mutateIndividualWeightProb,true);
    gProperties.checkAndGetPropertyValue("mutate_add_node_prob",&Helper::mutateAddNodeProb,true);
    gProperties.checkAndGetPropertyValue("mutate_add_link_prob",&Helper::mutateAddLinkProb,true);
    gProperties.checkAndGetPropertyValue("mate_only_prob",&Helper::mateOnlyProb,true);
    gProperties.checkAndGetPropertyValue("recur_only_prob",&Helper::recurOnlyProb,true);
    gProperties.checkAndGetPropertyValue("newstructure_tries",&Helper::newStructureTries,true);
    gProperties.checkAndGetPropertyValue("mutate_toggle_enable_prob",&Helper::mutateToggleEnableProb,true);
    gProperties.checkAndGetPropertyValue("allowMultisynapses",&Helper::allowMultisynapses,true);



    Helper::rangeW = OriginalEA2017SharedData::gNeuronWeightRange / 2;

    std::string logItemName = "items.log";
    gProperties.checkAndGetPropertyValue("logItemName",&logItemName,true);
    //logItemFile = std::ofstream("logs/items.log");// + logItemName);
    /*logItemFile = std::ofstream(//"logs/" +
                                logItemName);*/
    logItemFile.open(logItemName);

    std::string logItGatheredName = "itemsIter.log";
    gProperties.checkAndGetPropertyValue("logItGatheredName",&logItGatheredName,true);
    //logItGatheredFile = std::ofstream("logs/itemsIter.log");// + logItGatheredName);
    /*logItGatheredFile = std::ofstream(//"logs/"+
                                      logItGatheredName);*/
    logItGatheredFile.open(logItGatheredName);

    std::string logColorChangesName = "colorChanges.log";
    gProperties.checkAndGetPropertyValue("logColorChangesName",&logColorChangesName,true);
    //logChangesColorFile = std::ofstream("logs/colorChanges.log");// + logColorChangesName);
    /*logChangesColorFile = std::ofstream(//"logs/" +
                                        logColorChangesName);*/
    logChangesColorFile.open(logColorChangesName);

    std::string logGivenRewardName = "givenReward.log";
    gProperties.checkAndGetPropertyValue("logGivenRewardName",&logGivenRewardName,true);

    /*logGivenRewardFile = std::ofstream(//"logs/" +
                                        logGivenRewardName);*/
    logGivenRewardFile.open(logGivenRewardName);

    logItemFile << "c1 c2 c3 c4 c5 c6 c7 c8" << std::endl;

    int numberColors = 8;
    itemCounts = std::vector<int>(numberColors);
    for(int i = 0 ; i < numberColors; i++)
    {
        itemCounts[i] = 0;
    }

    //  iteration and generation counters
	_lifeIterationCount = -1;
	_generationCount = -1;

    //Given average individual reward
    _averageReward = 0.0;
    _nbRobotsCorrect = 0.0;
}

OriginalEA2017WorldObserver::~OriginalEA2017WorldObserver()
{
	// nothing to do.
}

void OriginalEA2017WorldObserver::reset()
{
	// nothing to do.
}

void OriginalEA2017WorldObserver::step()
{
    _lifeIterationCount++;
    
    updateMonitoring();

    if( _lifeIterationCount >= OriginalEA2017SharedData::gEvaluationTime ) // switch to next generation.
	{
        // update iterations and generations counters
        _lifeIterationCount = 0;
        _generationCount++;
    }
}


void OriginalEA2017WorldObserver::updateEnvironment()
{
}

void OriginalEA2017WorldObserver::updateMonitoring()
{
    // * Log at end of each generation
    //std::cout << gWorld->getIterations() << std::endl;
    if( _lifeIterationCount >= OriginalEA2017SharedData::gEvaluationTime ) // end of generation.
	{
		if ( gVerbose )
		{
            std::cout << "[gen:" << (gWorld->getIterations()/OriginalEA2017SharedData::gEvaluationTime)
                      << "]\n";
		}
        // Logging here
        double sumFitness = 0.0;
        int gatheredGenomes = 0;
        double forgetMeasure = 0.0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             sumFitness += (dynamic_cast<OriginalEA2017Controller*>(gWorld->getRobot(i)->getController()))
                     -> getFitness();
             OriginalEA2017Controller* ctrl = (dynamic_cast<OriginalEA2017Controller*>
                                         (gWorld->getRobot(i)->getController()));
             gatheredGenomes += (ctrl->_doEvoTopo? ctrl ->_genomesList.size():ctrl ->_genomesFList.size());

             if(ctrl -> _storedF.empty())
            {
                 forgetMeasure += -1.0;
            }
            else
            {
                forgetMeasure += ctrl -> forget();
            }
            //Forget measure for each robot ? what to do
            //Attention! sum of fitness is not always number of items
        }
        //std::cout << gWorld->getIterations() << " ";
        //<< (sumFitness  / gNumberOfRobots) / Collect2SharedData::gEvaluationTime
        //average forget measure. TODO other estimator/or all the data?


            std::cout << sumFitness
                      //<< " " << (forgetMeasure/ gNumberOfRobots)
                      <<  std::endl;

            if(fabs(sumFitness) < 0.0000000001)
            {
                logGivenRewardFile << 0.0 << " "
                               << 0.0 << " "
                               << std::endl;
            }
            else
            {
                logGivenRewardFile << _averageReward / sumFitness << " "
                               << _nbRobotsCorrect / sumFitness << " "
                               << std::endl;
            }
            _nbRobotsCorrect = 0.0;
            _averageReward = 0.0;

	}


    if (gWorld->getIterations() == (gMaxIt - 1))
    {
        double sumFitness = 0.0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             sumFitness += (dynamic_cast<OriginalEA2017Controller*>(gWorld->getRobot(i)
                                                    ->getController()))-> getFitness();

        }

        std::cout << "End fitness: " << sumFitness  / gNumberOfRobots
                  << " at it: " << gWorld->getIterations() << std::endl;
        if(OriginalEA2017SharedData::gSaveGenome)
        {
            for (int i = 0 ; i != gNumberOfRobots ; i++ )
            {
                OriginalEA2017Controller* ctrl = (dynamic_cast<OriginalEA2017Controller*>(gWorld->getRobot(i)->getController()));
                if(!ctrl->_doEvoTopo)
                        ctrl -> logGenomeF(OriginalEA2017SharedData::gOutGenomeFile + std::to_string(i) + ".log");
                else
                {
                    //TODO
                }
            }
        }
        logItemFile.close();
        logItGatheredFile.close();

    }
    int countItGathered = 0;
    int countPossible = 0;
    //HIC SUNT COLORIS
    //THINK OF THE SENSOR WE NEED!!
    switch (OriginalEA2017SharedData::gFitness)
    {
        case 2:
            for(int i = 0; i < gNbOfPhysicalObjects;i++)
            {                 
                double color = -2.0;
                //VARIANT with same color as item
                color = gPhysicalObjects[i]->getColorValue();
                //Value is 0.125 more to get symmetric color values between positive and negative values
                //color += 0.125; //TOTEST
                //For other VARIANT bool isSynchColor = false;
                bool isSynchColor = true;
                int nbRobotsCorrectColor = 0;
                if(listCollected[i].size() >= 2)
                {
                    //test if all agents (maybe more than 2) same color
                    for(auto it = listCollected[i].begin(); it != listCollected[i].end();it++)
                    {                        
                        /*if(color == -2.0)
                        {
                            //Unused now
                            color = it->second;
                            isSynchColor = true;
                        }
                        else
                        {*/
                        //TOCHECK on insertion in the list (maybe also controller output neuron) //Value is 0.125 higher!!
                        if(color == it->second)
                            nbRobotsCorrectColor++;
                            /*if(color != it->second)
                            {
                                isSynchColor = false;
                                break;
                            }*/
                        //}
                    }
                    isSynchColor = (nbRobotsCorrectColor >= 2);
                    if(isSynchColor)
                    {
                        gPhysicalObjects[i]->isWalked(0); //Default agent for callback (unused callback)
                        for(auto it = listCollected[i].begin(); it != listCollected[i].end();it++)
                        {
                            //Share reward of one item between the agents involved
                            //TOCHANGE for different reward
                            if(it->second == color)
                            {
                                dynamic_cast<OriginalEA2017Controller*>(gWorld->getRobot(it->first)
                                    ->getController())->updateFitness(1.0 / (double)nbRobotsCorrectColor);

                            }
                        }
                        _averageReward += 1.0 / (double)nbRobotsCorrectColor; //TOCHECK
                        _nbRobotsCorrect +=(double)nbRobotsCorrectColor;
                        int colorInt = int((color + 0.875) * 4.0); // int((color + 1.0) * 4.0);
                        countItGathered++;
                        itemCounts[colorInt]++;
                    }
                    countPossible++;
                }
                listCollected[i].clear();
            }
            for(int i = 0; i < gNumberOfRobots; i++)
            {
                Uint8 r, g, b;
                RobotWorldModel* wm = gWorld->getRobot(i)->getWorldModel();
                Uint32 pixel = getPixel32( gGroundSensorImage, wm->_xReal+0.5, wm->_yReal+0.5);
                SDL_GetRGB(pixel,gGroundSensorImage->format,&r,&g,&b);
                wm->_groundSensorValue[0] = r;
                wm->_groundSensorValue[1] = g;
                wm->_groundSensorValue[2] = b;
            }
            logItGatheredFile << countItGathered << " ";
            logItGatheredFile << countPossible << " ";
            logItGatheredFile << std::endl;
            break;
        default:
            break;
    }
    //TOERASE test seamless node mutation
    /*if(gWorld->getIterations() == 10000)
    {
        Helper::mutateLinkWeightsProb = 0.0;
        Helper::mutateAddNodeProb = 1.0;
        std::cout << "Now only node mutation" << std::endl;
    }*/
    if(gWorld->getIterations() ==
            OriginalEA2017SharedData::gTimeSeq[OriginalEA2017SharedData::gTaskIdx])
    {
        //std::cout << "----------------------------" << std::endl;
        //std::cout << "Task changed!" << std::endl;
        OriginalEA2017SharedData::gFitness =
                OriginalEA2017SharedData::gTaskSeq[OriginalEA2017SharedData::gTaskIdx];
        OriginalEA2017SharedData::gTaskIdx++;

        //Store representative for forget measure? For each robot? (previous task)
        for (int i = 0 ; i != gNumberOfRobots ; i++)
        {
            (dynamic_cast<OriginalEA2017Controller*>(gWorld->getRobot(i)->getController()))->storeRepresentative();
        }
    }
    for(int i = 0; i < gNumberOfRobots; i++)
    {
        OriginalEA2017Controller* ctrl = (dynamic_cast<OriginalEA2017Controller*>(gWorld->getRobot(i)->getController()));

        logChangesColorFile << ctrl->getNbColorChanges() << " ";
    }
    logChangesColorFile << std::endl;
    if( _lifeIterationCount >= OriginalEA2017SharedData::gEvaluationTime ) // end of generation.
    {

        int numberColors = 8;
        for(int i = 0 ; i < numberColors; i++)
        {
            logItemFile << itemCounts[i] << " ";
            itemCounts[i] = 0;
        }
        logItemFile << std::endl;

        /*for(int i = 0; i < gNumberOfRobots; i++)
        {
            OriginalEA2017Controller* ctrl = (dynamic_cast<OriginalEA2017Controller*>(gWorld->getRobot(i)->getController()));

            logChangesColorFile << ctrl->getNbColorChanges() << " ";
        }
        logChangesColorFile << std::endl;*/
    }
}

