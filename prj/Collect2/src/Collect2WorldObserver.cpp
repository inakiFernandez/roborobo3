/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Observers/AgentObserver.h"
#include "Observers/WorldObserver.h"
#include "Collect2/include/Collect2WorldObserver.h"
#include "Collect2/include/Collect2Controller.h"
#include "World/World.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <vector>

Collect2WorldObserver::Collect2WorldObserver( World* world ) : WorldObserver( world )
{
    _world = world;

	// ==== loading project-specific properties

	gProperties.checkAndGetPropertyValue("gSigmaRef",&Collect2SharedData::gSigmaRef,true);
    gProperties.checkAndGetPropertyValue("gPopSize",&Collect2SharedData::gPopulationSize,true);
    gProperties.checkAndGetPropertyValue("gClearPopulation",&Collect2SharedData::gClearPopulation,true);
    gProperties.checkAndGetPropertyValue("gStoreOwn",&Collect2SharedData::gStoreOwn,true);
    //gProperties.checkAndGetPropertyValue("gSelectionMethod",&Collect2SharedData::gSelectionMethod,true);

    gProperties.checkAndGetPropertyValue("gCommunicationBySensors",&Collect2SharedData::gCommunicationBySensors,true);
    gProperties.checkAndGetPropertyValue("gFitness",&Collect2SharedData::gFitness,true);

    std::string tasks;
    gProperties.checkAndGetPropertyValue("gTaskSeq",&tasks,true);
    std::vector<std::string> taskV;
    boost::algorithm::split(taskV, tasks, boost::algorithm::is_any_of(","));
    for(auto it = taskV.begin(); it != taskV.end(); it++)
    {
        int t = stoi((*it));
        Collect2SharedData::gTaskSeq.push_back(t);
    }
    std::cout << std::endl;
    if(Collect2SharedData::gTaskSeq.size() > 0)
    {
        Collect2SharedData::gFitness = Collect2SharedData::gTaskSeq[0];
        Collect2SharedData::gTaskIdx = 1;
    }
    listCollected.resize(gNbOfPhysicalObjects);

    std::string transitionTime;
    gProperties.checkAndGetPropertyValue("gTimeChange",&transitionTime,true);
    std::vector<std::string> timeV;
    boost::algorithm::split(timeV, transitionTime, boost::algorithm::is_any_of(","));
    for(auto it = timeV.begin(); it != timeV.end(); it++)
    {
        int t = stoi((*it));
        Collect2SharedData::gTimeSeq.push_back(t);
    }
    if(Collect2SharedData::gTimeSeq.size() != Collect2SharedData::gTaskSeq.size())
    {
        std::cerr << "Task sequence and time of transition sequence not same size." << std::endl;
        exit(-1);
    }

    gProperties.checkAndGetPropertyValue("gEvaluationTime",&Collect2SharedData::gEvaluationTime,true);

    gProperties.checkAndGetPropertyValue("gControllerType",&Collect2SharedData::gControllerType,true);

    gProperties.checkAndGetPropertyValue("gNbHiddenLayers",&Collect2SharedData::gNbHiddenLayers,true);
	gProperties.checkAndGetPropertyValue("gNbNeuronsPerHiddenLayer",&Collect2SharedData::gNbNeuronsPerHiddenLayer,true);
	gProperties.checkAndGetPropertyValue("gNeuronWeightRange",&Collect2SharedData::gNeuronWeightRange,true);
    gProperties.checkAndGetPropertyValue("gWithBias",&Collect2SharedData::gWithBias,true);

    gProperties.checkAndGetPropertyValue("gOutGenomeFile",&Collect2SharedData::gOutGenomeFile,true);

	// * iteration and generation counters

	_lifeIterationCount = -1;
	_generationCount = -1;

}

Collect2WorldObserver::~Collect2WorldObserver()
{
	// nothing to do.
}

void Collect2WorldObserver::reset()
{
	// nothing to do.
}

void Collect2WorldObserver::step()
{
    _lifeIterationCount++;
    
    updateMonitoring();

    if( _lifeIterationCount >= Collect2SharedData::gEvaluationTime ) // switch to next generation.
	{
        // update iterations and generations counters
        _lifeIterationCount = 0;
        _generationCount++;
    }
}


void Collect2WorldObserver::updateEnvironment()
{
}

void Collect2WorldObserver::updateMonitoring()
{
    // * Log at end of each generation

    if( _lifeIterationCount >= Collect2SharedData::gEvaluationTime ) // end of generation.
	{
		if ( gVerbose )
		{
            std::cout << "[gen:" << (gWorld->getIterations()/Collect2SharedData::gEvaluationTime)
                      << "]\n";
		}
        // Logging here
        double sumFitness = 0.0;
        double sumAvgLocalPopFitness = 0.0;
        int gatheredGenomes = 0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             sumFitness += (dynamic_cast<Collect2Controller*>(gWorld->getRobot(i)->getController()))
                     -> getFitness();
             sumAvgLocalPopFitness += (dynamic_cast<Collect2Controller*>
                                       (gWorld->getRobot(i)->getController())) -> getAvgPopFitness();
             gatheredGenomes += (dynamic_cast<Collect2Controller*>
                                 (gWorld->getRobot(i)->getController())) ->_genomesList.size();
        }
        std::cout << gWorld->getIterations() << " ";
        //<< (sumFitness  / gNumberOfRobots) / Collect2SharedData::gEvaluationTime

        // divided by two because each item gives 1 fitness point to both agents
        std::cout << sumFitness / 2 << " " << sumAvgLocalPopFitness / 2 << std::endl;
        //<< gatheredGenomes
	}

    if (gWorld->getIterations() == (gMaxIt - 1))
    {
        double sumFitness = 0.0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             sumFitness += (dynamic_cast<Collect2Controller*>(gWorld->getRobot(i)
                                                    ->getController()))-> getFitness();

        }

        std::cout << "End fitness: " << sumFitness  / gNumberOfRobots
                  << " at it: " << gWorld->getIterations() << std::endl;
        if(Collect2SharedData::gSaveGenome)
        {
            for (int i = 0 ; i != gNumberOfRobots ; i++ )
            {
                (dynamic_cast<Collect2Controller*>(gWorld->getRobot(i)->getController()))
                        -> logGenome(Collect2SharedData::gOutGenomeFile + std::to_string(i) + ".log");
            }
        }
    }
    switch (Collect2SharedData::gFitness)
    {
        case 2:
            for(int i = 0; i < gNbOfPhysicalObjects;i++)
            {
                if(listCollected[i].size() >= 2)
                {
                    gPhysicalObjects[i]->isWalked(0); //Default agent for callback (unused callback)
                    for(auto it = listCollected[i].begin(); it != listCollected[i].end();it++)
                    {
                        dynamic_cast<Collect2Controller*>(gWorld->getRobot((*it))
                            ->getController())->updateFitness(1.0);
                    }
                }
                listCollected[i].clear();
            }
            updateDisplay();
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
            break;
        default:
            break;
    }
    if(gWorld->getIterations() ==
            Collect2SharedData::gTimeSeq[Collect2SharedData::gTaskIdx])
    {
        std::cout << "Task changed!" << std::endl;
        Collect2SharedData::gFitness =
                Collect2SharedData::gTaskSeq[Collect2SharedData::gTaskIdx];
        Collect2SharedData::gTaskIdx++;
    }
}

