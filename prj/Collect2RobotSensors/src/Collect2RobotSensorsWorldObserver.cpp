/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Observers/AgentObserver.h"
#include "Observers/WorldObserver.h"
#include "Collect2RobotSensors/include/Collect2RobotSensorsWorldObserver.h"
#include "Collect2RobotSensors/include/Collect2RobotSensorsController.h"
#include "World/World.h"
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <vector>

Collect2RobotSensorsWorldObserver::Collect2RobotSensorsWorldObserver( World* world ) : WorldObserver( world )
{
    _world = world;

	// ==== loading project-specific properties

	gProperties.checkAndGetPropertyValue("gSigmaRef",&Collect2RobotSensorsSharedData::gSigmaRef,true);
    gProperties.checkAndGetPropertyValue("gPopSize",&Collect2RobotSensorsSharedData::gPopulationSize,true);
    gProperties.checkAndGetPropertyValue("gClearPopulation",&Collect2RobotSensorsSharedData::gClearPopulation,true);
    gProperties.checkAndGetPropertyValue("gStoreOwn",&Collect2RobotSensorsSharedData::gStoreOwn,true);
    //gProperties.checkAndGetPropertyValue("gSelectionMethod",&Collect2RobotSensorsSharedData::gSelectionMethod,true);

    gProperties.checkAndGetPropertyValue("gCommunicationBySensors",&Collect2RobotSensorsSharedData::gCommunicationBySensors,true);
    gProperties.checkAndGetPropertyValue("gFitness",&Collect2RobotSensorsSharedData::gFitness,true);

    std::string tasks;
    gProperties.checkAndGetPropertyValue("gTaskSeq",&tasks,true);
    std::vector<std::string> taskV;
    boost::algorithm::split(taskV, tasks, boost::algorithm::is_any_of(","));
    for(auto it = taskV.begin(); it != taskV.end(); it++)
    {
        int t = stoi((*it));
        Collect2RobotSensorsSharedData::gTaskSeq.push_back(t);
    }
    std::cout << std::endl;
    if(Collect2RobotSensorsSharedData::gTaskSeq.size() > 0)
    {
        Collect2RobotSensorsSharedData::gFitness = Collect2RobotSensorsSharedData::gTaskSeq[0];
        Collect2RobotSensorsSharedData::gTaskIdx = 1;
    }
    listCollected.resize(gNbOfPhysicalObjects);

    std::string transitionTime;
    gProperties.checkAndGetPropertyValue("gTimeChange",&transitionTime,true);
    std::vector<std::string> timeV;
    boost::algorithm::split(timeV, transitionTime, boost::algorithm::is_any_of(","));
    for(auto it = timeV.begin(); it != timeV.end(); it++)
    {
        int t = stoi((*it));
        Collect2RobotSensorsSharedData::gTimeSeq.push_back(t);
    }
    if(Collect2RobotSensorsSharedData::gTimeSeq.size() != Collect2RobotSensorsSharedData::gTaskSeq.size())
    {
        std::cerr << "Task sequence and time of transition sequence not same size." << std::endl;
        exit(-1);
    }

    gProperties.checkAndGetPropertyValue("gEvaluationTime",&Collect2RobotSensorsSharedData::gEvaluationTime,true);

    gProperties.checkAndGetPropertyValue("gControllerType",&Collect2RobotSensorsSharedData::gControllerType,true);

    gProperties.checkAndGetPropertyValue("gNbHiddenLayers",&Collect2RobotSensorsSharedData::gNbHiddenLayers,true);
	gProperties.checkAndGetPropertyValue("gNbNeuronsPerHiddenLayer",&Collect2RobotSensorsSharedData::gNbNeuronsPerHiddenLayer,true);
	gProperties.checkAndGetPropertyValue("gNeuronWeightRange",&Collect2RobotSensorsSharedData::gNeuronWeightRange,true);
    gProperties.checkAndGetPropertyValue("gWithBias",&Collect2RobotSensorsSharedData::gWithBias,true);

    gProperties.checkAndGetPropertyValue("gOutGenomeFile",&Collect2RobotSensorsSharedData::gOutGenomeFile,true);

	// * iteration and generation counters

	_lifeIterationCount = -1;
	_generationCount = -1;

}

Collect2RobotSensorsWorldObserver::~Collect2RobotSensorsWorldObserver()
{
	// nothing to do.
}

void Collect2RobotSensorsWorldObserver::reset()
{
	// nothing to do.
}

void Collect2RobotSensorsWorldObserver::step()
{
    _lifeIterationCount++;
    
    updateMonitoring();

    if( _lifeIterationCount >= Collect2RobotSensorsSharedData::gEvaluationTime ) // switch to next generation.
	{
        // update iterations and generations counters
        _lifeIterationCount = 0;
        _generationCount++;
    }
}


void Collect2RobotSensorsWorldObserver::updateEnvironment()
{
}

void Collect2RobotSensorsWorldObserver::updateMonitoring()
{
    // * Log at end of each generation
    //std::cout << gWorld->getIterations() << std::endl;
    if( _lifeIterationCount >= Collect2RobotSensorsSharedData::gEvaluationTime ) // end of generation.
	{
		if ( gVerbose )
		{
            std::cout << "[gen:" << (gWorld->getIterations()/Collect2RobotSensorsSharedData::gEvaluationTime)
                      << "]\n";
		}
        // Logging here
        double sumFitness = 0.0;
        double sumAvgLocalPopFitness = 0.0;
        int gatheredGenomes = 0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             sumFitness += (dynamic_cast<Collect2RobotSensorsController*>(gWorld->getRobot(i)->getController()))
                     -> getFitness();
             sumAvgLocalPopFitness += (dynamic_cast<Collect2RobotSensorsController*>
                                       (gWorld->getRobot(i)->getController())) -> getAvgPopFitness();
             gatheredGenomes += (dynamic_cast<Collect2RobotSensorsController*>
                                 (gWorld->getRobot(i)->getController())) ->_genomesList.size();
        }
        //std::cout << gWorld->getIterations() << " ";
        //<< (sumFitness  / gNumberOfRobots) / Collect2SharedData::gEvaluationTime

        // divided by two because each item gives 1 fitness point to both agents
        std::cout << sumFitness / 2 << std::endl;
	}

    if (gWorld->getIterations() == (gMaxIt - 1))
    {
        double sumFitness = 0.0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             sumFitness += (dynamic_cast<Collect2RobotSensorsController*>(gWorld->getRobot(i)
                                                    ->getController()))-> getFitness();

        }

        std::cout << "End fitness: " << sumFitness  / gNumberOfRobots
                  << " at it: " << gWorld->getIterations() << std::endl;
        if(Collect2RobotSensorsSharedData::gSaveGenome)
        {
            for (int i = 0 ; i != gNumberOfRobots ; i++ )
            {
                (dynamic_cast<Collect2RobotSensorsController*>(gWorld->getRobot(i)->getController()))
                        -> logGenome(Collect2RobotSensorsSharedData::gOutGenomeFile + std::to_string(i) + ".log");
            }
        }
    }
    switch (Collect2RobotSensorsSharedData::gFitness)
    {
        case 2:
            for(int i = 0; i < gNbOfPhysicalObjects;i++)
            {
                if(listCollected[i].size() >= 2)
                {
                    gPhysicalObjects[i]->isWalked(0); //Default agent for callback (unused callback)
                    for(auto it = listCollected[i].begin(); it != listCollected[i].end();it++)
                    {
                        dynamic_cast<Collect2RobotSensorsController*>(gWorld->getRobot((*it))
                            ->getController())->updateFitness(1.0);
                    }
                }
                listCollected[i].clear();
            }
            //updateDisplay();
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
            Collect2RobotSensorsSharedData::gTimeSeq[Collect2RobotSensorsSharedData::gTaskIdx])
    {
        //std::cout << "Task changed!" << std::endl;
        Collect2RobotSensorsSharedData::gFitness =
                Collect2RobotSensorsSharedData::gTaskSeq[Collect2RobotSensorsSharedData::gTaskIdx];
        Collect2RobotSensorsSharedData::gTaskIdx++;
    }
}

