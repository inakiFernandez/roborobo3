/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Observers/AgentObserver.h"
#include "Observers/WorldObserver.h"
#include "Increment/include/IncrementWorldObserver.h"
#include "Increment/include/IncrementController.h"
#include "World/World.h"


IncrementWorldObserver::IncrementWorldObserver( World* world ) : WorldObserver( world )
{
    _world = world;

	// ==== loading project-specific properties

	gProperties.checkAndGetPropertyValue("gSigmaRef",&IncrementSharedData::gSigmaRef,true);
    gProperties.checkAndGetPropertyValue("gPopSize",&IncrementSharedData::gPopulationSize,true);
    gProperties.checkAndGetPropertyValue("gClearPopulation",&IncrementSharedData::gClearPopulation,true);
    gProperties.checkAndGetPropertyValue("gStoreOwn",&IncrementSharedData::gStoreOwn,true);
    gProperties.checkAndGetPropertyValue("gSelectionMethod",&IncrementSharedData::gSelectionMethod,true);

    gProperties.checkAndGetPropertyValue("gFitness",&IncrementSharedData::gFitness,true);
    gProperties.checkAndGetPropertyValue("gEvaluationTime",&IncrementSharedData::gEvaluationTime,true);

    gProperties.checkAndGetPropertyValue("gNbHiddenLayers",&IncrementSharedData::gNbHiddenLayers,true);
	gProperties.checkAndGetPropertyValue("gNbNeuronsPerHiddenLayer",&IncrementSharedData::gNbNeuronsPerHiddenLayer,true);
	gProperties.checkAndGetPropertyValue("gNeuronWeightRange",&IncrementSharedData::gNeuronWeightRange,true);
    gProperties.checkAndGetPropertyValue("gWithBias",&IncrementSharedData::gWithBias,true);


	// * iteration and generation counters

	_lifeIterationCount = -1;
	_generationCount = -1;

}

IncrementWorldObserver::~IncrementWorldObserver()
{
	// nothing to do.
}

void IncrementWorldObserver::reset()
{
	// nothing to do.
}

void IncrementWorldObserver::step()
{
    _lifeIterationCount++;
    
    updateMonitoring();

    if( _lifeIterationCount >= IncrementSharedData::gEvaluationTime ) // switch to next generation.
	{
        // update iterations and generations counters
        _lifeIterationCount = 0;
        _generationCount++;
    }
}


void IncrementWorldObserver::updateEnvironment()
{
}

void IncrementWorldObserver::updateMonitoring()
{
    // * Log at end of each generation
    
    if( _lifeIterationCount >= IncrementSharedData::gEvaluationTime ) // end of generation.
	{
		// * monitoring: count number of active agents.
        
		int activeCount = 0;
		for ( int i = 0 ; i != gNumberOfRobots ; i++ )
		{
            if ( (dynamic_cast<IncrementController*>(gWorld->getRobot(i)->getController()))
                 ->getWorldModel()->isAlive() == true )
				activeCount++;
		}
        
		if ( gVerbose )
		{
            std::cout << "[gen:" << (gWorld->getIterations()/IncrementSharedData::gEvaluationTime)
                      << ";pop:" << activeCount << "]\n";
		}
        // Logging here TODO
        double totalFitness = 0.0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             totalFitness += (dynamic_cast<IncrementController*>(gWorld->getRobot(i)->getController()))
                     -> getFitness();
        }
        std::cout << totalFitness  / gNumberOfRobots << std::endl;
	}
    if (gWorld->getIterations() % 500 == 0)
    {
        double totalFitness = 0.0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             totalFitness += (dynamic_cast<IncrementController*>(gWorld->getRobot(i)->getController()))
                     -> getFitness();
        }
        //std::cout << totalFitness  / gNumberOfRobots << std::endl;
        //std::cout << "It: " << gWorld->getIterations() << std::endl;
    }
}

