/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Observers/AgentObserver.h"
#include "Observers/WorldObserver.h"
#include "TestNN/include/TestNNWorldObserver.h"
#include "TestNN/include/TestNNController.h"
#include "World/World.h"


TestNNWorldObserver::TestNNWorldObserver( World* world ) : WorldObserver( world )
{
    _world = world;

	// ==== loading project-specific properties


    gProperties.checkAndGetPropertyValue("gFitness",&TestNNSharedData::gFitness,true);
    gProperties.checkAndGetPropertyValue("gEvaluationTime",&TestNNSharedData::gEvaluationTime,true);

    gProperties.checkAndGetPropertyValue("gNbHiddenLayers",&TestNNSharedData::gNbHiddenLayers,true);
	gProperties.checkAndGetPropertyValue("gNbNeuronsPerHiddenLayer",&TestNNSharedData::gNbNeuronsPerHiddenLayer,true);
	gProperties.checkAndGetPropertyValue("gNeuronWeightRange",&TestNNSharedData::gNeuronWeightRange,true);
    gProperties.checkAndGetPropertyValue("gWithBias",&TestNNSharedData::gWithBias,true);

    gProperties.checkAndGetPropertyValue("gOutGenomeFile",&TestNNSharedData::gOutGenomeFile,true);

	// * iteration and generation counters

	_lifeIterationCount = -1;
	_generationCount = -1;

}

TestNNWorldObserver::~TestNNWorldObserver()
{
	// nothing to do.
}

void TestNNWorldObserver::reset()
{
	// nothing to do.
}

void TestNNWorldObserver::step()
{
    _lifeIterationCount++;
    
    updateMonitoring();

    if( _lifeIterationCount >= TestNNSharedData::gEvaluationTime ) // switch to next generation.
	{
        // update iterations and generations counters
        _lifeIterationCount = 0;
        _generationCount++;
    }
}


void TestNNWorldObserver::updateEnvironment()
{
}

void TestNNWorldObserver::updateMonitoring()
{
    // * Log at end of each generation
    
    if( _lifeIterationCount >= TestNNSharedData::gEvaluationTime ) // end of generation.
	{
		if ( gVerbose )
		{
            std::cout << "[gen:" << (gWorld->getIterations()/TestNNSharedData::gEvaluationTime)
                      << "]\n";
		}
        // Logging here
        double sumFitness = 0.0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             sumFitness += (dynamic_cast<TestNNController*>(gWorld->getRobot(i)->getController()))
                     -> getFitness();
        }
        std::cout << "It: " << gWorld->getIterations() << std::endl;
        std::cout << sumFitness  / gNumberOfRobots << std::endl;
	}
    //if (gWorld->getIterations() % 100 == 0)
        //std::cout << "It: " << gWorld->getIterations() << std::endl;
    if (gWorld->getIterations() == (gMaxIt - 1))
    {
        double sumFitness = 0.0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             sumFitness += (dynamic_cast<TestNNController*>(gWorld->getRobot(i)->getController()))
                     -> getFitness();

        }

        std::cout << "End fitness: " << sumFitness  / gNumberOfRobots
                  << " at it: " << gWorld->getIterations() << std::endl;

    }
}

