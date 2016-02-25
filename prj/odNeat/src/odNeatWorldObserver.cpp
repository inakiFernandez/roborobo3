/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */

#include "Observers/AgentObserver.h"
#include "Observers/WorldObserver.h"
#include "odNeat/include/odNeatWorldObserver.h"
#include "odNeat/include/odNeatController.h"
#include "World/World.h"


odNeatWorldObserver::odNeatWorldObserver( World* world ) : WorldObserver( world )
{
    _world = world;

	// ==== loading project-specific properties

	gProperties.checkAndGetPropertyValue("gSigmaRef",&odNeatSharedData::gSigmaRef,true);
    gProperties.checkAndGetPropertyValue("gPopSize",&odNeatSharedData::gPopulationSize,true);
    gProperties.checkAndGetPropertyValue("gClearPopulation",&odNeatSharedData::gClearPopulation,true);
    gProperties.checkAndGetPropertyValue("gStoreOwn",&odNeatSharedData::gStoreOwn,true);
    gProperties.checkAndGetPropertyValue("gSelectionMethod",&odNeatSharedData::gSelectionMethod,true);

    gProperties.checkAndGetPropertyValue("gCommunicationBySensors",&odNeatSharedData::gCommunicationBySensors,true);
    gProperties.checkAndGetPropertyValue("gFitness",&odNeatSharedData::gFitness,true);
    gProperties.checkAndGetPropertyValue("gEvaluationTime",&odNeatSharedData::gEvaluationTime,true);

    gProperties.checkAndGetPropertyValue("gOutGenomeFile",&odNeatSharedData::gOutGenomeFile,true);


    //OdNeat----------------------------------
    //Mutations
    gProperties.checkAndGetPropertyValue("mutate_only_prob",&Helper::mutateProb,true);
    gProperties.checkAndGetPropertyValue("mutate_link_weights_prob",&Helper::mutateLinkWeightsProb,true);
    gProperties.checkAndGetPropertyValue("mutate_add_node_prob",&Helper::mutateAddNodeProb,true);
    gProperties.checkAndGetPropertyValue("mutate_add_link_prob",&Helper::mutateAddLinkProb,true);
    gProperties.checkAndGetPropertyValue("mate_only_prob",&Helper::mateOnlyProb,true);
    gProperties.checkAndGetPropertyValue("recur_only_prob",&Helper::recurOnlyProb,true);
    gProperties.checkAndGetPropertyValue("newstructure_tries",&Helper::newStructureTries,true);

    //Evaluation, tabu and thresholds
    gProperties.checkAndGetPropertyValue("gMaturationPeriod",&odNeatSharedData::gMaturationPeriod,true);
    gProperties.checkAndGetPropertyValue("gCompatThreshold",&odNeatSharedData::gCompatThreshold,true);
    gProperties.checkAndGetPropertyValue("gTabuThreshold",&odNeatSharedData::gTabuThreshold,true);
    gProperties.checkAndGetPropertyValue("gTabuTimeout",&odNeatSharedData::gTabuTimeout,true);

    //Energy from items
    gProperties.checkAndGetPropertyValue("gEnergyItemValue",&odNeatSharedData::gEnergyItemValue,true);
    gProperties.checkAndGetPropertyValue("gEnergyConsumption",&odNeatSharedData::gEnergyConsumption,true);
    //----------------------------------
	// * iteration and generation counters

	_lifeIterationCount = -1;
	_generationCount = -1;

    _numberItems = 0;
}

odNeatWorldObserver::~odNeatWorldObserver()
{
	// nothing to do.
}

void odNeatWorldObserver::reset()
{
	// nothing to do.
}

void odNeatWorldObserver::step()
{
    _lifeIterationCount++;
    
    updateMonitoring();

    if( _lifeIterationCount >= odNeatSharedData::gEvaluationTime ) // switch to next generation.
	{
        // update iterations and generations counters
        _lifeIterationCount = 0;
        _generationCount++;
    }
}


void odNeatWorldObserver::updateEnvironment()
{
}

void odNeatWorldObserver::updateMonitoring()
{
    // * Log every N iterations
    
    if( _lifeIterationCount >= odNeatSharedData::gEvaluationTime )
	{
        // Logging here
        double sumFitness = 0.0; int sumLinks = 0; int numSpecies = 0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {
             sumFitness += (dynamic_cast<odNeatController*>(gWorld->getRobot(i)->getController()))
                     -> getFitness()/odNeatSharedData::gMaxEnergy;//Normalize by max energy
             sumLinks += (dynamic_cast<odNeatController*>(gWorld->getRobot(i)->getController()))
                     ->_genome->genes.size();
             numSpecies += (dynamic_cast<odNeatController*>(gWorld->getRobot(i)->getController()))
                     ->_pop.size();
        }
        std::cout << "It: " << gWorld->getIterations() << std::endl;
        std::cout << "Fitness: " <<  sumFitness  / gNumberOfRobots << std::endl;
        std::cout << "Links: " << (double)sumLinks/ gNumberOfRobots << std::endl;
        std::cout << "Species: " << (double)numSpecies/ gNumberOfRobots << std::endl;
        if (odNeatSharedData::gFitness == 1)
        {
            std::cout << "Number of items: " << _numberItems << std::endl;
        }
        _numberItems = 0;
	}
    if (gWorld->getIterations() == (gMaxIt - 1))
    {
        double sumFitness = 0.0; int sumItems = 0;
        for ( int i = 0 ; i != gNumberOfRobots ; i++ )
        {

             sumFitness += (dynamic_cast<odNeatController*>(gWorld->getRobot(i)->getController()))
                     -> getFitness();
             sumItems += (dynamic_cast<odNeatController*>(gWorld->getRobot(i)->getController()))
                     -> getItems();

        }

        std::cout << "End fitness: " << sumFitness/odNeatSharedData::gMaxEnergy / gNumberOfRobots
                  << " at it: " << gWorld->getIterations() << std::endl;
        if (odNeatSharedData::gFitness == 1)
        {
            //std::cout << "Number of items: " << sumItems << std::endl;
        }
        if(odNeatSharedData::gSaveGenome)
        {
            for (int i = 0 ; i != gNumberOfRobots ; i++ )
            {
                (dynamic_cast<odNeatController*>(gWorld->getRobot(i)->getController()))
                        -> logGenome(odNeatSharedData::gOutGenomeFile + std::to_string(i) + ".log");
            }
        }
    }
}
