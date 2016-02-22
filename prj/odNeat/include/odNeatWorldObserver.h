#ifndef ODNEATWORLDOBSERVER_H
#define ODNEATWORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "odNeat/include/odNeatSharedData.h"
#include "odneatgc/helper.h"
#include <set>

class odNeatWorldObserver : public WorldObserver
{
	private:
		void updateEnvironment();
        void updateMonitoring();

	protected:
		int _generationCount;
		int _lifeIterationCount;
        /*//Counters for collisions in OdNeat
        int _countGeneClockCollisions;
        int _countGenes;*/
	public:
		odNeatWorldObserver(World *world);
		~odNeatWorldObserver();
        int _numberItems;
        //Counters for OdNeat
        /*std::set<int> usedGeneCounters;
        for collisions in OdNeat
         * void incrementCollisions();
        void incrementCounterGenes();*/
        //----------
		void reset();
		void step();

		int getLifeIterationCount() { return _lifeIterationCount; }

};

#endif
