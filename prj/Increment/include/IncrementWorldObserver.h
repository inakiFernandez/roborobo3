/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */





#ifndef INCREMENTWORLDOBSERVER_H
#define INCREMENTWORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "Increment/include/IncrementSharedData.h"

//class World;

class IncrementWorldObserver : public WorldObserver
{
	private:
		void updateEnvironment();
        void updateMonitoring();

	protected:
		int _generationCount;
		int _lifeIterationCount;

	public:
		IncrementWorldObserver(World *world);
		~IncrementWorldObserver();

		void reset();
		void step();

		int getLifeIterationCount() { return _lifeIterationCount; }

};

#endif
