/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */





#ifndef TESTNNWORLDOBSERVER_H
#define TESTNNWORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "TestNN/include/TestNNSharedData.h"
#include "TestNN/include/TestNNSharedData.h"

//class World;

class TestNNWorldObserver : public WorldObserver
{
	private:
		void updateEnvironment();
        void updateMonitoring();

	protected:
		int _generationCount;
		int _lifeIterationCount;

	public:
		TestNNWorldObserver(World *world);
		~TestNNWorldObserver();

		void reset();
		void step();

		int getLifeIterationCount() { return _lifeIterationCount; }

};

#endif
