/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */





#ifndef ORIGINALWORLDOBSERVER_H
#define ORIGINALWORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "Original/include/OriginalSharedData.h"

//class World;

class OriginalWorldObserver : public WorldObserver
{
	private:
		void updateEnvironment();
        void updateMonitoring();

	protected:
		int _generationCount;
		int _lifeIterationCount;

	public:
		OriginalWorldObserver(World *world);
		~OriginalWorldObserver();
        //list for checking cooperative collecting
        std::vector< std::vector< int > > listCollected;

		void reset();
		void step();

		int getLifeIterationCount() { return _lifeIterationCount; }

};

#endif
