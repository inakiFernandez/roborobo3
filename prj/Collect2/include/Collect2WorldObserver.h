/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */





#ifndef COLLECT2WORLDOBSERVER_H
#define COLLECT2WORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "Collect2/include/Collect2SharedData.h"

//class World;

class Collect2WorldObserver : public WorldObserver
{
	private:
		void updateEnvironment();
        void updateMonitoring();

	protected:
		int _generationCount;
		int _lifeIterationCount;

	public:
		Collect2WorldObserver(World *world);
		~Collect2WorldObserver();
        //list for checking cooperative collecting
        std::vector< std::vector< int > > listCollected;

		void reset();
		void step();

		int getLifeIterationCount() { return _lifeIterationCount; }

};

#endif
