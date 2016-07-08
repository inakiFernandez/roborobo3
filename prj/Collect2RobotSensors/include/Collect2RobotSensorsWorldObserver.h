/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */





#ifndef COLLECT2ROBOTSENSORSWORLDOBSERVER_H
#define COLLECT2ROBOTSENSORSWORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "Collect2RobotSensors/include/Collect2RobotSensorsSharedData.h"

//class World;

class Collect2RobotSensorsWorldObserver : public WorldObserver
{
	private:
		void updateEnvironment();
        void updateMonitoring();

	protected:
		int _generationCount;
		int _lifeIterationCount;

	public:
		Collect2RobotSensorsWorldObserver(World *world);
		~Collect2RobotSensorsWorldObserver();
        //list for checking cooperative collecting
        std::vector< std::vector< int > > listCollected;

		void reset();
		void step();

		int getLifeIterationCount() { return _lifeIterationCount; }

};

#endif
