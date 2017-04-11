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
#include <fstream>
#include <sstream>

//class World;

class OriginalWorldObserver : public WorldObserver
{
	private:
		void updateEnvironment();
        void updateMonitoring();

        std::ofstream logItemFile;
        std::ofstream logItGatheredFile;
        std::ofstream logChangesColorFile;

        std::vector<int> itemCounts;


	protected:
		int _generationCount;
		int _lifeIterationCount;

	public:
		OriginalWorldObserver(World *world);
		~OriginalWorldObserver();
        //list for checking cooperative collecting
        std::vector< std::vector< std::pair<int,double > > > listCollected;

		void reset();
		void step();

		int getLifeIterationCount() { return _lifeIterationCount; }

};

#endif
