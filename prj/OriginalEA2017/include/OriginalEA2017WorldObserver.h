/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */





#ifndef ORIGINALEA2017WORLDOBSERVER_H
#define ORIGINALEA2017WORLDOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/Observer.h"
#include "Observers/WorldObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "OriginalEA2017/include/OriginalEA2017SharedData.h"
#include <fstream>
#include <sstream>

//class World;

class OriginalEA2017WorldObserver : public WorldObserver
{
	private:
		void updateEnvironment();
        void updateMonitoring();

        std::ofstream logItemFile;
        std::ofstream logItGatheredFile;
        std::ofstream logChangesColorFile;
        std::ofstream logGivenRewardFile;
        std::ofstream logRobotsPerItemFile;
        std::vector<int> itemCounts;
        double _averageReward;
        double _nbRobotsCorrect;
        int nbCollectedItems;


	protected:
		int _generationCount;
		int _lifeIterationCount;

	public:
		OriginalEA2017WorldObserver(World *world);
		~OriginalEA2017WorldObserver();
        //list for checking cooperative collecting
        std::vector< std::vector< std::pair<int,double > > > listCollected;

		void reset();
		void step();

		int getLifeIterationCount() { return _lifeIterationCount; }

};

#endif
