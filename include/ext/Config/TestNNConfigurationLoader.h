/*
 * MedeaConfigurationLoader.h
 */

#ifndef TESTNNCONFIGURATIONLOADER_H
#define TESTNNCONFIGURATIONLOADER_H

#include "Config/ConfigurationLoader.h"


class TestNNConfigurationLoader : public ConfigurationLoader
{
	private:

	public:
		TestNNConfigurationLoader();
		~TestNNConfigurationLoader();

		WorldObserver *make_WorldObserver(World* wm) ;
		RobotWorldModel *make_RobotWorldModel();
		AgentObserver *make_AgentObserver(RobotWorldModel* wm) ;
		Controller *make_Controller(RobotWorldModel* wm) ;
};



#endif
