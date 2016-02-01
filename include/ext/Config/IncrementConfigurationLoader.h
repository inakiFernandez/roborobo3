/*
 * MedeaConfigurationLoader.h
 */

#ifndef INCREMENTCONFIGURATIONLOADER_H
#define INCREMENTCONFIGURATIONLOADER_H

#include "Config/ConfigurationLoader.h"


class IncrementConfigurationLoader : public ConfigurationLoader
{
	private:

	public:
		IncrementConfigurationLoader();
		~IncrementConfigurationLoader();

		WorldObserver *make_WorldObserver(World* wm) ;
		RobotWorldModel *make_RobotWorldModel();
		AgentObserver *make_AgentObserver(RobotWorldModel* wm) ;
		Controller *make_Controller(RobotWorldModel* wm) ;
};



#endif
