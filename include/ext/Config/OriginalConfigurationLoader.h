/*
 * MedeaConfigurationLoader.h
 */

#ifndef ORIGINALCONFIGURATIONLOADER_H
#define ORIGINALCONFIGURATIONLOADER_H

#include "Config/ConfigurationLoader.h"


class OriginalConfigurationLoader : public ConfigurationLoader
{
	private:

	public:
		OriginalConfigurationLoader();
		~OriginalConfigurationLoader();

		WorldObserver *make_WorldObserver(World* wm) ;
		RobotWorldModel *make_RobotWorldModel();
		AgentObserver *make_AgentObserver(RobotWorldModel* wm) ;
		Controller *make_Controller(RobotWorldModel* wm) ;
};



#endif
