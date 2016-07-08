/*
 * MedeaConfigurationLoader.h
 */

#ifndef COLLECT2CONFIGURATIONLOADER_H
#define COLLECT2CONFIGURATIONLOADER_H

#include "Config/ConfigurationLoader.h"


class Collect2ConfigurationLoader : public ConfigurationLoader
{
	private:

	public:
		Collect2ConfigurationLoader();
		~Collect2ConfigurationLoader();

		WorldObserver *make_WorldObserver(World* wm) ;
		RobotWorldModel *make_RobotWorldModel();
		AgentObserver *make_AgentObserver(RobotWorldModel* wm) ;
		Controller *make_Controller(RobotWorldModel* wm) ;
};



#endif
