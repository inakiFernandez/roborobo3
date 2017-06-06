/*
 * MedeaConfigurationLoader.h
 */

#ifndef ORIGINALEA2017CONFIGURATIONLOADER_H
#define ORIGINALEA2017CONFIGURATIONLOADER_H

#include "Config/ConfigurationLoader.h"


class OriginalEA2017ConfigurationLoader : public ConfigurationLoader
{
	private:

	public:
		OriginalEA2017ConfigurationLoader();
		~OriginalEA2017ConfigurationLoader();

		WorldObserver *make_WorldObserver(World* wm) ;
		RobotWorldModel *make_RobotWorldModel();
		AgentObserver *make_AgentObserver(RobotWorldModel* wm) ;
		Controller *make_Controller(RobotWorldModel* wm) ;
};



#endif
