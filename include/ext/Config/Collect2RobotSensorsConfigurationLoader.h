/*
 * MedeaConfigurationLoader.h
 */

#ifndef COLLECT2ROBOTSENSORSCONFIGURATIONLOADER_H
#define COLLECT2ROBOTSENSORSCONFIGURATIONLOADER_H

#include "Config/ConfigurationLoader.h"


class Collect2RobotSensorsConfigurationLoader : public ConfigurationLoader
{
	private:

	public:
		Collect2RobotSensorsConfigurationLoader();
		~Collect2RobotSensorsConfigurationLoader();

		WorldObserver *make_WorldObserver(World* wm) ;
		RobotWorldModel *make_RobotWorldModel();
		AgentObserver *make_AgentObserver(RobotWorldModel* wm) ;
		Controller *make_Controller(RobotWorldModel* wm) ;
};



#endif
