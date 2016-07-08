#if defined PRJ_COLLECT2ROBOTSENSORS || !defined MODULAR

#include "Config/Collect2RobotSensorsConfigurationLoader.h"

#include "Collect2RobotSensors/include/Collect2RobotSensorsWorldObserver.h"
#include "Collect2RobotSensors/include/Collect2RobotSensorsAgentObserver.h"
#include "Collect2RobotSensors/include/Collect2RobotSensorsController.h"

#include "WorldModels/RobotWorldModel.h"

Collect2RobotSensorsConfigurationLoader::Collect2RobotSensorsConfigurationLoader()
{
}

Collect2RobotSensorsConfigurationLoader::~Collect2RobotSensorsConfigurationLoader()
{
	//nothing to do
}

WorldObserver* Collect2RobotSensorsConfigurationLoader::make_WorldObserver(World* wm)
{
	return new Collect2RobotSensorsWorldObserver(wm);
}

RobotWorldModel* Collect2RobotSensorsConfigurationLoader::make_RobotWorldModel()
{
	return new RobotWorldModel();
}

AgentObserver* Collect2RobotSensorsConfigurationLoader::make_AgentObserver(RobotWorldModel* wm)
{
	return new Collect2RobotSensorsAgentObserver(wm);
}

Controller* Collect2RobotSensorsConfigurationLoader::make_Controller(RobotWorldModel* wm)
{
	return new Collect2RobotSensorsController(wm);
}

#endif
