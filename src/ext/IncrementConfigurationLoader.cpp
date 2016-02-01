#if defined PRJ_INCREMENT || !defined MODULAR

#include "Config/IncrementConfigurationLoader.h"

#include "Increment/include/IncrementWorldObserver.h"
#include "Increment/include/IncrementAgentObserver.h"
#include "Increment/include/IncrementController.h"

#include "WorldModels/RobotWorldModel.h"

IncrementConfigurationLoader::IncrementConfigurationLoader()
{
}

IncrementConfigurationLoader::~IncrementConfigurationLoader()
{
	//nothing to do
}

WorldObserver* IncrementConfigurationLoader::make_WorldObserver(World* wm)
{
	return new IncrementWorldObserver(wm);
}

RobotWorldModel* IncrementConfigurationLoader::make_RobotWorldModel()
{
	return new RobotWorldModel();
}

AgentObserver* IncrementConfigurationLoader::make_AgentObserver(RobotWorldModel* wm)
{
	return new IncrementAgentObserver(wm);
}

Controller* IncrementConfigurationLoader::make_Controller(RobotWorldModel* wm)
{
	return new IncrementController(wm);
}

#endif
