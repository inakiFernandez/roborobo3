#if defined PRJ_ORIGINAL || !defined MODULAR

#include "Config/OriginalConfigurationLoader.h"

#include "Original/include/OriginalWorldObserver.h"
#include "Original/include/OriginalAgentObserver.h"
#include "Original/include/OriginalController.h"

#include "WorldModels/RobotWorldModel.h"

OriginalConfigurationLoader::OriginalConfigurationLoader()
{
}

OriginalConfigurationLoader::~OriginalConfigurationLoader()
{
	//nothing to do
}

WorldObserver* OriginalConfigurationLoader::make_WorldObserver(World* wm)
{
	return new OriginalWorldObserver(wm);
}

RobotWorldModel* OriginalConfigurationLoader::make_RobotWorldModel()
{
	return new RobotWorldModel();
}

AgentObserver* OriginalConfigurationLoader::make_AgentObserver(RobotWorldModel* wm)
{
	return new OriginalAgentObserver(wm);
}

Controller* OriginalConfigurationLoader::make_Controller(RobotWorldModel* wm)
{
	return new OriginalController(wm);
}

#endif
