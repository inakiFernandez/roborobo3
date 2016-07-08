#if defined PRJ_COLLECT2 || !defined MODULAR

#include "Config/Collect2ConfigurationLoader.h"

#include "Collect2/include/Collect2WorldObserver.h"
#include "Collect2/include/Collect2AgentObserver.h"
#include "Collect2/include/Collect2Controller.h"

#include "WorldModels/RobotWorldModel.h"

Collect2ConfigurationLoader::Collect2ConfigurationLoader()
{
}

Collect2ConfigurationLoader::~Collect2ConfigurationLoader()
{
	//nothing to do
}

WorldObserver* Collect2ConfigurationLoader::make_WorldObserver(World* wm)
{
	return new Collect2WorldObserver(wm);
}

RobotWorldModel* Collect2ConfigurationLoader::make_RobotWorldModel()
{
	return new RobotWorldModel();
}

AgentObserver* Collect2ConfigurationLoader::make_AgentObserver(RobotWorldModel* wm)
{
	return new Collect2AgentObserver(wm);
}

Controller* Collect2ConfigurationLoader::make_Controller(RobotWorldModel* wm)
{
	return new Collect2Controller(wm);
}

#endif
