#if defined PRJ_ORIGINALEA2017 || !defined MODULAR

#include "Config/OriginalEA2017ConfigurationLoader.h"

#include "OriginalEA2017/include/OriginalEA2017WorldObserver.h"
#include "OriginalEA2017/include/OriginalEA2017AgentObserver.h"
#include "OriginalEA2017/include/OriginalEA2017Controller.h"

#include "WorldModels/RobotWorldModel.h"

OriginalEA2017ConfigurationLoader::OriginalEA2017ConfigurationLoader()
{
}

OriginalEA2017ConfigurationLoader::~OriginalEA2017ConfigurationLoader()
{
	//nothing to do
}

WorldObserver* OriginalEA2017ConfigurationLoader::make_WorldObserver(World* wm)
{
	return new OriginalEA2017WorldObserver(wm);
}

RobotWorldModel* OriginalEA2017ConfigurationLoader::make_RobotWorldModel()
{
	return new RobotWorldModel();
}

AgentObserver* OriginalEA2017ConfigurationLoader::make_AgentObserver(RobotWorldModel* wm)
{
	return new OriginalEA2017AgentObserver(wm);
}

Controller* OriginalEA2017ConfigurationLoader::make_Controller(RobotWorldModel* wm)
{
	return new OriginalEA2017Controller(wm);
}

#endif
