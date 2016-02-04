#if defined PRJ_TESTNN || !defined MODULAR

#include "Config/TestNNConfigurationLoader.h"

#include "TestNN/include/TestNNWorldObserver.h"
#include "TestNN/include/TestNNAgentObserver.h"
#include "TestNN/include/TestNNController.h"

#include "WorldModels/RobotWorldModel.h"

TestNNConfigurationLoader::TestNNConfigurationLoader()
{
}

TestNNConfigurationLoader::~TestNNConfigurationLoader()
{
	//nothing to do
}

WorldObserver* TestNNConfigurationLoader::make_WorldObserver(World* wm)
{
	return new TestNNWorldObserver(wm);
}

RobotWorldModel* TestNNConfigurationLoader::make_RobotWorldModel()
{
	return new RobotWorldModel();
}

AgentObserver* TestNNConfigurationLoader::make_AgentObserver(RobotWorldModel* wm)
{
	return new TestNNAgentObserver(wm);
}

Controller* TestNNConfigurationLoader::make_Controller(RobotWorldModel* wm)
{
	return new TestNNController(wm);
}

#endif
