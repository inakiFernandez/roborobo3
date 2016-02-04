/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#ifndef TESTNNAGENTOBSERVER_H
#define TESTNNAGENTOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/AgentObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "TestNN/include/TestNNSharedData.h"

#include <iomanip>

class TestNNAgentObserver : public AgentObserver
{
	public:
		TestNNAgentObserver(RobotWorldModel *wm);
		~TestNNAgentObserver();

		void reset();
		void step();

};

#endif

