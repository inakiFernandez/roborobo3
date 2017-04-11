/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#ifndef ORIGINALAGENTOBSERVER_H
#define ORIGINALAGENTOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/AgentObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "Original/include/OriginalSharedData.h"

#include <iomanip>

class OriginalAgentObserver : public AgentObserver
{
	public:
		OriginalAgentObserver(RobotWorldModel *wm);
		~OriginalAgentObserver();

		void reset();
		void step();
        std::ofstream _agentLog;
        bool isRightColorValue(double v);
        bool isSameColor(double rCol, double objCol);
};

#endif

