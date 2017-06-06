/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#ifndef ORIGINALEA2017AGENTOBSERVER_H
#define ORIGINALEA2017AGENTOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/AgentObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "OriginalEA2017/include/OriginalEA2017SharedData.h"

#include <iomanip>

class OriginalEA2017AgentObserver : public AgentObserver
{
	public:
		OriginalEA2017AgentObserver(RobotWorldModel *wm);
		~OriginalEA2017AgentObserver();

		void reset();
		void step();
        std::ofstream _agentLog;
        bool isRightColorValue(double v);
        bool isSameColor(double rCol, double objCol);
};

#endif

