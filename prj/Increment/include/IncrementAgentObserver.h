/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#ifndef INCREMENTAGENTOBSERVER_H
#define INCREMENTAGENTOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/AgentObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "Increment/include/IncrementSharedData.h"

#include <iomanip>

class IncrementAgentObserver : public AgentObserver
{
	public:
		IncrementAgentObserver(RobotWorldModel *wm);
		~IncrementAgentObserver();

		void reset();
		void step();

};

#endif

