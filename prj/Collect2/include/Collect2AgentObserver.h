/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#ifndef COLLECT2AGENTOBSERVER_H
#define COLLECT2AGENTOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/AgentObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "Collect2/include/Collect2SharedData.h"

#include <iomanip>

class Collect2AgentObserver : public AgentObserver
{
	public:
		Collect2AgentObserver(RobotWorldModel *wm);
		~Collect2AgentObserver();

		void reset();
		void step();

};

#endif

