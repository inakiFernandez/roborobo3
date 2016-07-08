/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#ifndef COLLECT2ROBOTSENSORSAGENTOBSERVER_H
#define COLLECT2ROBOTSENSORSAGENTOBSERVER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Observers/AgentObserver.h"
#include "WorldModels/RobotWorldModel.h"
#include "Collect2RobotSensors/include/Collect2RobotSensorsSharedData.h"

#include <iomanip>

class Collect2RobotSensorsAgentObserver : public AgentObserver
{
	public:
		Collect2RobotSensorsAgentObserver(RobotWorldModel *wm);
		~Collect2RobotSensorsAgentObserver();

		void reset();
		void step();

};

#endif

