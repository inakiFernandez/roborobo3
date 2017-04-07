/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#include "Collect2RobotSensors/include/Collect2RobotSensorsAgentObserver.h"
#include "World/World.h"
#include "Utilities/Misc.h"
#include "RoboroboMain/roborobo.h"
#include "Collect2RobotSensors/include/Collect2RobotSensorsController.h"
#include <cmath>
#include "Collect2RobotSensors/include/Collect2RobotSensorsWorldObserver.h"
#include <string>


Collect2RobotSensorsAgentObserver::Collect2RobotSensorsAgentObserver( RobotWorldModel *wm )
{
    _wm = (RobotWorldModel*)wm;
}

Collect2RobotSensorsAgentObserver::~Collect2RobotSensorsAgentObserver()
{
}

void Collect2RobotSensorsAgentObserver::reset()
{
}

void Collect2RobotSensorsAgentObserver::step()
{
    // * send callback messages to objects touched or walked upon.
    
    // through distance sensors
    for( int i = 0 ; i < _wm->_cameraSensorsNb; i++)
    {
        int targetIndex = _wm->getObjectIdFromCameraSensor(i);

        // sensor ray bumped into a physical object
        if ( PhysicalObject::isInstanceOf(targetIndex))
        {
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            gPhysicalObjects[targetIndex]->isTouched(_wm->getId());
        }
    }
    
    // through floor sensor
    int targetIndex = _wm->getGroundSensorValue();
    // ground sensor is upon a physical object (OR: on a place marked with
    //this physical object footprint, cf. groundsensorvalues image)
    if(PhysicalObject::isInstanceOf(targetIndex))
    {
        switch (Collect2RobotSensorsSharedData::gFitness) {
        case 0:
            //Floreano's increment: nothing to do
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            gPhysicalObjects[targetIndex]->isWalked(_wm->getId());
            break;
        case 1:
            //Counting items
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            gPhysicalObjects[targetIndex]->isWalked(_wm->getId());
            dynamic_cast<Collect2RobotSensorsController*>(gWorld->getRobot(_wm -> _id)->getController())->updateFitness(1.0);
            break;
        case 2:
            //Recording the id of the robot for the item stepped upon
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            dynamic_cast<Collect2RobotSensorsWorldObserver*>(gWorld->getWorldObserver())
                                ->listCollected[targetIndex].push_back(_wm->getId());
            break;
        default:
            break;
        }
    }
}
