/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#include "Collect2/include/Collect2AgentObserver.h"
#include "World/World.h"
#include "Utilities/Misc.h"
#include "RoboroboMain/roborobo.h"
#include "Collect2/include/Collect2Controller.h"
#include <cmath>
#include "Collect2/include/Collect2WorldObserver.h"
#include <string>


Collect2AgentObserver::Collect2AgentObserver( RobotWorldModel *wm )
{
    _wm = (RobotWorldModel*)wm;
}

Collect2AgentObserver::~Collect2AgentObserver()
{
}

void Collect2AgentObserver::reset()
{
}

void Collect2AgentObserver::step()
{
    // * send callback messages to objects touched or walked upon.
    
    // through distance sensors
    for( int i = 0 ; i < _wm->_cameraSensorsNb; i++)
    {
        int targetIndex = _wm->getObjectIdFromCameraSensor(i);
        
        if ( PhysicalObject::isInstanceOf(targetIndex))   // sensor ray bumped into a physical object
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
        switch (Collect2SharedData::gFitness) {
        case 0:
            //Floreano's increment: nothing to do
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            gPhysicalObjects[targetIndex]->isWalked(_wm->getId());
            break;
        case 1:
            //Counting items
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            gPhysicalObjects[targetIndex]->isWalked(_wm->getId());
            dynamic_cast<Collect2Controller*>(gWorld->getRobot(_wm -> _id)->getController())->updateFitness(1.0);
            break;
        case 2:
            //Recording the id of the robot for the item stepped upon
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            dynamic_cast<Collect2WorldObserver*>(gWorld->getWorldObserver())
                                ->listCollected[targetIndex].push_back(_wm->getId());
            break;
        default:
            break;
        }
    }
}
