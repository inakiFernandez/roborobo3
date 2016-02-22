#include "odNeat/include/odNeatAgentObserver.h"
#include "World/World.h"
#include "Utilities/Misc.h"
#include "RoboroboMain/roborobo.h"
#include "odNeat/include/odNeatController.h"
#include <cmath>
#include "odNeat/include/odNeatWorldObserver.h"
#include <string>


odNeatAgentObserver::odNeatAgentObserver( RobotWorldModel *wm )
{
    _wm = (RobotWorldModel*)wm;
}

odNeatAgentObserver::~odNeatAgentObserver()
{
}

void odNeatAgentObserver::reset()
{
}

void odNeatAgentObserver::step()
{
    // * send callback messages to objects touched or walked upon.
    
    // through distance sensors
    for( int i = 0 ; i < _wm->_cameraSensorsNb; i++)
    {
        int targetIndex = _wm->getObjectIdFromCameraSensor(i);
        // sensor ray bumped into a physical object
        if ( PhysicalObject::isInstanceOf(targetIndex) )
        {
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            gPhysicalObjects[targetIndex]->isTouched(_wm->getId());
        }
    }
    
    // through floor sensor
    int targetIndex = _wm->getGroundSensorValue();
    // ground sensor is upon a physical object
    //(OR: on a place marked with this physical object footprint, cf. groundsensorvalues image)
    if ( PhysicalObject::isInstanceOf(targetIndex) )
    {
        targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
        gPhysicalObjects[targetIndex]->isWalked(_wm->getId());
        switch (odNeatSharedData::gFitness) {
        case 0:
            //Floreano's increment: nothing to do
            break;
        case 1:
            //Counting items
            dynamic_cast<odNeatController*>(gWorld->getRobot(_wm -> _id)->
                                            getController())->pickItem();
            dynamic_cast<odNeatWorldObserver*>(gWorld->getWorldObserver())->_numberItems++;
            break;
        default:
            break;
        }
    }
}
