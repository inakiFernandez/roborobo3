/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#include "TestNN/include/TestNNAgentObserver.h"
#include "World/World.h"
#include "Utilities/Misc.h"
#include "RoboroboMain/roborobo.h"
#include "TestNN/include/TestNNController.h"
#include <cmath>
#include "TestNN/include/TestNNWorldObserver.h"
#include <string>


TestNNAgentObserver::TestNNAgentObserver( RobotWorldModel *wm )
{
    _wm = (RobotWorldModel*)wm;
}

TestNNAgentObserver::~TestNNAgentObserver()
{
}

void TestNNAgentObserver::reset()
{
}

void TestNNAgentObserver::step()
{
    // * send callback messages to objects touched or walked upon.
    
    // through distance sensors
    for( int i = 0 ; i < _wm->_cameraSensorsNb; i++)
    {
        int targetIndex = _wm->getObjectIdFromCameraSensor(i);
        
        if ( PhysicalObject::isInstanceOf(targetIndex) )   // sensor ray bumped into a physical object
        {
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            gPhysicalObjects[targetIndex]->isTouched(_wm->getId());
        }
    }
    
    // through floor sensor
    int targetIndex = _wm->getGroundSensorValue();
    // ground sensor is upon a physical object (OR: on a place marked with this physical object footprint, cf. groundsensorvalues image)
    if ( PhysicalObject::isInstanceOf(targetIndex) )
    {
        targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
        gPhysicalObjects[targetIndex]->isWalked(_wm->getId());
        switch (TestNNSharedData::gFitness) {
        case 0:
            //Floreano's increment: nothing to do
            break;
        case 1:
            //Counting items
            dynamic_cast<TestNNController*>(gWorld->getRobot(_wm -> _id)->getController())->updateFitness(1.0);
            break;
        default:
            break;
        }
    }
}
