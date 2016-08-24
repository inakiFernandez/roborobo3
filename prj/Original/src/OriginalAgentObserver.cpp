/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */


#include "Original/include/OriginalAgentObserver.h"
#include "World/World.h"
#include "Utilities/Misc.h"
#include "RoboroboMain/roborobo.h"
#include "Original/include/OriginalController.h"
#include <cmath>
#include "Original/include/OriginalWorldObserver.h"
#include <string>


OriginalAgentObserver::OriginalAgentObserver( RobotWorldModel *wm )
{
    _wm = (RobotWorldModel*)wm;
    /*if(OriginalSharedData::gEvolutionLogFile.compare("") != 0)
    {
        std::string sId = std::to_string(_wm->_id);
        _agentLog.open(OriginalSharedData::gEvolutionLogFile + std::string("-")
                       + sId + std::string(".log"));
    }
    else
    {
        std::cerr << "[ERROR] Empty evolution filename" << std::endl;
        exit(-1);
    }*/
}

OriginalAgentObserver::~OriginalAgentObserver()
{    
}

void OriginalAgentObserver::reset()
{
}

void OriginalAgentObserver::step()
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
        switch (OriginalSharedData::gFitness) {
        case 0:
            //Floreano's increment: nothing to do
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            gPhysicalObjects[targetIndex]->isWalked(_wm->getId());
            break;
        case 1:
            //Counting items
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            if((OriginalSharedData::gWithCollectColorEffector &&
                    isRightColorValue(dynamic_cast<OriginalController*>(gWorld->getRobot(_wm -> _id)->getController())->getColorEffector()))
                || !OriginalSharedData::gWithCollectColorEffector)
            {
                    gPhysicalObjects[targetIndex]->isWalked(_wm->getId());
                    dynamic_cast<OriginalController*>(gWorld->getRobot(_wm -> _id)->getController())->updateFitness(1.0);
            }

            break;
        case 2:
            //Recording the id of the robot for the item stepped upon
            targetIndex = targetIndex - gPhysicalObjectIndexStartOffset;
            if((OriginalSharedData::gWithCollectColorEffector &&
                    isRightColorValue(
                    dynamic_cast<OriginalController*>(gWorld->getRobot(_wm -> _id)->getController())->getColorEffector()))
                || !OriginalSharedData::gWithCollectColorEffector)
            {
                double color = dynamic_cast<OriginalController*>(gWorld->getRobot(_wm->_id)
                                                         ->getController())->getColorEffector();
                std::pair<int, double> idRobAndColor = std::make_pair(_wm->_id,color);
                dynamic_cast<OriginalWorldObserver*>(gWorld->getWorldObserver())->listCollected[targetIndex].push_back(idRobAndColor);
            }
            break;
        default:
            break;
        }
    }
    //TODO measuring agent distance to item and fitness every step
    //out to log file
    /*double distToItem = 0;
    double distToAgent= 0;
    double itX, itY;
    double agentX,agentY;
    double otherX,otherY;
    if((gNbOfPhysicalObjects >= 1)
            && (gInitialNumberOfRobots >= 2))
    {
        itX = gPhysicalObjects[0]->getXCenter();//getPosition().x;//
        itY = gPhysicalObjects[0]->getYCenter(); //getPosition().y;//

        agentX = _wm->getXReal();
        agentY = _wm->getYReal();

        distToItem = sqrt(pow(itX - agentX, 2.0) + pow(itY - agentY, 2.0));
        int otherId;
        if(_wm->_id == 0)
        {
            otherId = 1;
        }
        else
        {
            otherId = 0;
        }
        otherX = gWorld->getRobot(otherId)->getWorldModel()->getXReal();
        otherY = gWorld->getRobot(otherId)->getWorldModel()->getYReal();
        //std::cout << agentX << " " << agentY << "; " << otherX << " " << otherY << std::endl;
        distToAgent = sqrt(pow(otherX - agentX, 2.0) + pow(otherY - agentY,2.0));

        //column format: fitness distToItem dist2OtherAgent populationSize
        _agentLog
             << dynamic_cast<OriginalController*>(
                    gWorld->getRobot(_wm->_id)->getController())->getFitness()
             << " " << distToItem << " "
             << distToAgent << " "
             //<< dynamic_cast<OriginalController*>(gWorld->getRobot(_wm->_id)->getController())
             //   ->_genomesList.size()
             << "\n";
    }
    if (gWorld->getIterations() == (gMaxIt - 1))
    {
        //_agentLog.close();
    }*/
}
bool OriginalAgentObserver::isRightColorValue(double v)
{
    bool result = false;
    switch (OriginalSharedData::gFitness) {
        case 1:
            if((v >= 0.75) && (v <=1.0))
            {
                result = true;
                //std::cout << v << std::endl;
            }

        break;
        case 2:
            if((v <= -0.75) && (v >= -1.0))
                result = true;
        break;
        default:
        break;
    }
    return result;
}
