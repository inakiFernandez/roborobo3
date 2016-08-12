/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#ifndef ORIGINALCONTROLLER_H
#define ORIGINALCONTROLLER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Utilities/Graphics.h"
#include "Controllers/Controller.h"
#include "WorldModels/RobotWorldModel.h"
#include "Original/include/OriginalAgentObserver.h"
#include <neuralnetworks/NeuralNetwork.h>
#include <iomanip>

#define RED 1
#define BLUE 2

using namespace Neural;
struct GC
{
    int robot_id;
    int gene_id;
    bool operator<(const GC &o)  const {
        if(gene_id == o.gene_id)
            return robot_id < o.robot_id;

        return gene_id < o.gene_id;
    }
    bool operator==(const GC &o)  const {
        return gene_id == o.gene_id && robot_id == o.robot_id;
    }
    friend std::ostream& operator<<(std::ostream& os, const GC& gene_clock);
};

class OriginalController : public Controller
{
private:
    int _iteration;
    int _birthdate; // iteration when this genome was initialized
    int _typeOfRobot; //RED or BLUE
    bool _withRobotType;

    std::string _nnType;
    std::vector<int> _nbHiddenNeuronsPerLayer;
    std::vector<int> _nbBiaisNeuronsPerLayer;
    NeuralNetwork* nn;

    //previous neural net todo
    //NeuralNetwork* previousNN;
    //forgetting measure todo
    //double forget();

    void createNN();
    
    void selectRandomGenome();
    void selectBestGenome();
    void selectRankBasedGenome();
    void selectTournament(double sp);
    void mutate(float sigma);

    void stepBehaviour();
    void stepEvolution();
    
    void broadcastGenome();
    void loadNewGenome();
    
    unsigned int computeRequiredNumberOfWeights();

    std::vector<double> _braitWeights;
    // evolutionary engine
    std::vector<double> _genome; // current genome in evaluation
    GC _genomeId;

    double _currentFitness;
    float _currentSigma;
    int _lifetime;
    

    // ANN
    double _minValue;
    double _maxValue;
    unsigned int _nbInputs;
    unsigned int _nbOutputs;
    unsigned int _nbHiddenLayers;
    std::vector<unsigned int>* _nbNeuronsPerHiddenLayer;
    
    void storeGenome(std::vector<double> genome, GC senderId, double fitness);
    void storeOwnGenome();
    void resetRobot();
    
public:

    OriginalController(RobotWorldModel *wm);
    ~OriginalController();

    void reset();
    void step();
    void updateFitness(double delta);
    int getBirthdate() { return _birthdate; }
    double getFitness(){ return _currentFitness;}
    int getRobotType()
    {
        return _typeOfRobot;
    }

    void readGenome(std::string s);
    void logGenome(std::string s);
    std::map<GC, std::vector<double> > _genomesList;
    std::map<GC, double > _fitnessList;

};


#endif

