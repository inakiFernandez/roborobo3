/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#ifndef TESTNNCONTROLLER_H
#define TESTNNCONTROLLER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Utilities/Graphics.h"
#include "Controllers/Controller.h"
#include "WorldModels/RobotWorldModel.h"
#include "TestNN/include/TestNNAgentObserver.h"
#include <neuralnetworks/NeuralNetwork.h>

#include <iomanip>

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

class TestNNController : public Controller
{
private:
    int _iteration;
    int _birthdate; // evaluation when this controller was initialized.
    double _currentFitness;
    double _recordedFitness;

    std::vector<double> _parameters;
    std::string _nnType;
    std::vector<int> _nbHiddenNeuronsPerLayer;
    std::vector<int> _nbBiaisNeuronsPerLayer;
    NeuralNetwork* nn;

    void createNN();
    
    void mutate(float sigma);

    void stepBehaviour();
    
    unsigned int computeRequiredNumberOfWeights();

    std::vector<double> _genome; // current genome in evaluation
    GC _genomeId;
    std::vector<double> _currentGenome;
    int _lifetime;
    
    // ANN
    double _minValue;
    double _maxValue;
    unsigned int _nbInputs;
    unsigned int _nbOutputs;
    unsigned int _nbHiddenLayers;
    std::vector<unsigned int>* _nbNeuronsPerHiddenLayer;
    
    void resetRobot();
    
public:

    TestNNController(RobotWorldModel *wm);
    ~TestNNController();

    void reset();
    void step();
    void updateFitness(double delta);
    int getBirthdate() { return _birthdate; }
    double getFitness(){ return _currentFitness;}

    void loadTestGenome(std::string filename);

};


#endif

