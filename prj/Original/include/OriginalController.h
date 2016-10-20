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

#include <iomanip>
#include <set>
#include <tuple>
#include <map>

#include <neuralnetworks/NeuralNetwork.h>
#include <odneatgc/network.h>
#include <odneatgc/genome.h>

using namespace Neural;
using namespace ODNEATGC;

/*struct GC
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
};*/

class OriginalController : public Controller
{
private:
    int _iteration;
    int _birthdate; // iteration when this genome was initialized

    bool _withCollectColorEffector;


    void selectRandomGenome();
    void selectBestGenome();
    void selectRankBasedGenome();
    void selectTournament(double sp);
    void mutateFixed(float sigma);
    void mutateEvoTopo(float sigma);
    void mutate(float sigma);

    void stepBehaviour();
    void stepEvolution();
    
    void broadcastGenome();
    void loadNewGenome();
    
    unsigned int computeRequiredNumberOfWeights();

    std::vector<double> _braitWeights;
    // evolutionary engine
    std::vector<double> _genomeF; // current genome in evaluation (fixed topo)
    Genome *_genome; // current genome in evaluation (evo topo)

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
    
    std::string _nnType;
    std::vector<int> _nbHiddenNeuronsPerLayer;
    std::vector<int> _nbBiaisNeuronsPerLayer;

    NeuralNetwork* nnF;
    Network *nn;

    //previous neural net todo
    //NeuralNetwork* previousNN;
    //forgetting measure todo
    //double forget();

    void createNN();

    //gene (link) and neuron local counters
    int _g_count; int _n_count;

    void storeGenome(Genome* genome, GC senderId, double fitness, int nCount, int gCount);
    void storeGenomeF(std::vector<double> genome, GC senderId, double fitness);

    void storeOwnGenome();
    void storeOwnGenomeF();

    void resetRobot();


    int roundDown(int numToRound, int multiple);
    std::vector<double> getBraitenberg();

public:

    OriginalController(RobotWorldModel *wm);
    ~OriginalController();

    void reset();
    void step();
    void updateFitness(double delta);
    int getBirthdate() { return _birthdate; }
    double getFitness(){ return _currentFitness;}

    double getColorEffector()
    {
        return 2.0 * ((double)_wm->getRobotLED_redValue()/256.0) - 1.0;
    }

    //Forgetting measures
    std::vector<double> _storedF;
    Genome* _storedG;
    void storeRepresentative()
    {
        if(!_doEvoTopo)
        {
            _storedF = _genomeF;
        }
        else
        {
            _storedG = _genome->duplicate();
        }
    }
    double forget();
    double forgetStructural();

    void readGenomeF(std::string s);
    //TODO read evo topo
    //void readGenome(std::string s)
    void logGenomeF(std::string s);
    //TODO log evo topo
    //void logGenome(std::string s)

    bool _doEvoTopo;
    std::map<GC, std::vector<double> > _genomesFList;
    std::map<GC, Genome* > _genomesList;

    std::map<GC, double > _fitnessList;

};


#endif

