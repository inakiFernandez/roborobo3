/**
 * @author Nicolas Bredeche <nicolas.bredeche@upmc.fr>
 * NNlib: Leo Cazenille <leo.cazenille@upmc.fr>
 */



#ifndef INCREMENTCONTROLLER_H
#define INCREMENTCONTROLLER_H

#include "RoboroboMain/common.h"
#include "RoboroboMain/roborobo.h"
#include "Utilities/Graphics.h"
#include "Controllers/Controller.h"
#include "WorldModels/RobotWorldModel.h"
#include "Increment/include/IncrementAgentObserver.h"
#include <neuralnetworks/NeuralNetwork.h>

#include <iomanip>

using namespace Neural;


class IncrementController : public Controller
{
	private:
		int _iteration;
        int _birthdate; // evaluation when this controller was initialized.
        double _currentFitness;

		std::vector<double> _parameters;
		std::string _nnType;
		std::vector<int> _nbHiddenNeuronsPerLayer;
		std::vector<int> _nbBiaisNeuronsPerLayer;
		NeuralNetwork* nn;

		void createNN();
    
        //bool _isAlive; // agent stand still if not.
        bool _isNewGenome;
    
        void selectRandomGenome();
        void mutate(float sigma);

        void stepBehaviour();
        void stepEvolution();
    
        void broadcastGenome();
        void loadNewGenome();
    
        unsigned int computeRequiredNumberOfWeights();
    
        bool getNewGenomeStatus() { return _isNewGenome; }
        void setNewGenomeStatus( bool __status ) { _isNewGenome = __status; }
    
        // evolutionary engine
        std::vector<double> _genome; // current genome in evaluation
        std::vector<double> _previousGenome; // 1+1-online-ES surviving genome
        std::map<int, std::vector<double> > _genomesList;
        std::map<int, double > _fitnessList;
        int _populationSize;
        std::map<int,int> _birthdateList; // store the birthdate of the received controllers (useful for monitoring).
        std::vector<double> _currentGenome;
        float _currentSigma;
    
        // ANN
        double _minValue;
        double _maxValue;
        unsigned int _nbInputs;
        unsigned int _nbOutputs;
        unsigned int _nbHiddenLayers;
        std::vector<unsigned int>* _nbNeuronsPerHiddenLayer;
    
        void storeGenome(std::vector<double> genome, int senderId, int senderBirthdate);
        void resetRobot();
    
	public:

        IncrementController(RobotWorldModel *wm);
		~IncrementController();

		void reset();
		void step();
        void updateFitness(double delta);
        int getBirthdate() { return _birthdate; }

};


#endif

