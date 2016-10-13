#include <odneatgc/network.h>
#include <odneatgc/genome.h>
#include <stdlib.h>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <fstream>

using namespace ODNEATGC;
std::vector<double> generateRandomInputs(int dim)
{
    std::vector<double> result;
    //std::cout << "Inputs" << std::endl;
    for(unsigned int i=0; i < dim; i++)
    {
      double rdm = (double)(rand() % 10000)/5000.0 - 1.0;
      result.push_back(rdm);
      //std::cout << rdm << ", ";
    }
    //std::cout << std::endl;

    return result;
}

std::vector<double> activateNN(Network* nAct, std::vector<double> inputs)
{
     nAct->load_sensors(&(inputs[0]));
    if (!(nAct->activate ()))
        {
            std::cerr << "[ERROR] Activation of ANN not correct" << std::endl;
            exit (-1);
        }
    std::vector<double> outputs;
    //std::cout << "Outputs" << std::endl;
    for (auto out_iter  = nAct->outputs.begin();
             out_iter != nAct->outputs.end();
             out_iter++)
      {
	double outVal = (*out_iter)->activation;
	//			std::cout << outVal << ", ";
	  outputs.push_back(outVal);
      }
    //std::cout << std::endl;
     
    return outputs;
}

double functionError(Network* nTest, std::vector<std::vector<double>> inputBase,
		       std::vector<std::vector<double>> outReference)
{
  //Average euclidean distance between outputs, over all inputs samples in the base
  double  result = 0.0;
  for(auto itI = inputBase.begin(), itO = outReference.begin(); 
      itI != inputBase.end(); itI++, itO++)
    {
      std::vector<double> outs = activateNN(nTest, (*itI));
      double errorOneSample = 0.0;
      for(auto itDimO = (*itO).begin(), itTestO = outs.begin();
	  itDimO != (*itO).end(); itDimO++, itTestO++)
	{   	  
	  double squaredDifference = ((*itDimO) - (*itTestO)) * ((*itDimO) - (*itTestO));
	    errorOneSample+= squaredDifference;
	}
      result += sqrt(errorOneSample);
    }
  return result / (double)inputBase.size();
}

int main(int argc, char* argv[])
{
  //Create a random ANN (fully connected perceptron)
  srand (time(NULL));
  if(argc!=2)
  {
     std::cerr << "[ERROR] Call as ./evaluateRandomMultisynapse {0|1} for {noMultiSynapse,Multisynapse}" << std::endl;
     exit (-1);
  }
  int allowMulti = atoi(argv[1]);
  Helper::allowMultisynapses = allowMulti==1; 
  Helper::mutateToggleEnableProb=0.10;
  Helper::mutateLinkWeightsProb=1.0;
  Helper::mutateAddNodeProb=0.05;
  Helper::mutateAddLinkProb=0.15;
  std::string strAllow = (Helper::allowMultisynapses?"-Multi":"-NoMulti");
  int nIn = 20;
  int nO = 5; 
  GC id;
   
  std::vector<std::vector<double>> inputSet;
  		     
  int numberSamples = 2000;
  //Store reference output values
  for(unsigned int i = 0; i < numberSamples;i++)
  {
    std::vector<double> inputSample = generateRandomInputs(nIn);
    inputSet.push_back(inputSample);
  }	

  unsigned int nbRunsSameSamples = 30;
  for(unsigned int j=0; j < nbRunsSameSamples; j++)
  {
    std::ofstream oFile("logsRobustnessMut/logsMultiNoMultiSame/sameSample-"+strAllow+"-I"+std::to_string(nIn) +"-O"+std::to_string(nO)+"-Run"+ std::to_string(j)+".log");
        
    id.robot_id = -1;
    id.gene_id = 1;
    Genome* g = new Genome(id,nIn,nO);
    
    g->initialize_link_weights();
    Network* n = g->genesis();
    //Add structure at random and mutate its weights
    unsigned int numberNodes = 50;
    double probNode = 0.8;
    int tries = 100;
    int idR = -1;
    int nodeId = nIn + nO + 1;
    int geneId = nIn * nO + 1;
    //Do some random node mutations
    for(unsigned int i = 0; i< numberNodes; i++)
      {
	if(((double)(rand() % 10000)/5000.0) < probNode)
	  {
	    g -> mutate_add_node(tries,idR, nodeId,geneId);
	    n = g->genesis();
	  }
      }
    unsigned int numberLinks = 100;
    double probLink = 0.9;
    //Do some random link mutations
    for(unsigned int i = 0; i< numberLinks; i++)
      {
	if(((double)(rand() % 10000)/5000.0) < probLink)
	  {
	    g -> mutate_add_link(tries,idR, geneId);
	    n = g->genesis();
	  }
      }
    
    double sigma = 0.1;
    double probWeight = 0.8;
    unsigned int numberWeightMutations = 100;
    
    //Do some random weight mutations
    for(unsigned int i = 0; i< numberWeightMutations; i++)
      {
	if(((double)(rand() % 10000)/5000.0) < probWeight)
	  {
	    g -> mutate_link_weights(sigma);
	    n = g->genesis();
	  }
      }
    
    
    std::vector<std::vector<double>> outputReference;
    
    for(unsigned int i = 0; i < numberSamples;i++)
      {
	outputReference.push_back(activateNN(n, inputSet[i]));
      }
    //Perform random perturbations (weights ?and/or structure?)
    //Then measure the distance w.r.t the original function
    unsigned int numberPerturbations = 500;
    double sigmaPerturbations = 0.1;
    for(unsigned int i = 0; i< numberPerturbations; i++)
      {
	id.gene_id++;
	g = g->mutate(sigmaPerturbations,idR,id, nodeId,geneId);
	
	n = g->genesis();
	std::cout 
	//oFile
	  << functionError(n, inputSet, outputReference) << std::endl;
	  std::stringstream os;
	  os << "logsTest/" << i << ".nn";
	  g -> print_to_filename(os.str().c_str());
      }
    oFile.close();
  }
}
