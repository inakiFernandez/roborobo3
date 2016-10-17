#include <odneatgc/network.h>
#include <odneatgc/genome.h>
#include <stdlib.h>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <fstream>
//#include <>

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
  
//Testing neural network mutations (visualization after topological mutations)
int main()
{
  srand (time(NULL));
  int nIn = 5;
  int nO = 3;
  GC id;

  /* std::vector<double> inTest;
       inTest.push_back(1.0);
       inTest.push_back(1.0);
	  
      std::vector<double> outTest = activateNN(n, inTest);
      for (auto it = outTest.begin();it!=outTest.end();it++)
      {
          std::cout << (*it) << ", ";
      }
      std::cout << std::endl;*/
  std::vector<std::vector<double>> inputSet;
 		     
  int numberSamples = 2000;
 
  for(unsigned int i = 0; i < numberSamples;i++)
  {
    std::vector<double> inputSample = generateRandomInputs(nIn);
    inputSet.push_back(inputSample);
    // outputReference.push_back(activateNN(n, inputSample));
  }
  /*unsigned int nbRunsSameSamples = 30;
  for(unsigned int j=0; j < nbRunsSameSamples; j++)
  {
    std::ofstream oFile("logsSame/sameSampleNodes-Post-I"+std::to_string(nIn) +"-O"+std::to_string(nO)+"-Run"+ std::to_string(j)+".log");
        
    id.robot_id = -1;
    id.gene_id = 1;
    Genome* g = new Genome(id,nIn,nO);

    g->initialize_link_weights();
    unsigned int numberInitialMutates = 20;
    double sigmaInitialMut = 0.5;
    
    for(unsigned int i = 0; i< numberInitialMutates; i++)
      {
	g->mutate_link_weights(sigmaInitialMut);
      }
    Network* n = g->genesis();
     std::vector<std::vector<double>> outputReference;

    for(unsigned int i = 0; i < numberSamples;i++)
      {
	outputReference.push_back(activateNN(n, inputSet[i]));
      }
    //std::stringstream os;
    //os << "logsNode/" << -1 << ".nn";
    //g -> print_to_filename(os.str().c_str());
    unsigned int numberNodes = 500;
    
    int tries = 100;
    int idR = -1;
    int nodeId = nIn + nO + 1;
    int geneId = nIn * nO + 1;
    
    //Do some random node mutations
    for(unsigned int i = 0; i< numberNodes; i++)
      {
	g -> mutate_add_node(tries,idR, nodeId,geneId);
	n = g->genesis();
	//std::cout 
	oFile
	  << functionError(n, inputSet, outputReference)  << std::endl;
	/*std::stringstream os;
	os << "logsNode/" << i << ".nn";
	g -> print_to_filename(os.str().c_str());*/
  /*  }
   oFile.close();								    
  }*/


    id.robot_id = -1;
    id.gene_id = 1;    
    Genome* g = new Genome(id,nIn,nO);
    int tries = 100;
    Network* n = g->genesis();
    int idR = -1;
    int nodeId = nIn + nO + 1;
    int geneId = nIn * nO + 1;
    

    int numberLinks = 10;
    std::vector<std::vector<double>> outputReference;
    
    for(unsigned int i = 0; i < numberSamples;i++)
      {
	outputReference.push_back(activateNN(n, inputSet[i]));
      }	 
    Helper::allowMultisynapses = true; //false; //
    for(unsigned int i = 0; i < numberLinks; i++)
      {
	if(!g -> mutate_add_link(tries,idR,geneId))
	  {  std::cout << "Not link added" << std::endl;
	  }
	    n = g->genesis();
	    std::stringstream os;
	    os << "logsTestCountLinks/link" << i << ".nn";
	    g -> print_to_filename(os.str().c_str());


	    std::cout 
	      << functionError(n, inputSet, outputReference)  
	      << " (";
	       std::vector<Gene*>::iterator curgene;
	      for(curgene=g->genes.begin();curgene!=g->genes.end();curgene++)
	      {
		std::cout << g->getNumberSynapses((*curgene)) << ", ";
	      }
	      std::cout << std::endl;
       }


    // std::cout << "Error with itself: " 
    //<< functionError(n, inputSet, outputReference) << std::endl;
    /* unsigned int numberToggle = 2;
  for(unsigned int i = 0; i < numberToggle; i++)
  {
    g -> mutate_toggle_enable(1);
    n = g->genesis();      
	  
      std::vector<double> outTest = activateNN(n, inTest);
      for (auto it = outTest.begin();it!=outTest.end();it++)
      {
          std::cout << (*it) << ", ";
      }
      std::cout << std::endl;
	//std::cout << functionError(n, inputSet, outputReference) << std::endl;
    std::stringstream os;
    os << "logsTog/" << i << ".nn";
    g -> print_to_filename(os.str().c_str());
  }*/

    //  std::cout << functionError(n, inputSet, outputReference) << std::endl;
    /*   int numberLinks = 1;
	 /*		   
     Helper::allowMultisynapses = true;
    for(unsigned int i = 0; i < numberLinks; i++)
      {
	if(!g -> mutate_add_link(tries,idR,geneId))
	  {  std::cout << "Not link added" << std::endl;
	  }
	    n = g->genesis();
   std::stringstream os;
    os << "logs/link" << i << ".nn";
    g -> print_to_filename(os.str().c_str());

	        //std::cout << "Error with net" << i << ": ";
     std::cout 
       //<<"ERR:" 
	       << functionError(n, inputSet, outputReference) 
       //<<"\n\n\n\n\n"
       << std::endl;
       }
       std::vector<double> inTest;
       inTest.push_back(1.0);
       inTest.push_back(-1.0);
	  
      std::vector<double> outTest = activateNN(n, inTest);
      for (auto it = outTest.begin();it!=outTest.end();it++)
      {
          std::cout << (*it) << ", ";
      }
      std::cout << std::endl;*/

}
