//#include <odneatgc/network.h>
//#include <odneatgc/genome.h>
//Generate with NEAT ANNs
#include <network.h>
#include <genome.h>
#include <stdlib.h>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <fstream>
//#include <learnFunction.h>
#include "auxlearnfunction.h"

NEAT::Genome* createRandomGenome(int nIn,int nO, int id) //GC id)
{
    NEAT::Genome* randomG = NEAT::Genome::makeGenome(id,nIn,nO);//new
    int neuronId = nIn + nO + 1;
    double linkId = nIn * nO + 1; //int
    std::vector<NEAT::Innovation*> innovs;
    int nbNodeMut = nIn + nO + 5; //TOERASE A LOT OF MUTATIONS FOR FUNCTION TO LEARN
    NEAT::Network* n;
    for(int i = 0 ; i < nbNodeMut;i++)
    {
          randomG->mutate_add_node(innovs, neuronId, linkId);//NEAT::newlink_tries,-1,neuronId,linkId);//Helper, newStructureTries
          n = randomG->genesis(id); //noparam
    }
    int nbLinkMut = nbNodeMut * nbNodeMut / 4;//TOERASE A LOT OF MUTATIONS FOR FUNCTION TO LEARN
    for(int i = 0 ; i< nbLinkMut;i++)
    {
        randomG->mutate_add_link(innovs, linkId, NEAT::newlink_tries);//NEAT::newlink_tries,-1,linkId);//Helper, newStructureTries
        n = randomG->genesis(id); //noparam
    }
    double sigma = 0.1;
    int nbWeightMut = 20 * 4;//TOERASE A LOT OF MUTATIONS FOR FUNCTION TO LEARN
    for(int i = 0 ; i< nbWeightMut;i++)
    {
       randomG->mutate_link_weights(sigma,1.0,NEAT::GAUSSIAN);//(sigma);
       n = randomG->genesis(id); //noparam
    }
    return randomG;
}

//Generate dataset
int main(int argc, char *argv[])
{
  srand (time(NULL));
  int nIn; int nO;

  if (argc != 3)
  {
    std::cerr << "A parameter file and a basename are required to run the experiments!"
         << std::endl;
    return -1;
  }

  //Load in the params

  //GC id;
  //id.robot_id = -1;
  //id.gene_id = 0;
  int id = 0;

  //Mutation parameters
  NEAT::newlink_tries= 100;//Helper::newStructureTries = 100;
  NEAT::mutate_link_weights_prob=1.0;//Helper::mutateLinkWeightsProb=1.0;
  NEAT::mutateIndividualWeightProb = 1.0;//Helper::mutateIndividualWeightProb = 1.0;
  NEAT::allowMultisynapses = false;//Helper::allowMultisynapses = false;


  int numberSamples = 2000;

  std::string paramfile = argv[1];
  if(!loadProperties(paramfile))
  {
      std::cerr << "[ERROR] Wrong parameter filename" << std::endl;
      exit(-1);
  }

  NEAT::gProperties.checkAndGetPropertyValue("nIn",&nIn,true);
  NEAT::gProperties.checkAndGetPropertyValue("nO",&nO,true);
  NEAT::gProperties.checkAndGetPropertyValue("numberSamples",&numberSamples,true);
  NEAT::gProperties.checkAndGetPropertyValue("newstructure_tries",&NEAT::newlink_tries,true); //&Helper::newStructureTries,true);
  NEAT::gProperties.checkAndGetPropertyValue("mutate_individual_weight_prob",&NEAT::mutateIndividualWeightProb,true); //&Helper::mutateIndividualWeightProb,true);
  NEAT::gProperties.checkAndGetPropertyValue("allowMultisynapses",&NEAT::allowMultisynapses,true);//&Helper::allowMultisynapses,true);
  NEAT::gProperties.checkAndGetPropertyValue("bias",&NEAT::withBias,true);//&Helper::withBias,true);

  if(NEAT::withBias)//if(Helper::withBias)
      nIn++;
  NEAT::Genome* g1 = createRandomGenome(nIn,nO,id);
  //id.gene_id++;
  NEAT::Network* n1 = g1->genesis(id); //noparam
  id++;
  NEAT::Genome* g2 = createRandomGenome(nIn,nO,id);
  NEAT::Network* n2 = g2->genesis(id); //noparam

  std::vector<std::vector<double>> inputSet, inputSet2;
  std::vector<std::vector<std::vector<double>>> outputSet(2, std::vector<std::vector<double>>());


  std::stringstream ssFname;
  ssFname << "../datasets/o1-" << argv[2] << ".nn";
  std::ofstream nn1File(ssFname.str().c_str());
  g1 -> print_to_file(nn1File);//to_filename() ssFname.str().c_str());
  nn1File.close();

  ssFname.str(std::string());
  ssFname << "../datasets/o2-" << argv[2] << ".nn";
  std::ofstream nn2File(ssFname.str().c_str());
  g2 -> print_to_file(nn2File);//to_filename() ssFname.str().c_str());
  nn2File.close();

  double inputBound = 1; // / sqrt(nIn); for sphere function

  if(NEAT::withBias)//if(Helper::withBias)
      nIn--;
  for(unsigned int i = 0; i < numberSamples;i++)
  {
    std::vector<double> inputSample = generateRandomInputs(nIn, inputBound);
    inputSet.push_back(inputSample);
    std::vector<double> inputSample2 = generateRandomInputs(nIn, inputBound);
    inputSet2.push_back(inputSample2);

    // Check if input sample are the same?
    outputSet[0].push_back(activateNN(n1, inputSample));

    outputSet[1].push_back(activateNN(n2, inputSample2));

  }


    ssFname.str(std::string());
    ssFname << "../datasets/in-" << argv[2] << ".dat";

  std::ofstream inputOFile(ssFname.str().c_str());
  for(auto &input:inputSet)
  {
      for(auto &elt:input)
      {
          inputOFile << elt << " ";
      }
      inputOFile << "\n";
  }
  inputOFile.close();
  ssFname.str(std::string());
  ssFname << "../datasets/in2-" << argv[2] << ".dat";

  std::ofstream inputOFile2(ssFname.str().c_str());
  for(auto &input:inputSet2)
  {
      for(auto &elt:input)
      {
          inputOFile2 << elt << " ";
      }
      inputOFile2 << "\n";
  }
  inputOFile2.close();
  ssFname.str(std::string());
  ssFname << "../datasets/o1-" << argv[2] << ".dat";

  std::ofstream output1OFile(ssFname.str().c_str());
  for(auto &output:outputSet[0])
  {
      for(auto &elt:output)
      {
          output1OFile << elt << " ";
      }
      output1OFile << "\n";
  }
  output1OFile.close();
  ssFname.str(std::string());
  ssFname << "../datasets/o2-" << argv[2] << ".dat";

  std::ofstream output2OFile(ssFname.str().c_str());
  for(auto &output:outputSet[1])
  {
      for(auto &elt:output)
      {
          output2OFile << elt << " ";
      }
      output2OFile << "\n";
  }
  output2OFile.close();

}
