#include <odneatgc/network.h>
#include <odneatgc/genome.h>
#include <stdlib.h>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <fstream>
#include <learnFunction.h>
#include "auxlearnfunction.h"

Genome* createRandomGenome(int nIn,int nO, GC id)
{

    Genome* randomG = new Genome(id,nIn,nO);
    int neuronId = nIn + nO + 1;
    int linkId = nIn * nO + 1;

    int nbNodeMut = nIn + nO + 5; //TOERASE A LOT OF MUTATIONS FOR FUNCTION TO LEARN
    Network* n;
    for(int i = 0 ; i < nbNodeMut;i++)
    {
          randomG->mutate_add_node(Helper::newStructureTries,-1,neuronId,linkId);
          n = randomG->genesis();
    }
    int nbLinkMut = nbNodeMut * nbNodeMut / 4;//TOERASE A LOT OF MUTATIONS FOR FUNCTION TO LEARN
    for(int i = 0 ; i< nbLinkMut;i++)
    {
        randomG->mutate_add_link(Helper::newStructureTries,-1,linkId);
        n = randomG->genesis();
    }
    double sigma = 0.1;
    int nbWeightMut = 20 * 4;//TOERASE A LOT OF MUTATIONS FOR FUNCTION TO LEARN
    for(int i = 0 ; i< nbWeightMut;i++)
    {
       randomG->mutate_link_weights(sigma);
       n = randomG->genesis();
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

  GC id;
  id.robot_id = -1;
  id.gene_id = 0;

  //Mutation parameters
  Helper::newStructureTries = 100;
  Helper::mutateLinkWeightsProb=1.0;
  Helper::mutateIndividualWeightProb = 1.0;
  Helper::allowMultisynapses = false;


  int numberSamples = 2000;

  std::string paramfile = argv[1];//"config/template-params";
  if(!loadProperties(paramfile))
  {
      std::cerr << "[ERROR] Wrong parameter filename" << std::endl;
      exit(-1);
  }

  gProperties.checkAndGetPropertyValue("nIn",&nIn,true);
  gProperties.checkAndGetPropertyValue("nO",&nO,true);
  gProperties.checkAndGetPropertyValue("numberSamples",&numberSamples,true);
  gProperties.checkAndGetPropertyValue("newstructure_tries",&Helper::newStructureTries,true);
  gProperties.checkAndGetPropertyValue("mutate_individual_weight_prob",&Helper::mutateIndividualWeightProb,true);
  gProperties.checkAndGetPropertyValue("allowMultisynapses",&Helper::allowMultisynapses,true);
  gProperties.checkAndGetPropertyValue("bias",&Helper::withBias,true);

  if(Helper::withBias)
      nIn++;
  Genome* g1 = createRandomGenome(nIn,nO,id);
  id.gene_id++;
  Genome* g2 = createRandomGenome(nIn,nO,id);

  Network* n1 = g1->genesis();
  Network* n2 = g2->genesis();

  std::vector<std::vector<double>> inputSet, inputSet2;
  std::vector<std::vector<std::vector<double>>> outputSet(2, std::vector<std::vector<double>>());

  double inputBound = 1; // / sqrt(nIn); for sphere function

  if(Helper::withBias)
      nIn--;
  for(unsigned int i = 0; i < numberSamples;i++)
  {
    //generate outputs for f1 and f2 from SAME inputs?
    std::vector<double> inputSample = generateRandomInputs(nIn, inputBound);
    inputSet.push_back(inputSample);
    std::vector<double> inputSample2 = generateRandomInputs(nIn, inputBound);
    inputSet2.push_back(inputSample2);

    // Check if input sample are the same?

    outputSet[0].push_back(activateNN(n1, inputSample));
    n1->flush();
    outputSet[1].push_back(activateNN(n2, inputSample2));
    n2->flush();
  }
//paramfile.substr(paramfile.find("/")+1,paramfile.size()
    std::stringstream ssFname;
    ssFname << "../datasets/o1-" << argv[2] << ".nn";
    //g1 -> print_to_filename(("datasets/o1-" + argv[2] +".nn").c_str());

    g1 -> print_to_filename(ssFname.str().c_str());
    ssFname.str(std::string());
    ssFname << "../datasets/o2-" << argv[2] << ".nn";

    //("datasets/o2-" + argv[2] +".nn").c_str()
    g2 -> print_to_filename(ssFname.str().c_str());
    ssFname.str(std::string());
    ssFname << "../datasets/in-" << argv[2] << ".dat";
    //"datasets/in-" + argv[2] +".dat"
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
  //"datasets/in2-" + argv[2] +".dat"
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
  //"datasets/o1-" + argv[2] +".dat"
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
  //"datasets/o2-" + argv[2] +".dat"
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

    //functionError(n, inputSet, outputReference)
}

//Normalize outputs (by maxabsolute value)
/*double max = maxAbs(outputSet[0]);
for(int i = 0; i < outputSet[0].size(); i++)
{
    for(int j = 0; j < outputSet[0][i].size(); j++)
    {
        outputSet[0][i][j] = outputSet[0][i][j]/max;
    }
}
max = maxAbs(outputSet[1]);
for(int i = 0; i < outputSet[1].size(); i++)
{
    for(int j = 0; j < outputSet[1][i].size(); j++)
      outputSet[1][i][j] = outputSet[1][i][j]/max;
}*/
