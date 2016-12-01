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

using namespace ODNEATGC;

std::vector<Genome*> initPopulation(unsigned int mu, int nI, int nO)
{
    std::vector<Genome*> result;
    for(unsigned int j=0; j < mu; j++)
    {
        GC id;
        id.robot_id = -1;
        id.gene_id = j+1;
        Genome* g = new Genome(id,nI,nO);

        g->initialize_link_weights();

        result.push_back(g);
    }

    return result;
}
std::vector<double> evaluatePopulation(const std::vector<Genome*>& pop,const std::vector<std::vector<double>>&inBase,
                                       const std::vector<std::vector<double>>& outBase)
{
    std::vector<double> result;
    for(auto &g:pop)
    {
        result.push_back(functionError(g->genesis(), inBase,outBase));
    }
    return result;
}
bool compare(std::pair<Genome*,double> i,std::pair<Genome*,double> j)
{
    return (i.second<j.second);
}
std::pair<std::vector<Genome*>,std::vector<double>> select(const std::vector<Genome*>& pop,
                                                           const std::vector<double>& fitness,
                                                           unsigned int number)
{
    //returns the "number" best genomes in the population
    std::vector<Genome*> resultG;
    std::vector<double> resultF;

    std::vector<std::pair<Genome*,double>> popAndFitness;
    for(unsigned int i=0; i < pop.size();i++)
    {
        popAndFitness.push_back(std::make_pair(pop[i],fitness[i]));
    }

    std::sort(popAndFitness.begin(),popAndFitness.end(),compare);
    int numberStored = 0;
    for(auto &individual:popAndFitness)
    {
        resultG.push_back(individual.first);
        resultF.push_back(individual.second);
        numberStored++;
        if(numberStored == number)
            break;
    }
    return std::make_pair(resultG,resultF);
}

std::vector<Genome*> mutate(const std::vector<Genome*>& par, const GC& idGenome, int &nodeId, int &gc, double stdDev)
{
    std::vector<Genome*> result;
    GC idNewGenome;
    idNewGenome.robot_id = -1;
    idNewGenome.gene_id = idGenome.gene_id + 1;

    for(auto &parent:par)
    {
        result.push_back(parent->mutate(stdDev, -1, idNewGenome,nodeId, gc));
        idNewGenome.gene_id++;
    }
    return result;
}


std::vector<std::vector<double>> readDoubleFile(std::string filename)
{
    std::ifstream inputFile(filename);
    std::string line;
    std::vector<std::vector<double>> result;
    if (inputFile)
    {
         while(std::getline(inputFile,line))
         {
             std::istringstream ss(line);
             double value;
             std::vector<double> instance;
             while ( ss >> value )
             {
                 instance.push_back(value);
             }
             result.push_back(instance);
         }
    }
    return result;
}
//Evolutionary algorithm
int main()
{
  srand (time(NULL));
  int nIn = 5;
  int nO = 1; // same dimension for identity function

  GC id;
  id.robot_id = -1;
  id.gene_id = 0;
  int G = 1;
  int durationPerFunction = 100;
  int mu = 1;
  int lambda = 1;

  //Mutation parameters
  Helper::newStructureTries = 100;

  Helper::mutateLinkWeightsProb=1.0;
  Helper::mutateIndividualWeightProb = 1.0;
  Helper::mutateAddNodeProb=0.25;
  Helper::mutateAddLinkProb=0.35;
  Helper::mutateToggleEnableProb = 0.05;
  Helper::allowMultisynapses = false;
  double sigma = 0.1;

  int numberSamples = 2000;

  std::string paramfile = "config/template-params";
  if(!loadProperties(paramfile))
  {
      std::cerr << "[ERROR] Wrong parameter filename" << std::endl;
      exit(-1);
  }
  gProperties.checkAndGetPropertyValue("newstructure_tries",&Helper::newStructureTries,true);
  gProperties.checkAndGetPropertyValue("mutate_link_weights_prob",&Helper::mutateLinkWeightsProb,true);
  gProperties.checkAndGetPropertyValue("mutate_individual_weight_prob",&Helper::mutateIndividualWeightProb,true);
  gProperties.checkAndGetPropertyValue("mutate_add_node_prob",&Helper::mutateAddNodeProb,true);
  gProperties.checkAndGetPropertyValue("mutate_add_link_prob",&Helper::mutateAddLinkProb,true);
  gProperties.checkAndGetPropertyValue("mutate_toggle_enable_prob",&Helper::mutateToggleEnableProb,true);
  gProperties.checkAndGetPropertyValue("allowMultisynapses",&Helper::allowMultisynapses,true);
  gProperties.checkAndGetPropertyValue("sigma",&sigma,true);
  gProperties.checkAndGetPropertyValue("nIn",&nIn,true);
  gProperties.checkAndGetPropertyValue("nO",&nO,true);
  gProperties.checkAndGetPropertyValue("G",&G,true);
  gProperties.checkAndGetPropertyValue("durationPerFunction",&durationPerFunction,true);
  gProperties.checkAndGetPropertyValue("mu",&mu,true);
  gProperties.checkAndGetPropertyValue("lambda",&lambda,true);
  gProperties.checkAndGetPropertyValue("numberSamples",&numberSamples,true);

  //Population structures
  std::vector<double> fitness;
  std::vector<Genome*> population;
  std::pair<std::vector<Genome*>,std::vector<double>> parents;
  std::vector<Genome*> children;
  std::vector<double> fitnessChildren;
  int neuronId, linkId;


  std::vector<std::vector<double>> inputSet;
  std::vector<std::vector<std::vector<double>>> outputSet(2, std::vector<std::vector<double>>());

  inputSet = readDoubleFile("datasets/in-" + paramfile.substr(paramfile.find("/")+1,paramfile.size()) +".dat");
  outputSet[0] = readDoubleFile("datasets/o1-" + paramfile.substr(paramfile.find("/")+1,paramfile.size()) +".dat");
  outputSet[1] = readDoubleFile("datasets/o2-" + paramfile.substr(paramfile.find("/")+1,paramfile.size()) +".dat");

  //Runs
  unsigned int nbRunsSameSamples = 1;//30;
  for(unsigned int i=0; i < nbRunsSameSamples; i++)
  {
      int f = 0;
      //EA
      //InitPopulation
      population = initPopulation(mu, nIn, nO);
      neuronId = nIn + nO + 1;
      linkId = nIn * nO + 1;
      //Evaluate initial population (mu)
      fitness = evaluatePopulation(population,inputSet, outputSet[f]);
      //Evo Loop
      for(unsigned int j=0; j < G; j++)
      {
        parents = select(population,fitness, lambda);
        children = mutate(parents.first,id,neuronId,linkId,sigma);
        fitnessChildren = evaluatePopulation(children,inputSet, outputSet[f]);

        std::vector<Genome*> fullPop = joinVectors(population, children);
        std::vector<double> fullPopFitness = joinVectors(fitness, fitnessChildren);

        std::pair<std::vector<Genome*>,std::vector<double>> p = select(fullPop,fullPopFitness, mu);

        population = p.first;
        fitness  = p.second;

        //log?
        std::cout << j << " "
                  << f << " "
                  << average(fitness) << " "
                  << best(fitness) << " "
                  //<< worst(fitness) << " "
                  << avgNumberNodes(population) << " "
                  //<< avgNumberLinks(population)
                  << "\n";
        //change function
        if(((j % durationPerFunction) == 0) && (j != 0))
        {
            f = (f + 1) % 2;
            fitness = evaluatePopulation(population,inputSet, outputSet[f]);
        }
      }
  }

  /*
    std::ofstream oFile("logsSame/sameSampleNodes-Post-I"+std::to_string(nIn) +"-O"+std::to_string(nO)+"-Run"+ std::to_string(j)+".log");

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



}
