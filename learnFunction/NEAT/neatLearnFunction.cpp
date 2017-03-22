
/*
 Copyright 2001 The University of Texas at Austin

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/

#include <iostream>
#include <vector>
#include "neat.h"
#include "population.h"
#include "experiments.h"
using namespace std;

//To avoid circular dependency
NEAT::Network* _nnModel1; //first task model net ODNEATGC
NEAT::Network* _nnModel2; //second task model net ODNEATGC
int main(int argc, char *argv[])
{
  NEAT::Population *p=0;

  //RANDOM SETUP, to the clock tick
  unsigned int seed =(unsigned)time(NULL) + (unsigned)clock();
  //seed = 1481293865; // !!!!
  //std::cerr << seed << std::endl;
  srand(seed);

  if (argc != 5)
  {
    cerr << "A NEAT parameters file (.ne file), a number of generations, idRun and a expBasename"
            << " are required to run the experiments!"
         << endl;
    return -1;
  }

  //Load in the params
  NEAT::load_neat_params(argv[1],false);
  int gens = atoi(argv[2]);
  //third parameter is folder Run for out logfiles, fourth is nameExperiment for getting datasets
  p = rdmNNFunction_test(gens,argv[1], argv[3],argv[4]);


  /*
  //Test a genome file
  Genome *g; Network *n; CartPole *thecart; thecart=new CartPole(true,0);
  g=Genome::new_Genome_load("tested"); n=g->genesis(1); Organism *org= new Organism(0, g, 1);
  thecart->nmarkov_long=true; thecart->generalization_test=false;
  pole2_evaluate(org,0,thecart); cout<<"made score "<<org->fitness<<endl;
  */
  if (p)
    delete p;

  return(0);

}

