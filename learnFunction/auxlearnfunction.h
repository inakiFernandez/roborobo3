#ifndef AUXLEARNFUNCTION_H
#define AUXLEARNFUNCTION_H
#include <stdlib.h>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <vector>
//#define ODNEAT_FUNCTIONS
//#ifdef ODNEAT_FUNCTIONS
//#include <odneatgc/network.h>
//#include <odneatgc/genome.h>
//#endif //ODNEAT_FUNCTIONS

#include <ExtendedProperties.h>
#include <fstream>
#include<network.h>
#include<genome.h>
#include<neat.h>

//#ifdef ODNEAT_FUNCTIONS
//using namespace ODNEATGC;
//#endif //ODNEAT_FUNCTIONS

//using namespace NEAT;//no NEAT

template <typename T> std::vector<T> joinVectors(const std::vector<T>& v1, const std::vector<T>& v2);

std::vector<double> generateRandomInputs(int dim,double inBound);

std::vector<double> f1(std::vector<double> inputs, int outDim);
std::vector<double> f2(std::vector<double> inputs, int outDim);


double maxAbs(std::vector<std::vector<double> >vv);

double average(std::vector<double> v);
double best(std::vector<double> v);
double worst(std::vector<double> v);

std::vector<double> targetFunction(std::vector<double> inputInstance,int functionSelector, int dimOut);


//#ifdef ODNEAT_FUNCTIONS
std::vector<double> activateNN(NEAT::Network* nAct, std::vector<double> inputs); //noNEAT

double functionError(NEAT::Network* nTest, std::vector<std::vector<double> > inputBase, //noNEAT
               std::vector<std::vector<double>> outReference);

//static ExtendedProperties gProperties;
bool loadProperties(std::string filename);
//#endif //ODNEAT_FUNCTIONS


double avgNumberNodes(std::vector<NEAT::Genome*> vG);
double avgNumberLinks(std::vector<NEAT::Genome*> vG);

#endif // AUXLEARNFUNCTION_H
