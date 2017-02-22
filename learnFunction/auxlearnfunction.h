#ifndef AUXLEARNFUNCTION_H
#define AUXLEARNFUNCTION_H
#include <stdlib.h>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <stdio.h>
#include <vector>
#define ODNEAT_FUNCTIONS
#ifdef ODNEAT_FUNCTIONS
#include <odneatgc/network.h>
#include <odneatgc/genome.h>
#endif //ODNEAT_FUNCTIONS

#include <ExtendedProperties.h>
#include <fstream>

#ifdef ODNEAT_FUNCTIONS
using namespace ODNEATGC;
#endif //ODNEAT_FUNCTIONS


template <typename T> std::vector<T> joinVectors(const std::vector<T>& v1, const std::vector<T>& v2);
//{
//    std::vector<T> result;
//    for(auto &element:v1)
//    {
//        result.push_back(element);
//    }
//    for(auto &element:v2)
//    {
//        result.push_back(element);
//    }
//    return result;
//}
std::vector<double> generateRandomInputs(int dim,double inBound);

std::vector<double> f1(std::vector<double> inputs, int outDim);
std::vector<double> f2(std::vector<double> inputs, int outDim);


double maxAbs(std::vector<std::vector<double> >vv);

double average(std::vector<double> v);
double best(std::vector<double> v);
double worst(std::vector<double> v);

std::vector<double> targetFunction(std::vector<double> inputInstance,int functionSelector, int dimOut);

#ifdef ODNEAT_FUNCTIONS
std::vector<double> activateNN(Network* nAct, std::vector<double> inputs);

double functionError(Network* nTest, std::vector<std::vector<double> > inputBase,
               std::vector<std::vector<double>> outReference);

extern ExtendedProperties gProperties;
bool loadProperties(std::string filename);
#endif //ODNEAT_FUNCTIONS
double avgNumberNodes(std::vector<Genome*> vG);
double avgNumberLinks(std::vector<Genome*> vG);

#endif // AUXLEARNFUNCTION_H
