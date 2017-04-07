#include "auxlearnfunction.h"


#include <NEAT/neat.h>


std::vector<double> generateRandomInputs(int dim,double inBound)
{
    //Returns a vector of size dim with uniformly distributed values between -1 and 1
    std::vector<double> result;
    for(unsigned int i=0; i < dim; i++)
    {
      double rdm = ((double)(rand() % 10000)/5000.0 - 1.0) * inBound;
      result.push_back(rdm);
    }
    if(NEAT::withBias)//if(Helper::withBias)
        result.push_back(1.0);
    return result;
}

std::vector<double> f1(std::vector<double> inputs, int outDim)
{
    std::vector<double> result;
    //Sphere function
    for(auto i=0;i< outDim ; i++)
    {
        double partialResult = 0.0;
        for(auto &elt:inputs)
        {
            partialResult += elt * elt;
        }
        if(partialResult > 1.0)
        {
            std::cout << "Output out of bounds" << std::endl;
        }
        result.push_back(partialResult);
    }

    return result;
}
std::vector<double> f2(std::vector<double> inputs, int outDim)
{
    std::vector<double> result;
    for(auto i=0;i< outDim ; i++)
        {
            double partialResult = 1.0;
            for(auto &elt:inputs)
            {
                partialResult *=  sin(elt);
            }
            result.push_back(partialResult);
        }
    return result;
}
double avgNumberNodes(std::vector<NEAT::Genome*> vG)
{
    double result = 0.0;

    for(auto &g:vG)
    {
        result += g->nodes.size();
    }

    return result/vG.size();
}
double avgNumberLinks(std::vector<NEAT::Genome*> vG)
{
    double result = 0.0;

    for(auto &g:vG)
    {
        result += g->genes.size();
    }

    return result/vG.size();
}
double maxAbs(std::vector<std::vector<double>>vv)
{
    double result = 0.0;
    for(auto &v : vv)
    {
        for(auto &elt:v)
        {
            if (fabs(elt)> result)
               result = fabs(elt);
        }
    }
    return result;
}

double average(std::vector<double> v)
{
    double result = 0.0;
    for(auto &elt:v)
    {
        result += elt;
    }
    return result/(double)v.size();
}
double best(std::vector<double> v)
{
    double result = 1000.0;
    for(auto &elt:v)
    {
        if(elt < result)
            result = elt;
    }
    return result;
}
double worst(std::vector<double> v)
{
    double result = -1000.0;
    for(auto &elt:v)
    {
        if(elt > result)
            result = elt;
    }
    return result;
}
std::vector<double> targetFunction(std::vector<double> inputInstance,int functionSelector, int dimOut)
{
    switch(functionSelector)
    {
        case 0:
            return f1(inputInstance, dimOut);
            break;
        case 1:
            return f2(inputInstance, dimOut);
            break;
        default:
            return std::vector<double>();
    }
}


std::vector<double> activateNN(NEAT::Network* nAct, std::vector<double> inputs)
{

    nAct->load_sensors(&(inputs[0]));
    if (!(nAct->activate ()))
        {
            std::cerr << "[ERROR] Activation of ANN not correct"
                      << std::endl;
            exit (-1);
        }
    int nnDepth = 30;
    for(int i=0; i < nnDepth; i++)
        nAct->activate ();
    std::vector<double> outputs;

    for (auto out_iter  = nAct->outputs.begin();
             out_iter != nAct->outputs.end();
             out_iter++)
      {
        double outVal = (*out_iter)->activation;
        outputs.push_back(outVal);
      }

    nAct->flush();
    return outputs;
}

double functionError(NEAT::Network* nTest, std::vector<std::vector<double> > inputBase,
               std::vector<std::vector<double>> outReference)
{
  //Average quadratic error between outputs, over all inputs samples in the base
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
//ExtendedProperties gProperties;
bool loadProperties(std::string filename)
{
    std::ifstream in(filename.c_str());

    if ( !in.is_open() ) // WAS: if ( in == NULL )
        return false;
    NEAT::gProperties.load(in);
    in.close();
    return true;
}

template <typename T> std::vector<T> joinVectors(const std::vector<T>& v1, const std::vector<T>& v2)
{
    std::vector<T> result;
    for(auto &element:v1)
    {
        result.push_back(element);
    }
    for(auto &element:v2)
    {
        result.push_back(element);
    }
    return result;
}
