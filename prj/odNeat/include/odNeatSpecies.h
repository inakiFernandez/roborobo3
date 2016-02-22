#ifndef ODNEATSPECIES_H
#define ODNEATSPECIES_H

#include "odneatgc/genome.h"
#include "odneatgc/helper.h"
#include "odNeat/include/odNeatSharedData.h"
#include <map>

using namespace ODNEATGC;
struct genomeFitness
{
    Genome* g;
    double f;
    genomeFitness()
    {
        g=NULL;
        f = -1.0;
    }
    ~genomeFitness()
    {
        g=NULL;
        f = -2.0;
    }
};

class odNeatSpecies
{
    public:
        int _id;
        double _speciesFitness;
        //map from Genome ID to (Genome and fitness)
        std::map<GC, genomeFitness> _lGenomes;

        void add(Genome* g, double f); int remove(GC idG);
        bool has(GC gId);

        void computeSpeciesFitness();

        odNeatSpecies(int id);
        ~odNeatSpecies();
};

#endif // ODNEATSPECIES_H
