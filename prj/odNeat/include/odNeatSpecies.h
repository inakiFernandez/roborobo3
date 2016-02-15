#ifndef ODNEATSPECIES_H
#define ODNEATSPECIES_H

#include "odneatgc/genome.h"
#include "odneatgc/helper.h"
#include "odNeat/include/odNeatSharedData.h"
#include <map>

using namespace ODNEATGC;

class odNeatSpecies
{
    public:
        int _id;
        std::map<GC, message> _lGenomes;
        double _speciesFitness;

        void add(message m);
        void remove(GC idG);
        void computeSpeciesFitness();

        bool has(Genome* g);

        odNeatSpecies(int id);
        ~odNeatSpecies();
};

#endif // ODNEATSPECIES_H
