
#include "odNeat/include/odNeatSpecies.h"
#include <map>


odNeatSpecies::odNeatSpecies(int id)
{
    _id = id;
    _speciesFitness = 0.0;
}

odNeatSpecies::~odNeatSpecies()
{
}
bool odNeatSpecies::has(Genome* g)
{
    for(auto it = _lGenomes.begin(); it != _lGenomes.end(); it++)
    {
        if(it->first == g->genome_id)
            return true;
    }
    return false;
}

void odNeatSpecies::computeSpeciesFitness()
{
    double avgAdjFitness=0.0;
    //TODO
    _speciesFitness = avgAdjFitness;
}

void odNeatSpecies::add(message m)
{
    _lGenomes.insert(std::make_pair(std::get<0>(m)->genome_id, m));
}
