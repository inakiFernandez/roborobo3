
#include "odNeat/include/odNeatSpecies.h"
#include <map>


odNeatSpecies::odNeatSpecies(int id)
{
    _id = id;
    _speciesFitness = 0.0;
}

odNeatSpecies::~odNeatSpecies()
{
    std::cerr << "TODO destructor species\n";
    //TODO => do not destroy genomes. Needed when recomputing species
}
bool odNeatSpecies::has(GC gId)
{
    for(auto it = _lGenomes.begin(); it != _lGenomes.end(); it++)
    {
        if(it->first == gId) return true;
    }
    return false;
}

void odNeatSpecies::computeSpeciesFitness()
{
    double totalFit = 0.0;

    for (auto itS = _lGenomes.begin(); itS != _lGenomes.end(); itS++)
    {
        //Adjust each genome's fitness
        totalFit += itS->second.f / _lGenomes.size();
    }
    //Average adjusted fitness
    _speciesFitness = totalFit / _lGenomes.size();
}

void odNeatSpecies::add(Genome* g, double f)
{
    genomeFitness gF; gF.g = g; gF.f = f;

    _lGenomes.insert(std::make_pair(g->genome_id, gF));
}

int odNeatSpecies::remove(GC idG)
{
    int result = _lGenomes.erase(idG);
    if(result == 0)
    {
        std::cerr << "Tried to erase unexisting genome ("
                  << idG.robot_id << ";" << idG.gene_id << ")"
                  << " in species " << _id << std::endl;
        exit(-1);
    }
    computeSpeciesFitness();
    return result;
}
/*
double odNeatSpecies::worseAdjFitness()
{
    double result = 10000.0;
    for (auto itS = _lGenomes.begin(); itS != _lGenomes.end(); itS++)
    {
        double adjF = std::get<1>(itS->second) / _lGenomes.size();
        if ( adjF < result)
            result = adjF;
    }
    return result;
}*/

