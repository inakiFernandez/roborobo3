#include "odneatgc/helper.h"
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sstream>
using namespace ODNEATGC;

    //Global variables: definition
    // Prob. of mutating the selected genome or the result of x-over
    double     Helper::mutateProb=0;
    double     Helper::mutateLinkWeightsProb=0;
    double     Helper::mutateIndividualWeightProb = 0;
    double     Helper::mutateAddNodeProb=0;
    double     Helper::mutateAddLinkProb=0;
    double     Helper::mutateToggleEnableProb = 0;
    // Prob. of mating
    double     Helper::mateOnlyProb=0;
    // Probability of forcing choice of ONLY links that are naturally recurrent
    double     Helper::recurOnlyProb=0;
    // Number of tries mutate_add_link or mutate_add_node will attempt to find an open link
    int        Helper::newStructureTries=0;

    double     Helper::coefE= 5.0; //1.0;
    double     Helper::coefD= 10.0; //1.5;
    double     Helper::coefW= 0.4;

    double     Helper::rangeW=10.0;

    bool Helper::allowMultisynapses = false;

    double Helper::randFloat()
    {
        return rand() / (double) RAND_MAX;
    }

    int Helper::randPosNeg()
    {
        if (rand()%2)
            return 1;
        else
            return -1;
    }
    int Helper::randInt(int x, int y)
    //Random int between x and y (inclusive)
    {
        return rand()%(y-x+1)+x;
    }
    int Helper::getUnitCounts(const char *string,const char *set)
    {
        int count = 0;
        short last = 0;
        while(*string)
        {
            last = *string++;

            for(int i =0; set[i]; i++)
            {
                if(last == set[i])
                {
                    count++;
                    last = 0;
                    break;
                }
            }
        }
        if(last) count++;
        return count;
    }

    double Helper::fSigmoid(double activesum,double slope)
    {
        //NON-SHIFTED STEEPENED
        return (1/(1+(exp(-(slope*activesum)))));
    }

    double Helper::gaussRand()
    {
        static int iset=0;
        static double gset;
        double fac,rsq,v1,v2;

        if (iset==0)
        {
            do
            {
                v1=2.0*(randFloat())-1.0;
                v2=2.0*(randFloat())-1.0;
                rsq=v1*v1+v2*v2;
            } while (rsq>=1.0 || rsq==0.0);
            fac=sqrt(-2.0*log(rsq)/rsq);
            gset=v1*fac;
            iset=1;
            return v2*fac;
        }
        else
        {
            iset=0;
            return gset;
        }
    }

    bool Helper::load_odneat_params(const char *filename, bool output)
    {
        std::ifstream paramFile(filename);

        if(!paramFile) {
            return false;
        }
        char curword[128];

        // **********LOAD IN PARAMETERS************** //
        if(output)
            printf("NEAT READING IN %s", filename);

        paramFile>>curword; paramFile>>mutateProb;
        paramFile>>curword; paramFile>>mutateLinkWeightsProb;
        paramFile>>curword; paramFile>>mutateAddNodeProb;
        paramFile>>curword; paramFile>>mutateAddLinkProb;
        paramFile>>curword; paramFile>>mateOnlyProb;
        paramFile>>curword; paramFile>>recurOnlyProb;
        paramFile>>curword; paramFile>>newStructureTries;

        if(output)
        {
            printf("mutate_only_prob=%f\n",mutateProb);
            printf("mutate_link_weights_prob=%f\n",mutateLinkWeightsProb);
            //printf("mutate_toggle_enable_prob=%f\n",mutateToggleEnableProb);
            //printf("mutate_gene_reenable_prob=%f\n",mutateGeneReenableProb);
            printf("mutate_add_node_prob=%f\n",mutateAddNodeProb);
            printf("mutate_add_link_prob=%f\n",mutateAddLinkProb);
            printf("mate_only_prob=%f\n",mateOnlyProb);
            printf("recur_only_prob=%f\n",recurOnlyProb);
            printf("newstructure_tries=%d\n",newStructureTries);
        }
        paramFile.close();
        return true;
    }

    /*void print_Genome_tofile(Genome *g,const char *filename)
    {
        std::string file = "";
        file += filename;

        std::ofstream oFile(file.c_str());

        g->print_to_file(oFile);

        oFile.close();
    }*/
