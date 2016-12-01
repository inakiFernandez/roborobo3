#ifndef HELPER_H
#define HELPER_H
#include <iomanip>
namespace ODNEATGC
{
    struct GC
    {
        int robot_id;
        int gene_id;
        bool operator<(const GC &o)  const {
            if(gene_id == o.gene_id)
                return robot_id < o.robot_id;

            return gene_id < o.gene_id;
        }
        bool operator==(const GC &o)  const {
            return gene_id == o.gene_id && robot_id == o.robot_id;
        }
        //friend std::ostream& operator<<(std::ostream& os, const GC& gene_clock);
    };
    class Helper
    {
        public:
            //Global variables: definition
            static double     mutateProb; // Prob. of a mutation
            static double     mutateLinkWeightsProb;
            static double     mutateIndividualWeightProb;
            static double     mutateAddNodeProb;
            static double     mutateAddLinkProb;
            static double     mutateToggleEnableProb;
            // Prob. of mating
            static double     mateOnlyProb;
            // Probability of forcing choice of ONLY links that are naturally recurrent
            static double     recurOnlyProb;
            // Number of tries mutate_add_link or mutate_add_node will attempt to find an open link
            static int        newStructureTries;

            static double     coefE;
            static double     coefD;
            static double     coefW;

            static double     rangeW;

            static bool withBias;

            static bool allowMultisynapses;

            static double randFloat();
            static int randPosNeg();
            static int randInt(int x, int y);
            static int getUnitCounts(const char *string,const char *set);

            static double fSigmoid(double activesum,double slope);
            static double gaussRand();
            static bool load_odneat_params(const char *filename, bool output);
    };
    //#################AUXILIARY#############################
    /*std::ostream& ODNEATGC::operator<<(std::ostream& os, const ODNEATGC::GC& gene_clock)
    {
        os << "(R" << gene_clock.robot_id << "; G" << gene_clock.gene_id << ")";
        return os;
    }*/
} // namespace ODNEATGC

#endif // HELPER_H
