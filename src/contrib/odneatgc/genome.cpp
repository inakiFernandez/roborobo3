#include "odneatgc/genome.h"
#include "odneatgc/network.h"
#include "odneatgc/gene.h"
#include "odneatgc/helper.h"

#include "Utilities/Misc.h"

#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <cstring>
#include <map>
#include <set>
#include <limits>
#include <boost/range/adaptor/reversed.hpp>

using namespace ODNEATGC;

Genome::Genome(GC id, std::ifstream &iFile)
{
    //max word size of 128 characters
    char curword[128];
    //max line size of 1024 characters
    char curline[1024];
    char delimiters[] = " \n";

    int done=0;
    genome_id=id;

    iFile.getline(curline, sizeof(curline));
    int wordcount = Helper::getUnitCounts(curline, delimiters);
    int curwordnum = 0;

    //Loop until file is finished, parsing each line
    while (!done)
    {
        if (curwordnum > wordcount || wordcount == 0)
        {
            iFile.getline(curline, sizeof(curline));
            wordcount = Helper::getUnitCounts(curline, delimiters);
            curwordnum = 0;
        }

        std::stringstream ss(curline);

        ss >> curword;

        //Check for end of Genome
        if (strcmp(curword,"genomeend")==0)
        {
            GC id;
            ss >> curword; id.robot_id = atoi(curword);
            ss >> curword; id.gene_id = atoi(curword);

            if (!(id==genome_id)) printf("ERROR: id mismatch in genome");
            done=1;
        }

        //Ignore genomestart if it hasn't been gobbled yet
        else if (strcmp(curword,"genomestart")==0)
            ++curwordnum;
        //Ignore comments surrounded by - they get printed to screen
        else if (strcmp(curword,"/*")==0)
        {
            ss >> curword;
            while (strcmp(curword,"*/")!=0)
                ss >> curword;
        }
        //Read in a node
        else if (strcmp(curword,"node")==0)
        {
            NNode *newnode;
            char argline[1024];
            curwordnum = wordcount + 1;

            ss.getline(argline, 1024);
            //Allocate the new node
            newnode=new NNode(argline);
            //Add the node to the list of nodes
            nodes.push_back(newnode);
        }
        //Read in a Gene
        else if (strcmp(curword,"gene")==0)
        {
            Gene *newgene;
            char argline[1024];
            curwordnum = wordcount + 1;

            ss.getline(argline, 1024);
            //Allocate the new Gene
            newgene=new Gene(argline,nodes);

            //Add the gene to the genome
            genes.push_back(newgene);
        }

    }
    species = -1;
    phenotype =NULL;
    nbFitnessUpdates = 0;
}

Genome::Genome(GC id,int num_in,int num_out) {

    //Temporary lists of nodes
    std::vector<NNode*> inputs; std::vector<NNode*> outputs;

    std::vector<NNode*>::iterator curnode1; //Node iterator1
    std::vector<NNode*>::iterator curnode2; //Node iterator2

    //For creating the new genes
    NNode *newnode; Gene *newgene; int ncount, count;

    genome_id=id; phenotype = NULL;

    //Create the inputs and outputs
    //Build the input nodes. Last one is bias
    for(ncount=1;ncount<=num_in;ncount++)
    {
        //Set common initial Gene Clock
        innov innovClock; innovClock.idR = -1; innovClock.gc = ncount;

        if (ncount<num_in)
            newnode=new NNode(SENSOR,innovClock,INPUT);
        else
            newnode=new NNode(SENSOR,innovClock,BIAS);

        //Add the node to the lists of nodes
        nodes.push_back(newnode); inputs.push_back(newnode);
    }
    //Build the output nodes
    for(ncount=num_in+1;ncount<=num_in+num_out;ncount++)
    {
        innov innovClock; innovClock.idR = -1; innovClock.gc = ncount;
        newnode=new NNode(NEURON,innovClock,OUTPUT);
        //Add the node to the lists of nodes
        nodes.push_back(newnode); outputs.push_back(newnode);
    }

    //Create the links: just connect inputs straight to outputs
    //Initially fully-connected single-layer perceptron
    count = 0;
    //Loop over the outputs
    for(curnode1=outputs.begin();curnode1!=outputs.end();++curnode1)
    {
        //Loop over the inputs
        for(curnode2=inputs.begin();curnode2!=inputs.end();++curnode2)
        {
            count++; innov innovClock; innovClock.idR = -1; innovClock.gc = count;
            //Connect each input to each output with 0 weight, to be mutated later
            newgene=new Gene(0, (*curnode2), (*curnode1),false,innovClock);
            //Add the gene to the genome
            genes.push_back(newgene);
        }
    }
    GC m = {id.robot_id, -1}; GC d = {id.robot_id, -1};
    mom_id = m; dad_id = d; species = -1;
    nbFitnessUpdates = 0;
}

Genome *Genome::duplicate()
{
    //Collections for the new Genome
    std::vector<NNode*> nodes_dup; std::vector<Gene*> genes_dup;
    //Iterators for the old Genome
    std::vector<NNode*>::iterator curnode; std::vector<Gene*>::iterator curgene;

    //New item pointers
    NNode *newnode; Gene *newgene;
    //For forming a gene
    NNode *inode; NNode *onode;
    Genome *newgenome;

    //Duplicate NNodes
    for(curnode=nodes.begin();curnode!=nodes.end();++curnode)
    {
        newnode=new NNode(*curnode);
        //Remember this node's old copy
        (*curnode)->dup=newnode;
        nodes_dup.push_back(newnode);
    }
    //Duplicate Genes
    for(curgene=genes.begin();curgene!=genes.end();++curgene)
    {
        //First find the nodes connected by the gene's link
        inode=(((*curgene)->lnk)->in_node)->dup;
        onode=(((*curgene)->lnk)->out_node)->dup;

        newgene=new Gene(*curgene,inode,onode);
        genes_dup.push_back(newgene);
    }

    //Finally, return the genome
    newgenome=new Genome(this->genome_id,nodes_dup,genes_dup);
    newgenome->mom_id = this->mom_id; newgenome->dad_id = this->dad_id;
    newgenome->nbFitnessUpdates = 0;
    newgenome->species = -1;
    //newgenome->genesis();

    return newgenome;
}
Genome::Genome(GC id, std::vector<NNode*> n, std::vector<Gene*> g)
{
    GC m = {id.robot_id, -1}; GC d = {id.robot_id, -1};
    genome_id=id; nodes=n; genes=g; mom_id = m; dad_id = d; species = -1; phenotype =NULL; nbFitnessUpdates = 0;
}
Genome::~Genome() 
{    
    std::vector<NNode*>::iterator curnode;
    std::vector<Gene*>::iterator curgene;

    for(curnode=nodes.begin();curnode!=nodes.end();++curnode)
        delete (*curnode);
    nodes.clear();
    for(curgene=genes.begin();curgene!=genes.end();++curgene)
        delete (*curgene);
    genes.clear();

    //delete (phenotype);
}

Network *Genome::genesis()
{
    std::vector<NNode*>::iterator curnode;
    std::vector<Gene*>::iterator curgene;
    NNode *newnode;
    Link *curlink; Link *newlink;

    //Inputs and outputs will be collected here for the network
    //All nodes are collected in an all_list, used for later safe destruction of the net
    std::vector<NNode*> inlist; std::vector<NNode*> outlist; std::vector<NNode*> all_list;

    //Gene translation variables
    NNode *inode; NNode *onode;
    //New network
    Network *newnet;

    //Create the nodes
    for(curnode=nodes.begin();curnode!=nodes.end();++curnode)
    {
        newnode=new NNode((*curnode)->type,(*curnode)->node_id,(*curnode)->gen_node_label);

        //Check for input or output designation of node
        if ((((*curnode)->gen_node_label)==INPUT) || (((*curnode)->gen_node_label)==BIAS))
            inlist.push_back(newnode);
        else if (((*curnode)->gen_node_label)==OUTPUT)
            outlist.push_back(newnode);
        //Keep track of all nodes, not just input and output
        all_list.push_back(newnode);
        //Have the node specifier point to the node it generated
        (*curnode)->analogue=newnode;
    }

    //Create the links by iterating through the genes
    for(curgene=genes.begin();curgene!=genes.end();++curgene)
    {
        //Only create the link if the gene is enabled
        if (((*curgene)->enable)==true)
        {
            curlink=(*curgene)->lnk;
            inode=(curlink->in_node)->analogue;
            onode=(curlink->out_node)->analogue;

            //NOTE: This line could be run through a recurrency check if desired
            newlink=new Link(curlink->weight,inode,onode,curlink->is_recurrent);
            (onode->incoming).push_back(newlink); (inode->outgoing).push_back(newlink);
        }
    }

    //Create the new network
    newnet=new Network(inlist,outlist,all_list,genome_id);

    //Attach genotype and phenotype together
    newnet->genotype=this;

    if(phenotype != NULL)
        delete phenotype;
    phenotype=newnet;

    return newnet;
}


Genome *Genome::mutate(float sigma, int idR,GC idNewGenome,int &nodeId,int &gc)
{    
    Genome *new_genome;
    new_genome = this -> duplicate();

    //If the genome's ID is the same as the idNewGenome
    //mutation => just after crossover: no need to change the ID's
    if(!(this->genome_id == idNewGenome))
    {
        new_genome->genome_id = idNewGenome;
        new_genome->mom_id = this->genome_id;
        GC d = {idNewGenome.robot_id, -1};
        new_genome->dad_id = d;
    }
    //Choose the mutation depending on probabilities
    if (Helper::randFloat () < Helper::mutateAddNodeProb)
    {
        //Innovation numbers as gene clocks
        if(!(new_genome->mutate_add_node(Helper::newStructureTries,idR,nodeId,gc)))
        {
            //No node added, no connection found to split. Maybe try again?            
        }
    }
    else
    {   //links CANNOT be added directly after a node  because the phenotype
        // will not be appropriately altered to reflect the change
        if (Helper::randFloat() < Helper::mutateAddLinkProb)
        {
            if (!(new_genome->mutate_add_link(Helper::newStructureTries,idR,gc)))
            {
                //No link was added. Maybe all links already present
            }
        }
        else
        {
            if (Helper::randFloat() < Helper::mutateToggleEnableProb) {
                new_genome->mutate_toggle_enable(1);

            }
            else
            {
                //If no structural mutation, mutate weights
                if (Helper::randFloat() < Helper::mutateLinkWeightsProb){
                    new_genome->mutate_link_weights(sigma);
                }
            }
            /*if (Helper::randfloat()<Helper::mutate_gene_reenable_prob) {
                new_genome->mutate_gene_reenable();
            }*/
        }
    }
    return new_genome;
}

void Genome::mutate_link_weights(double power)
{    
    bool multisynapseAdaptedSigma = true;//false; //

    if(!Helper::allowMultisynapses || !multisynapseAdaptedSigma)
    {
        for(auto &curgene : genes)
        {
            if(curgene->enable)
            {
                if (Helper::randFloat() < Helper::mutateIndividualWeightProb)
                {
                    curgene->lnk->weight += getGaussianRand(0, power);
                    //curgene->lnk->weight = capWeights(curgene->lnk->weight);
                }
            }
        } //for
    }
    else //multisynapses (divide sigma std deviation by sqrt(number of synapses btw the same nodes)         
    {
         //TOTEST If multysynapses, then mutate only one of them (? the newest?).With full sigma
        /*std::set<std::pair<innov,innov>> alreadyMut;

        for(auto &curgene :  boost::adaptors::reverse(genes))
        {
            if(curgene -> enable)
            {
                if (Helper::randFloat() < Helper::mutateIndividualWeightProb)
                {
                    auto inNodeID = curgene->lnk->in_node->node_id;
                    auto outNodeID = curgene->lnk->out_node->node_id;

                    auto pair = std::make_pair(inNodeID,outNodeID);

                    if(alreadyMut.find(pair) == alreadyMut.end())
                    {
                        curgene->lnk-> weight += getGaussianRand(0, power);
                                    //curgene->lnk->weight = capWeights(curgene->lnk->weight);
                        alreadyMut.insert(pair);
                    }
                }
            }
        } //for*/

        for(auto &curgene : genes)
        {
            if(curgene -> enable)
            {
                int nb = getNumberSynapses(curgene);
                if (Helper::randFloat() < 1.0/(double)nb) //Helper::mutateIndividualWeightProb)
                {

                    curgene->lnk-> weight += getGaussianRand(0, power );// sqrt((double)nb) );
                              //curgene->lnk-> weight = capWeights(curgene->lnk-> weight);
                }
            }
        } //for
    }
}
int Genome::getNumberSynapses(Gene* gene)
{
    int result = 0;
    std::vector<Gene*>::iterator curgene;
    for(curgene=genes.begin();curgene!=genes.end();curgene++)
    {
        if((*curgene) -> enable)
        {
            if(gene->lnk->in_node == (*curgene)->lnk->in_node
                    && gene->lnk->out_node == (*curgene)->lnk->out_node)
            {
                result++;
            }

        }
    }
    return result;
}

void Genome::initialize_link_weights()
{
    std::vector<Gene*>::iterator curgene;
    //Loop on all genes
    for(curgene=genes.begin();curgene!=genes.end();curgene++)
    {
        if((*curgene) -> enable)
        {
            // weights: random init between -1 and +1
            ((*curgene)-> lnk) -> weight = (double)(rand() % 10000)/5000.0 - 1.0;
        }
    } //end for loop
}

double Genome::capWeights(double w)
{
    double result = w;
    double range = 2 * Helper::rangeW;
    // bouncing upper/lower bounds
    if ( result < -Helper::rangeW )
    {
        double overflow = - (result - -Helper::rangeW );
        overflow = overflow - 2*range * (int)( overflow / (2*range) );
        if ( overflow < range )
            result = -Helper::rangeW + overflow;
        else // overflow btw range and range*2
            result = -Helper::rangeW + range - (overflow-range);
    }
    else if ( result > Helper::rangeW )
    {
        double overflow = result - Helper::rangeW;
        overflow = overflow - 2 * range * (int)(overflow / (2*range));
        if ( overflow < range )
            result = Helper::rangeW - overflow;
        else // overflow btw range and range*2
            result = Helper::rangeW - range + (overflow-range);
    }
    /*double result = w;
    if(result < -Helper::rangeW)
        result = -Helper::rangeW;
    if(result > +Helper::rangeW)
        result = +Helper::rangeW;*/

    return result;
}
void Genome::mutate_toggle_enable(unsigned int times)
{   //TOTEST
    int genenum;
    std::vector<Gene*>::iterator thegene;  //Gene to toggle
    std::vector<Gene*>::iterator checkgene;  //Gene to check
    int genecount;

    for (unsigned int i=0; i < times; i++)
    {
        //Choose a random gene
        genenum=Helper::randInt(0,genes.size()-1);
        thegene=genes.begin();
        for(genecount=0;genecount<genenum;genecount++)
            ++thegene;

        //Toggle enable on this gene
        if (((*thegene)->enable)==true)
        {
            //make sure that another gene connects out of the in-node and out-node
            //if not: a section of network will become isolated
            checkgene=genes.begin();
            while( checkgene!=genes.end()
                   &&
                !(
                  ( ((*checkgene)->lnk)->in_node == ((*thegene)->lnk)->in_node )
                 && (((*checkgene)->enable))
                 && !((*checkgene)->innovation_num == (*thegene)->innovation_num))
                  &&
                   !(
                     ( ((*checkgene)->lnk)->out_node == ((*thegene)->lnk)->out_node )
                    && (((*checkgene)->enable))
                    && !((*checkgene)->innovation_num == (*thegene)->innovation_num)))
                ++checkgene;

            //Disable the gene if it's safe to do so
            if (checkgene!=genes.end())
                (*thegene)->enable=false;
        }
        else
            (*thegene)->enable=true;
    }
}
bool Genome::mutate_add_node(int tries,int idR,int &nodeId, int &gc)
{
    std::vector<Gene*>::iterator thegene;  //random gene containing the original link
    int genenum;  //Random gene number

    //Nodes connected by the gene
    NNode *in_node; NNode *out_node;
    //Link of random gene
    Link *thelink;

    Gene *newgene1; Gene *newgene2; NNode *newnode;

    //Weight of original link
    double oldweight;

    //Take a few tries to find an open node
    int trycount; bool found;

    //First, find a random gene already in the genome
    trycount=0; found=false;

    //Random uniform choice of genes
    while ((trycount<tries)&&(!found))
    {
        //This old totally random selection is bad- splitting
        //inside something recently splitted adds little power
        //to the system (should use a gaussian if doing it this way)
        genenum=Helper::randInt(0,genes.size()-1);

        //find the gene
        thegene=genes.begin();
        for(int genecount=0;genecount<genenum;genecount++)
            ++thegene;

        //If either the gene is disabled,
        //REMOVED FOR TEST or it has a bias input, try again***************************************************************************************************
        if (!(  ( (*thegene)->enable==false )
                //|| ((*thegene)->lnk->in_node->gen_node_label==BIAS)
                ))
            found=true;

        ++trycount;
    }

    //If we couldn't find anything so say goodbye
    if (!found)
        return false;
    //Disable the gene
    (*thegene)->enable=false;

    //Extract the link
    thelink=(*thegene)->lnk; oldweight=(*thegene)->lnk->weight;

    //Extract the nodes
    in_node=thelink->in_node; out_node=thelink->out_node;

    innov innovClock; innovClock.idR = idR; innovClock.gc = nodeId;

    //Create the new NNode (Hidden neuron)
    newnode=new NNode(NEURON,innovClock,HIDDEN);
    nodeId++;
    double preCoeff = 0.20311035;//res_250It_lim5 optimized with CMA-ES//1.0 otherwise
    double postCoeff = 5.0; // res_250It_lim5 //1.0 otherwise
    //Create the new Genes
    if (thelink->is_recurrent)
    {
        innovClock.idR = idR; innovClock.gc = gc;
        newgene1=new Gene(preCoeff,in_node,newnode,true,innovClock);
        //TODO test coeffs to keep mutation seamless
        innovClock.idR = idR; innovClock.gc = gc + 1;
        newgene2=new Gene(postCoeff*oldweight,newnode,out_node,false,innovClock);
        gc = gc + 2;
    }
    else
    {
        innovClock.idR = idR; innovClock.gc = gc;

        //? maybe oldweight to pre-connection!! Unlike NEAT. Needs to be applied before non-linearity
        //Attention: possible problem with postCoeff*oldweight not in boundaries
        newgene1=new Gene(preCoeff,in_node,newnode,false,innovClock);
        //coeffs to keep mutation seamless [0.203,5.0] SEAMLESS SPLIT NODE MUTATION
        innovClock.idR = idR; innovClock.gc = gc + 1;
        newgene2=new Gene(postCoeff*oldweight,newnode,out_node,false,innovClock);
        gc = gc + 2;
    }

    //Now add the new NNode and new Genes to the Genome
    //Add genes in correct order
    add_gene(genes,newgene1); add_gene(genes,newgene2);
    node_insert(nodes,newnode);
    return true;
} 

bool Genome::mutate_add_link(int tries,int idR,int &gc)
{
    //Random node numbers and iterators
    int nodenum1,nodenum2;
    std::vector<NNode*>::iterator thenode1,thenode2;

    int nodecount;  //Counter for finding nodes
    int trycount; //attempts to find an unconnected pair of nodes
    NNode* nodep1; NNode* nodep2; //Pointers to the nodes
    std::vector<Gene*>::iterator thegene; //Searches for existing link
    bool found=false;  //was an open pair was found?

    int recurflag; //is proposed link recurrent?
    Gene *newgene;  //new Gene

    double newweight;  //new weight for the new link

    //bool do_recur; bool loop_recur;
    int first_nonsensor;

    //used to avoid getting stuck in an infinite loop checking
    //for recursion
    //Note that we check for recursion to control the frequency of
    //adding recurrent links rather than to prevent any particular
    //kind of error
    int thresh=nodes.size()*nodes.size();
    int count=0;

    //Make attempts to find an unconnected pair
    trycount=0;

    //Decide whether to make this recurrent
    //TOCHECK is this really useful/necessary?
    /*if (Helper::randFloat() < Helper::recurOnlyProb)
        do_recur=true;
    else do_recur=false;*/

    //Find the first non-sensor so that the to-node won't look at sensors as
    //possible destinations
    first_nonsensor=0;
    thenode1=nodes.begin();
    while(((*thenode1)->get_type())==SENSOR)
    {
        first_nonsensor++;
        ++thenode1;
    }

    //Here is the recurrent finder loop- it is done separately
   /* if (do_recur)
    {
        while(trycount<tries)
        {
            //Some of the time try to make an auto recur loop
            //TOCHECK is this really useful/necessary?
            if (Helper::randFloat()<0.5)
            {
                loop_recur=true;
            }
            else
                loop_recur=false;

            if (loop_recur)
            {
                nodenum1=Helper::randInt(first_nonsensor,nodes.size()-1);
                nodenum2=nodenum1;
            }
            else
            {
                //Choose random nodenums
                nodenum1=Helper::randInt(0,nodes.size()-1);
                nodenum2=Helper::randInt(first_nonsensor,nodes.size()-1);
            }

            //Find the first node
            thenode1=nodes.begin();
            for(nodecount=0;nodecount<nodenum1;nodecount++)
                ++thenode1;

            //Find the second node
            thenode2=nodes.begin();
            for(nodecount=0;nodecount<nodenum2;nodecount++)
                ++thenode2;

            nodep1=(*thenode1);
            nodep2=(*thenode2);

            //See if an active link already exists  ALSO STOP AT END OF GENES!!
            //Don't allow SENSORS to get input
            thegene=genes.begin();
            while ( (thegene!=genes.end())
                    && ((nodep2->type)!=SENSOR)
                    &&
                    !(
                        (((*thegene)->lnk)->in_node==nodep1)&&
                        (((*thegene)->lnk)->out_node==nodep2) &&
                        (*thegene)->enable
                      )
                  )
                ++thegene;


            if (thegene!=genes.end())
                trycount++;
            else
            {
                count=0;
                recurflag=phenotype->is_recur(nodep1, nodep2, count, thresh);
                //recurflag=phenotype->is_recur(nodep1->analogue,nodep2->analogue,
                 //                             count,thresh);

                //ADDED: CONSIDER connections out of outputs recurrent
                if ((nodep1->gen_node_label)==OUTPUT)
                    recurflag=true;

                //Make sure it finds the right kind of link (recur)
                if (!(recurflag))
                    trycount++;
                else
                {
                    trycount=tries; //to exit loop
                    found=true;
                }
            }

        }
    }
    else
    {*/
        //Loop to find a link
        while(trycount<tries)
        {
            //Choose random nodenums
            nodenum1=Helper::randInt(0,nodes.size()-1);
            nodenum2=Helper::randInt(first_nonsensor,nodes.size()-1);

            //Find the first node
            thenode1=nodes.begin();
            for(nodecount=0;nodecount<nodenum1;nodecount++)
                ++thenode1;

            //Find the second node
            thenode2=nodes.begin();
            for(nodecount=0;nodecount<nodenum2;nodecount++)
                ++thenode2;

            nodep1=(*thenode1); nodep2=(*thenode2);

            //See if a link already exists  ALSO STOP AT END OF GENES!!!!
            //Don't allow SENSORS to get input

            thegene=genes.begin();
            bool isExistingConnection = ((*thegene)->lnk->in_node==nodep1)
                                          && ((*thegene)->lnk->out_node==nodep2)
                                          && (*thegene)->enable;
            while ( (thegene!=genes.end())
                   && (nodep2->type!=SENSOR)
                   && (!isExistingConnection || Helper::allowMultisynapses)
                  )
            {
                ++thegene;
                if(thegene!=genes.end())
                    isExistingConnection = ((*thegene)->lnk->in_node==nodep1)
                                              && ((*thegene)->lnk->out_node==nodep2)
                                              && (*thegene)->enable;
            }


            if (thegene!=genes.end())
                trycount++;
            else
            {
                count=0;
                //recurflag=phenotype->is_recur(nodep1,nodep2, count,thresh);
                recurflag=phenotype->is_recur(nodep1->analogue,nodep2->analogue,
                                              count,thresh);

                //ADDED: CONSIDER connections out of outputs recurrent
                if ((nodep1->gen_node_label)==OUTPUT)
                    recurflag=true;

                //TOERASE NO RECURRENT LINKS
                if (recurflag)
                    trycount++;
                else
                {
                    trycount=tries; //to exit loop
                    found=true;
                }
            }

        } //End of normal link finding loop
    //}
    //Continue only if an open link was found
    if (found)
    {
        //If it was supposed to be recurrent, make sure it gets labeled that way
        //?TOCHECK is this necessary? if (do_recur) recurflag=1;

        newweight= 0.0; //Helper::randPosNeg()*Helper::randFloat()*1.0;

        innov innovClock; innovClock.idR = idR; innovClock.gc = gc;

        //Create the new gene
        newgene=new Gene(newweight,nodep1,nodep2,recurflag,innovClock);
        gc++;

        //Now add the new Genes to the Genom in correct order
        add_gene(genes,newgene);
        return true;
    }
    else
        return false;
}


//Adds a new gene that has been created through a mutation in the
//*correct order* into the list of genes in the genome
void Genome::add_gene(std::vector<Gene*> &glist,Gene *g) 
{
    std::vector<Gene*>::iterator curgene;

    innov inum=g->innovation_num;

    curgene=glist.begin();

    while ((curgene!=glist.end())&&
           (((*curgene)->innovation_num)<inum))
        ++curgene;

    glist.insert(curgene,g);

}


void Genome::node_insert(std::vector<NNode*> &nlist,NNode *n) 
{
    std::vector<NNode*>::iterator curnode;
    innov id=n->node_id;

    curnode=nlist.begin();
    while ( (curnode!=nlist.end()) &&
            (((*curnode)->node_id) < id))
        ++curnode;

    nlist.insert(curnode,n);

}

Genome *Genome::mate(Genome *g,GC genomeid,double fitness1,double fitness2)
{
    //The baby Genome will contain these new NNodes and Genes
    std::vector<NNode*> newnodes; std::vector<Gene*> newgenes;
    Genome *new_genome;

    //Checks for link duplication
    std::vector<Gene*>::iterator curgene2;

    //Iterators for moving through the two parents' genes
    std::vector<Gene*>::iterator p1gene; std::vector<Gene*>::iterator p2gene;

    //Innovation numbers for genes inside parents' Genomes
    innov p1innov; innov p2innov;

    Gene *chosengene = NULL;  //Gene chosen for baby to inherit

    //NNodes connected to the chosen Gene
    NNode *inode; NNode *onode;
    NNode *new_inode; NNode *new_onode;

    //For checking if NNodes exist already
    std::vector<NNode*>::iterator curnode;

    /*//Set to true if we want to disabled a chosen gene
    bool disable;
    disable=false;*/

    Gene *newgene;

    //Tells if the first genome (this one) has better fitness or not
    bool p1better; bool skip;

    //Which genome is better?
    //Worse genome should not be allowed to add extra structural baggage
    //If they are the same, the smaller one is considered better
    if (fitness1>fitness2)
        p1better=true;
    else if (fitness1==fitness2)
    {
        if (genes.size()<(g->genes.size()))
            p1better=true;
        else p1better=false;
    }
    else p1better=false;

    //Make sure all sensors and outputs are included
    for(curnode=(g->nodes).begin();curnode!=(g->nodes).end();++curnode)
    {
        if ((((*curnode)->gen_node_label)==INPUT)||
                (((*curnode)->gen_node_label)==BIAS)||
                (((*curnode)->gen_node_label)==OUTPUT))
        {
            //Create a new node off the sensor or output
            new_onode=new NNode((*curnode));
            //Add the new node
            node_insert(newnodes,new_onode);
        }
    }
    //Now move through the Genes of each parent until both genomes end
    p1gene=genes.begin();
    p2gene=(g->genes).begin();
    while(!((p1gene==genes.end()) && (p2gene==(g->genes).end())))
    {
        //Default to not skipping a chosen gene
        skip=false;
        if (p1gene==genes.end())
        {
            chosengene=*p2gene; ++p2gene;
            //Skip excess from the worse genome
            if (p1better) skip=true;
        }
        else if (p2gene==(g->genes).end())
        {
            chosengene=*p1gene; ++p1gene;
            //Skip excess from the worse genome
            if (!p1better) skip=true;
        }
        else
        {
            //Extract current innovation numbers
            p1innov=(*p1gene)->innovation_num;
            p2innov=(*p2gene)->innovation_num;

            if (p1innov==p2innov)
            {
                if (Helper::randFloat()<0.5)
                    chosengene=*p1gene;
                else
                    chosengene=*p2gene;
                /*
                //If one is disabled, the corresponding gene in the offspring
                //will likely be disabled
                if ((((*p1gene)->enable)==false)||
                        (((*p2gene)->enable)==false))
                    if (Helper::randFloat()<0.75) disable=true;
                */
                ++p1gene; ++p2gene;

            }
            else if (p1innov<p2innov)
            {
                chosengene=*p1gene; ++p1gene;

                if (!p1better) skip=true;
            }
            else if (p2innov<p1innov)
            {
                chosengene=*p2gene; ++p2gene;
                if (p1better) skip=true;
            }
        }

        //Check to see if the chosengene conflicts with an already chosen gene
        //i.e. do they represent the same link
        curgene2=newgenes.begin();
        while ((curgene2!=newgenes.end())&&
               (!((((((*curgene2)->lnk)->in_node)->node_id)==((((chosengene)->lnk)->in_node)->node_id))&&
                  (((((*curgene2)->lnk)->out_node)->node_id)==((((chosengene)->lnk)->out_node)->node_id))&&
                  ((((*curgene2)->lnk)->is_recurrent)== (((chosengene)->lnk)->is_recurrent)) ))&&
               (!((((((*curgene2)->lnk)->in_node)->node_id)==((((chosengene)->lnk)->out_node)->node_id))&&
                  (((((*curgene2)->lnk)->out_node)->node_id)==((((chosengene)->lnk)->in_node)->node_id))&&
                  (!((((*curgene2)->lnk)->is_recurrent)))&&
                  (!((((chosengene)->lnk)->is_recurrent))) )))
            ++curgene2;
        //Links conflicts, abort adding
        if (curgene2!=newgenes.end()) skip=true;

        if (!skip)
        {
            //Now add the chosengene to the baby, next check for the nodes, add them
            //if not in the baby Genome already
            inode=(chosengene->lnk)->in_node;
            onode=(chosengene->lnk)->out_node;

            //Check for inode in the newnodes list
            if (inode->node_id<onode->node_id)
            {
                //inode before onode
                //Checking for inode's existence
                curnode=newnodes.begin();
                while((curnode!=newnodes.end())
                      && (!((*curnode)->node_id==inode->node_id)))
                    ++curnode;

                if (curnode==newnodes.end())
                {
                    //Here we know the node doesn't exist so we have to add it
                    new_inode=new NNode(inode);
                    node_insert(newnodes,new_inode);
                }
                else
                {
                    new_inode=(*curnode);
                }

                //Checking for onode's existence
                curnode=newnodes.begin();
                while((curnode!=newnodes.end())&&
                      (!((*curnode)->node_id==onode->node_id)))
                    ++curnode;
                if (curnode==newnodes.end())
                {
                    //Here we know the node doesn't exist so we have to add it
                    new_onode=new NNode(onode);

                    node_insert(newnodes,new_onode);
                }
                else
                {
                    new_onode=(*curnode);
                }

            }
            //If the onode has a higher id than the inode we want to add it first
            else
            {
                //Checking for onode's existence
                curnode=newnodes.begin();
                while((curnode!=newnodes.end())&&
                      (!((*curnode)->node_id==onode->node_id)))
                    ++curnode;
                if (curnode==newnodes.end())
                {
                    //Here we know the node doesn't exist so we have to add it
                    new_onode=new NNode(onode);

                    node_insert(newnodes,new_onode);
                }
                else
                {
                    new_onode=(*curnode);

                }

                //Checking for inode's existence
                curnode=newnodes.begin();
                while((curnode!=newnodes.end())&&
                      (!((*curnode)->node_id==inode->node_id)))
                    ++curnode;
                if (curnode==newnodes.end())
                {
                    //Here we know the node doesn't exist so we have to add it
                    new_inode=new NNode(inode);
                    node_insert(newnodes,new_inode);
                }
                else
                {
                    new_inode=(*curnode);

                }
            } //End NNode checking section- NNodes are now in new Genome

            //Add the Gene
            newgene=new Gene(chosengene,new_inode,new_onode);
            //if (disable){ newgene->enable=false; disable=false;} //TOCHECK if recover toggle_enable
            newgenes.push_back(newgene);
        }
    }

    new_genome=new Genome(genomeid,newnodes,newgenes);

    new_genome->mom_id = genome_id;
    new_genome->dad_id = g->genome_id;
    //Return the baby Genome
    return (new_genome);
}

double Genome::dissimilarity(Genome *g)
{
    double result = 0.0;
    double exc = 0.0 , disj = 0.0, weight = 0.0;

    std::vector<Gene*>::iterator it1 = genes.begin();
    std::vector<Gene*>::iterator it2 = g->genes.begin();

    int maxLength; int common = 0;
    if(genes.size() < g->genes.size())
        maxLength = g->genes.size();
    else
        maxLength = genes.size();

    bool done = (it1 == genes.end()) && (it2 == g -> genes.end());

    while(!done)
    {
        if( (it2 == g -> genes.end()) && !(it1 == genes.end()))
        {
            //There are still genes of the first genome
            exc += 1.0;
            it1++;
        }
        else if( (it1 == genes.end()) && !(it2 == g -> genes.end()))
        {
            //There are still genes of the second genome
            exc += 1.0;
            it2++;
        }
        else //There are still genes of both genomes
        {
            Gene* g1 = *it1; Gene* g2 = *it2;

            if((g1 -> innovation_num) == (g2 -> innovation_num))
            {
                if(g1 -> enable && g2 ->enable)
                {
                    weight += fabs( (g1->lnk)->weight
                                    - (g2->lnk)->weight );
                    common++;
                }
                it1++;
                it2++;
            }
            else if((g1 -> innovation_num) < (g2 -> innovation_num))
            {
                disj += 1.0; it1++;
            }
            else if((g1 -> innovation_num) > (g2 -> innovation_num))
            {
                disj += 1.0; it2++;
            }
            else
            {
                std::cerr << "[ERROR] Innovation numbers not comparable (?)"
                          << std::endl;
                exit(-1);
            }
        }

        done = (it1 == genes.end()) && (it2 == g->genes.end());
    }

    result = (Helper::coefE * exc/maxLength)
            + (Helper::coefD * disj/maxLength)
            + (Helper::coefW * weight /common);

    //std::cout << "E: " << (Helper::coefE * exc/maxLength)
    //          << ", D: " << (Helper::coefD * disj/maxLength)
    //            << ", W: " << (Helper::coefW * weight /common) << std::endl;
    //std::cout << "E: " << (exc/maxLength)
    //          << ", D: " << (disj/maxLength)
    //            << ", W: " << (weight /common) << std::endl;
    /*std::cout << "E: " << exc
                  << ", D: " << disj
                    << ", W: " << weight/common << std::endl;*/
    //std::cout << result << std::endl;
    return result;
}

void Genome::print_to_file(std::ostream &outFile)
{
    std::vector<NNode*>::iterator curnode;
    std::vector<Gene*>::iterator curgene;
    outFile <<"genomestart " << "R" <<genome_id.robot_id << ", G"
            << genome_id.gene_id <<std::endl;

    //Output the nodes
    for(curnode=nodes.begin();curnode!=nodes.end();++curnode) {
        (*curnode)->print_to_file(outFile);
    }
    //Output the genes
    for(curgene=genes.begin();curgene!=genes.end();++curgene) {
        (*curgene)->print_to_file(outFile);
    }
    outFile << "genomeend " << "R" <<genome_id.robot_id << ", G" << genome_id.gene_id << std::endl;
}

void Genome::print_to_filename(const char *filename)
{
    std::ofstream oFile(filename);
    print_to_file(oFile);
    oFile.close();
}
