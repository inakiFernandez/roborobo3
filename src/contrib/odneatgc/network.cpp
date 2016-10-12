#include "odneatgc/network.h"
#include "odneatgc/link.h"
#include "odneatgc/helper.h"

#include <iostream>
#include <sstream>
#include <iomanip>

using namespace ODNEATGC;

Network::Network(std::vector<NNode*> in,std::vector<NNode*> out,std::vector<NNode*> all,GC netid)
{
    inputs=in; outputs=out; all_nodes=all; numnodes=-1; numlinks=-1; net_id=netid;
}

Network::Network(GC netid)
{
    numnodes=-1; numlinks=-1; net_id=netid;
}

Network::~Network() 
{
    // Kill off all the nodes and links
    destroy();
    //Nodes and links already freed
}

// Puts the network back into an initial state
void Network::flush() 
{
    for(std::vector<NNode*>::iterator curnode=outputs.begin();curnode!=outputs.end();++curnode)
        (*curnode)->flushback();
}

// Debugger: Checks network state
void Network::flush_check() 
{
    std::vector<NNode*>::iterator curnode;
    std::vector<NNode*>::iterator location;
    std::vector<NNode*> seenlist;  //List of nodes not to doublecount

    for(curnode=outputs.begin();curnode!=outputs.end();++curnode)
    {
        location= std::find(seenlist.begin(),seenlist.end(),(*curnode));
        if (location==seenlist.end())
        {
            seenlist.push_back(*curnode);
            (*curnode)->flushback_check(seenlist);
        }
    }
}

// If all outputs are not active then return true
bool Network::outputsoff() 
{
    for(std::vector<NNode*>::iterator curnode=outputs.begin();curnode!=outputs.end();++curnode)
        if (((*curnode)->activation_count)==0)
            return true;
    return false;
}

bool Network::hiddenoff()
{
    for(auto curnode=all_nodes.begin();curnode!=all_nodes.end();curnode++)
    {
        //Ignore SENSORS
        if ((*curnode)->type!=SENSOR)
        {
            if (!(*curnode)->active_flag)
                return true;
        }
    }
    return false;
}
// Activates the net such that all outputs are active. Returns true on success
bool Network::activate() 
{
    std::vector<NNode*>::iterator curnode;
    std::vector<Link*>::iterator curlink;
    double add_amount;
    bool onetime; //ensure at least activate once
    int abortcount=0;  //in case the output is truncated from the network

    //Keep activating until all the outputs have become active
    //(This only happens on the first activation, because after that, they
    // are always active)
    onetime=false;
    for(curnode=all_nodes.begin();curnode!=all_nodes.end();curnode++)
    {
        if ((*curnode)->type!=SENSOR)
             (*curnode)->active_flag=false;
    }
    while(outputsoff()||hiddenoff()||!onetime)
    {
        //std::cout << "ABORTCOUNT " << abortcount << std::endl;
        ++abortcount;
        if (abortcount==200)
            return false;
        //TOCHECK recurrent connections??
        // For each node, compute the sum of its incoming activation
        for(curnode=all_nodes.begin();curnode!=all_nodes.end();curnode++)
        {
            //Ignore SENSORS
            if ((*curnode)->type!=SENSOR)
            {
                //std::cout << "Cur."<<(*curnode)->gen_node_label << ", nodeId: " << (*curnode)->node_id.gc << std::endl;
                (*curnode)->activesum=0.0;
                //does it have any active inputs?

                (*curnode)->active_flag=true; //false; //

                //For each incoming connection, add its activity to the activesum
                for(curlink=((*curnode)->incoming).begin();
                    curlink!=((*curnode)->incoming).end();
                    ++curlink)
                {
                    if(!(*curlink)->is_recurrent)
                    {
                        if ((*curlink)->in_node->active_flag||(*curlink)->in_node->type==SENSOR)
                        {
                            //std::cout << "In."<<(*curlink)->in_node->gen_node_label <<
                              //           ", nodeId: " << (*curlink)->in_node->node_id.gc << std::endl;
                            double w = (*curlink)->weight;
                            double preActivation = (*curlink)->in_node->get_active_out();

                            add_amount=w * preActivation;

                            //std::cout <<"activation: " << preActivation << std::endl;

                            (*curnode)->activesum+=add_amount;
                        }
                        else
                        {
                            (*curnode)->active_flag=false;
                            //std::cout << "Bad. Missing activation: " << (*curlink)->in_node->node_id.gc << std::endl;
                        }
                    }
                    else
                    {
                        double w = (*curlink)->weight;
                        double preActivation = (*curlink)->in_node->get_active_out();

                        add_amount=w * preActivation;

                        (*curnode)->activesum+=add_amount;
                    }

                } //End for incoming links
                //Only activate if some active input came in
                if ((*curnode)->active_flag)
                {
                    //Now run the net activation through an activation function
                    if ((*curnode)->ftype==SIGMOID)
                        (*curnode)->activation=tanh((*curnode)->activesum);
                    else
                    {
                        std::cerr << "[ERROR] No valid activation function defined"
                                  << std::endl; exit(-1);
                    }
                    //Increment the activation_count, First activation cannot be from nothing!
                    (*curnode)->activation_count++;
                    //std::cout << "final: " << (*curnode)->activation << std::endl;

                }
            } //End if !=SENSOR

        } //End for all nodes


        // Activation function for all the non-sensor nodes
        /*for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode)
        {
            if ((*curnode)->type!=SENSOR)
            {
                //Only activate if some active input came in
                if ((*curnode)->active_flag)
                {
                    //Now run the net activation through an activation function
                    if ((*curnode)->ftype==SIGMOID)                        
                        (*curnode)->activation=tanh((*curnode)->activesum);
                    else
                    {
                        std::cerr << "[ERROR] No valid activation function defined"
                                  << std::endl; exit(-1);
                    }
                    //Increment the activation_count, First activation cannot be from nothing!
                    (*curnode)->activation_count++;

                }
            }
        }*/
        onetime=true;
        //std::cout <<"\n\n\n\n\n\n" <<std::endl;
    }
    return true;
}

// Add an input
void Network::add_input(NNode *in_node)
{
    inputs.push_back(in_node);
}

// Add an output
void Network::add_output(NNode *out_node)
{
    outputs.push_back(out_node);
}

// Takes an array of sensor values and loads it into SENSOR inputs ONLY
void Network::load_sensors(double *sensvals) 
{
    std::vector<NNode*>::iterator sensPtr;
    for(sensPtr=inputs.begin();sensPtr!=inputs.end();++sensPtr)
    {
        //only load values into SENSORS (not BIASes)
        if (((*sensPtr)->type)==SENSOR)
        {
            (*sensPtr)->sensor_load(*sensvals);
            sensvals++;
        }
    }
}

void Network::load_sensors(const std::vector<float> &sensvals) 
{

    std::vector<NNode*>::iterator sensPtr;
    std::vector<float>::const_iterator valPtr;

    for(valPtr = sensvals.begin(), sensPtr = inputs.begin();
        sensPtr != inputs.end() && valPtr != sensvals.end();
        ++sensPtr, ++valPtr)
    {
        //only load values into SENSORS (not BIASes)
        if (((*sensPtr)->type)==SENSOR)
        {
            (*sensPtr)->sensor_load(*valPtr);
        }
    }
}

// The following two methods recurse through a network from outputs
// down in order to count the number of nodes and links in the network.
// This can be useful for debugging genotype->phenotype spawning 
// (to make sure their counts correspond)

int Network::nodecount()
{
    int counter=0;
    std::vector<NNode*>::iterator curnode; std::vector<NNode*>::iterator location;
    std::vector<NNode*> seenlist;  //List of nodes not to doublecount

    for(curnode=outputs.begin();curnode!=outputs.end();++curnode)
    {
        location = std::find(seenlist.begin(),seenlist.end(),(*curnode));
        if (location==seenlist.end())
        {
            counter++; seenlist.push_back(*curnode);
            nodecounthelper((*curnode),counter,seenlist);
        }
    }
    numnodes=counter;
    return counter;
}

void Network::nodecounthelper(NNode *curnode,int &counter,std::vector<NNode*> &seenlist)
{
    std::vector<Link*> innodes=curnode->incoming;
    std::vector<Link*>::iterator curlink; std::vector<NNode*>::iterator location;

    if (!((curnode->type)==SENSOR))
    {
        for(curlink=innodes.begin();curlink!=innodes.end();++curlink)
        {
            location = std::find(seenlist.begin(),seenlist.end(),((*curlink)->in_node));
            if (location==seenlist.end())
            {
                counter++;
                seenlist.push_back((*curlink)->in_node);
                nodecounthelper((*curlink)->in_node,counter,seenlist);
            }
        }
    }
}

int Network::linkcount()
{
    int counter=0;
    std::vector<NNode*>::iterator curnode;
    std::vector<NNode*> seenlist;  //List of nodes not to doublecount

    for(curnode=outputs.begin();curnode!=outputs.end();++curnode)
        linkcounthelper((*curnode),counter,seenlist);

    numlinks=counter;
    return counter;
}

void Network::linkcounthelper(NNode *curnode,
                              int &counter,std::vector<NNode*> &seenlist)
{
    std::vector<Link*> inlinks=curnode->incoming;
    std::vector<Link*>::iterator curlink;
    std::vector<NNode*>::iterator location;

    location = std::find(seenlist.begin(),seenlist.end(),curnode);
    if ((!((curnode->type)==SENSOR))&&(location==seenlist.end()))
    {
        seenlist.push_back(curnode);

        for(curlink=inlinks.begin();curlink!=inlinks.end();++curlink)
        {
            counter++;
            linkcounthelper((*curlink)->in_node,counter,seenlist);
        }
    }
}

// Destroy will find every node in the network and 
// delete them one by one.  Since deleting a node deletes its incoming
// links, all nodes and links associated with a network will be destructed
void Network::destroy() 
{
    std::vector<NNode*>::iterator curnode;
    // Erase all nodes from all_nodes list
    for(curnode=all_nodes.begin();curnode!=all_nodes.end();++curnode)
        delete (*curnode);
}

// This checks a POTENTIAL link between a potential in_node and 
//potential out_node to see if it must be recurrent 
bool Network::is_recur(NNode *potin_node,NNode *potout_node,int &count,int thresh) 
{
    std::vector<Link*>::iterator curlink;
    ++count;  //Count the node as visited

    if (count>thresh)
        return false; //Short out the whole thing- loop detected

    if (potin_node==potout_node)
        return true;
    else
    {
        //Check back on all links...
        for(curlink=(potin_node->incoming).begin();curlink!=(potin_node->incoming).end();curlink++)
        {
            //But skip links that are already recurrent
            //(We want to check back through the forward flow of signals only
            if (!((*curlink)->is_recurrent))
            {
                if (is_recur((*curlink)->in_node,potout_node,count,thresh))
                    return true;
            }
        }
        return false;
    }
}

//Find the maximum number of neurons between an ouput and an input
int Network::max_depth() 
{
    std::vector<NNode*>::iterator curoutput; //The current output we are looking at
    int cur_depth; //The depth of the current node
    int max=0; //The max depth

    for(curoutput=outputs.begin();curoutput!=outputs.end();curoutput++)
    {
        cur_depth=(*curoutput)->depth(0,this);
        if (cur_depth>max)
            max=cur_depth;
    }
    return max;
}
