bool doTest = true;
          if(doTest)
          {
              start_genome = Genome::makeGenome(id,2,1);
              Network* n = start_genome->genesis(1);
              start_genome->genes[0]->lnk->weight = 0.0;
              start_genome->genes[1]->lnk->weight = 1.0;

              int nodeid = 4;
              double geneid = 3;
              if(!addNode(start_genome,1,3,nodeid,geneid))
              {
                  std::cerr << "[ERROR] Wrong in out neuron ids" << std::endl;
                  exit(-1);
