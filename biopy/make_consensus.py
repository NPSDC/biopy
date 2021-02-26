from ITrees import Tree, NodeData
import sys, random, os
import Nodes
import biopy.parseNewick

def read_bparts(bpart_file, write_file, min_group_len=10):
    #groups_tree=dict()
    if not os.path.exists(bpart_file):
        sys.exit("Invalid partition split file")
    f=open(write_file,"w")
    groups_criteria = 0
    total_groups = 0
    with open(bpart_file) as ifile:
        while(1):
            group_info = ifile.readline().rstrip().split("\t")
            if(group_info==[""]):
                break
            if len(group_info) != 2:
                print(group_info)
                sys.exit("Invalid name in groups")
            group = group_info[0].rstrip()
            alltaxa = set(group.split("_"))
            
            n_part = int(group_info[1].rstrip())
            clades = dict()
            #print(group_info)
            total_groups += 1                
            for i in range(n_part):
                bpart = ifile.readline().rstrip().split("\t")
                txps = bpart[0].rstrip().split("_")
                clades[str(txps)] = int(bpart[1].rstrip())
            
            if len(alltaxa) >= min_group_len:
                if(len(alltaxa) != len(group.split("_"))):
                    print(group,alltaxa)
                groups_criteria += 1
                
                f.write(get_newick_consensus(clades, NodeData, alltaxa).replace(":0.0","")+";\n")
    print("Total number of groups that meet criteria", groups_criteria)
    print("Total number of groups", total_groups)
    f.close()

def get_newick_consensus(clades, dataclass, alltaxa, threshold = 0.5, outgroup=None):
    consensus=Tree(name='consensus_%2.1f' % float(threshold),data=dataclass)
    for (c,s) in clades.items():
        node = Nodes.Node(data=dataclass())
        node.data.support=s
        node.data.taxon=set(eval(c))
        consensus.add(node)
    consensus.node(consensus.root).data.support=None
    consensus.node(consensus.root).data.taxon=alltaxa
    # we sort the nodes by no. of taxa in the clade, so root will be the last
    consensus_ids=consensus.all_ids()
    consensus_ids.sort(lambda x,y:len(consensus.node(x).data.taxon)-len(consensus.node(y).data.taxon))
    # now we just have to hook each node to the next smallest node that includes all taxa of the current 
    for i,current in enumerate(consensus_ids[:-1]): # skip the last one which is the root
        #print '----'
        #print 'current: ',consensus.node(current).data.taxon
        # search remaining nodes
        for parent in consensus_ids[i+1:]:
            #print 'parent: ',consensus.node(parent).data.taxon
            if consensus.node(parent).data.taxon.issuperset(consensus.node(current).data.taxon):
                break
        else:
            sys.exit('corrupt tree structure?')
        # internal nodes don't have taxa
        if len(consensus.node(current).data.taxon)==1:
            consensus.node(current).data.taxon=consensus.node(current).data.taxon.pop()
            # reset the support for terminal nodes to maximum
            #consensus.node(current).data.support=max_support
        else:
            consensus.node(current).data.taxon=None
        consensus.link(parent,current)
    # eliminate root taxon name
    consensus.node(consensus_ids[-1]).data.taxon=None 
    return consensus.toNewick()

def main():
    bpart_file="../../terminus_ase_noor/data/term/Sample1/cluster_bipart_splits.txt"
    write_file="../../terminus_ase_noor/data/term/Sample1/consensus_splits.txt"
    read_bparts(bpart_file, write_file)

if __name__ == "__main__":
    main()