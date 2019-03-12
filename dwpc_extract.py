# -*- coding: utf-8 -*-
"""
@author: C Heiser
Mar 2019

usage: dwpc_extract.py [-h] -n NODES [NODES ...] outfile

Generate DWPC data between input nodes in Hetionet V1.0

positional arguments:
  outfile               name of output .csv file to write to cwd

optional arguments:
  -h, --help            show this help message and exit
  -n NODES [NODES ...], --nodes NODES [NODES ...]
                        list of node names to compare

"""
import numpy as np
import pandas as pd
import argparse
# get all the hetio functions
import random
from hetio.readwrite import *
from hetio.pathtools import *
from hetio.stats import *


def get_nodes_by_name(graph, nodes):
    out_nodes = [] # init list of nodes of interest
    for n in graph.get_nodes():
        if n.name in nodes:
            out_nodes.append(n)

    return out_nodes

def get_all_DWPCs(graph, node_list, possible_metapaths):
    '''
    return pd.DataFrame with DWPC between all combinations of nodes in node_list for paths of len<=max_length
    '''
    paths = {'pair':[], 'source':[], 'target':[], 'metapath':[], 'paths':[], 'DWPC':[]} # init output dictionary
    for pair in list(itertools.combinations(node_list, 2)): # iterate through possible pairs of nodes of interest
        print('Analyzing paths between nodes {} and {}'.format(pair[0].get_id()[1], pair[1].get_id()[1])) # status update
        for meta in possible_metapaths: # iterate through metapaths of len<=max_length between source and target genes
            path = paths_between(graph=graph, source=pair[0], target=pair[1], metapath=meta) # get all paths between source and target of metapath type
            if len(path)!=0: # if that metapath exists between the source and target genes, append results to dictionary
                paths['pair'].append(pair)
                paths['source'].append(pair[0])
                paths['target'].append(pair[1])
                paths['metapath'].append(meta)
                paths['paths'].append(path)
                paths['DWPC'].append(DWPC(path, damping_exponent=0.4))

    return pd.DataFrame(paths)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generate DWPC data between input nodes in Hetionet V1.0')
    parser.add_argument('outfile', help='name of output .csv file to write to cwd')
    parser.add_argument('-n', '--nodes', help='list of node names to compare', required=True, nargs='+', type=str)
    args = parser.parse_args()

    # Read Hetionet v1.0
    url = 'https://github.com/dhimmel/hetionet/raw/{}/{}'.format(
        '00bf0b6f8886821d91cfdf00eadad145a7a1b6da',
        'hetnet/json/hetionet-v1.0.json.bz2',
    )
    my_graph = read_graph(url) # read actual graph into graph object
    my_metagraph = my_graph.metagraph # read metagraph into metagraph object

    # get list of possible Gene-to-Gene metapaths up to len==max_length for downstream calcs
    pos_metapaths = my_metagraph.extract_metapaths(source='Gene', target='Gene', max_length=3)

    # 'query 1' from Erin:
    #genes1 = ["TFE3","TYK2","CPT1A","NUCB1","ENTPD4","DDX17","KLC1","JPH4","FAM214B","WDR48","CPT2"]
    # genes not in HetioNet (that I know of) = ["H2A","H2B1"]

    # get list of nodes from args
    in_nodes = get_nodes_by_name(graph=my_graph, nodes=args.nodes)

    # calculate DWPCs
    out = get_all_DWPCs(graph=my_graph, node_list=in_nodes, possible_metapaths=pos_metapaths)

    # write results to .csv file
    out.to_csv(args.outfile)
