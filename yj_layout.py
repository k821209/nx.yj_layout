import math
import networkx as nx
import numpy as np
from tqdm import tqdm
from random import randint

def YJlayout(gnx_in): 
    #= description
    #= select largest connected Vertexcluster and draw layout using defined rule
	#= In principle, the rule tries to strech out high degree nodes from each other while low degree nodes stay near high degree hubs. 
	#= The edge length is proportional to min(node1.degree, node2.degree), and inverse proportional to abs(node1.degree - node2.degree)
    #= rule description can be found in document. 
	#= genes,pos,sizes = yj_layout.YJlayout(gnx), where gnx is graph from networkx
	
	
    #+ equivalent igraph.clusters()
    clusters = []
    for node in gnx_in.node.keys():
        try:
            if node in set(clusters[-1]):
                pass
            else:
                if len(nx.shortest_path(gnx_in,node).keys()) > 2:
                    clusters.append(sorted(nx.shortest_path(gnx_in,node).keys()))
        except IndexError:
            clusters.append(sorted(nx.shortest_path(gnx_in,node).keys()))
    clusters.sort(key=lambda x : len(x),reverse=True)
    #-

    each_cluster_genename_list = clusters[0]
    each_cluster_genename_list.sort(key=lambda x : gnx_in.degree(x),reverse=True)

    
    def cos_modify(x):
        return 1-np.cos(x*np.pi)
    def sin_modify(x,c):
        return np.sin(c*x*np.pi/2)
    def sigmoid(x):
        return 1 / (1 + math.exp(-x))

    max_d         = 50
    min_d         = 10
    rad_size      = 1
    start_x       = 0
    start_y       = 0

    genename        =  each_cluster_genename_list[0] #'AT2G01570'
    dicG2pos        = {}
    def G2pos(genename,option=(start_x,start_y,rad_size)):
        xaxis, yaxis, irad_size  = option
        gravity                  = gnx_in.degree(genename) 

        interactors_list         = gnx_in.neighbors(genename)
        interactors_gravity_list = [gnx_in.degree(x) for x in interactors_list]
        edge_l_list              = ['%0.3f'%(float(min(gravity,gravity_i)) /((abs(gravity - gravity_i)+1))) if min(gravity,gravity_i) > 1 else 0 for gravity_i in interactors_gravity_list]
        dicG2L = dict(zip(interactors_list,edge_l_list))
        #if genename == 'AT4G10930':
        #    print(dicG2L['AT1G73790'])
        #    print(interactors_list)
        #    print(interactors_gravity_list)
        #    print(edge_l_list)
        dicEdge2theta_interval   = dict(zip(list(set(edge_l_list)),[2*np.pi/edge_l_list.count(x) for x in set(edge_l_list)])) 
        dicEdge2interactors      = {}
        for i, edge_l in enumerate(edge_l_list):
            try:
                dicEdge2interactors[edge_l].append(interactors_list[i])
            except KeyError:
                dicEdge2interactors[edge_l] = [interactors_list[i]]
        dicEdge2interactors_keys =  dicEdge2interactors.keys()
        dicEdge2interactors_keys.sort()
        for j,edge_l in enumerate(dicEdge2interactors_keys):
            theta_interval = dicEdge2theta_interval[edge_l]
            interactors    = dicEdge2interactors[edge_l]
            for i, interactor in enumerate(interactors):
                theta = (theta_interval*j)/3.2 + (theta_interval*i)
                try:
                    (x1,y1,size) = dicG2pos[interactor]
                except KeyError:
                    l                    = float(edge_l)*max_d + min_d
                    x1                   = xaxis + l*np.cos(theta)
                    y1                   = yaxis + l*np.sin(theta)
                    size                 = float(irad_size + 0.01*theta_interval + 2*gnx_in.degree(interactor))#irad_size*(sigmoid(theta_interval)) 
                    dicG2pos[interactor] = (x1,y1,size)


        if option == (start_x,start_y,rad_size):
            dicG2pos[genename] = (xaxis,yaxis,rad_size*2)
        return 1
    G2pos(genename)
    todolist = each_cluster_genename_list
    done     = []
    while len(dicG2pos.keys()) != len(each_cluster_genename_list):
        b = 0
        for genename in tqdm(todolist):
            if genename in done:
                continue
            try:
                b += G2pos(genename, dicG2pos[genename])
                done.append(genename)
            except KeyError:
                pass
        if b == 0:
            break


    genes = clusters[0]
    pos   = {k:(v[0],v[1]) for k,v in dicG2pos.iteritems()}
    sizes  = [v[2] for k,v in dicG2pos.iteritems()]
    return genes,pos,sizes