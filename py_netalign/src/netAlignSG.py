##### Arif Khan #####
### arif.khan@pnnl.gov ###

import warnings
import igraph as ig
import networkx as nx
import collections
import scipy as sp
import numpy as np
import bMatching as bm
import  Evaluation as ev
import importlib
import matplotlib.pylab as plt
import time
import os
import subprocess
import pandas as pd
import math

def read_graphml(f1,f2,fl=None,input_dir=None):
    
    if input_dir is not None:
        f1=input_dir+f1
        f2=input_dir+f2
        
    
    G1 = ig.read(f1,format="graphml")
    G2 = ig.read(f2,format="graphml")
    
    L=None
    if fl is not None:
        if input_dir is not None:
            fl=input_dir+fl
        L=ig.read(fl,format="graphml")
    
    return G1,G2,L


def read_networks_edgelist(f1,f2,fl=None,input_dir=None):
    
    if input_dir is not None:
        f1=input_dir+f1
        f2=input_dir+f2
        
    
    G1 = ig.Graph.Read_Edgelist(f1)   
    G2 = ig.Graph.Read_Edgelist(f2)
    
    L=None
    if fl is not None:
        if input_dir is not None:
            fl=input_dir+fl
        L=ig.read(fl,format="graphml")
    
    return G1,G2,L

def read_networks(f1,f2,fl=None,input_dir=None):
    
    if input_dir is not None:
        f1=input_dir+f1
        f2=input_dir+f2
        
    T=nx.read_leda(f1)
    nx.write_graphml(T,'graph.graphml')
    G1 = ig.read('graph.graphml',format="graphml")
    
    T=nx.read_leda(f2)
    nx.write_graphml(T,'graph.graphml')
    G2 = ig.read('graph.graphml',format="graphml")
    
    L=None
    if fl is not None:
        if input_dir is not None:
            fl=input_dir+fl
        L=ig.read(fl,format="graphml")
    
    return G1,G2,L


def get_alignment(G1,G2,M):
    
    n1=G1.vcount()
    n2=G2.vcount()
    Align=set([])
   
    for (i,j) in M:    
        #print(i,j)
        u=G1.vs[i]['id']
        v=G2.vs[j]['id'] 
        #print(u,v)
        Align.add((u,v))
   
    return Align 


def save_alignment(L,filename):
    
    fl=open(filename,"w")
    for (u,v) in L:
        t=u+" "+v+"\n"
        fl.write(t)
    
    fl.close()
    

def evaluate(AQ,G1,G2,AL):

    AQ.mapping_set=AL
    if AQ.true_mapping_set == None:
        qual=AQ.evaluate(False,False,False,False,False,True,False,False,False,False)
        score=qual["NCV-GS3"]
    else:
        qual=AQ.evaluate(False,False,True,False,False,True,False,False,False,False)
        score=float(qual["F-NC"])

    return score, qual


def igraph2mtx(G,n,fname='graph.mtx',bipartite=False):
    c1='%%MatrixMarket matrix coordinate real general'
    c2='% Generated'
    
    edgeG=G.get_edgelist()
    #print(edgeG)
    m=int(len(edgeG))
    
    try:
        weights=G.es['weight']
    except:
        weights=[1]*m

    
    if bipartite:
        nl=n
        nr=G.vcount()-nl
        h=str(nl)+" "+str(nr)+" "+str(m)
        
    else:
        nl=n
        nr=n
        h=str(nl)+" "+str(nr)+" "+str(m*2)
    
    f=open(fname,'w')
    
    f.write(c1+"\n")
    f.write(c2+"\n")
    f.write(h+"\n")
    #print(c1)
    #print(c2)
    #print(h)
    for i in range(m):
        (u,v)=edgeG[i]
        w=weights[i]
        
        if bipartite:
            v=v-nl
            #print(u,v,w)
            f.write(str(u+1)+" "+str(v+1)+" "+str(w)+"\n")
        else:
            #print(u,v,w)
            #print(v,u,w)
            f.write(str(u+1)+" "+str(v+1)+" "+str(w)+"\n")
            f.write(str(v+1)+" "+str(u+1)+" "+str(w)+"\n")
    
    f.close()
    

def netAlignSC(f1,f2,fl,alpha=1,beta=1,ft=None,input_dir=None):

    G1,G2,L=read_networks_edgelist(f1,f2,fl,input_dir)

    igraph2mtx(G1,G1.vcount(),fname='graph-A.mtx')
    igraph2mtx(G2,G2.vcount(),fname='graph-B.mtx')
    igraph2mtx(L,G1.vcount(),fname='graph-L.mtx',bipartite=True)

    '''
    #### Initialize evaluation class
    if ft is not None:
        AQ=ev.AlignmentQuality(input_dir+f1, input_dir+f2, None, input_dir+ft, None, None)
    else:
        AQ=ev.AlignmentQuality(input_dir+f1, input_dir+f2, None, None, None, None)
    '''

    command=['/home/khan242/netAlign/netalignmc/netalign/netAlign', 'graph', '-osc.aln']
    alpha='-a '+str(alpha)
    beta='-b '+str(beta)
    nofinalize ='-f'
    verbose = '-v'
    command.append(alpha)
    command.append(beta)
    command.append(verbose)
    command.append(nofinalize)
    print(command)
    result = subprocess.run(command, stdout=subprocess.PIPE)
    print(result.stdout.decode('utf-8'))

    f=open('sc.aln','r')
    M=[]
    for l in f:
        l=l.strip('\n')
        l=l.split(' ')
        M.append((int(l[0])-1,int(l[1])-1))
    f.close()
   
    #AL=get_alignment(G1,G2,M)
    #s,q=evaluate(AQ,G1,G2,AL)
    
    command=['rm','graph-A.mtx','graph-B.mtx','graph-L.mtx']
    result = subprocess.run(command, stdout=subprocess.PIPE)  
    
    #return AL,s,q

if __name__ == '__main__':
    
    warnings.filterwarnings("ignore")
    
    input_dir="/home/khan242/netAlign/data/synthetic_networks/"
    res_dir="/home/khan242/netAlign/results/synthetic_networks/star/"
    
    
    f1="star.g"
    f2="email-Enron.txt"
    ft=None
    
    data=[]
    rows=[]

    for alpha in  [0]:
        beta=1-alpha
        rows=[]
        a=alpha
        print("\n******************************\n")
        fl="L_1_star_email-Enron.graphml"
        resf="star_email-Enron_"+str(a)+".aln"
        AL,s,q=netAlignSC(f1,f2,fl,alpha,beta,ft,input_dir=input_dir)
        print('Aligned')
        save_alignment(AL,res_dir+resf)
        print('saved')    
#         rows.append(round(s,3))
        
#         print("\n******************************\n")
#         fl="L_0.1_yeast0_yeast5.graphml"
#         resf="yeast0_yeast5_Y2H1_0.1_"+str(a)+".aln"
#         AL,s,q=netAlignSC(f1,f2,fl,alpha,beta,ft,input_dir=input_dir)
#         save_alignment(AL,res_dir+resf)
#         print('Score: ',s,q) 
#         rows.append(round(s,3))
        
#         print("\n******************************\n")
#         fl="L_0.25_yeast0_yeast5.graphml"
#         resf="yeast0_yeast5_Y2H1_0.25_"+str(a)+".aln"
#         AL,s,q=netAlignSC(f1,f2,fl,alpha,beta,ft,input_dir=input_dir)
#         save_alignment(AL,res_dir+resf)
#         print('Score: ',s,q) 
#         rows.append(round(s,3))
        
#         print("\n******************************\n")
#         fl="L_0.5_yeast0_yeast5.graphml"
#         resf="yeast0_yeast5_Y2H1_0.5_"+str(a)+".aln"
#         AL,s,q=netAlignSC(f1,f2,fl,alpha,beta,ft,input_dir=input_dir)
#         save_alignment(AL,res_dir+resf)
#         print('Score: ',s,q) 
#         rows.append(round(s,3))
        
#         print("\n******************************\n")
#         fl="L_1_yeast0_yeast5.graphml"
#         resf="yeast0_yeast5_Y2H1_1_"+str(a)+".aln"
#         AL,s,q=netAlignSC(f1,f2,fl,alpha,beta,ft,input_dir=input_dir)
#         save_alignment(AL,res_dir+resf)
#         print('Score: ',s,q) 
#         rows.append(round(s,3))
        
        data.append(rows)
#     rows=[-1,-1,-1,-1, 0.581]
    '''
    print("\n******************************\n")
    fl="L_1_yeast0_yeast5.graphml"
    resf="yeast0_yeast5_Y2H1_all.aln"
    AL,s,q=netAlignSC(f1,f2,fl,alpha=0,beta=1,ft=ft,input_dir=input_dir)
    save_alignment(AL,res_dir+resf)
    print('Score: ',s,q) 
    rows.append(s)
    '''
    
#     data.append(rows)
    
#     df = pd.DataFrame(np.array(data),columns=['K=100%'], index=['a=1', 'a=0.75', 'a=0.5', 'a=0.25','a=0'])
#     df.to_pickle("results_star.pkl")
#     #'''





