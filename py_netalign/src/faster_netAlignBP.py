#!/usr/bin/env python
# coding: utf-8

# In[11]:


import warnings
import igraph as ig
import networkx as nx
import collections
import scipy as sp
import numpy as np
import bMatching as bm
import  Evaluation as ev
import copy
import importlib
import matplotlib.pylab as plt
import time
import pickle


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
    for i in range(len(M)):
        if i < n1:
            j=M[i][0][1]-n1
            
            u=G1.vs[i]['id']
            v=G2.vs[j]['id']
            
            Align.add((u,v))
        else:
            break
    
    return Align 



def save_alignment(L,filename):
    
    fl=open(filename,"w")
    for (u,v) in L:
        t=u+" "+v+"\n"
        fl.write(t)
    
    fl.close()



def round_heuristic(L,perm,w):
    
    wk=np.concatenate((np.ravel(w),np.ravel(w)))
    
    for i in range(int(L.nnz/2),L.nnz):
        wk[i]=wk[perm[i]]
    
    L.data=wk
    
    return bm.bSuitor(L,1)
    


def evaluate(AQ,Mbase,bestM,bestScore):
     
    bestk=-1
    score=-1
    for k,M in Mbase.items():
        
        AL=get_alignment(G1,G2,M)
        AQ.mapping_set=AL
        
        if AQ.true_mapping_set == None:
            qual=AQ.evaluate(False,False,False,False,False,True,False,False,False,False)
            score=qual["NCV-GS3"]
        else:
            qual=AQ.evaluate(False,False,True,False,False,True,False,False,False,False)
            score=qual["F-NC"]
        
        if score > bestScore:
            bestScore=score
            bestk=k
    
    if bestk > -1:        ### new best
        bestM=copy.deepcopy(Mbase[k])
    
    return bestM,bestScore



def othermax(L,perm,y,z):
       
    w=np.concatenate((np.ravel(y),np.ravel(y)))
    
    for i in range(L.nnz):
        if i < int(L.nnz/2):
            w[i]=z[i]
        else:
            w[i]=w[perm[i]]
    
 
    L.data=w    ##### So the left is loaded with z and right is y
    
    L=sp.sparse.csr_matrix(L)
    row,col=L.nonzero()
    (nrow,ncol)=L.shape
    
    indptr=L.indptr
    indices=L.indices
    
    maxr=0
    second=0
    for i in range(nrow):
        t=np.array(L.getrow(i).data)
        temp=list(np.sort(t[(-t).argsort()[:2]]))
        #print('t:',t,temp)
        if len(temp) == 2:
            second =temp[0]
            maxr=temp[1]
        else:
            second=0
            maxr=temp[0]
        
        if second < 0:
            second=0
        if maxr < 0:
            maxr=0
        
        for j in range(indptr[i],indptr[i+1]):
            #print('data: ',L.data[j],maxr,second)
            if L.data[j] == maxr:
                L.data[j]=second
            else:
                L.data[j]=maxr
            #print('dataU: ',L.data[j],maxr,second)
    
    ll=int(L.nnz/2)
    omax_y=np.array(L.data[0:ll])
    omax_z=np.array(L.data[ll:L.nnz])
    L=sp.sparse.coo_matrix(L)
    
    return omax_y,omax_z


def create_overlaps_matrix(G1,G2,L=None):
    
    ### G1: network 1
    ### G2: network 2
    ### L : Bipartite graph where 
    ### left side is G1 vertices and 
    ### right side is G2 vertices
    
    n1=G1.vcount()
    n2=G2.vcount()
    
    ### for iterating in L 
    G1_vid=list(range(n1))
    G2_vid=list(range(n1,n1+n2))
    
    ###
    nS=0
    if L is None:  
        ### L is not avaiable so assume complete bipartite
        ### and assume edges are indexed accordingly
        nS=n1*n2 
        
    else:
        nS=L.ecount()
    
    #S=ig.Graph()  
    #S.add_vertices(list(range(nS)))
    edgeS=[]
    
    for u in G1_vid:
        
        nbor1=G1.neighbors(u)
        
        for v in G2_vid:
            uv=-1
            if L is not None:
                uv=L.get_eid(u,v,directed=False,error=False)
                
                #if  uv == -1:
                    #print("Sparse L!!")
                    #continue
            else:
                uv= u*n2+(v-n1) #### Be careful
            
            if uv== -1:
                continue ## Just sanity check
            
            nbor2=G2.neighbors(v-n1)
            

            for i in nbor1:
                for jj in nbor2:
                    j=jj+n1 ### shifting the vertex id for bipartite graph

                    ij=-1
                    if L is not None:
                        ### Now check whether the neighbors has cross edge
                        ij=L.get_eid(i,j,directed=False,error=False)

                        #if ij == -1:
                            #print("L Sparse 2 !!")
                            #continue
                    else:
                        ij=(i*n2)+jj
                    
                    if ij==-1:
                        continue ### Just sanity check

                    
                    edgeS.append((uv,ij))

    #S.add_edges(edgeS)
    
    #edges = self.get_edgelist() 
    weights = [1] * len(edgeS)
    S = sp.sparse.csr_matrix((weights, zip(*edgeS)), shape=(nS, nS))
    S = S + sp.sparse.triu(S, 1).T + sp.sparse.tril(S, -1).T
    S=S/2
    
    ###### Permutation S is different from L permutation
    ###### Perm S is for transposing and perm L is for concatenating
    indptr=S.indptr
    indices=S.indices
    permS=np.array(S.indices)
    newptr=np.array(indptr)
    size=S.shape
    for i in range(int(size[0])):
        for j in range(indptr[i],indptr[i+1]):
            c=indices[j]
            permS[newptr[c]]=j
            newptr[c]+=1
    
    if L is None:
        row=[]
        col=[]
        rowS=[]
        colS=[]
        data=[1]*(n1*n2*2) #### Factor 2 is for symmetry
        for i in range(n1):
            for j in range(n2):
                row.append(i)
                col.append(j+n1)
                rowS.append(j+n1)
                colS.append(i)
                
        row = row + rowS
        col = col + colS
        row  = np.array(row)
        col  = np.array(col)
        data = np.array(data)
        L = sp.sparse.csr_matrix((data, (row, col)), shape=(n1+n2, n1+n2))
        
        ###### Permutation
        ll=int(L.nnz/2)
        perm=list(range(L.nnz))
        for i in range(ll,L.nnz):
            perm[i]=perm[i]-ll
        
        
    else:
        edgeL=L.get_edgelist()
        nL=L.vcount()
        weights=L.es['weight']
        L = sp.sparse.csr_matrix((weights, zip(*edgeL)), shape=(nL, nL))
        L = L + sp.sparse.triu(L, 1).T + sp.sparse.tril(L, -1).T
        
        ###### Permutation
        ll=int(L.nnz/2)
        perm=list(range(L.nnz))
        row,col=L.nonzero()
        indptr=L.indptr
        indices=L.indices
        for i in range(ll,L.nnz):
            c=row[i]
            r=col[i]           
            for j in range(indptr[r],indptr[r+1]):
                if c == indices[j]:
                    perm[i]=j
                    break
                else:
                    perm[i]=-1          
    
    L=sp.sparse.coo_matrix(L)
    w=np.array(L.data[0:ll])
    perm=np.array(perm)
    permS=np.array(permS)
    
    pickle.dump( S, open( "S.p", "wb" ) )
    pickle.dump( L, open( "L.p", "wb" ) )
    pickle.dump( w, open( "w.p", "wb" ) )
    pickle.dump( permS, open( "permS.p", "wb" ) )
    pickle.dump( perm, open( "perm.p", "wb" ) )
    
    return S,L,w,permS,perm



def netAlign(f1,f2,fl=None,true_mapping=None,resf=None,alpha=1,beta=2,gamma=0.99,maxiter=100,nbatch=20,input_dir=None,res_dir=None):
    
    if input_dir is None:
        input_dir="/home/khan242/netAlign/data/synthetic_networks/"
    if res_dir is None:
        res_dir="/home/khan242/netAlign/results/"  
     
    if resf is None:
        t=f1.split(".")
        resf=t[0]
    
        t=f2.split(".")
        resf=resf+"_"+t[0]+".aln"
    
        resp=resf.split(".")
    
    #### Output
    bestM=[]
    bestScore=-1
    Mbase={}
    
    #### Reading networks
    t1=time.perf_counter()
    #G1,G2,L=read_graphml(f1,f2,fl,input_dir)
    G1,G2,L=read_networks(f1,f2,fl,input_dir)
    t2=time.perf_counter()
    
    print('Network Reading Done: ',round((t2-t1),2))
    S,L,w,permS,perm=create_overlaps_matrix(G1,G2,L)
    t3=time.perf_counter()
    print('Computing Overlap Matrix (S.shape, L.shape, L.nnz) :', S.shape, S.nnz, L.shape, L.nnz, ', Time: ',round((t3-t2),2))
    
    ts=time.perf_counter()
    
    #### Initialize evaluation class
    if true_mapping:
        AQ=ev.AlignmentQuality(input_dir+f1, input_dir+f2, None, input_dir+true_mapping, None, None)
    else:
        AQ=ev.AlignmentQuality(input_dir+f1, input_dir+f2, None, None, None, None)

    
    Lsize=int(L.nnz/2)
    indptr=S.indptr
    indices=S.indices
    Sval=np.array(S.data)
    Sk=np.array(S.data)
    
    l3=time.perf_counter()
    print('Initialization Done ', round((l3-ts),2))
      
    ###### First iteration

    F=np.full((len(Sval)), beta)       ### Line 3
    d=(alpha*w)                        ### Line 4
    for i in range(Lsize):             ### Line 4
        tot=0
        for j in range(indptr[i],indptr[i+1]):
            tot=tot+F[j]
        d[i]=tot
        
    l4=time.perf_counter()
    print('Line 3/4: ', round((l4-l3),2))
    
    yk=np.array(d)                     ### Line 5
    zk=np.array(d)                     ### Line 6
    
    l5=time.perf_counter()
    print('Line othermax: ', round((l5-l4),2))
    
    tval=yk+zk-d                       ### Line 7
    for i in range(Lsize):
        for j in range(indptr[i],indptr[i+1]):
            Sk[j]=tval[i]-F[j] 
    
    
    yk1=gamma*yk                       ### Line 8
    zk1=gamma*zk
    Sk1=gamma*Sk
    
    l6=time.perf_counter()
    print('Line updates: ', round((l6-l5),2))
    
    Mbase[0]=round_heuristic(L,perm,yk) ### Line 10
    Mbase[1]=round_heuristic(L,perm,zk) ### Line 11
    
    te=time.perf_counter()
    print('Iteration 1 : ',round((te-ts),2),round((l4-l3),2),round((l5-l4),2),round((l6-l5),2),round((te-l6),2))
    
    batch_id=2
    for t in range(500):
        
        if t >= maxiter:
            break
        ts=time.perf_counter()
        l3=time.perf_counter()
        
        F=beta*Sval+Sk[permS]             ### Line 3
        np.clip(F,0,beta,out=F)
        
        for i in range(Lsize):            ### Line 4
            tot=0
            for j in range(indptr[i],indptr[i+1]):
                tot=tot+F[j]+alpha*w[i]
            d[i]=tot
          
        l4=time.perf_counter()
        
        omax_y,omax_z=othermax(L,perm,yk1,zk1)
        yk=d - omax_y                         ### Line 5
        zk=d - omax_z                         ### Line 6
        
        l5=time.perf_counter()
        
        tval=yk+zk-d
        for i in range(Lsize):                ### Line 7 
            tv=tval[i]
            for j in range(indptr[i],indptr[i+1]):
                Sk[j]=tv-F[j]
        
        yk1=gamma*yk + (1-gamma)*yk1          ### Line 8
        zk1=gamma*zk + (1-gamma)*zk1
        Sk1=gamma*Sk + (1-gamma)*Sk1
        
        l6=time.perf_counter()

        Mbase[batch_id]=round_heuristic(L,perm,yk)   ### Line 10
        Mbase[batch_id+1]=round_heuristic(L,perm,zk) ### Line 11
        
        te=time.perf_counter()
        print('Iteration',(t+2),': ',round((te-ts),2),round((l4-l3),2),round((l5-l4),2),round((l6-l5),2),round((te-l6),2))
        
        batch_id=batch_id+2
        if batch_id == nbatch-1:
            batch_id=0
            bestM,bestScore=evaluate(AQ,Mbase,bestM,bestScore)
            print('Batch score: ',bestScore)
            
    
    AL=get_alignment(G1,G2,bestM)
    save_alignment(AL,res_dir+resf)
    return AL


if __name__ == '__main__':
    
    warnings.filterwarnings("ignore")
    
    input_dir="/home/khan242/netAlign/data/synthetic_networks/"
    res_dir="/home/khan242/netAlign/results/"
    
    
    f1="yeast0_Y2H1.gw"
    f2="yeast5_Y2H1.gw"

    resf="yeast0_yeast5_Y2H1.aln"
    netAlign(f1,f2,None,None,resf,maxiter=1,nbatch=4,res_dir=res_dir)
    
    '''
    f1="yeast0_Y2H1.gw"
    f2="yeast15_Y2H1.gw"
    resf="yeast0_yeast15_Y2H1.aln"
    
    netAlign(f1,f2,resf,None,"ALL",input_dir,res_dir)
    
    f1="yeast0_Y2H1.gw"
    f2="yeast20_Y2H1.gw"
    resf="yeast0_yeast20_Y2H1.aln"
    
    netAlign(f1,f2,resf,None,"ALL",input_dir,res_dir)
    
    f1="yeast0_Y2H1.gw"
    f2="yeast25_Y2H1.gw"
    resf="yeast0_yeast25_Y2H1.aln"
    
    netAlign(f1,f2,resf,None,"ALL",input_dir,res_dir)
    '''







