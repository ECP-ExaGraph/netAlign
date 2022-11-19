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
import copy
import importlib
import matplotlib.pylab as plt
import time
import pickle
import pandas as pd


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
            if M[i][0][1] > -1:
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
    
    wk=wk[perm]
    #for i in range(int(L.nnz/2),L.nnz):
        #wk[i]=wk[perm[i]]
    
    L.data=wk
    
    #r,c=L.nonzeros()
    #for i
    
    return bm.bSuitor(L,1)
    


def evaluate(AQ,G1,G2,Mbase,bestM,bestScore):
     
    bestk=-1
    score=-1
    #alll={}
    for k,M in Mbase.items():
        
        AL=get_alignment(G1,G2,M)
        AQ.mapping_set=AL
        #for (x,y) in AL:
            #if x == y:
                #print(x,y)
        #alll[k]=AL
        if AQ.true_mapping_set == None:
            qual=AQ.evaluate(False,False,False,False,False,True,False,False,False,False)
            score=qual["NCV-GS3"]
        else:
            qual=AQ.evaluate(False,False,True,False,False,True,False,False,False,False)
            score=float(qual["F-NC"])
        
        print('Score: ',score, qual)
        if score > bestScore:
            bestScore=score
            bestk=k
    
    #for i in range(4):
        #for j in range(i+1,4):
            #print(i,j,len(alll[i]),len(alll[j]),len(alll[i]-alll[j]))
    if bestk > -1:        ### new best
        bestM=copy.deepcopy(Mbase[k])
    
    return bestM,bestScore



def othermax(L,perm,y,z):
       
    w=np.concatenate((np.ravel(z),np.ravel(z)))
    w=w[perm]
    for i in range(int(L.nnz/2)):
        w[i]=y[i]
 
    L.data=w    ##### So the left is loaded with z and right is y
    
    #print(y)
    #print(z)
    #print(w)
    #print('LLLLL')
    #print(L)
    
    L=sp.sparse.csr_matrix(L)
    row,col=L.nonzero()
    (nrow,ncol)=L.shape
    
    indptr=L.indptr
    indices=L.indices
    
    maxr=0
    second=0
    for i in range(nrow):
        t=np.array(L.getrow(i).data)
        maxr=0
        second=0
        if len(t) <= 0:
            maxr=0
            second=0
        else:
            temp=list(-np.sort(-t))
            maxr=temp[0]
            second=maxr
            
            for r in range(len(temp)):
                if temp[r] < second:
                    second=temp[r]
                    break
            if maxr == second:
                second=0
            #print('t:',t,temp,maxr,second)
            if  maxr < 0:
                maxr=0
            if second < 0:
                second =0
        
        for j in range(indptr[i],indptr[i+1]):
            #print('data: ',L.data[j],maxr,second)
            if L.data[j] == maxr:
                L.data[j]=second
            else:
                L.data[j]=maxr
            #print('dataU: ',L.data[j],maxr,second)
    
    #print('zzzzz')
    #print(L)
    ll=int(L.nnz/2)
    omax_z=np.array(L.data[0:ll])
    omax_y=np.array(L.data[0:ll]) ## dummy
    
    i=0
    for j in range(ll,L.nnz):
        omax_y[perm[j]]=L.data[j]
    L=sp.sparse.coo_matrix(L)
    
    #print('omax')
    #print(omax_y)
    #print(omax_z)
    return omax_y,omax_z


def create_overlaps_matrix(G1,G2,L=None,saved=False,tok=''):
    
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
    
    if saved == False:    

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

        f=tok

        pickle.dump( S, open( "S_"+f+".p", "wb" ) )
        pickle.dump( L, open( "L_"+f+".p", "wb" ) )
        pickle.dump( w, open( "w_"+f+".p", "wb" ) )
        pickle.dump( permS, open( "permS_"+f+".p", "wb" ) )
        pickle.dump( perm, open( "perm_"+f+".p", "wb" ) )
    
    else:
        
        f=tok
        S = pickle.load( open( "S_"+f+".p", "rb" ) )
        L = pickle.load( open( "L_"+f+".p", "rb" ) )
        w = pickle.load( open( "w_"+f+".p", "rb" ) )
        permS = pickle.load( open( "permS_"+f+".p", "rb" ) )
        perm = pickle.load( open( "perm_"+f+".p", "rb" ) )
    
    return S,L,w,permS,perm



def netAlign(f1,f2,fl=None,true_mapping=None,resf=None,alpha=1,beta=1,gamma=0.99,maxiter=200,nbatch=40,input_dir=None,res_dir=None,saved=False,tok=''):
    
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
    prev_score=-1
    Mbase={}
    
    #### Reading networks
    t1=time.perf_counter()
    #G1,G2,L=read_graphml(f1,f2,fl,input_dir)
    G1,G2,L=read_networks(f1,f2,fl,input_dir)
    t2=time.perf_counter()
    
    print('Network Reading Done: ',round((t2-t1),2))
    S,L,w,permS,perm=create_overlaps_matrix(G1,G2,L,saved,tok)
    
    t3=time.perf_counter()
    print('Computing Overlap Matrix (S.shape, L.shape, L.nnz) :', S.shape, S.nnz, L.shape, L.nnz, ', Time: ',round((t3-t2),2))
    
    ts=time.perf_counter()
    
    #### Initialize evaluation class
    if true_mapping is not None:
        AQ=ev.AlignmentQuality(input_dir+f1, input_dir+f2, None, input_dir+true_mapping, None, None)
    else:
        AQ=ev.AlignmentQuality(input_dir+f1, input_dir+f2, None, None, None, None)

     
    Lsize=int(L.nnz/2)
    indptr=S.indptr
    indices=S.indices
    Sval=np.full(S.nnz, 1)
    Sk=np.full(S.nnz, 1)
    
    l3=time.perf_counter()
    print('Initialization Done ', round((l3-ts),2))
      
    ###### First iteration

    F=np.full((len(Sval)), beta)           ### Line 3
    d=(alpha*w)+beta*np.ravel(S.sum(1))    ### Line 4
        
    l4=time.perf_counter()
    print('Line 3/4: ', round((l4-l3),2))
    
    yk=np.array(d)                     ### Line 5
    zk=np.array(d)                     ### Line 6
    
    l5=time.perf_counter()
    print('Line othermax: ', round((l5-l4),2))
    
    Sk=(sp.sparse.diags(yk+zk-d)*S).data-F  ### Line 7
    
    yk1=gamma*yk                            ### Line 8
    zk1=gamma*zk
    Sk1=gamma*Sk
    gamma=gamma*gamma
    l6=time.perf_counter()
    print('Line updates: ', round((l6-l5),2))
    
    Mbase[0]=round_heuristic(L,perm,yk) ### Line 10
    Mbase[1]=round_heuristic(L,perm,zk) ### Line 11
    
    te=time.perf_counter()
    print('Iteration 1 : ',round((te-ts),2),round((l4-l3),2),round((l5-l4),2),round((l6-l5),2),round((te-l6),2))
    
    batch_id=2
    for t in range(1000):
        
        if t >= maxiter:
            break
        ts=time.perf_counter()
        l3=time.perf_counter()
        
        F=beta*Sval+Sk1[permS]             ### Line 3
        np.clip(F,0,beta,out=F)
        
        S.data=F
        d=(alpha*w)+np.ravel(S.sum(1))    ### Line 4
        S.data=Sval
          
        l4=time.perf_counter()
        
        omax_y,omax_z=othermax(L,perm,yk1,zk1)
     
        yk=d - omax_y                         ### Line 5
        zk=d - omax_z                         ### Line 6
        
        l5=time.perf_counter()
        
        diag=(yk+zk-d)+0.0001 #### Avoiding zero scaling
        Sk=(sp.sparse.diags(diag)*S).data-F  ### Line 7
       
        
        yk1=gamma*yk + (1-gamma)*yk1            ### Line 8
        zk1=gamma*zk + (1-gamma)*zk1
        Sk1=gamma*Sk + (1-gamma)*Sk1
        gamma=gamma*gamma
        
        l6=time.perf_counter()

        Mbase[batch_id]=round_heuristic(L,perm,yk)   ### Line 10
        Mbase[batch_id+1]=round_heuristic(L,perm,zk) ### Line 11
        
        te=time.perf_counter()
        print('Iteration',(t+2),': ',round((te-ts),2),round((l4-l3),2),round((l5-l4),2),round((l6-l5),2),round((te-l6),2))
        
        batch_id=batch_id+2
        breaker=False
        if batch_id == nbatch:
            batch_id=0
            prev_score=bestScore
            bestM,bestScore=evaluate(AQ,G1,G2,Mbase,bestM,bestScore)
            print('Best score: ',bestScore)
            if bestScore == prev_score:
                breaker = True
        
        if breaker and t >= 100:
            break
    
    AL=get_alignment(G1,G2,bestM)
    save_alignment(AL,res_dir+resf)
    return bestScore


if __name__ == '__main__':
    
    warnings.filterwarnings("ignore")
    
    input_dir="/home/khan242/netAlign/data/synthetic_networks/"
    res_dir="/home/khan242/netAlign/results/synthetic_networks/"
    
    
    f1="yeast0_Y2H1.gw"
    f2="yeast5_Y2H1.gw"
    ft="true_node_mapping.txt"
    
    data=[]
    rows=[]
    #'''
    print("\n******************************\n")
    fl="L_0.05_yeast0_yeast5.graphml"
    resf="yeast0_yeast5_Y2H1_0.05_1.aln"
    s=netAlign(f1,f2,fl,ft,resf,alpha=1,beta=0,res_dir=res_dir,tok='0.05')
    rows.append(s)
    
    print("\n******************************\n")
    fl="L_0.1_yeast0_yeast5.graphml"
    resf="yeast0_yeast5_Y2H1_0.1_1.aln"
    s=netAlign(f1,f2,fl,ft,resf,alpha=1,beta=0,res_dir=res_dir,tok='0.1')
    rows.append(s)
    
    print("\n******************************\n")
    fl="L_0.25_yeast0_yeast5.graphml"
    resf="yeast0_yeast5_Y2H1_0.25_1.aln"
    s=netAlign(f1,f2,fl,ft,resf,alpha=1,beta=0,res_dir=res_dir,tok='0.25')
    rows.append(s)
    
    print("\n******************************\n")
    fl="L_0.5_yeast0_yeast5.graphml"
    resf="yeast0_yeast5_Y2H1_0.5_1.aln"
    s=netAlign(f1,f2,fl,ft,resf,alpha=1,beta=0,res_dir=res_dir,tok='0.5')
    rows.append(s)
    
    print("\n******************************\n")
    fl="L_1_yeast0_yeast5.graphml"
    resf="yeast0_yeast5_Y2H1_1_1.aln"
    s=netAlign(f1,f2,fl,ft,resf,alpha=1,beta=0,res_dir=res_dir,tok='1')
    rows.append(s)
    
    data.append(rows)
    
    for a in  [0.9, 0.75, 0.50, 0.25, 0.1]:
        b=1-a
        rows=[]
        
        print("\n******************************\n")
        fl="L_0.05_yeast0_yeast5.graphml"
        resf="yeast0_yeast5_Y2H1_0.05_"+str(a)+".aln"
        s=netAlign(f1,f2,fl,ft,resf,alpha=a,beta=b,res_dir=res_dir,saved=True,tok='0.05')
        rows.append(s)
        
        print("\n******************************\n")
        fl="L_0.1_yeast0_yeast5.graphml"
        resf="yeast0_yeast5_Y2H1_0.1_"+str(a)+".aln"
        s=netAlign(f1,f2,fl,ft,resf,alpha=a,beta=b,res_dir=res_dir,saved=True,tok='0.1')
        rows.append(s)
        
        print("\n******************************\n")
        fl="L_0.25_yeast0_yeast5.graphml"
        resf="yeast0_yeast5_Y2H1_0.25_"+str(a)+".aln"
        s=netAlign(f1,f2,fl,ft,resf,alpha=a,beta=b,res_dir=res_dir,saved=True,tok='0.25')
        rows.append(s)
        
        print("\n******************************\n")
        fl="L_0.5_yeast0_yeast5.graphml"
        resf="yeast0_yeast5_Y2H1_0.5_"+str(a)+".aln"
        s=netAlign(f1,f2,fl,ft,resf,alpha=a,beta=b,res_dir=res_dir,saved=True,tok='0.5')
        rows.append(s)
        
        print("\n******************************\n")
        fl="L_1_yeast0_yeast5.graphml"
        resf="yeast0_yeast5_Y2H1_1_"+str(a)+".aln"
        s=netAlign(f1,f2,fl,ft,resf,alpha=a,beta=b,res_dir=res_dir,saved=True,tok='1')
        rows.append(s)
        
        data.append(rows)
    
    rows=[-1,-1,-1,-1,0.581]
    '''
    print("\n******************************\n")
    fl="L_1_yeast0_yeast5.graphml"
    resf="yeast0_yeast5_Y2H1_all.aln"
    s=netAlign(f1,f2,None,ft,resf,alpha=1,beta=1,res_dir=res_dir,saved=True,tok='all',maxiter=3,nbatch=8)
    rows.append(s)
    '''
    
    data.append(rows)
    
    df = pd.DataFrame(np.array(data),columns=['K=0.05', 'K=0.1', 'K=0.25', 'K=0.50', 'K=1'], index=['a=1', 'a=0.9', 'a=0.75', 'a=0.5', 'a=0.25', 'a=0.1', 'a=0'])
    df.to_pickle("results_danai.pkl")
    #'''





