#import sys
import os
import scipy.io as spio
import scipy as sp
import numpy as np
#import networkx as nx
#import community
#import matplotlib.pylab as pl
from functools import cmp_to_key
import heapq
import queue
import time


#def edge_induce(G,M):
    #return A

def comparator(i,j):
    
    #print(i,j)
    
    if G.data[i] > G.data[j]:
        return -1
    if G.data[i] == G.data[j]:
        if G.indices[i] > G.indices[j]:
            return -1
    return 1

def neighborlist_sort(G):
    
    permutation=[]
    for i in range(len(G.indptr)-1):
        start=G.indptr[i]
        end=G.indptr[i+1]
        permutation.extend(sorted(list(range(start,end)),key=cmp_to_key(comparator)))
    
    permutation=np.array(permutation)
    G.indices=G.indices[permutation]
    G.data=G.data[permutation]
    

def verify_matching(M,b):
    
    
    for i in range(len(M)):
        for j in range(b):
            if M[i][j][1]>=0:
                k=M[i][j][1]
                flag=-1
                for (x,y) in M[k]:
                    if y == i:
                        flag=1
                        #print(i,M[i][j])
                        break
                if flag == -1 :
                    print(i,k,M[k])
                    return False                            
    return True

def bSuitorCSR(G,b):
    
    #m=G.nnz
    n=len(G.indptr)-1
    
    Q=queue.Queue(maxsize=n)
    Q1=queue.Queue(maxsize=n)
    tQ=Q
    tQ1=Q1
    
    bval1=[0]*n
    
    cur=[0]*n
    alive=[1]*n
    
    #if type(b)=="int":
        #print("Here")
    bval=[b]*n
    #else:
        #bval=b
    
    
    M=[]
    for i in range(n):
        tlist=[]
        for j in range(bval[i]):
            tlist.append((0.0,-1))
        
        M.append(tlist)
        Q.put(i)
        cur[i]=G.indptr[i]
    
    print("Memory Allocation Done")
    
    neighborlist_sort(G)
    
    print("Neighborlist Sorting Done")   
    
    while(True):
        print("Iteration: ",tQ.qsize())
        while(True):
            if tQ.empty():
                break
            i=tQ.get(False)
            for k in range(bval[i]):
                for p in range(cur[i],G.indptr[i+1]):
                    j=G.indices[p]
                    v=G.data[p]
                    if v<=0.0 or i==j:
                        continue
                    
                    partner=(v,i)
                    if partner > M[j][0]:
                        kicked=heapq.heapreplace(M[j],partner)
                        
                        if kicked[0] > 0.0:
                            looser=kicked[1]
                            if bval1[looser] == 0 and alive[looser]==1:
                                tQ1.put(looser,False)
                            
                            bval1[looser]=bval1[looser]+1
                        
                        if p == G.indptr[i+1]-1:
                            alive[i]=0
                        
                        cur[i]=p+1
                        
                        break
                    
                
                if p>=G.indptr[i+1]:
                    cur[i]=p
                    alive[i]=0
        
        if tQ1.empty():
            break
        tempQ=tQ1
        tQ1=tQ
        tQ=tempQ
        for i in range(n):
            bval[i]=bval1[i]
            bval1[i]=0
    
    
    return M


 
def bSuitor(S,b,verbose=False):
    
    (nVer,t)=S.shape
    G=[]
    
    #if verbose:
        #S1=sp.sparse.coo_matrix(S)
        #H=S1.getH()
        #print(type(H),type(S1))
        #print(H-S1)
    
    for i in range(nVer):
        r=S.getrow(i)
        t=list(zip(r.data,r.indices))
        G.append(t)
    
    if verbose:
        print("# of vertices: ",nVer)
        print("# of edges: ",S.nnz)
    #m=G.nnz
    n=len(G) 
   
    Q=queue.Queue(maxsize=n)
    Q1=queue.Queue(maxsize=n)
    tQ=Q
    tQ1=Q1
    
    bval1=[0]*n
    
    cur=[0]*n
    alive=[1]*n
    
    #if type(b)=="int":
        #print("Here")
    bval=[b]*n
    #else:
        #bval=b
    
    
    M=[]
    for i in range(n):
        tlist=[]
        for j in range(bval[i]):
            tlist.append((0.0,-1))
        
        M.append(tlist)
        Q.put(i)
        cur[i]=0
    
    if verbose:
        print("Memory Allocation Done")
    
    for i in range(n):
        G[i].sort(reverse=True)
    
    
    if verbose:
        print("Neighborlist Sorting Done")   
    
    while(True):
        #if verbose:
            #print("Iteration: ",tQ.qsize())
        while(True):
            if tQ.empty():
                break
            i=tQ.get(False)
            for k in range(bval[i]):
                for p in range(cur[i],len(G[i])): #G
                    suitor=G[i][p]
                    j=suitor[1]  #G
                    v=suitor[0]     #G
                    if v<=0.0 or i==j:
                        continue
                    
                    partner=(v,i)
                    if partner > M[j][0]:
                        kicked=heapq.heapreplace(M[j],partner)
                        
                        if kicked[0] > 0.0:
                            looser=kicked[1]
                            if bval1[looser] == 0 and alive[looser]==1:
                                tQ1.put(looser,False)
                            
                            bval1[looser]=bval1[looser]+1
                        
                        if p >= len(G[i])-1:  #G
                            alive[i]=0
                        
                        cur[i]=p+1
                        
                        break
                    
                
                if p>=len(G[i]): #G
                    cur[i]=p
                    alive[i]=0
        
        if tQ1.empty():
            break
        tempQ=tQ1
        tQ1=tQ
        tQ=tempQ
        for i in range(n):
            bval[i]=bval1[i]
            bval1[i]=0
    
    
    if verbose:
        if verify_matching(M,b):
            print("Matching is Verified !!")
        else:
            print("Matching Bad..!")

        card=0
        weight=0.0
        for i in range(nVer):
            for j in range(b):
                weight=weight+M[i][j][0]
                if M[i][j][1]>-1:
                    card+=1
        print("Matching Weight: ", weight/2.0, "Cardinality: ",card/2)
            
    return M


###################### MAIN #########################
if __name__ == "__main__":
    
    dataFolder="../data/"
    resFolder="../results/"

    for fileName in os.listdir(dataFolder):
        if fileName.endswith(".btx"):
            print("Processing: "+fileName)
            fileName="cond-mat.btx"
            b=1

            outFile=fileName
            outFile=outFile.split(".")
            outFile=outFile[0]
            problemName=outFile

            fileName=dataFolder+fileName

            header=spio.mminfo(fileName)
            nVer=header[0]
            nEdge=header[2]
            print("Number of vertices: ",nVer)
            print("Number of edges: ",nEdge)


            S=spio.mmread(fileName)
           

    #        G=S.tocsr()
    #        del S
    #        M=bSuitorCSR(G,b)
            
            t0=time.time()
            M=bSuitor(S,b,True)
            t1=time.time()
            for i in range(10):
                print(M[i])
            print("Total Time: ", t1-t0)
        
        
