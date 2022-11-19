'''
Created on May 2, 2015

@author: Lei Meng
'''

import sys
import math
import scipy
import scipy.stats

class AlignmentQuality:
    
    mapping_set = None
    true_mapping_set = None
    goterm_dict1 = None 
    goterm_dict2 = None
    node_set1 = None
    node_set2 = None
    edge_set1 = None
    edge_set2 = None
    
    qlty_NC_P = None
    qlty_NC_R = None
    qlty_NC_F = None
    qlty_GS3 = None
    qlty_NCV = None
    qlty_NCV_GS3 = None
    qlty_GC = None
    qlty_PF_P = None
    qlty_PF_R = None
    qlty_PF_F = None

    def __init__(self, graph1, graph2, mapping_name, true_mapping_name, goterm_name1, goterm_name2):
        """Set the file names used to evalute the alignment quality"""
        if mapping_name != "" and mapping_name != None:
            self.mapping_set = self.read_file(mapping_name)

        if true_mapping_name != "" and true_mapping_name != None:
            self.true_mapping_set = self.read_file(true_mapping_name)
            
        if goterm_name1 != "" and goterm_name1 != None:
            self.goterm_dict1 = self.read_go_term(goterm_name1)
            
        if goterm_name2 != "" and goterm_name2 != None:
            self.goterm_dict2 = self.read_go_term(goterm_name2) 
            
        if graph1 != "" and graph1 != None:
            (self.node_set1, self.edge_set1) = self.read_graph(graph1)

        if graph2 != "" and graph2 != None:
            (self.node_set2, self.edge_set2) = self.read_graph(graph2)
        
    def evaluate(self, eval_NC_P, eval_NC_R, eval_NC_F, eval_GS3, eval_NCV, eval_NCV_GS3, eval_GC, eval_PF_P, eval_PF_R, eval_PF_F):      
        
        res = dict()
        if eval_NC_P == True or eval_NC_R == True or eval_NC_F == True:
            self.NC()
            if eval_NC_P == True:
                res["P-NC"] = self.qlty_NC_P
            if eval_NC_R == True:
                res["R-NC"] = self.qlty_NC_R
            if eval_NC_F == True:
                res["F-NC"] = self.qlty_NC_F
                
        if eval_GS3 == True:
            self.GS3()
            res["GS3"] = self.qlty_GS3
            
        if eval_NCV == True:
            self.NCV()
            res["NCV"] = self.qlty_NCV
            
        if eval_NCV_GS3 == True:
            self.NCV_GS3()
            res["NCV-GS3"] = self.qlty_NCV_GS3
            
        if eval_GC == True:
            self.GC()
            res["GC"] = self.qlty_GC
        
        if eval_PF_P == True or eval_PF_R == True or eval_PF_F == True:
            self.PF()
            if eval_PF_P == True:
                res["P-PF"] = self.qlty_PF_P
            if eval_PF_R == True:
                res["R-PF"] = self.qlty_PF_R
            if eval_PF_F == True:
                res["F-PF"] = self.qlty_PF_F
        return res
                
    def read_graph(self, name):
        gw_file = open(name, 'r')
        gw_file.readline()
        gw_file.readline()
        gw_file.readline()
        gw_file.readline()

        n_nodes1 = int(gw_file.readline().strip())
        nodes1 = []
        for i in range(n_nodes1):
                node = gw_file.readline().strip().strip("{}|")
                nodes1.append(node)

        edge_set1 = set()
        n_edges1 = int(gw_file.readline().strip())
        for i in range(n_edges1):
                tokens = tuple(gw_file.readline().strip().split())
                v1 = nodes1[int(tokens[0])-1]
                v2 = nodes1[int(tokens[1])-1]
                edge_set1.add((v1, v2))

        return (set(nodes1), edge_set1)

    def read_file(self, name):
        pair_set = set([])
        #print('Reading file: ',name)
        fi = open(name, 'r')
        for row in fi:
            (first, second) = row.rstrip().split()
            pair_set.add((first, second))
        return pair_set
                
                
    def read_go_term(self, name):
        go_term_dict = {}
        fi = open(name, 'rb')
        for row in fi:
            protein, goterm = row.rstrip().split('\t')
    
            if protein in go_term_dict:
                go_term_dict[protein].add(goterm)
            else:
                go_term_dict[protein] = set([goterm])
        return go_term_dict


    def NC(self):
        """Node correctness (require true node mapping)""" 
        if self.true_mapping_set == None:
            print("Please provide true node mapping")
            return

        (self.qlty_NC_P, self.qlty_NC_R, self.qlty_NC_F) = self.Accuray(self.mapping_set, self.true_mapping_set)
    
    def Accuray(self, pred_set, true_set):
        """Compute precision, recall, and F-score of the prediction"""
        len_overlap = len(pred_set & true_set)
        len_pred_set = len(pred_set)
        len_true_set = len(true_set)
        
        #universe_p=set([u for (u,v) in pred_set])#+[v for (u,v) in pred_set])
        #universe_t=set([u for (u,v) in true_set])#+[v for (u,v) in true_set])
        #print('Check: ',universe_t-universe_p)              
        
        if len_pred_set == 0:
            precision = '0'
        else:
            precision = len_overlap * 1.0/len_pred_set
            
        if len_true_set == 0:
            recall = '0'
        else:
            recall = len_overlap * 1.0/len_true_set
            
        if precision == '0' or recall == '0' or len_overlap == 0:
            F1 = '0'
        else:
            F1 = precision*recall*2.0/(precision+recall)
        
        return (precision, recall, F1)
        
    def GS3_edge_conservation(self, alignment, edges1, edges2):
        """Compute the number of edges from graph2 aligned with each edge from graph1"""
        conserved_set={}
        for (n1, n2) in edges1:
            total = 0
            conserved = 0
            if n1 in alignment and n2 in alignment:
                for u in alignment[n1]:
                    for v in alignment[n2]:
                        if u != v:
                            total = total + 1
    
                        if (u, v) in edges2 or (v, u) in edges2:
                            conserved = conserved + 1
            if total > 0:
                conserved_set[(n1, n2)] = (conserved, total)
        return (conserved_set)


    def GS3(self):
        """NCV, GS3, NCV-GS3"""
        if self.mapping_set == None or self.edge_set1 == None or self.edge_set2 == None:
            print("Cannot compute GS3 score. You need to provide two graphs and a mapping file.")
            return
        
        alignment1={}
        alignment2={}
        for (n1, n2) in self.mapping_set:
            if n1 in alignment1:
                alignment1[n1].add(n2)
            else:
                alignment1[n1] = set([n2])

            if n2 in alignment2:
                alignment2[n2].add(n1)
            else:
                alignment2[n2] = set([n1])


        #number of possible edges between f(n1) and f(n2)
        conserved_set1=self.GS3_edge_conservation(alignment1, self.edge_set1, self.edge_set2)
        conserved_set2=self.GS3_edge_conservation(alignment2, self.edge_set2, self.edge_set1)
    
        overlap = 0#number of conserved edge pairs
        conserved_graph1 = 0
        conserved_graph2 = 0
        for (n1, n2) in conserved_set1:
            (conserved, total) = conserved_set1[(n1, n2)]
            overlap = overlap + conserved
            conserved_graph1 = conserved_graph1 + total
    
        for (n1, n2) in conserved_set2:
            (conserved, total) = conserved_set2[(n1, n2)]
            conserved_graph2 = conserved_graph2 + total
    
        if conserved_graph1+conserved_graph2-overlap == 0:
            GS3 = '0'
        else:
            GS3 = 1.0 * overlap / (conserved_graph1+conserved_graph2-overlap)
        #print('GS3: ',GS3)
        self.qlty_GS3 = GS3

    def NCV(self):
        """Compute node coverage"""
        node_set1 = set([])
        node_set2 = set([])
        for (n1, n2) in self.mapping_set:
            node_set1.add(n1)
            node_set2.add(n2)
            
        if (len(self.node_set1) + len(self.node_set2)) == 0:
            coverage = '0'
        else:
            coverage = 1.0 * (len(node_set1) + len(node_set2)) / (len(self.node_set1) + len(self.node_set2))
        
        #print('coverage: ',coverage)
        self.qlty_NCV = coverage

    def NCV_GS3(self):
        """Combined NCV and GS3"""
        #print('call: ',self.qlty_NCV,self.qlty_GS3)
        if self.qlty_NCV == None:
            self.NCV()
         
        
        if self.qlty_GS3 == None:
            self.GS3()
            
        self.qlty_NCV=float(self.qlty_NCV)
        self.qlty_GS3=float(self.qlty_GS3)
        if self.qlty_NCV > 0 and self.qlty_GS3 > 0:
            self.qlty_NCV_GS3 = (math.sqrt (1.0 * self.qlty_GS3 * self.qlty_NCV))
        else:
            self.qlty_NCV_GS3=0.0
        self.qlty_NCV = None
        self.qlty_GS3 = None
        
    def GC(self, k=1):
        """GO ontology correctness (require Go term)"""
        if self.goterm_dict1 == None or self.goterm_dict2 == None:
            print("Please provide GO term for both species")
            return

        overlap = 0
        total = 0
        for (n1, n2) in self.mapping_set:    
            if n1 in self.goterm_dict1 and n2 in self.goterm_dict2:
                goSet1 = self.goterm_dict1[n1]
                goSet2 = self.goterm_dict2[n2]
                
                if len(goSet1) >=k and len(goSet2) >= k:
                    total = total + 1
                    inter = goSet1.intersection(goSet2)
                    if len(inter) >= k:
                        overlap = overlap + 1
        if total == 0:
            gc = 0
        else:
            gc = overlap * 1.0 / total
            self.qlty_GC = gc


    def GO_related_genes(self, nodes, go_dict):
        res = {}
        for node in nodes:
            if node in go_dict:
                go_set = go_dict[node]
                for go in go_set:
                    if go in res:
                        res[go].add(node)
                    else:
                        res[go] = set([node])
        return res

    def PF(self):
        """Protein function prediction (require Go term)"""
        if self.goterm_dict1 == None or self.goterm_dict2 == None:
            print("Please provide GO term for both species")
            return
        
        ncount1 = 0
        ncount2 = 0
        alncount = 0

        #ground truth go term set
        go_true_set = set([])
        for node in self.node_set1:
            if node in self.goterm_dict1:
                ncount1 = ncount1 + 1

        for node in self.node_set2:
            if node in self.goterm_dict2:
                ncount2 = ncount2 + 1

        for (n1, n2) in self.mapping_set:
            if n1 in self.goterm_dict1 and n2 in self.goterm_dict2:
                alncount = alncount + 1

        #the total number of node pairs from different networks, in which both nodes are annotated with at least one go term
        M = ncount1 * ncount2
        #the total number of aligned node pairs, in which both nodes are annotated with at least one go term
        n = alncount

        go_nodes1 = self.GO_related_genes(self.node_set1, self.goterm_dict1)#dictionary, where keyed on GO terms and the values are the genes annotated with the key
        go_nodes2 = self.GO_related_genes(self.node_set2, self.goterm_dict2)
        go_set = set(go_nodes1.keys()).union(set(go_nodes2.keys()))

        go_pred_set = set([])



        for go in go_set:

            #find all the nodes from network 1 associated with go
            gn1 = {}
            if go in go_nodes1:
                gn1 = go_nodes1[go]

            #find all the nodes from network 2 assocaited with go
            gn2 = {}
            if go in go_nodes2:
                gn2 = go_nodes2[go]

            #number of nodes from network 1 associated with go
            len1 = len(gn1)
            #number of nodes from network 2 associated with go
            len2 = len(gn2)

            if len1 + len2 < 3:
                continue

            for n1 in gn1:
                go_true_set.add((n1+'_y1', go))

            for n2 in gn2:
                go_true_set.add((n2+'_y2', go))

            #find the number of aligned nodes that are both annotated with g
            k = 0
            alignedNodes1 = {}
            alignedNodes2 = {}
            for (n1, n2) in self.mapping_set:
                if n1 in gn1 and n2 in gn2:
                    k = k + 1
                    if n1 in alignedNodes1:
                        alignedNodes1[n1].add(n2)
                    else:
                        alignedNodes1[n1] = set([n2])

                    if n2 in alignedNodes2:
                        alignedNodes2[n2].add(n1)
                    else:
                        alignedNodes2[n2] = set([n1])

            #now we start to hide node and transfer labels between aligned nodes
            for (n1, n2) in self.mapping_set:

                #we do not transfer go term go if neither of the aligned nodes are associated with it
                if n1 not in gn1 and n2 not in gn2:
                    continue

                #we transfer go from n1 to n2 if n1 is associated with go
                #we hide the label of n2, and more importantly the label of
                #n2 does not contribute to the significant test
                if n1 in gn1 and n2 in self.goterm_dict2:
                    #but we first need to compute statistic significance
                    #we decrease 1 if n2 is also associated with go
                    if n2 in gn2:
                        len2 = len2 - 1

                    #the number of node pairs between gn1 and gn2
                    #also, the maximum possible pairs that could be aligned between
                    #nodes from the two networks that are associated with go
                    N = len(gn1) * len(gn2)

                    #if only less than 3 nodes are assocaited with go term, we continue
                    if N == 0:
                        continue

                    #count the number of aligned node pairs, in which both nodes are associated with go
                    #n2 cannot be involved
                    if n2 in gn2:
                        k0 = k - len(alignedNodes2[n2])
                    else:
                        k0 = k
                    p = 1 - scipy.stats.hypergeom.cdf(k0-1, M, n, N)

                    #we predict n2 to be assocaited with go
                    if p < 0.05:
                        go_pred_set.add((n2+"_y2", go))

                if n2 in gn2 and n1 in self.goterm_dict1:
                    if n1 in gn1:
                        len1 = len1 - 1
                    N = len(gn1) * len(gn2)
                    if N == 0:
                        continue

                    if n1 in gn1:
                        k0 = k - len(alignedNodes1[n1])
                    else:
                        k0 = k

                    p = 1 - scipy.stats.hypergeom.cdf(k0-1, M, n, N)
                    if p < 0.05:
                        go_pred_set.add((n1+"_y1", go))

        (self.qlty_PF_P, self.qlty_PF_R, self.qlty_PF_F) = self.Accuray(go_pred_set, go_true_set)
