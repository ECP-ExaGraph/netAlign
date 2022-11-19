# encoding: utf-8

import sys
import Evaluation

def main(argv=None):

   if len(argv) != 3:
       print('usage: main.py config_file output_file')
       sys.exit()
        
   config_file = argv[1].rstrip()
   output_file = argv[2].rstrip()
   graph1 = None
   graph2 = None
   mapping_name = None
   true_mapping_name = None
   goterm_name1 = None
   goterm_name2 = None
   measures = None
    
   eval_NC_P = None
   eval_NC_R = None
   eval_NC_F = None
   eval_GS3 = None
   eval_NCV = None
   eval_NCV_GS3 = None
   eval_GC = None
   eval_PF_P = None
   eval_PF_R = None
   eval_PF_F = None
    
   fi = open(config_file, 'r')
   for row in fi:
       (key, value) = row.rstrip().split(":")
       key = key.strip()
       value = value.strip()
       if key == "Network1":
           graph1 = value
       elif key == "Network2":
           graph2 = value
       elif key == "Aligned node pairs":
           mapping_name = value
       elif key == "Ground truth node mapping":
           true_mapping_name = value        
       elif key == "GO data for Network1":
           goterm_name1 = value
       elif key == "GO data for Network2":
           goterm_name2 = value
       elif key == "Measures":
           measures = value.split(",")      
 
   for measure in measures:
       measure = measure.strip()
       if measure == 'P-NC':
           eval_NC_P = True
       elif measure == 'R-NC':
           eval_NC_R = True
       elif measure == 'F-NC':
           eval_NC_F = True
       elif measure == 'NCV':
           eval_NCV = True
       elif measure == 'GS3':
           eval_GS3 = True
       elif measure == 'NCV-GS3':
           eval_NCV_GS3 = True
       elif measure == 'GC':
           eval_GC = True
       elif measure == 'P-PF':
           eval_PF_P = True     
       elif measure == 'R-PF':
           eval_PF_R = True 
       elif measure == 'F-PF':
           eval_PF_F = True            
    
    #print graph1,graph2,mapping_name, true_mapping_name, goterm_name1, goterm_name2
    #print eval_NC_P, eval_NC_R, eval_NC_F, eval_GS3, eval_NCV, eval_NCV_GS3, eval_GC, eval_PF_P, eval_PF_R, eval_PF_F
                                                                  
   quality = Evaluation.AlignmentQuality(graph1, graph2, mapping_name, true_mapping_name, goterm_name1, goterm_name2);
   result = quality.evaluate(eval_NC_P, eval_NC_R, eval_NC_F, eval_GS3, eval_NCV, eval_NCV_GS3, eval_GC, eval_PF_P, eval_PF_R, eval_PF_F)
   
   fo = open(output_file, 'w')
   for key in measures:
       if key in result:
           value = result[key]
           fo.write(key+"\t"+str(value)+"\n")
   fo.close()
   return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv))
