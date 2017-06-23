'''
Created on 2010-07-26

@author: Shobhit Jain

@contact: shobhit@cs.toronto.edu
'''

import urllib2

class Parser(object):
    '''
    Parser class implement functions for parsing input data files.
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        super(Parser, self).__init__()
        
        
     
    def _obo_parser(self, obo_file):    
        '''
        Parser for Gene Ontology obo files. go_annotations variable is 
        updated.
        '''
#        file = urllib2.urlopen(obo_file)
        file = open(obo_file)
        flag = 0
        obsolete_flag = 0
        for line in file:
            line = line.strip()
            if line.startswith('[Term]'): 
                flag = 1
                node = ''
                name = ''
                domain = ''
                parent = set()
                obsolete_flag = 0
            elif flag == 1 and line == '' and obsolete_flag != 1: 
                flag = 0
                self._add_node(node)
                for term in parent:
                    self._add_node(term)
                    self._add_edge(term, node)
                self.go_annotations[node] = {'name':name, 
                                             'domain':domain, 
                                             'gene':set(), 
                                             'ancestors':set(), 
                                             'cluster':{}
                                             }
            elif flag == 1 and line.startswith("id"):
                node = line.split('id:')[1].strip()   
            elif flag == 1 and line.startswith("namespace"):
                domain = line.split('namespace:')[1].strip()
            elif flag == 1 and line.startswith("name"):
                name = line.split('name:')[1].strip()
            elif flag == 1 and line.startswith("is_a"):
                parent.add(line.split(' ')[1])
            elif flag == 1 and line.startswith("relationship"):
                parent.add(line.split(' ')[2])
            elif flag == 1 and line.startswith("is_obsolete: true"):
                obsolete_flag = 1
        file.close()

                
                
                
    def _go_annotations(self, taxid, cd):
        '''
        Parser for gene annotation file (SGD/human). go_annotations variable is updated.
        '''
        file = urllib2.urlopen("http://www.ebi.ac.uk/QuickGO/GAnnotation?tax=%s&format=tsv&limit=1000000000&evidence=IDA,IPI,EXP" % taxid)
        for line in file:
            line = line.strip()
            if line != "" and not line.startswith('!'):
                line = line.split('\t')
                term = line[6].strip()
                gene = line[1].strip()
                code = line[9].strip()
                if term in self.go_annotations and code != cd:
                    self.go_annotations[term]['gene'].add(gene)
                    
                
                
                

                        
                        
