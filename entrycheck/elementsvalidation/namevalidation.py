'''
Created on Mar 23, 2013
Generate an object NAME, responsible to check, counts and do whathever it might be necessary on that field.
such as: get the hits on Google, and many more awesomeness. 
@author: Andrea Cadeddu andrea.cadeddu@northwestern.edu
'''
import abstractvalidation
class Name(abstractvalidation.AbstractElement):
    
    def __init__(self, name =""):
        self.me=name
        self.name=name
        
    def isValid(self):
        if self.name=="" or self.name is None:
            return False
        if self.NO_NO_CHARS in self.name: #this shall be done with regexp
            return False
        else:
            return True
        
    def getGoogleHits(self):
        ''' 
        searches the big G and return the number of hits associated to self.name
        temporarily hijacked to return only 0.
        UNDER CONSTRUCTION
        '''
        
        hits=0 
        return str(hits)
    
