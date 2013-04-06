'''
Created on Mar 23, 2013

@author: Andrea Cadeddu andrea.cadeddu@northwestern.edu
'''

class AbstractElement(object):
    '''
    takes care of the default ops
    '''
    NO_NO_CHARS="!'.\!@#$%^&*()"

    def __init__(self, params):
        '''
        Constructor
        '''
        self.me=params
        
    def get(self):
        return self.me
    def __repr__(self):
        return str(self.me)
    def __str__(self):
        return str(self.me)
    def __add__(self, whatelse):
        return str(self.me) + str(whatelse)
    def __radd__(self, whatelse):
        return str(self.me) + str(whatelse)