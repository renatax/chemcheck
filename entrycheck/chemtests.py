'''
Created on Mar 22, 2013

@author: ernia

perform a series of tests to verify the inputs.
should probably become a real object one day.
'''
#import the modified kyle scripts
from elementsvalidation.stereofix import KyleFix
#import string


def testem(formarguments):
    '''
    an ugly case handler for the various subtests: each function returns True/False;
    '''
    results={}
    
    
    results['name'],out= checkName(formarguments['name'])
    results['products']= checkProducts(formarguments['products'])
    results['conditions']=checkConditions(formarguments['conditions'])
    results['synthons']= checkSynthons(formarguments['synthons'])
    results['doi']= checkDoi(formarguments['doi'])
    results['reaction'],out, stuff= checkReaction(formarguments['reaction'], formarguments['products'])       
    print formarguments['products']        
    if False in results.values():
        return False
    else:
        return True


def checkName(name):
    import elementsvalidation.namevalidation
    out=""
    if name == "":
        out =out + "name field empty"
        return False
    else:
        if name =="notthisname": #test only! remove me!
            return False
        name=elementsvalidation.namevalidation.Name(name)
        if not name.isValid():
            return False, "Invalid name <br>"
        hits=name.getGoogleHits()
        
        return True,out
    
def checkProducts(prods):
    import elementsvalidation.abstractvalidation
    
    if prods==""or prods is None:
        print "product fails"
        return False
    else: 
        validator=elementsvalidation.abstractvalidation.AbstractElement(prods)
        return validator.checkbraces(prods)
    
    
def checkConditions(conds):
    if conds=="" or conds is None:
        print "conditions fails"
        return False
    else:
        #conds.split('>')
        
        return True
    

def checkDoi(name):
    if name == "":
        print "doi fails"
        return False
    else:
        return True

def checkSynthons(name):
    if name == "" or name=='0':
        print "synthons fails"
        return False
    else:
        return True

def checkReaction(reaction, prods, addH=False, env=True ):
    '''
    passes the reaction and the products to kyle test script, returns 
    the boolean "isOk" and the full output of the script "out";
    
    '''
    
    
    
    out="fails"
    stuff={}
    if reaction == "" or reaction is None:
        #print "reaction fails"
        return False, out, stuff
    else:
        #kyleok, out, stuff = rdkit_fix.kylescript(reaction,prods)
        #print prods
        #generate a kylefixobject
        kylefix=KyleFix(reaction,prods, addH,env)
        kyleok= kylefix.isOk()
        out=kylefix.getOut() + "<br>"
        stuff=kylefix.getStuff()
        #
        return kyleok, out, stuff #comment this line when ready
#         kyleok,out=rdkit_fix.kylescript(reaction, prods) 
#         print out
#         return kyleok
