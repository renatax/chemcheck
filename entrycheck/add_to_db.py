#!/usr/bin/env python

import pymongo
from rdkit import Chem
from rdkit.Chem import AllChem

db = pymongo.MongoClient('localhost').datatest

fields = ['name','reaction_smarts','products','condition','synthon_count']
named_fields = ['Name','Reaction SMARTS', 'Retron SMILES', 'Conditions',
                'Synthon Count']
                
def start_query():
    full_string = raw_input('Please enter full reaction string:\t')
    if full_string:
        upload(full_string)
    else:
        print "Whatever, I didn't want to do anything anyway."
        return
    if input_another():
       start_query()
        
def upload(full_string):
    info = full_string.strip().split('\t')
    if len(info) > 5:
	document = {'reference':info.pop()}
    else:
	document = {}
    while(info):
        document[fields[len(info)]] = info.pop()
    print document
    document['products'] = document['products'].split('.')
    document['synthon_count'] = int(document['synthon_count'])
    print db.retro.count()
    if present_document(document):
        document['_id'] = db.retro.count()+1000
	print document
    	db.retro.save(document)
	print db.retro.count()
        return
    else:
        upload_full()
        return
        
def upload_full():
    print 'Please fill in the following fields'
    print 'For fields with many values, please separate values with "."'
    document = {}
    for item, field in zip(fields, named_fields):
        document[item] = raw_input('%s:\t' %field)
    document['products'] = document['products'].split('.')
    document['synthon_count'] = int(document['synthon_count'])
    doi = raw_input('Enter Reference?\t')
    if doi in ['Yes','yes','Y','y']:
	document['reference'] = raw_input('Enter DOI:\t')
    if present_document(document):
        document['_id'] = db.retro.count() + 1000 
        print document
        #db.retro.save(document)
        return
    else:
        print 'Please check input and try again'
        return
        
def present_document(document):
    for key, value in zip(named_fields, fields):
        print '%s:\t%s' %(key, str(document[value]))
    if not check_rdkit(document):
        print 'Input fails through RDKit, please recheck and try again'
        return False
    good = raw_input('Is this okay?\t')
    if good in ['yes','Yes','y','Y']:
        return True
    else:
        return False

def check_rdkit(document):
    rxn = AllChem.ReactionFromSmarts(document['reaction_smarts'])
    reactants = [Chem.MolFromSmarts(mol) for mol in document['products']]
    reactants = tuple(reactants)
    pset = rxn.RunReactants(reactants)
    if not pset:
        return False
    print 'Resulting synthons:'
    print '.'.join([Chem.MolToSmiles(mol) for mol in pset[0]])
    good = raw_input('Are these synthons correct?\t')
    if good in ['yes','Yes','y','Y']:
        return True
    else:
        return False

def input_another():
    another = raw_input('Put in another entry?\t')
    if another in ['Y','y','yes','Yes']:
        return True
    return False
    
if __name__ == '__main__':
    start_query()
