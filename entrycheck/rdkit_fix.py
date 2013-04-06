#!/usr/bin/env python2
import sys
import itertools
import copy
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

##########################################################################


def mol2smiles(mols):
    """Generate a SMILES for each molecule in the list"""
    return [Chem.MolToSmiles(m, isomericSmiles=True) for m in mols]


##########################################################################
def smiles2mol(smiles):
    """Convert each SMILES from the list to molecule"""
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    return [m for m in mols if m is not None]


##########################################################################
def run_reactants(rxn, reactants):
    # Number atoms in reactants
    number_atoms(reactants)

    # Create reactant atom map
    reactantAtomMap = create_reactant_map(reactants)

    # Run reaction transform
    try:
        productList = rxn.RunReactants(tuple(reactants))
    except:
        sys.stderr.write('Error: Applying transform failed.\n')
        sys.exit(1)

    # Create template atom maps
    [rTemplateList, pTemplateList, rTemplateAtomMaps, pTemplateAtomMaps,
        IsMatchChiral] = create_template_map(reactants, rxn)
    if len(IsMatchChiral) != len(productList):
        sys.stderr.write('Error: Template maps does not match products.\n')
        sys.exit()

    # Remove non-chiral matches AND Check for Ring-opening reactions
    newProductList = []
    for idx, products in enumerate(productList):
        if IsMatchChiral[idx]:
            newProducts = list(products)
            possibleOverlaps = True
            while possibleOverlaps:
                [possibleOverlaps, newProducts] = merge(newProducts)
            newProductList.append(tuple(newProducts))
    productList = tuple(newProductList)

    # Parse products to correct chirality
    for matchIdx, products in enumerate(productList):
        for product in products:
            for productAtom in product.GetAtoms():
                if not productAtom.HasProp('Core'):
                    # Unique atom ID
                    atomIdx = productAtom.GetIsotope()
                    # Map Atoms
                    tmp = reactantAtomMap[
                        atomIdx] if atomIdx in reactantAtomMap else None
                    reactantAtom = reactants[tmp[0]
                                             ].GetAtomWithIdx(tmp[1]) if not tmp == None else None
                    tmp = rTemplateAtomMaps[matchIdx][atomIdx] if atomIdx in rTemplateAtomMaps[matchIdx] else None
                    rTemplateAtom = rTemplateList[matchIdx][tmp[
                        0]].GetAtomWithIdx(tmp[1]) if tmp != None else None
                    tmp = pTemplateAtomMaps[matchIdx][atomIdx] if atomIdx in pTemplateAtomMaps[matchIdx] else None
                    pTemplateAtom = pTemplateList[matchIdx][tmp[
                        0]].GetAtomWithIdx(tmp[1]) if tmp != None else None
                    # Correct Chiraliry of Product Atom
                    fix_chirality(reactantAtom,
                                  productAtom, rTemplateAtom, pTemplateAtom)

    # Remove atom mappings
    for products in productList:
        for product in products:
            for atom in product.GetAtoms():  # remove atom mappings
                atom.SetIsotope(0)
                atom.ClearProp('Core')

    return productList


##########################################################################
def number_atoms(reactants):
    """Number atoms using the isotope field"""
    atomIdx = 1
    for reactant in reactants:
        for atom in reactant.GetAtoms():
            atom.SetIsotope(atomIdx)  # use Isotope field to store atom maps
            atom.SetProp('Core', '')
            atomIdx += 1


##########################################################################
def create_reactant_map(mols):
    """Creates a dictionary atomMap where the keys are the unique atom
    number stored in the Isotope field and values are 2-element lists
    containing the molecule number and the atom index"""
    molIdx = 0
    atomMap = {}
    for mol in mols:
        for atom in mol.GetAtoms():
            atomMap[atom.GetIsotope()] = [molIdx, atom.GetIdx()]
        molIdx += 1
    return atomMap


##########################################################################
def create_template_map(reactants, rxn):
    """Generates atom maps for reactant templates and product templates. Atom maps are stored
    in a list of dictionaries. Each element of the list corresponds to single match between
    reactants and reactant templates. The dictionary keys are the unique atom indices stored
    in the Isotope field. The dictionary values contain lists of two elements: 1) the template
    index (0...number of templates) and 2) the atom index (0...number of atoms in the template).
    Given a unique atom index, the atom map data structure allows one to rapidly locate that
    atom within the reactant and product templates.  NOTE: In the case of product templates
    only "mapped atoms" (atoms with atom map numbers in the templates) are included"""

    # Extract templates
    rTemplates = [rxn.GetReactantTemplate(
        x) for x in range(rxn.GetNumReactantTemplates())]
    pTemplates = [rxn.GetProductTemplate(
        x) for x in range(rxn.GetNumProductTemplates())]

    # Generate substructure matches between reactants and reactant templates
    rMatches = []
    rMatchesChiral = []
    for rTemplate, reactant in itertools.izip(rTemplates, reactants):
        matches = reactant.GetSubstructMatches(rTemplate, uniquify=0)
        rMatches.append(matches)
        rMatchesChiral.append([is_chiral_match(
            reactant, rTemplate, match) for match in matches])

    # Generate all permutations of matches
    IsMatchChiral = []
    for matches in itertools.product(*rMatchesChiral):
        IsMatchChiral.append(False if False in matches else True)

    # Generate atom map for reactant templates
    rTemplateAtomMaps = []
    rTemplateList = []
    matchIdx = 0
    for idx, matches in enumerate(itertools.product(*rMatches)):
        if IsMatchChiral[idx]:
            rTemplateAtomMaps.append({})
            rTemplateList.append(copy.deepcopy(
                rTemplates))  # make a copy (NOTE: DOESN'T COPY CHIRALITY)
            for rTemplateIdx, match in enumerate(matches):
                for rTemplateAtomIdx, rAtomIdx in enumerate(match):
                    rAtom = reactants[rTemplateIdx].GetAtomWithIdx(rAtomIdx)
                    rTemplateAtom = rTemplateList[matchIdx][
                        rTemplateIdx].GetAtomWithIdx(rTemplateAtomIdx)
                    rTemplateAtom.SetIsotope(rAtom.GetIsotope())
                    rTemplateAtom.SetChiralTag(rTemplates[rTemplateIdx].GetAtomWithIdx(rTemplateAtomIdx).GetChiralTag())
                    rTemplateAtomMaps[matchIdx][rAtom.GetIsotope(
                    )] = [rTemplateIdx, rTemplateAtomIdx]
            matchIdx += 1

    # Generate atom map for product templates
    pTemplateAtomMaps = []
    pTemplateList = []
    for matchIdx, rTemplates in enumerate(rTemplateList):
        pTemplateAtomMaps.append({})
        pTemplateList.append(list(pTemplates))  # make a copy
        for rTemplate in rTemplates:
            for rTemplateAtom in rTemplate.GetAtoms():
                if rTemplateAtom.HasProp('molAtomMapNumber'):
                    for pTemplateIdx, pTemplate in enumerate(pTemplateList[matchIdx]):
                        for pTemplateAtomIdx, pTemplateAtom in enumerate(pTemplate.GetAtoms()):
                            if (pTemplateAtom.HasProp('molAtomMapNumber') and
                                    (pTemplateAtom.GetProp('molAtomMapNumber') == rTemplateAtom.GetProp('molAtomMapNumber'))):
                                pTemplateAtom.SetIsotope(
                                    rTemplateAtom.GetIsotope())
                                pTemplateAtomMaps[matchIdx][rTemplateAtom.GetIsotope()] = [pTemplateIdx, pTemplateAtomIdx]

    return [rTemplateList, pTemplateList, rTemplateAtomMaps, pTemplateAtomMaps, IsMatchChiral]


##########################################################################
def is_chiral_match(molecule, template, match):
    '''
    Given a molecule, a template, and a match between them (as generated by
    GetSubstructureMatch()), this function checks whether or not the template
    has the same chirality as the molecule.
    '''
    for templateAtomIdx, templateAtom in enumerate(template.GetAtoms()):
        if ((templateAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW) or
                (templateAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)):
                moleculeAtom = molecule.GetAtomWithIdx(match[templateAtomIdx])
                if ((moleculeAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW) or
                        (moleculeAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)):

                    # Compare Chirality
                    templateNeigbors = [match[atom.GetIdx(
                    )] for atom in templateAtom.GetNeighbors()]
                    moleculeNeighbors = [atom.GetIdx(
                    ) for atom in moleculeAtom.GetNeighbors()]

                    # Check Parity of permutation
                    if are_perms_equal_parity(templateNeigbors, moleculeNeighbors):
                        if templateAtom.GetChiralTag() == moleculeAtom.GetChiralTag():
                            continue
                        else:
                            return False
                    else:
                        if templateAtom.GetChiralTag() != moleculeAtom.GetChiralTag():
                            continue
                        else:
                            return False
    return True


##########################################################################
def merge(products):
    '''
    Given a list of products, this function looks through all pairs of products to
    see if they share common atoms (i.e., atoms with the same unique atom labels
    stored in the Isotope field).  If so, it merges these products into a single
    product.
    '''
    for i1 in range(len(products)):
        for i2 in range(i1 + 1, len(products)):
            idx1 = [atom.GetIsotope() for atom in products[i1].GetAtoms()]
            idx2 = [atom.GetIsotope() for atom in products[i2].GetAtoms()]
            if len(set(idx1).intersection(set(idx2))) > 0:
                prod = combine_fragments(products[i1], products[i2])
                del products[i2]
                del products[i1]
                products.insert(i1, prod)
                return [True, products]
    return [False, products]


##########################################################################
def combine_fragments(m1, m2):
    '''
    Given two molecules that share common atoms (i.e., atoms with the same unique
    atom labels stored in the Isotope field), this function stitches them back
    together using common atoms and common bonds.
    '''
    # build dictionary of atomidx and atommaps
    atomidx = {}
    for atom in m1.GetAtoms():
        atomidx[atom.GetIsotope()] = atom.GetIdx()

    # make new editable molecule from m1
    newmol = Chem.EditableMol(m1)

    # add atoms from m2
    for atom in m2.GetAtoms():
        if atom.GetIsotope() not in atomidx:
            # add the atom
            idx = newmol.AddAtom(atom)
            atomidx[atom.GetIsotope()] = idx

    # add bonds from m2
    m = newmol.GetMol()
    for bond in m2.GetBonds():
        atom1 = atomidx[bond.GetBeginAtom().GetIsotope()]
        atom2 = atomidx[bond.GetEndAtom().GetIsotope()]
        if m.GetBondBetweenAtoms(atom1, atom2) is None:
            newmol.AddBond(atom1, atom2, bond.GetBondType())

    m = newmol.GetMol()
    Chem.SanitizeMol(m)
    return m


##########################################################################
def fix_chirality(rAtom, pAtom, rTAtom, pTAtom):
    '''
    Assign the appropriate tetrahedral chirality to the product atom, given the
    information present in the reaction atom, the reactant template atom, and the
    product template atom.
    '''
    if(rAtom == None):
        if(rTAtom == None):
            if(pTAtom == None):
                # This case corresponds to an atom in the product that was introduced
                # by the template.  RDkit handles this appropriately, so we do
                # nothing.
                return
            else:
                sys.stderr.write(
                    'Error: Should not be here (1) in fix_chirality().\n')
                sys.exit()
        else:
            sys.stderr.write(
                'Error: Should not be here (2) in fix_chirality().\n')
            sys.exit()
    elif(rAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED):
        if(rTAtom == None):
            # rAtom-Achiral, rTAtom-None
            sys.stderr.write(
                'Error: Should not be here (3) in fix_chirality().\n')
            sys.exit()
        elif(rTAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED):
            if(pTAtom == None):
                # rAtom-Achiral, rTAtom-Achiral, pTAtom-None
                sys.stderr.write(
                    'Error: Should not be here (4) in fix_chirality().\n')
                sys.exit()
            elif(pTAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED):
                # rAtom-Achiral, rTAtom-Achiral, pTAtom-Achiral --> UNSPECIFIED
                pAtom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
                return
            else:
                # rAtom-Achiral, rTAtom-Achiral, pTAtom-Chiral --> FORCE pTAtom
                # CHIRALITY
                pAtom.SetChiralTag(pTAtom.GetChiralTag())
                return
        else:  # rTAtom is chiral
            # rAtom-Achiral, rTAtom-Chiral --> UNSPECIFIED
            pAtom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
            return
    else:  # rAtom is chiral
        if(rTAtom == None):
            if(pTAtom == None):  # 19
                # PRESERVE CHIRALITY
                copy_chirality(rAtom, pAtom)
                return
            elif(pTAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED):
                # PRESERVE CHIRALITY
                copy_chirality(rAtom, pAtom)
                return
            else:  # 21
                # ADD CHIRALITY
                pAtom.SetChiralTag(pTAtom.GetChiralTag())
                return
        elif(rTAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED):
            if(pTAtom == None):
                # PRESERVE CHIRALITY
                copy_chirality(rAtom, pAtom)
                return
            elif(pTAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED):
                # rAtom-Chiral, rTAtom-Achiral, pTAtom-Achiral --> PRESERVE
                # rAtom CHIRALITY
                copy_chirality(rAtom, pAtom)
                return
            else:
                # rAtom-Chiral, rTAtom-Achiral, pTAtom-Chiral --> FORCE pTAtom
                # CHIRALITY
                pAtom.SetChiralTag(pTAtom.GetChiralTag())
                return
        else:  # rTAtom is chiral
            if(pTAtom == None):
                # PRESERVE CHIRALITY
                copy_chirality(rAtom, pAtom)
                return
            elif(pTAtom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED):
                # rAtom-Chiral, rTAtom-Chiral, pTAtom-Achiral --> REMOVE
                # CHIRALITY
                pAtom.SetChiralTag(Chem.rdchem.ChiralType.CHI_UNSPECIFIED)
                return
            elif(pTAtom.GetChiralTag() == rTAtom.GetChiralTag()):  # same chirality
                # PRESERVE CHIRALITY
                copy_chirality(rAtom, pAtom)
                return
            else:  # different chirality
                # INVERT CHIRALITY
                copy_chirality(rAtom, pAtom)
                pAtom.InvertChirality()
                return


##########################################################################
def remove_Hs(mol):
    '''
    RDkit's Chem.RemoveHs() function screws up cis/trans stereochemistry by erasing
    the directional bonds involving hydrogen.  The workaround is to move all directional
    bonds off of hydrogen-containing bonds.  Then we can use RDkit's RemoveHs() to
    remove all hydrogens.
    '''
    # Create a new molcule from scratch
    rdmol = Chem.Mol()
    rdemol = Chem.EditableMol(rdmol)

    # add atoms
    for atom in mol.GetAtoms():
        rdemol.AddAtom(atom)

    # add bonds
    for bond in mol.GetBonds():
        startAtomIdx = bond.GetBeginAtom().GetIdx()
        endAtomIdx = bond.GetEndAtom().GetIdx()
        rdemol.AddBond(startAtomIdx, endAtomIdx, bond.GetBondType())

    # convert editable molecule to normal molecule
    rdmol = rdemol.GetMol()

    # add stereobonds (but not to hydrogen atoms)
    for bondIdx, bond in enumerate(mol.GetBonds()):
        if bond.GetBondDir() != Chem.BondDir.NONE:
            # check for hydrogen
            if bond.GetBeginAtom().GetAtomicNum() == 1:
                otherAtom = bond.GetEndAtom()
            elif bond.GetEndAtom().GetAtomicNum() == 1:
                otherAtom = bond.GetBeginAtom()
            else:
                rdmol.GetBondWithIdx(bondIdx).SetBondDir(bond.GetBondDir())
                continue

            # Move bond directionality to the non-hydrogen atom
            for otherBond in otherAtom.GetBonds():
                if ((otherBond.GetBondType() == Chem.rdchem.BondType.SINGLE)
                        and (otherBond.GetOtherAtom(otherAtom).GetAtomicNum() != 1)):
                        if bond.GetBondDir() == Chem.rdchem.BondDir.ENDUPRIGHT:
                            rdmol.GetBondWithIdx(otherBond.GetIdx()).SetBondDir(Chem.rdchem.BondDir.ENDDOWNRIGHT)
                        elif bond.GetBondDir() == Chem.rdchem.BondDir.ENDDOWNRIGHT:
                            rdmol.GetBondWithIdx(otherBond.GetIdx(
                            )).SetBondDir(Chem.rdchem.BondDir.ENDUPRIGHT)
                        break

    return Chem.RemoveHs(rdmol)


##########################################################################
def copy_chirality(a1, a2):
    a2.SetChiralTag(a1.GetChiralTag())

    neighbors1 = [atom.GetIsotope() for atom in a1.GetNeighbors()]
    neighbors2 = [atom.GetIsotope() for atom in a2.GetNeighbors()]
    print neighbors1, neighbors2
    # Check if chirality was lost
    if len(neighbors2) < 3:
        return

    # Add Hydrogens
    while len(neighbors1) < 4:
        neighbors1.append(-1)
    while len(neighbors2) < 4:
        neighbors2.append(-1)

    # Check for ambiguous chirality
    d1 = list(set(neighbors1) - set(neighbors2))
    if len(d1) > 1:
        return  # ambiguous chirality
    elif len(d1) == 1:
        d2 = list(set(neighbors2) - set(neighbors1))
        neighbors2 = [d1[0] if x == d2[0] else x for x in neighbors2]

    # Check Parity of Permutation
    if not are_perms_equal_parity(neighbors1, neighbors2):
        a2.InvertChirality()


##########################################################################
def are_perms_equal_parity(perm0, perm1):
    """Check if 2 permutations are of equal parity.

    Assume that both permutation lists are of equal length
    and have the same elements. No need to check for these
    conditions.
    """
    perm1 = perm1[:]  # copy this list so we don't mutate the original

    transCount = 0
    for loc in range(len(perm0) - 1):                         # Do (len - 1) transpositions
        p0 = perm0[loc]
        p1 = perm1[loc]
        if p0 != p1:
            sloc = perm1[loc:].index(p0) + \
                loc          # Find position in perm1
            perm1[loc], perm1[sloc] = p0, p1          # Swap in perm1
            transCount += 1

    # Even number of transpositions means equal parity
    if (transCount % 2) == 0:
        return True
    else:
        return False

##########################################################################


def unique_products(transformSMARTS, reactantSMILES):
    # Initialize reaction
    rxn = AllChem.ReactionFromSmarts(transformSMARTS)
    rxn.Initialize()

    # Add explicit hydrogens
    reactants = [Chem.AddHs(
        Chem.MolFromSmiles(smiles)) for smiles in reactantSMILES.split('.')]

    # Run Reaction
    ps = run_reactants(rxn, reactants)

    # Get unique products
    uniqueProductsList = []
    for p in ps:
        uniqueProductsList.append('.'.join([Chem.MolToSmiles(
            remove_Hs(x), isomericSmiles=True) for x in p]))
    uniqueProductsList = list(set(uniqueProductsList))
    return uniqueProductsList


####################################################################################
# # Process the command line arguments.
# parser = argparse.ArgumentParser()
# parser.add_argument('smarts', type=str, help='transform SMARTS')
# parser.add_argument('smiles', type=str, nargs='+',
#                     help='space separated list of products SMILES')
# parser.add_argument('--explicitHs', type=bool, default=False,
#                     help='use explicit hydrogens?')
# args = parser.parse_args()
###################################################################################

'''hijacking kyle script for the webpage'''
def kylescript (smarts, smiles, explicitHs=False):
    
# Simple sanity checks. Never trust a user! ;-)
    reactant_smarts, product_smarts = [s.split('.') for s in smarts.split(
        '>>')]
    if len(reactant_smarts) != len(smiles):
        sys.stderr.write('Error: transform cannot be applied: '
                         'invalid number of reactant SMILES.\n')
        sys.exit(1)
    
    
    # Read the transform from SMARTS.
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        rxn.Initialize()
    except:
        sys.stderr.write('Error: Invalid transform SMARTS.\n')
        sys.exit(1)
    
    if explicitHs:
        # Add explicit hydrogens
        reactants = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in smiles]
    
        # Run Reaction
        ps = run_reactants(rxn, reactants)
    
        # Get unique products
        uniqueProductsList = []
        for p in ps:
            uniqueProductsList.append('.'.join([Chem.MolToSmiles(
                remove_Hs(x), isomericSmiles=True) for x in p]))
        uniqueProductsList = list(set(uniqueProductsList))
    else:
        reactants = [Chem.MolFromSmiles(s) for s in smiles]
        ps = run_reactants(rxn, reactants)
        uniqueProductsList = []
        for p in ps:
            uniqueProductsList.append(
                '.'.join([Chem.MolToSmiles(x, isomericSmiles=True) for x in p]))
        uniqueProductsList = list(set(uniqueProductsList))
    
    
    print 'Transform:  ' + smarts
    print 'Reactants:  ' + '.'.join(smiles)
    print 'Products:   ' + '   '.join(uniqueProductsList)
    out= 'Transform:  ' + smarts + '\n Reactants:  ' + '.'.join(smiles) + '\n Products:   ' + '   '.join(uniqueProductsList)
    return [True, out]
