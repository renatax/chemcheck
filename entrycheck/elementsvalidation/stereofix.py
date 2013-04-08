import copy
import itertools
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import BondDir, BondType, ChiralType
from environments import compact, contains_any
import logging
logging.basicConfig(level=logging.DEBUG)

class KyleFix():    
    '''
    INPUT:
        smarts: the daylight smarts of the reaction to fix
        smilesstr: the daylight smiles of compounds to process with the reaction
        explicitHs:  - defaults to False - 
    OUTPUT: 
        out: contains the error or a string containing the output of the reaction
        stuff: a dictionary containing reaction, products and reactants. 
        
        modification done by adding the class keyword, and calling runreactants in init with the parameters specified. This is allegedly temporary.
        so :
        added class declaration, and "self" as first argument in every function declaration and preposed to any function call; redirected all the print statements, and exits to
        to append to self.out.
    
    '''
    
    out=[]
    stuff={}
    isok=False
    
    
    def __init__(self, smarts, smilesstr, needsExplicitHs=False, needsEnvironments=True):
        self.out=[]
        self.stuff={}
        self.isok=False
        
        if smarts is None or smarts is "" or smilesstr is None or smilesstr is "":
            self.out.append('fill the Reaction and Test fields, and check the format\n')
            return  # False, self.out, {}
        logging.debug('\nbeginning\n input:\n smarts:  ' + smarts + 'smilesstr: ' + smilesstr)    
        smarts=smarts.encode('ascii')
        smilesstr=smilesstr.encode('ascii')
        smiles=[s for s in smilesstr.split()]
        # Simple sanity checks. Never trust a user! ;-)
        try:
            reactant_smarts, product_smarts = [s.split('.') for s in smarts.split('>>')]
        except:
            self.out.append("Error: the reaction is not correct: are the arrows there?")
            return
        
        if len(reactant_smarts) != len(smiles):
            self.out.append('<br> Error: transform cannot be applied: invalid number of reactant SMILES.<br>\n input smiles:' + str(smiles))
            self.out.append("<br> lenReactantSmarts= " + str(len(reactant_smarts)) +" len smiles = " + str(len(smiles)) +smarts)
            return# False, self.out, {}
        logging.debug('feeding runreactants with: smarts:  ' + smarts + ' smiles: ' + ".".join(smiles))
        try:
            #the join is ugly but necessary until we converge on a single format
            ps = self.run_reactants(smarts, ".".join(smiles), needsExplicitHs, needsEnvironments)
        except:
            self.out.append('<br> Error: operation failed - check the syntax of the various fields.\n ')
            return
        self.out.append('Transform:  ' + smarts + '<br>\n   ')
        self.out.append('Reactants: ' +  ".".join(smiles)) 
        self.out.append('Products:   ' + ' '.join(ps))
        self.stuff = {'transform':smarts, 'reactants':smiles, 'products':ps} 
        self.isok= True
    
    def getOut(self):
        return "".join(self.out)
    def getListOut(self):
        return self.out
    def getOutHtml(self):
        htmlout="<br>".join(self.out)
        return htmlout
    def getStuff(self):
        return self.stuff
    def isOk(self):
        return self.isok
    
    
    def mol2smiles(self, mols):
        """Generate a SMILES for each molecule in the list"""
        return [Chem.MolToSmiles(m, isomericSmiles=True) for m in mols]
    
    def smiles2mol(self,smiles):
        """Convert each SMILES from the list to molecule"""
        mols = [Chem.MolFromSmiles(s) for s in smiles]
        return [m for m in mols if m is not None]
        
    
    def run_reactants(self, rxnSmarts, reactantSmiles,
                      needsExplicitHs=False, needsEnvironments=True):
        """Return a list of unique reaction products.
    
        Function is a wrapper around the RDKit built-in RunReaction() which attempts to
        address its limitations related to handling stereochemistry in reaction cores.
        """
        # By default, rewrite the transform using SMARTS environments.
        logging.debug("processing run_reactants")
        if needsEnvironments:
            logging.debug('environment needed: compacting: '+ rxnSmarts)
            rxnSmarts = compact(rxnSmarts)
            logging.debug('environmental compacting done: '+ rxnSmarts)
        # Generate reaction from input SMARTS.
        try:
            logging.debug('trying rdkit reaction from smarts with ' + rxnSmarts)
            rxn = AllChem.ReactionFromSmarts(rxnSmarts)
            rxn.Initialize()
        except:
            logging.error('invalid SMARTS: {0}'.format(rxnSmarts))
            self.out.append('Error: invalid SMARTS: {0}.\n'.format(rxnSmarts))
            sys.exit(1)
    
        # Convert reactants SMILES to molecules and add explicit H atoms need if
        # requested.
        logging.debug('reaction initialized')
        reactants = [Chem.MolFromSmiles(s) for s in reactantSmiles.split('.')]
        if needsExplicitHs:
            reactants = [self.add_Hs(m) for m in reactants]
    
        if None not in reactants:
            [m.UpdatePropertyCache() for m in reactants]
        else:
            self.out.append('Error: invalid input SMILES.\n')
            sys.exit(1)
    
        # Number atoms in reactants.
        self.number_atoms(reactants)
    
        # Create reactant atom map.
        reactantAtomMap = self.create_reactant_map(reactants)
    
        # Run reaction transform.
        try:
            productList = rxn.RunReactants(tuple(reactants))
        except:
            self.out.append('Error: Applying transform failed.\n')
            sys.exit(1)
    
        # Create template atom maps.
        [rTemplateList, pTemplateList, rTemplateAtomMaps, pTemplateAtomMaps, IsMatchChiral] = self.create_template_map(reactants, rxn)
        if len(IsMatchChiral) != len(productList):
            self.out.append('Error: Template maps does not match products.\n')
            sys.exit(1)
    
        # Remove non-chiral matches AND check for ring-opening reactions.
        newProductList = []
        for idx, products in enumerate(productList):
            if IsMatchChiral[idx]:
                newProducts = list(products)
                possibleOverlaps = True
                while possibleOverlaps:
                    [possibleOverlaps, newProducts] = self.merge(newProducts)
                newProductList.append(tuple(newProducts))
        productList = tuple(newProductList)
    
        # Parse products to correct chirality.
        for matchIdx, products in enumerate(productList):
            for product in products:
                for productAtom in product.GetAtoms():
                    if not productAtom.HasProp('Core'):
                        # Unique atom ID
                        atomIdx = productAtom.GetIsotope()
    
                        # Map Atoms
                        tmp = reactantAtomMap[atomIdx] if atomIdx in reactantAtomMap.keys() else None
                        reactantAtom = reactants[tmp[0]].GetAtomWithIdx(tmp[1]) if tmp is not None else None
                        tmp = rTemplateAtomMaps[matchIdx][atomIdx] if atomIdx in rTemplateAtomMaps[matchIdx].keys() else None
                        rTemplateAtom = rTemplateList[matchIdx][tmp[0]].GetAtomWithIdx(tmp[1]) if tmp is not None else None
                        tmp = pTemplateAtomMaps[matchIdx][atomIdx] if atomIdx in pTemplateAtomMaps[matchIdx].keys() else None
                        pTemplateAtom = pTemplateList[matchIdx][tmp[0]].GetAtomWithIdx(tmp[1]) if tmp is not None else None
    
                        # Correct Chirality of Product Atom
                        self.fix_chirality(reactantAtom, productAtom, rTemplateAtom, pTemplateAtom)
  
        # Remove Hs if necessary
        if needsExplicitHs:
            tempProductList = []
            for products in productList:
                tempProductList.append([self.remove_Hs(p) for p in products])
            productList = tempProductList
        
        # Remove atom mappings.
        for products in productList:
            for product in products:
                for atom in product.GetAtoms():
                    atom.SetIsotope(0)
                    atom.ClearProp('Core')
    
        # Get unique product sets.
        uniqueProductSets = []
        for products in productList:
            uniqueProductSets.append('.'.join(set([Chem.MolToSmiles(p, isomericSmiles=True) for p in products])))
    
        return self.uniquify(uniqueProductSets)
    
    
    def uniquify(self, seq):
        """ remove duplicates from list """
        seen = set()
        seen_add = seen.add
        return [x for x in seq if x not in seen and not seen_add(x)]
    
    
    def are_perms_equal_parity(self, perm0, perm1):
        """Check if 2 permutations are of equal parity.
    
        Assume that both permutation lists are of equal length
        and have the same elements. No need to check for these
        conditions.
        """
        perm1 = perm1[:]  # copy this list so we don't mutate the original
    
        transCount = 0
        for loc in range(len(perm0) - 1):           # Do (len - 1) transpositions
            p0 = perm0[loc]
            p1 = perm1[loc]
            if p0 != p1:
                sloc = perm1[loc:].index(p0) + loc  # Find position in perm1
                perm1[loc], perm1[sloc] = p0, p1    # Swap in perm1
                transCount += 1
    
        # Even number of transpositions means equal parity
        if (transCount % 2) == 0:
            return True
        else:
            return False
    
    
    def number_atoms(self, reactants):
        """Number atoms using the Isotope field."""
        atomIdx = 1
        for reactant in reactants:
            for atom in reactant.GetAtoms():
                atom.SetIsotope(atomIdx)
                atom.SetProp('Core', '')
                atomIdx += 1
    
    
    def create_reactant_map(self, mols):
        """Map
    
        Creates a dictionary atomMap where the keys are the unique atom number
        stored in the Isotope field and values are 2-element lists containing the
        molecule number and the atom index.
        """
        molIdx = 0
        atomMap = {}
        for mol in mols:
            for atom in mol.GetAtoms():
                atomMap[atom.GetIsotope()] = [molIdx, atom.GetIdx()]
            molIdx += 1
        return atomMap
    
    
    def create_template_map(self, reactants, rxn):
        """Generate atom  maps for reactant and product templates.
    
        Atom maps are stored in a list of dictionaries. Each element of the list
        corresponds to single match between reactants and reactant templates. The
        dictionary keys are the unique atom indices stored in the Isotope field.
        The dictionary values contain lists of two elements: 1) the template index
        (0...number of templates) and 2) the atom index (0...number of atoms in the
        template).  Given a unique atom index, the atom map data structure allows
        one to rapidly locate that atom within the reactant and product templates.
        NOTE: In the case of product templates only "mapped atoms" (atoms with atom
        map numbers in the templates) are included.
        """
    
        # Extract templates.
        rTemplates = [rxn.GetReactantTemplate(x)
                      for x in range(rxn.GetNumReactantTemplates())]
        pTemplates = [rxn.GetProductTemplate(x)
                      for x in range(rxn.GetNumProductTemplates())]
    
        # Generate substructure matches between reactants and reactant templates
        rMatches = []
        rMatchesChiral = []
        for rTemplate, reactant in itertools.izip(rTemplates, reactants):
            matches = reactant.GetSubstructMatches(rTemplate, uniquify=0)
            rMatches.append(matches)
            rMatchesChiral.append([self.is_chiral_match(reactant, rTemplate, match) for match in matches])
    
        # Generate all permutations of matches.
        IsMatchChiral = []
        for matches in itertools.product(*rMatchesChiral):
            IsMatchChiral.append(False if False in matches else True)
    
        # Generate atom map for reactant templates.
        rTemplateAtomMaps = []
        rTemplateList = []
        matchIdx = 0
        for idx, matches in enumerate(itertools.product(*rMatches)):
            if IsMatchChiral[idx]:
                rTemplateAtomMaps.append({})
                rTemplateList.append(copy.deepcopy(rTemplates))  # make a copy (NOTE: DOESN'T COPY CHIRALITY)
                for rTemplateIdx, match in enumerate(matches):
                    for rTemplateAtomIdx, rAtomIdx in enumerate(match):
                        rAtom = reactants[rTemplateIdx].GetAtomWithIdx(rAtomIdx)
                        rTemplateAtom = rTemplateList[matchIdx][rTemplateIdx].GetAtomWithIdx(rTemplateAtomIdx)
                        rTemplateAtom.SetIsotope(rAtom.GetIsotope())
                        rTemplateAtom.SetChiralTag(rTemplates[rTemplateIdx].GetAtomWithIdx(rTemplateAtomIdx).GetChiralTag())
                        rTemplateAtomMaps[matchIdx][rAtom.GetIsotope()] = [rTemplateIdx, rTemplateAtomIdx]
                matchIdx += 1
    
        # Generate atom map for product templates.
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
                                    (pTemplateAtom.GetProp('molAtomMapNumber') ==
                                     rTemplateAtom.GetProp('molAtomMapNumber'))):
                                    pTemplateAtom.SetIsotope(rTemplateAtom.GetIsotope())
                                    pTemplateAtomMaps[matchIdx][rTemplateAtom.GetIsotope()] = [pTemplateIdx, pTemplateAtomIdx]
    
        return [rTemplateList, pTemplateList, rTemplateAtomMaps, pTemplateAtomMaps, IsMatchChiral]
    
    
    def is_chiral_match(self, molecule, template, match):
        """Return True if a molecule and a template has the same chirality.
    
        Given a molecule, a template, and a match between them (as generated by
        GetSubstructureMatch()), this function checks whether or not the template
        has the same chirality as the molecule.
    
        Note
        This version works only with "development" (i.e. from github or svn repo) version
        of RDKit. Current stable release (2012_12_1) does not support stereochemistry in
        substructure matching.
        """
    
        #### Extract the match as a separate molecule
        # label the match atoms
        for atomIdx in match:
            molecule.GetAtomWithIdx(atomIdx).SetProp('match', '')
        # create an editable molecule
        fragmentMol = Chem.EditableMol(molecule)
        # delete all the unlabeled atoms
        atomIdxList = range(molecule.GetNumAtoms())
        for a in molecule.GetAtoms():
            if not a.HasProp('match'):
                Idx2 = next(i2 for i2, i1 in enumerate(atomIdxList) if i1 == a.GetIdx())
                fragmentMol.RemoveAtom(Idx2)
                del atomIdxList[Idx2]
        # convert to normal molecule
        fragmentMol = fragmentMol.GetMol()
        Chem.SanitizeMol(fragmentMol)
        # generate smiles
        fragmentSmiles = Chem.MolToSmiles(fragmentMol, isomericSmiles=True)
        # remove labels from molecule
        molecule.ClearProp('match')
    
        #### Generate a molecule from the template
        templateMol = Chem.EditableMol(Chem.Mol())
        # Copy atoms and bonds
        for atomIdx, atom in enumerate(template.GetAtoms()):
            templateMol.AddAtom(molecule.GetAtomWithIdx(match[atomIdx]))
        for bondIdx, bond in enumerate(template.GetBonds()):
            templateMol.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
        # Convert to molecule
        templateMol = templateMol.GetMol()
        Chem.SanitizeMol(templateMol)
        # Copy stereochemistry, both for atoms and bonds
        for atomIdx, atom in enumerate(templateMol.GetAtoms()):
            atom.SetChiralTag(template.GetAtomWithIdx(atomIdx).GetChiralTag())
        for bondIdx, bond in enumerate(templateMol.GetBonds()):
            bond.SetBondType(template.GetBondWithIdx(bondIdx).GetBondType())
            bond.SetBondDir(template.GetBondWithIdx(bondIdx).GetBondDir())
        # Generate smiles
        templateSmiles = Chem.MolToSmiles(templateMol, isomericSmiles=True)
    
        # if one has chirality and the other doesn't, we want a match.  must look for the non-chiral pattern
        # in the chiral molecule
        chiralSet = set(['@', '@@', '\\', '/'])
        if contains_any(templateSmiles, chiralSet):
            return Chem.MolFromSmiles(templateSmiles).HasSubstructMatch(Chem.MolFromSmiles(fragmentSmiles), useChirality=True)
        else:
            return Chem.MolFromSmiles(fragmentSmiles).HasSubstructMatch(Chem.MolFromSmiles(templateSmiles), useChirality=True)
    
    #DUPLICATE
#     def contains_any(str, set):
#         for c in set:
#             if c in str:
#                 return True
#         return False
#     
#     
    def merge(self, products):
        """Merge products sharing the same atoms.
    
        Given a list of products, this function looks through all pairs of products to
        see if they share common atoms (i.e., atoms with the same unique atom labels
        stored in the Isotope field).  If so, it merges these products into a single
        product.
        """
        for i1 in range(len(products)):
            for i2 in range(i1 + 1, len(products)):
                idx1 = [atom.GetIsotope() for atom in products[i1].GetAtoms()]
                idx2 = [atom.GetIsotope() for atom in products[i2].GetAtoms()]
                if len(set(idx1).intersection(set(idx2))) > 0:
                    prod = self.combine_fragments(products[i1], products[i2])
                    del products[i2]
                    del products[i1]
                    products.insert(i1, prod)
                    return [True, products]
        return [False, products]
    
    
    def combine_fragments(self, m1, m2):
        """Stich two molecules together.
    
        Given two molecules that share common atoms (i.e., atoms with the same unique
        atom labels stored in the Isotope field), this function stitches them back
        together using common atoms and common bonds.
        """
    
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
    
    
    def fix_chirality(self, rAtom, pAtom, rTAtom, pTAtom):
        """Correct chirality.
    
        Assign the appropriate tetrahedral chirality to the product atom, given the
        information present in the reaction atom, the reactant template atom, and the
        product template atom.
        """
        if rAtom is None:
            if rTAtom is None:
                if pTAtom is None:
                    # This case corresponds to an atom in the product that was introduced
                    # by the template.  RDkit handles this appropriately, so we do nothing.
                    return
                else:
                    self.out.append('Error: Should not be here (1) in fix_chirality().\n')
                    sys.exit()
            else:
                self.out.append('Error: Should not be here (2) in fix_chirality().\n')
                sys.exit()
        elif rAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
            if rTAtom is None:
                # rAtom-Achiral, rTAtom-None
                self.out.append('Error: Should not be here (3) in fix_chirality().\n')
                sys.exit()
            elif rTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                if pTAtom is None:
                    # rAtom-Achiral, rTAtom-Achiral, pTAtom-None
                    self.out.append('Error: Should not be here (4) in fix_chirality().\n')
                    sys.exit()
                elif(pTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED):
                    # rAtom-Achiral, rTAtom-Achiral, pTAtom-Achiral --> UNSPECIFIED
                    pAtom.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
                    return
                else:
                    # rAtom-Achiral, rTAtom-Achiral, pTAtom-Chiral --> FORCE pTAtom CHIRALITY
                    pAtom.SetChiralTag(pTAtom.GetChiralTag())
                    return
            else:  # rTAtom is chiral
                # rAtom-Achiral, rTAtom-Chiral --> UNSPECIFIED
                pAtom.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
                return
        else:  # rAtom is chiral
            if rTAtom is None:
                if pTAtom is None:  # 19
                    # Preserve chirality.
                    self.copy_chiralilty(rAtom, pAtom)
                    return
                elif pTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                    # Preserve chirality.
                    self.copy_chiralilty(rAtom, pAtom)
                    return
                else:  # 21
                    # Add chirality.
                    pAtom.SetChiralTag(pTAtom.GetChiralTag())
                    return
            elif rTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                if pTAtom is None:
                    # Preserve chirality.
                    self.copy_chiralilty(rAtom, pAtom)
                    return
                elif(pTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED):
                    # rAtom-Chiral, rTAtom-Achiral, pTAtom-Achiral --> PRESERVE rAtom CHIRALITY
                    self.copy_chiralilty(rAtom, pAtom)
                    return
                else:
                    # rAtom-Chiral, rTAtom-Achiral, pTAtom-Chiral --> FORCE pTAtom CHIRALITY
                    pAtom.SetChiralTag(pTAtom.GetChiralTag())
                    return
            else:  # rTAtom is chiral
                if pTAtom is None:
                    # Preserve chirality.
                    self.copy_chiralilty(rAtom, pAtom)
                    return
                elif pTAtom.GetChiralTag() == ChiralType.CHI_UNSPECIFIED:
                    # rAtom-Chiral, rTAtom-Chiral, pTAtom-Achiral --> REMOVE CHIRALITY
                    pAtom.SetChiralTag(ChiralType.CHI_UNSPECIFIED)
                    return
                elif pTAtom.GetChiralTag() == rTAtom.GetChiralTag():  # same chirality
                    # Preserve chirality.
                    self.copy_chiralilty(rAtom, pAtom)
                    return
                else:
                    # Different chiralitiy, invert it.
                    self.copy_chiralilty(rAtom, pAtom)
                    pAtom.InvertChirality()
                    return
    
    
    def copy_chirality(self, a1, a2):
        a2.SetChiralTag(a1.GetChiralTag())
    
        neighbors1 = [atom.GetIsotope() for atom in a1.GetNeighbors()]
        neighbors2 = [atom.GetIsotope() for atom in a2.GetNeighbors()]
    
        # Check if chirality was lost
        if len(neighbors2) < 3:
            return
    
        # Add Hydrogens
        while len(neighbors1) < 4:
            neighbors1.append(-1)
        while len(neighbors2) < 4:
            neighbors2.append(-1)
    
        # Check for ambiguous chirality.
        d1 = list(set(neighbors1) - set(neighbors2))
        if len(d1) > 1:
            return
        elif len(d1) == 1:
            d2 = list(set(neighbors2) - set(neighbors1))
            neighbors2 = [d1[0] if x == d2[0] else x for x in neighbors2]
    
        # Check Parity of Permutation
        if not self.are_perms_equal_parity(neighbors1, neighbors2):
            a2.InvertChirality()
    
    
    def copy_chirality_Hs(self, a1, a2):
        a2.SetChiralTag(a1.GetChiralTag())
    
        neighbors1 = [(atom.GetIsotope() if ((atom.GetAtomicNum() != 1) and (atom.GetIsotope() != 0)) else atom.GetSymbol()) for atom in a1.GetNeighbors()]
        neighbors2 = [(atom.GetIsotope() if ((atom.GetAtomicNum() != 1) and (atom.GetIsotope() != 0)) else atom.GetSymbol()) for atom in a2.GetNeighbors()]
    
        # Add Hydrogens
        while len(neighbors1) < 4:
            neighbors1.append('H')
        while len(neighbors2) < 4:
            neighbors2.append('H')
    
        # Check for ambiguous chirality.
        d1 = list(set(neighbors1) - set(neighbors2))
        if len(d1) > 1:
            return
        elif len(d1) == 1:
            d2 = list(set(neighbors2) - set(neighbors1))
            neighbors2 = [d1[0] if x == d2[0] else x for x in neighbors2]
    
        # Check Parity of Permutation
        if not self.are_perms_equal_parity(neighbors1, neighbors2):
            a2.InvertChirality()
    

    
    def add_Hs(self, mol):
        """Add explicit hydrogens.
    
        Returns a new molecule with H atoms added with atoms with implicit valence.
        """
        # Get a list of implicit H atoms to be added for each atom in the molecule.
        implicitH = []
        for atom in mol.GetAtoms():
            implicitH.append(atom.GetImplicitValence())
    
        # Check to see if we need to do anything
        if implicitH == []:
            return mol
        elif all([v == 0 for v in implicitH]):
            return mol
    
        #Define H atom.
        mol_H = Chem.MolFromSmiles('[H]')
        newatom = mol_H.GetAtoms()[0]
    
        # Define an editable molecule similar to given molecule (Don't just copy
        # because then we are unable to change stereochemistry information)
        rdemol = Chem.EditableMol(Chem.Mol())
        for atom in mol.GetAtoms():
            rdemol.AddAtom(atom)
        for bond in mol.GetBonds():
            rdemol.AddBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    
        # Add H atom in place of implicit H atoms in the editable molecule.
        for idx, i in enumerate(implicitH):
            if i != 0:
                for j in range(i):
                    newatomidx = rdemol.AddAtom(newatom)
                    rdemol.AddBond(idx, newatomidx)
    
        # Convert editable molecule to normal molecule.
        rdmol = rdemol.GetMol()
        for atom in rdmol.GetAtoms():
            if atom.GetAtomicNum() == 1:
                atom.GetBonds()[0].SetBondType(BondType.SINGLE)
        # Add Chirality information
        for atomIdx, atom in enumerate(mol.GetAtoms()):
            rdmol.GetAtomWithIdx(atomIdx).SetChiralTag(atom.GetChiralTag())
        # Add Bond information
        for bond in mol.GetBonds():
            newBond = rdmol.GetBondBetweenAtoms(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            newBond.SetBondType(bond.GetBondType())
            newBond.SetBondDir(bond.GetBondDir())
        try:
            Chem.SanitizeMol(rdmol)

            return rdmol
        except:
            self.out.append('Error: output molecule of add_Hs() not valid.\n')
    
    
    def remove_Hs(self, mol):
        """Remove explicit hydrogens from a molecule.
    
        Returns a new molecule after removing all the explicit H atoms present in
        the molecule and keeping stereo-chemistry of molecule preserved.
        """
    
        # Add stereobonds, but not to hydrogen atoms.
        for bondIdx, bond in enumerate(mol.GetBonds()):
            if bond.GetBondDir() != Chem.BondDir.NONE:
                # check for hydrogen
                if bond.GetBeginAtom().GetAtomicNum() == 1:
                    otherAtom = bond.GetEndAtom()
                elif bond.GetEndAtom().GetAtomicNum() == 1:
                    otherAtom = bond.GetBeginAtom()
                else:
                    mol.GetBondWithIdx(bondIdx).SetBondDir(bond.GetBondDir())
                    continue
    
                # Move bond directionality to the non-hydrogen atom
                for otherBond in otherAtom.GetBonds():
                    if ((otherBond.GetBondType() == BondType.SINGLE)
                        and (otherBond.GetOtherAtom(otherAtom).GetAtomicNum() != 1)):
                            if bond.GetBondDir() == BondDir.ENDUPRIGHT:
                                mol.GetBondWithIdx(otherBond.GetIdx()).SetBondDir(BondDir.ENDDOWNRIGHT)
                            elif bond.GetBondDir() == BondDir.ENDDOWNRIGHT:
                                mol.GetBondWithIdx(otherBond.GetIdx()).SetBondDir(BondDir.ENDUPRIGHT)
                            break
    

        # Create a new molecule from scrach
        rdemol = Chem.EditableMol(Chem.Mol())
    
    # add all non-hydrogen atoms
        atomMap = {}
        for atomIdx, atom in enumerate(mol.GetAtoms()):
            if atom.GetAtomicNum() != 1:
                atomMap[atomIdx] = rdemol.AddAtom(atom)
        for bond in mol.GetBonds():
            if ((bond.GetBeginAtom().GetAtomicNum() != 1) and
                (bond.GetEndAtom().GetAtomicNum() != 1)):
                rdemol.AddBond(atomMap[bond.GetBeginAtomIdx()], atomMap[bond.GetEndAtomIdx()])
        
        # Convert editable molecule to normal molecule.
        rdmol = rdemol.GetMol()
    # Add Chirality information
        for atomIdx, atom in enumerate(mol.GetAtoms()):
            if atom.GetAtomicNum() != 1:
                if atom.GetChiralTag() != ChiralType.CHI_UNSPECIFIED:
                    self.copy_chirality_Hs(atom, rdmol.GetAtomWithIdx(atomMap[atomIdx]))
        # Add Bond information
        for bond in mol.GetBonds():
            if ((bond.GetBeginAtom().GetAtomicNum() != 1) and
                (bond.GetEndAtom().GetAtomicNum() != 1)):
                newBond = rdmol.GetBondBetweenAtoms(atomMap[bond.GetBeginAtomIdx()], atomMap[bond.GetEndAtomIdx()])
                newBond.SetBondType(bond.GetBondType())
                newBond.SetBondDir(bond.GetBondDir())        
        
        try:
            Chem.SanitizeMol(rdmol)
            return rdmol
        except:
            self.out.append('invalid molecule: sanitization failed')
