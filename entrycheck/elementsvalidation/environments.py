import re
import Queue
from rdkit import Chem
from rdkit.Chem import AllChem
import logging
logging.basicConfig(level=logging.DEBUG)
## DUPLICATE!!
def contains_any(str, set):
    for c in set:
        if c in str:
            logging.debug("str contains one from set")
            return True
    logging.debug("str don't contains set")
    return False


def compact(transform):
    """Compact transform.

    Identifies core atoms in a transform and rewrites it using SMARTS
    environments.
    """
    chiralSet = set('@\\/')
    if any(c in chiralSet for c in transform.split('>>')[0]):
        return transform

    # Generate a reaction.
    try:
        logging.debug("trying with rdkit reaction from smarts" + transform )
        rxn = AllChem.ReactionFromSmarts(transform)
        logging.debug("done "+ transform +"" + str(rxn))
    except:
        logging.debug("failed with rdkit reaction from smarts" + transform+ "type: "+ str(type(transform)))
        pass
        #return
    logging.debug("rxn generated" + str(rxn))
    # Find its reactant and product templates.
    rTemplates = [rxn.GetReactantTemplate(i) for i in range(rxn.GetNumReactantTemplates())]
    pTemplates = [rxn.GetProductTemplate(i) for i in range(rxn.GetNumProductTemplates())]
    logging.debug("rTemplates: " + str(rTemplates))
    logging.debug("pTemplates: " + str(pTemplates))
    # Generate reactant and product atom maps (key = molAtomMapNumber, value =
    # [molecule #, atom #] )
    rAtomMap = make_atom_map(rTemplates)
    pAtomMap = make_atom_map(pTemplates)

    # Find core atoms and mark with property 'Core'.
    rCoreAtoms = find_cores(rTemplates, rAtomMap, pTemplates, pAtomMap)

    # Extract environment smarts.
    rAtomEnvironment = {}
    for rAtom in rCoreAtoms:
        molAtomMapNumber = rAtom.GetProp('molAtomMapNumber')
        temp = rAtomMap[molAtomMapNumber]
        smarts = get_environment(rAtom, rTemplates[temp[0]])
        if smarts != '*':
            rAtomEnvironment[molAtomMapNumber] = smarts

    # Strip both the reactant and product templates down to their cores.
    rCoreTemplateString = strip(rTemplates)
    pCoreTemplateString = strip(pTemplates)

    # Add environments to reactant template string.
    for molAtomMapNumber in rAtomEnvironment:
        patt = ':%d]' % int(molAtomMapNumber)
        repl = '$(' + rAtomEnvironment[molAtomMapNumber] + '):%d]' % int(molAtomMapNumber)
        rCoreTemplateString = re.sub(patt, repl, rCoreTemplateString)

    # Concatenate final transform
    newTransform = '>>'.join([rCoreTemplateString, pCoreTemplateString])
    return re.sub('-,:', '', newTransform)


def find_cores(rTemplates, rAtomMap, pTemplates, pAtomMap):
    """Identify and mark core atoms."""

    rCoreAtoms = []
    for idx, rTemplate in enumerate(rTemplates):
        for rAtom in rTemplate.GetAtoms():
            if rAtom.HasProp('molAtomMapNumber'):
                molAtomMapNumber = rAtom.GetProp('molAtomMapNumber')
                if molAtomMapNumber in pAtomMap:
                    temp = pAtomMap[molAtomMapNumber]
                    pAtom = pTemplates[temp[0]].GetAtomWithIdx(temp[1])
                    if not are_atoms_equal(rAtom, pAtom):
                        rCoreAtoms.append(rAtom)
                        rAtom.SetProp('Core', '')
                        pAtom.SetProp('Core', '')
    return rCoreAtoms


def make_atom_map(templates):
    """Make atom map."""
    atomMap = {}
    for idx, tmpl in enumerate(templates):
        for atom in tmpl.GetAtoms():
            if atom.HasProp('molAtomMapNumber'):
                atomMap[atom.GetProp('molAtomMapNumber')] = [idx, atom.GetIdx()]
    return atomMap


def strip(templates):
    """Strip templates to their cores.

    From each supplied template, function removes atoms which are not part of the
    template core.
    """

    cores = []
    for tmpl in templates:
        core = Chem.EditableMol(tmpl)

        atomIdxList = range(tmpl.GetNumAtoms())
        for atom in tmpl.GetAtoms():
            if not atom.HasProp('Core') and atom.HasProp('molAtomMapNumber'):
                idx = next(i2 for i2, i1 in enumerate(atomIdxList) if i1 == atom.GetIdx())
                core.RemoveAtom(idx)
                del atomIdxList[idx]

        core = core.GetMol()
        Chem.SanitizeMol(core)

        cores.append(Chem.MolToSmarts(core, isomericSmiles=True))

    return '.'.join(cores)


def are_atoms_equal(rAtom, pAtom):
    """Return True if atoms have the same environment, False otherwise."""

    # First compare neighboring atoms.
    rNeighbors = [(atom.GetProp('molAtomMapNumber') if atom.HasProp('molAtomMapNumber') else '#%d' % atom.GetAtomicNum()) for atom in rAtom.GetNeighbors()]
    pNeighbors = [(atom.GetProp('molAtomMapNumber') if atom.HasProp('molAtomMapNumber') else '#%d' % atom.GetAtomicNum()) for atom in pAtom.GetNeighbors()]

    if sorted(rNeighbors) != sorted(pNeighbors):
        return False

    # Get indices from sorted neighbor lists.
    rNeighborIdxs = [i[0] for i in sorted(enumerate(rNeighbors), key=lambda x: x[1])]
    pNeighborIdxs = [i[0] for i in sorted(enumerate(pNeighbors), key=lambda x: x[1])]

    # Now compare sorted bonds.
    rBonds = rAtom.GetBonds()
    pBonds = pAtom.GetBonds()

    for idx, rNeighborIdx in enumerate(rNeighborIdxs):
        rBond = rBonds[rNeighborIdx]
        pBond = pBonds[pNeighborIdxs[idx]]
        if rBond.GetBondType() != pBond.GetBondType():
            return False
        elif rBond.GetBondDir() != pBond.GetBondDir():
            return False
    return True


def get_environment(atom, mol):
    if not atom.HasProp('molAtomMapNumber'):
        print 'Error: get_environment only works for mapped atoms.'
        return None
    else:
        molAtomMapNumber = atom.GetProp('molAtomMapNumber')

    # Graph traversal to build a new molecule.
    atomQueue = Queue.Queue()
    atom.SetProp(molAtomMapNumber, '')
    atomQueue.put(atom)

    fragToMolMap = {}
    molToFragMap = {}
    atomCount = 0
    fragment = Chem.EditableMol(Chem.Mol())
    fragment.AddAtom(atom)
    molToFragMap[atom.GetIdx()] = atomCount
    fragToMolMap[atomCount] = atom.GetIdx()
    atomCount += 1

    while not atomQueue.empty():
        current = atomQueue.get()
        for a in current.GetNeighbors():
            if (not a.HasProp('Core')) and (not a.HasProp(molAtomMapNumber)):
                a.SetProp(molAtomMapNumber, '')
                atomQueue.put(a)

                # add atom to molecule
                fragment.AddAtom(a)
                molToFragMap[a.GetIdx()] = atomCount
                fragToMolMap[atomCount] = a.GetIdx()
                atomCount += 1

                # add bond to molecule
                fragment.AddBond(molToFragMap[current.GetIdx()], molToFragMap[a.GetIdx()])

    fragment = fragment.GetMol()
    Chem.SanitizeMol(fragment)

    # Add bond information.
    for fragBond in fragment.GetBonds():
        begin = fragToMolMap[fragBond.GetBeginAtom().GetIdx()]
        end = fragToMolMap[fragBond.GetEndAtom().GetIdx()]
        molBond = mol.GetBondBetweenAtoms(begin, end)

        if molBond.GetBondType() != Chem.BondType.SINGLE:
            fragBond.SetBondType(molBond.GetBondType())
        fragBond.SetBondDir(molBond.GetBondDir())

    for a in fragment.GetAtoms():
        if a.HasProp('molAtomMapNumber') and a.GetProp('molAtomMapNumber') != molAtomMapNumber:
            a.ClearProp('molAtomMapNumber')

    # Generate smarts.
    smarts = Chem.MolToSmarts(fragment, isomericSmiles=True)
    smarts = re.sub('\[[^\[\]]*:%d\]' % int(molAtomMapNumber), '*', smarts)
    return smarts


if __name__ == '__main__':
    insmarts = '[#6,#7,#8,#16:7][C:5](=[O:6])[O:4][C:2]([#6,#7,#8,#16:1])=[O:3]>>[#6,#7,#8,#16:1][C:2]([O:4])=[O:3].[*:7][C:5]([O])=[O:6]'
    expsmarts = '[C$(*([#6,#7,#8,#16])=O):5][O$(*C([#6,#7,#8,#16])=O):4]>>[O:4].[C:5]O'
    print "Testing... "
    if compact(insmarts) == expsmarts:
        print "Success"
    else:
        print "Failure"
