#!/usr/bin/env python
"""
Process a flat file to produce a changed file, and a file of changes using the syngenta tautomer processing
rules
"""
from rdkit import Chem
from rdkit.Chem import AllChem

import os, sys, copy
import CBR18_als1 as CBR18
rxns, exclusions = CBR18.GetRxnsAndExclusions(CBR18.reactions, CBR18.exclusions)

def mol_with_atom_index(mol):
    '''take a molecule and label it for drawing.  return molecule'''
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx())
    return mol
def mol_with_crank_index(mol):
    '''take a molecule, set the int molAtomMapNumber property to the canonical rank, return the molecule'''
    rank = list(Chem.CanonicalRankAtoms(mol, breakTies=False))
    for i, j in enumerate(rank):
        mol.GetAtomWithIdx(i).SetAtomMapNum(j + 1)
        mol.GetAtomWithIdx(i).SetIntProp('molAtomMapNumber', j + 1)
    return mol
def OrderMatch(mol, matches):
    '''take a molecule and the list of indices of each match.  Return an array of the priority order of each match'''
    mol = mol_with_crank_index(mol)
    MatchScore = []
    for m in matches:
        Sum = 0
        for j in m:
            Sum = Sum + mol.GetAtomWithIdx(j).GetIntProp('molAtomMapNumber')
        MatchScore.append(Sum)
    return [MatchScore.index(x) for x in sorted(MatchScore, reverse=False)] #want low to high
def CheckOverlap(A,B):
    '''A is the match to the reactant, B is the list of matches to the exclusion smarts.  Return True if they overlap'''
#    print(A, ' : ', B)
    Q = set(A)
    for R in B:
        ref = set(R)
        if ref.issubset(Q) or Q.issubset(ref):
            return True
    return False
def RemoveIndices(mol):
    '''remove the atommapnum, so that they do not affect the canonical smiles; return clean mol'''
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    return mol
def FixValence(mol, msg):
    '''fix a valence error caused by a persistent H.  Take a mol and error msg, return clean mol or None'''
#    print('FV: ',msg, Chem.MolToSmiles(mol))
    pt = Chem.GetPeriodicTable()
    #msg = 'Explicit valence for atom # 5 C, 5, is greater than permitted'
    words = msg.split()
    rwmol = Chem.RWMol(mol)
    A = rwmol.GetAtomWithIdx(int(words[5]))
    ReduceBy = int(words[7][0]) - pt.GetDefaultValence(words[6][:-1])
    for  atom  in  A.GetNeighbors(): 
        if atom.GetAtomicNum() == 1 and ReduceBy > 0:
            rwmol.RemoveAtom(atom.GetIdx())
            ReduceBy = ReduceBy - 1
    if ReduceBy == 0:
        try:
            romol = rwmol.GetMol()
            Chem.SanitizeMol(romol)
            print('FV success')
            return romol
        except ValueError as e:
            print ('FixValence product rejected due to %s'%e, Chem.MolToSmiles(mol), Chem.MolToSmiles(romol))
            return None
    else:
        print('FV fail ')
        return None
 
def process_tautomer(Name, MOL, EX, verbose=False):
    '''take a molecule, and as appropriate apply rxns to it, subject to exclusions
    rxns is an ordered array of reactions, exclusions is a dict indexed by the reaction name
    return possible altered molecule, Change flag, rxns with updated counts.  Name is just for debugging
    Richard Lewis Oct 2020'''
    Changed = False # record if a mol is changed at all by this routine
    GeneratedSmiles = [] # this is used to check for cycling
    mol = RemoveIndices(copy.deepcopy(MOL)) # as the indices affect the smiles strings
    Canon = Chem.MolToSmiles(mol) # for testing uniqueness
    GeneratedSmiles.append(Canon)
    molH = Chem.AddHs(mol) # need this for some of the substructure smarts with explicit H's in them
    mol = mol_with_crank_index(mol) # add the rank property
    for i,rxn in enumerate(rxns):
        ''' Apply rules in order as above
        if match to a reagent, try further
        if matches exclusion for the reaction, continue
        for the number of matches to the reagent query
           find the order of matches using CanonicalRank of the matched atoms
           perform the reaction to create a new mol
           If after all the changes, the molecule cannot be aromatised, revert to the last valid molecule'''
        Matches = mol.GetSubstructMatches(rxn[0])
        if Matches:
            if verbose: print('rxn: ', Name, rxn[2], Canon)
            Order = OrderMatch(mol, Matches) #find the canonical order of matches
            for Ind in Order: 
                Match = Matches[Ind]
#                print('h1', rxn[2], Matches, Match)
                skip = False
                if rxn[2] in exclusions: # does a reaction have associated exclusions?
                    for j, esmarts in enumerate(exclusions[rxn[2]]): 
                      #print('excl1: ', esmarts[2], esmarts[1])
                      if molH.HasSubstructMatch(esmarts[0]):
                        ematches = molH.GetSubstructMatches(esmarts[0])
                        if verbose: print('excl2: ', rxn[2], esmarts[2], ematches)
                        skip = CheckOverlap(Match, ematches) #is the exclusion for the same part of the molecule?
                        if skip:
                            exclusions[rxn[2]][j][3] = exclusions[rxn[2]][j][3] + 1
                            EX = EX + 1
                            break #out of the exclusions loop
                if skip:
                    Matches = None
                    break # to the next reaction as the match has overlap with an exclusion
#                print('Past excl')
                ps = rxn[1].RunReactant(mol,0) #or RunReactants([mol])
                if verbose: print('h2 ', len(ps), rxn[2], Name, Canon)
                if len(ps) == 0:
                    if verbose: print('no reaction')
                    continue
                #mol = mol
                else: # a reaction has occurred
                    for j,p in enumerate(ps):
                        if p[0] is None: continue
                        if verbose:
                            print ('product ', Name, rxn[2], Chem.MolToSmiles(mol), Chem.MolToSmiles(p[0]))
                            Chem.SanitizeMol(p[0]) # hard crash
                        try:
                            Chem.SanitizeMol(p[0])
                            nmol = RemoveIndices(p[0]) 
                            smi = Chem.MolToSmiles(nmol)
                        except ValueError as e:
                            print ('product rejected due to %s'%e, Name, rxn[2], Chem.MolToSmiles(mol), Chem.MolToSmiles(p[0]))
                            msg = str(e)
                            if 'Explicit valence' in msg:
                                nmol = FixValence(p[0], msg)
                                if nmol is not None: smi = Chem.MolToSmiles(nmol)
                            else:
                                continue # try the next product
                        if nmol is None: continue
                        if smi in GeneratedSmiles: 
                            continue
                        # nmol is now a new distinct valid molecule
                        NMatches = nmol.GetSubstructMatches(rxn[0]) # check that the reaction was applied to the right match
                        if NMatches and CheckOverlap(Match, NMatches):
                            # wrong reaction site
                            continue # and get the next product
                        mol = copy.deepcopy(nmol)
                        Canon = Chem.MolToSmiles(mol) # for testing uniqueness
                        GeneratedSmiles.append(Canon)
                        molH = Chem.AddHs(mol) # need this for some of the substructure smarts with explicit H's in them
                        rxns[i][3] = rxns[i][3]+1 # increment the reaction count
                        if verbose: print('CH: ', Name, rxn[2], Canon )
                        Changed = True
                        break # out of ps loop
          #now back to the next match in the structure.
        #now back to the next reaction
    if Changed: #don't bother if it is the same molecule
            Canon = GeneratedSmiles.pop() # get last valid molecule
            mol = Chem.MolFromSmiles(Canon)
            if verbose: print ('changed structure ', Name, Canon)
    return RemoveIndices(mol), rxns, Changed, exclusions, EX

def CopyProps(newmol, mol, OrigSmi, Name):
  for Prop in mol.GetPropNames():
    newmol.SetProp(Prop, mol.GetProp(Prop))
  newmol.SetProp('Original Smiles', OrigSmi)
  newmol.SetProp('_Name', Name)
  return newmol
def main():
  import argparse
  parser = argparse.ArgumentParser(description='Tautomer model',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.set_defaults()
  parser.add_argument('-i', '--input', default=None, help='input SDF file', required=True)
  parser.add_argument('-o','--output', default=None, help='SDF output file', required=True)
  parser.add_argument('-c','--changes', default=None, help='SDF change file', required=True)
  args = parser.parse_args()
  print(args)
  rxns, exclusions = GetRxnsAndExclusions(CRB18.reactions, CRB18.exclusions)
  nchange = 0
  EX = 0
  ID = 'Small Molecule Concept Preferred Key'
  SdfIn = Chem.SDMolSupplier(args.input)
  SdfOut = Chem.SDWriter(args.output)
  SdfC = Chem.SDWriter(args.changes)
  SdfU = Chem.SDWriter('Unchanged.sdf')
  nc = 0
  for mol in SdfIn:
    if mol is None:
        continue
    Name = mol.GetProp(ID)
    OrigSmi = Chem.MolToSmiles(mol)
    newmol, rxns, Changed, exclusions, EX = process_tautomer(Name, mol, rxns, exclusions, False, EX)
    
    if Changed:
        nchange = nchange + 1
        newmol = CopyProps(newmol, mol, OrigSmi, Name)
        SdfOut.write(newmol)
        SdfC.write(newmol)
        SdfU.write(mol)
    else:
        SdfOut.write(mol)
    nc = nc + 1
    if nc % 10000 == 0: print('%d compounds processed, %d changes made, %d exclusions fired'%(nc, nchange,EX))
  print('%d compounds processed, %d changes made'%(nc, nchange))
  print('%d exclusions' %EX)
if __name__ == "__main__":
  main()
