#!/usr/bin/env python
doc = """
Process a flat file to produce a changed file, and a file of changes using the syngenta tautomer processing rules
rules
"""

import os, sys, argparse

from Fix_Tautomer import process_tautomer, CopyProps
from rdkit import Chem

def FileHandle(s):
    stem =  s[-3:]
    assert stem in ['smi', 'sdf', 'txt']
    if stem == 'smi' or stem == 'txt':
      FH = Chem.SmilesWriter(s, includeHeader=False)
      #FH.SetProps(['_Name'])
    else:
      FH = Chem.SDWriter(s)
    return FH
def main():
  parser = argparse.ArgumentParser(description=doc,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
  parser.set_defaults(same='Unchanged.sdf', before=None)
  parser.add_argument('-i', '--input', dest='input', default=None, help='input SDF (.sdf) or smiles (.smi) file or .txt file (eg ChEMBL)', required=True)
  parser.add_argument('-o','--output', dest='output', default=None, help='output file', required=True)
  parser.add_argument('-c','--changes', dest='changes',default=None, help='change file', required=True)
  parser.add_argument('-b','--before', dest='before',default=None, help='before change file')
  parser.add_argument('-u','--same', dest='same', default=None, help='unchanged file')
  args = parser.parse_args()
  try:
    assert os.path.exists(args.input)
    stem =  args.input[-3:]
    assert stem in ['smi', 'sdf', '.gz', 'txt']
    if stem == 'smi':
      Suppl = Chem.SmilesMolSupplier(args.input)#, delimiter='\t')#, titleLine=False)
    elif stem == 'sdf':
      Suppl = Chem.SDMolSupplier(args.input)
    elif stem == 'txt':
      Suppl = open(args.input, 'r')
      line = next(Suppl)
      if 'chembl_id' not in line or 'canonical_smiles' not in line:
        print(f'{args.input} not in Chembl format: header is {line}')
        sys.exit(2)
    else:
      import gzip
      Inf = gzip.open(args.input)
      Suppl = Chem.ForwardSDMolSupplier(Inf)
  except:
    parser.print_help()
    sys.exit(1)
  nchange = 0
  EX = 0
  SdfOut = FileHandle(args.output)
  SdfC = FileHandle(args.changes)
  SdfU = None
  if args.same != None:
    SdfU = FileHandle(args.same)
  SdfB = None
  if args.before != None:
    SdfB = FileHandle(args.before)
  nc = 0
  for mol in Suppl:
    if stem == 'txt':
        words = mol.split()
        mol = Chem.MolFromSmiles(words[1])
        if mol is None:
           print(f'failed to process {words[0]} {words[1]}')
           continue
        mol.SetProp('_Name', words[0])
    if mol is None:
        continue
    Name = mol.GetProp('_Name')
    OrigSmi = Chem.MolToSmiles(mol)
    #print('O:',OrigSmi)
    newmol, rxns, Changed, exclusions, EX = process_tautomer(Name, mol, EX, verbose=False)
    #print('O1:', Chem.MolToSmiles(mol))
    if Changed:
        #print('N:',Chem.MolToSmiles(newmol))
        nchange = nchange + 1
        newmol = CopyProps(newmol, mol, OrigSmi, Name)
        SdfOut.write(newmol)
        SdfC.write(newmol)
        if SdfB is not None: SdfB.write(mol)
    else:
        SdfOut.write(mol)
        if SdfU is not None: SdfU.write(mol)
    nc = nc + 1
#    if nchange == 1: sys.exit(1)
    if nc % 10000 == 0: print('%d compounds processed, %d changes made, %d exclusions fired'%(nc, nchange,EX))
  print('%d compounds processed, %d changes made'%(nc, nchange))
  #print('%d exclusions' %EX)
if __name__ == "__main__":
  main()
