# Tautomer assignment using Smarts

Tautomerisation occurs when a proton migrates from one site to another. The
classic example is between 2-hydroxypyridine and pyridone. Tautomerism is an
equilibrium between the different states, and depends on factors such as solvent,
temperature and overall structure. Ideally, this would be represented as a
Boltzmann population, but for 2D databases this is not generally possible. For
model-building exercises using 2D descriptors, business rules are needed to
specify the most likely state. Baker and coworkers at Syngenta published a series
of 473 rules<sup>1</sup> based on QM calculations, to determine the most stable state of two
tautomeric forms. The rules were published as RXN files in v2000 and v3000
formats for use in Pipeline Pilot. This repo describes the conversion and curation
of these rules into Reaction SMARTS for use in RDKit.

In model building, the form of the 2D structure used to generate descriptors should be
consistent, even if it is not ‘right’. These 473 tautomer business rules at least try to
enforce that, and are more sophisticated with respect to exceptions. There are
other commercial programs that can create tautomer populations (Unicon from
Rarey group, U. Hamburg, Tauthor from Molecular Discovery), but they handle only
the more obvious cases. There is also some existing code
(https://rdkit.blogspot.com/2020/01/trying-out-new-tautomer.html). This approach
offers that possibility to identify those compounds that can exist as tautomers, using
the more sophisticated set of rules from Syngenta. If desired, the program could be
used to create input to a more rigorous scoring approach such as QM, using fast
methods like xtb with solvent corrections.

## Installation
The code requires python with rdkit installed.  The code has been tested under python 3.13, rdkit 2025.03.2 but will probably work with earlier versions.

## Usage
python tautomer.py [-h] -i INPUT -o OUTPUT -c CHANGES [-b BEFORE] [-u SAME]

Process a flat file to produce a changed file, and a file of changes using the
syngenta tautomer processing rules rules

options:
  -h, --help            show this help message and exit
  -i, --input INPUT     input SDF (.sdf) or smiles (.smi) file or .txt file
                        (eg ChEMBL) (default: None)
  -o, --output OUTPUT   output file (default: None)
  -c, --changes CHANGES
                        change file (default: None)
  -b, --before BEFORE   before change file (default: None)
  -u, --same SAME       unchanged file (default: None)
## Test
The code was run on the file chembl_35_chemreps.txt downloaded from the chembl website.
2474576 compounds were processed, 61352 changes were made.  Errors were seen in the processing of entries with [P, As]hexafluorides.  An earlier version (chembl_30) was tested and Chris Baker kindly compared the results to his set-up using pipeline pilot.  The results were the same.

## Conclusions
The Syngenta rules for assigning tautomer structures have been implemented in
RDKit and are available for general use.

## Hat tip to Paolo Tosco for much help.

## References
1. Tautomer Standardization in Chemical Databases: Deriving Business Rules from Quantum Chemistry Christopher M. Baker et al. J. Chem. Inf. Model. 2020, 60, 8, 3781–3791