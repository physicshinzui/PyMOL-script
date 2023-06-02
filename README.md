# pyMolBuilder
Written by Shinji Iida

## Tested environment
- PyMOL v2.5.4
- Python 3.8.5 

## How to use 

### Case 1. Fabricate a peptide with terminal capping.
```
python molBuilder --seq aa_sequence --mode fab
```
`aa_sequence` is written in terms of one-letter amino acid: e.g., `BAAAAAZ`.

**Terminal types**:
- 'B' : Acetyl group 
- 'Z' : N-methyl group
- 'X' : NH2 (amin group)

### Case 2. Acetyl and N-methyl capping at the terminus
```
python molBuilder --pdb pdb_file_name --mode cap
```