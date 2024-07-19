# pyMolBuilder
Written by Shinji Iida

## Tested environment
- PyMOL v2.5.4
- Python 3.8.5 

> [!NOTE]
> You may need to add the path of pymol module to PYTHONPATH so that you can execute `import pymol`.
> e.g. `export PYTHONPATH=/opt/homebrew/Cellar/pymol/3.0.0/libexec/lib/python3.12/site-packages:$PYTHONPATH`

## How to use 

### Case 1. Fabricate a peptide with terminal capping.
```
python molBuilder.py --seq aa_sequence --mode fab
```
`aa_sequence` is written in terms of one-letter amino acid: e.g., `BAAAAAZ`.

**Terminal types**:
- 'B' : Acetyl group 
- 'Z' : N-methyl group
- 'X' : NH2 (amin group)

### Case 2. Acetyl and N-methyl capping at the terminus
```
python molBuilder.py --pdb pdb_file_name --mode cap
```
>[!NOTE]
> - If the C atom of the backbone at the C-terminus does not exist in the iput PDB file, this program does not add NME to the end. 
> - Likewise, the N atom of the backbone at the N-terminus does not exist, it does not add ACE.
> - This process will not raise an error.
