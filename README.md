# AutoTerminusCap
Written by Shinji Iida

## Tested environment
- PyMOL v2.5.4
- Python 3.8.5 

## How to use 

### Case 1. Polymerize an amino-acid residues.
```
python molBuilder --seq "JAAAO" --mode polymerize
```

### Case 2. Acetyl and N-methyl capping at the terminus
```
python molBuilder --pdb pdb_file_name --mode cap
```