from pymol import editor
from pymol import cmd
import pymol
import sys
import argparse 

# Written by Shinji iida
class Builder():
    def __init__(self) -> None:
        self.aa_dict = {
        'A' : 'ala',
        'C' : 'cys',
        'D' : 'asp',
        'E' : 'glu',
        'F' : 'phe',
        'G' : 'gly',
        'H' : 'his',
        'I' : 'ile',
        'K' : 'lys',
        'L' : 'leu',
        'M' : 'met',
        'N' : 'asn',
        'P' : 'pro',
        'Q' : 'gln',
        'R' : 'arg',
        'S' : 'ser',
        'T' : 'thr',
        'V' : 'val',
        'W' : 'trp',
        'Y' : 'tyr',
        'J' : 'ace', # Acetyl group 
        'O' : 'nme', # N-methyl group
        'U' : 'nhh'  # NH2 (amin group)
        }

        self.check_pymol_version()

    @staticmethod
    def check_pymol_version() -> None:
        current_version = pymol.cmd.get_version()[0]
        tested_version = "2.5.4"
        print(f"Current version: {current_version}")
        if current_version != tested_version:
            print(f"NOTE: The current version ({current_version}) is not the tested version ({tested_version})")

    def polymerize(self, sequence, object_name, first_residue="1", outpdb='peptide.pdb') -> None:
        """
        Args: 
            object_name: created object name
            sequence: 1-letter string of amino-acid residues you want to create
            first_residue: integer that specify where is the first residue in the string you give.
        """
        if len(sequence):
            code = sequence[0]

        # Create a fragment of code (3 letter) named `object_name`
        cmd.fragment(self.aa_dict[code], object_name)

        # Alter the residue number of the first residue: `resi=1`
        cmd.alter(object_name,f"resi={first_residue}")

        # Pick C atom of `object_name` for editing
        cmd.edit(object_name+" and name C")

        # loop from the 2nd letter
        # i.e., XXXXXXX
        #        ^from here
        for code in sequence[1:]:
            # Attach subsequent amino-acid residues
            # e.g. 
            #     sequence=AKG; sequence[1:] is "KG"
            #     In this loop, code would be K first, and then G. 
            #     The letter is fed into the dictionary that stores 3 letter amino-acid names
            #     , giving the corresponding 3-letter name.
            editor.attach_amino_acid("pk1", self.aa_dict[code]) 
                                            #^ 3-letter name
            # NOTE:
            # After adding an amino-acid residues, 
            # pk1 changes s.t. it specify C atom in the added residue.

        cmd.save(f"{outpdb}")

    def cap(self, in_pdb_file, nter_cap="ace", cter_cap="nme", outpdb='capped') -> None:
        """
        Args: 
            pdb_file: PDB file name (String)
            nter_cap: "ace"
            cter_cap: "nme" or "nhh"
            outpdb: 'capped'
        """
        # NOTE: Single chain is ASSUMED in this function. 
        
        object_name = "pdb_obj"
        cmd.load(in_pdb_file, object_name)
        
        def get_first_last_residue_info(object_name):
            atoms = cmd.get_model(object_name)
            residues = [(atm.resn, atm.resi) for atm in atoms.atom]
            print(residues[0], residues[-1])
            first_resn, first_resi = residues[0]
            last_resn, last_resi = residues[-1]

            if first_resn.upper() == "ACE": 
                sys.exit("The N-terminus has already been capped by ACE.")
            if last_resn.upper() in ["NME", "NHH", "NHE", "NH2"]: 
                sys.exit("The N-terminus has already been capped by NME.")
            return (first_resn, first_resi), (last_resn, last_resi)
        
        # Check if the N-terminus has already been capped
        get_first_last_residue_info(object_name)

        # Clean the data before adding caps.
        # Note that this removes the caps in the object.
        cmd.remove("not polymer") 

        # Their info is used to add caps to them in subsequent lines.  
        (first_resn, first_resi), (last_resn, last_resi) = get_first_last_residue_info(object_name)

        # TODO: I wanna put this condition: If the selection below returns zero atom, then raise an error. 
        # nter_cap is conneted to the object pk1 (Nitrogen)         
        cmd.edit(object_name+f" and resi {first_resi} and name N") # pk1 is made here 
        editor.attach_amino_acid("pk1", nter_cap)

        # cter_cap is conneted to the object pk1 (Carbon)         
        cmd.edit(object_name+f" and resi {last_resi} and name C") # pk1 is updated here 
        editor.attach_amino_acid("pk1", cter_cap)

        cmd.save(f"{outpdb}.pdb")

def main():
    # TODO: The arguments are in a mess at the moment, so I need to consider what arguments should come.
    p = argparse.ArgumentParser()
    grp = p.add_mutually_exclusive_group(required=True)
    grp.add_argument("--seq")
    grp.add_argument("--pdb")
    #p.add_argument("-o", "--outpdb")
    p.add_argument("--mode", required=True, choices=["polymerize","cap"])
    args = p.parse_args()       

    builder = Builder()
    
    if args.mode == "polymerize":
        builder.polymerize(sequence=args.seq, object_name="poly", outpdb="peptide.pdb")
    
    elif args.mode == "cap":
        builder.cap(in_pdb_file=args.pdb, outpdb="capped")

main()
