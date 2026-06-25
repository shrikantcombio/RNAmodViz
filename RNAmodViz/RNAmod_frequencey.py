# Load basic libs
import os
import sys
import pandas as pd
import logging
# Gemmi for CIF mol parsing
import gemmi
from collections import Counter
from pathlib import Path
import glob

# 1. Create a master logger
logger = logging.getLogger()
logger.setLevel(logging.DEBUG) # Capture everything at the root level

# 2. Create the File Handler (Logs everything from INFO and above)
file_handler = logging.FileHandler('RNAmod_frequency.log')
file_handler.setLevel(logging.INFO)
file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(file_formatter)

# 3. Create the Console Handler (Logs ONLY WARNING and above to stdout)
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.WARNING)
console_formatter = logging.Formatter('%(levelname)s: %(message)s')
console_handler.setFormatter(console_formatter)

# 4. Add both handlers to the logger
logger.addHandler(file_handler)
logger.addHandler(console_handler)

def get_rna_modifications(cif_file):
    # Read the CIF block
    doc = gemmi.cif.read_file(cif_file)
    block = doc.sole_block()
    
    # Generate Gemmi's high-level Structure object from the block
    st = gemmi.make_structure_from_block(block)
    
    # Standard Canonical Nucleic acids
    canonical_nucleic = {"A", "G", "C", "U", "DA", "DG", "DC", "DT"}
    
    # Standard Amino Acids (20 standard + common variants/hetero-atoms)
    amino_acids = {
        "ARG", "LYS", "ASP", "GLU", "ASN", "GLN", "HIS", "SER", "THR", "TYR", 
        "TRP", "PHE", "CYS", "MET", "LEU", "ILE", "PRO", "VAL", "ALA", "GLY",
        "MSE" # Selenomethionine (common in PDB structures)
    }
    
    # Solvents and crystallization agents to ignore
    common_solvents = {"HOH", "DOD", "WAT", "GOL", "EDO", "SO4", "PO4", "CL", "MG"}

    modifications = set()  # Use a set to prevent duplicate names

    # Loop through models, chains, and residues elegantly
    for model in st:
        for chain in model:
            for res in chain:
                resname = res.name.strip()
                
                # Filter: Must NOT be standard nucleic, NOT an amino acid, and NOT solvent
                if (resname not in canonical_nucleic and 
                    resname not in amino_acids and 
                    resname not in common_solvents):
                    
                    # Ensure it's part of an RNA/DNA polymer or flagged as HETERO
                    modifications.add(resname)
                    
    return list(modifications)

# Update your working directory
work_dir = "/home/labuser/Projects/rcsb_NAdb/structures/"
data_subdir = "CCD_Modified_Nucleotides"

def main():

    freq = Counter()

    for cif_file in glob.glob(
            f"{work_dir}/*.cif"):
        
        mods = get_rna_modifications(cif_file)
        
        freq.update(mods)

    freq_df = pd.DataFrame(
        freq.items(),
        columns=[
            "id",
            "frequency"
        ]
    )

    freq_df.to_csv("RNAmod_frequency.csv", index=False)
    logging.info(f"Wrote RNA modified freuency to a CSV file.")

if __name__ == "__main__":
    main()