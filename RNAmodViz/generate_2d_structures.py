# Load basic libs
import os
import sys
from pathlib import Path
import pandas as pd
import logging

# 1. Create a master logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

if logger.hasHandlers():
    logger.handlers.clear()

# 2. Create the File Handler (Logs everything from INFO and above)
file_handler = logging.FileHandler('modres_RNA.log')
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

#---------------------------------------------------#
# Gemmi for CIF mol parsing
import gemmi
# RDKit for chemical drawing stuffs
from rdkit import Chem
from rdkit.Chem import Draw
#---------------------------------------------------#

#---------------------------------------------------#
#Basic functions to downlaod and generate 2D and 3D
#structure figures using gemmi and rdkit 
#---------------------------------------------------#

# Module 1: RNA canonical nucleotides
def classify_parent_base(metadata):

    parent = metadata["parent_comp_id"]

    mapping = {
        "A": "Adenosine",
        "G": "Guanosine",
        "C": "Cytidine",
        "U": "Uridine"
    }

    return mapping.get(
        parent,
        "Unknown"
    )

# Module 2: Parse CCD CIF
def parse_ccd_cif(cif_file):
    """
    Parse CCD CIF file and return metadata dictionary.
    Parse a CCD CIF file and return:
        - residue ID
        - residue name
        - canonical smiles (if available)
        - atom table
        - bond table
    """

    doc = gemmi.cif.read_file(cif_file)
    block = doc.sole_block()

    residue_id = block.find_value("_chem_comp.id")
    residue_name = block.find_value("_chem_comp.name")

    smiles = None
    
    metadata = {

        # Basic identifiers
        "id":
            block.find_value("_chem_comp.id"),

        "name":
            block.find_value("_chem_comp.name"),

        "type":
            block.find_value("_chem_comp.type"),

        "formula":
            block.find_value("_chem_comp.formula"),

        "formula_weight":
            block.find_value("_chem_comp.formula_weight"),

        "one_letter_code":
            block.find_value(
                "_chem_comp.one_letter_code"
            ),

        "parent_comp_id":
            block.find_value(
                "_chem_comp.mon_nstd_parent_comp_id"
            ),

        "synonyms":
            block.find_value(
                "_chem_comp.pdbx_synonyms"
            ),

        # descriptors
        "smiles":
            None,

        "inchi":
            None,

        "inchikey":
            None,

        # statistics
        "atom_count":
            0,

        "bond_count":
            0
    }
        
    try:

        descriptor_table = block.find(
            "_pdbx_chem_comp_descriptor.",
            [
                "type",
                "descriptor"
            ]
        )

        for row in descriptor_table:

            descriptor_type = row[0].strip()

            descriptor = (
                row[1]
                .strip()
                .strip('"')
                .strip("'")
            )

            if (
                descriptor_type ==
                "SMILES_CANONICAL"
            ):
                metadata["smiles"] = descriptor

            elif (
                descriptor_type ==
                "InChI"
            ):
                metadata["inchi"] = descriptor

            elif (
                descriptor_type ==
                "InChIKey"
            ):
                metadata["inchikey"] = descriptor

    except Exception:
        pass

    return metadata

# Module 2: Extract Atom Table
def extract_atoms(block):

    atoms = {}

    try:
        atom_table = block.find(
            "_chem_comp_atom.",
            [
                "atom_id",
                "type_symbol"
                ]
                )
    
    except Exception:
        raise ValueError("Unable to parse atom count!")

    atom_count = 0

    for row in atom_table:

        atom_name = row[0].strip()
        atom_type = row[1].strip().capitalize()

        atoms[atom_name] = atom_type

        atom_count += 1

    return atoms, atom_count

# Module 3: Extract Bond Table
def extract_bonds(block):

    bonds = []

    try:
        
        bond_table = block.find(
            "_chem_comp_bond.",
            [
                "atom_id_1",
                "atom_id_2",
                "value_order"
                ]
                )
    
    except Exception:
        raise ValueError("Unable to parse bond records!")

    bond_count = 0

    for row in bond_table:

        bonds.append(
            (
                row[0].strip(),
                row[1].strip(),
                row[2].strip()
            )
        )

        bond_count += 1

    return bonds, bond_count

#Module 4: Build RDKit Molecule from SMILES
def mol_from_smiles(smiles):

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Unable to parse SMILES!")

    return mol

# Module 5: Build RDKit Molecule from CCD Atoms/Bonds
BOND_TYPES = {
    "SING": Chem.BondType.SINGLE,
    "DOUB": Chem.BondType.DOUBLE,
    "TRIP": Chem.BondType.TRIPLE,
    "AROM": Chem.BondType.AROMATIC,
    "DELO": Chem.BondType.AROMATIC
}

def mol_from_ccd(cif_file):

    # Get metadata
    metadata = parse_ccd_cif(cif_file)
    
    doc = gemmi.cif.read_file(cif_file)
    block = doc.sole_block()
    
    atoms, atom_count = extract_atoms(block)
    bonds, bond_count = extract_bonds(block)
    
    metadata["atom_count"] = atom_count
    metadata["bond_count"] = bond_count
    
    rwmol = Chem.RWMol()

    atom_map = {}

    for atom_name, element in atoms.items():

        idx = rwmol.AddAtom(
            Chem.Atom(element)
        )

        atom_map[atom_name] = idx

    for atom1, atom2, bond_type in bonds:

        rdkit_bond = BOND_TYPES.get(
            bond_type,
            Chem.BondType.SINGLE
        )

        rwmol.AddBond(
            atom_map[atom1],
            atom_map[atom2],
            rdkit_bond
        )

    mol = rwmol.GetMol()
    
    try:
        Chem.SanitizeMol(mol)
    
    except Exception as e:
        logger.warning(
            f"{os.path.basename(cif_file)} "
            f"sanitization failed: {e}"
            )

    return mol, metadata

# Module 6: Draw Molecule
def draw_molecule(
        mol,
        outfile,
        size=(800, 800)
    ):
    logger.info(f"Writing figure: {outfile}")
    Draw.MolToFile(
        mol,
        outfile,
        size=size
    )

def process_ccd_file(
        work_dir,
        cif_file,
        fig_dir,
        outfile
        ):
    
    #---------------------------------------------------#
    #METADATA RECORDS
    #{
    #'id': 'PSU',
    #'name': 'PSEUDOURIDINE',
    #'type': 'RNA LINKING',
    #'formula': 'C9 H12 N2 O6 P',
    #'formula_weight': '304.17',
    #'one_letter_code': '?',
    #'parent_comp_id': 'U',
    #'synonyms': '?',
    #'smiles': 'O[C@H]1[C@@H](O)...',
    #'inchi': 'InChI=1S/...',
    #'inchikey': '...',
    #'atom_count': 29,
    #'bond_count': 31
    #}
    #---------------------------------------------------#

    metadata = parse_ccd_cif(cif_file)
    
    mol, metadata = mol_from_ccd(cif_file)

    metadata["parent_base"] = (
        classify_parent_base(metadata)
    )

    mol, build_method = build_molecule(
        metadata,
        cif_file
    )

    metadata["build_method"] = (
        build_method
    )

    draw_molecule(
        mol,
        outfile
    )

    metadata["png_file"] = outfile

    return metadata

def build_molecule(metadata, cif_file):
    """
    Build RDKit molecule.

    Priority:
        1. Canonical SMILES
        2. CCD atom/bond reconstruction
    """

    # -------------------------
    # Method 1: SMILES
    # -------------------------

    smiles = metadata.get("smiles")

    if smiles:

        try:

            mol = Chem.MolFromSmiles(smiles)

            if mol is not None:

                logger.info(
                    f"{metadata['id']}: "
                    "built from SMILES"
                )

                return mol, "SMILES"

        except Exception as e:

            logger.warning(
                f"{metadata['id']}: "
                f"SMILES failed ({e})"
            )

    # -------------------------
    # Method 2: CCD
    # -------------------------

    try:

        mol = mol_from_ccd(cif_file)

        logger.info(
            f"{metadata['id']}: "
            "built from CCD"
        )

        return mol, "CCD"

    except Exception as e:

        logger.warning(
            f"{metadata['id']}: "
            f"CCD failed ({e})"
        )

    raise RuntimeError(
        f"Unable to build molecule "
        f"for {metadata['id']}"
    )

#---------------------------------------------------#
#Main function
#---------------------------------------------------#
def main():
    logger.info(f"Configuring workspace and file systems.....")
    
    # Update your working directory
    work_dir = "/home/labuser/Projects/shrikant_libs/modres_helper/"

    # Batch Processing for All Modified Nucleotides
    import glob

    ccd_dir = os.path.join(work_dir, "CCD_Modified_Nucleotides")
    fig_dir = os.path.join(work_dir, "figures")
    os.makedirs(fig_dir, exist_ok=True)
    
    logging.info(f"\n")
    logger.info(f"Work directory name:{work_dir}")
    logger.info(f"Modified nucleotides directory name: {ccd_dir}")
    logger.info(f"Output 2D figures directory name: {fig_dir} \n")

    success = 0
    failed = 0
    metadata_records = []
    metadatatocsv = os.path.join(work_dir, "modified_nucleotide_catalog.csv")

    for cif_file in glob.glob(
            os.path.join(ccd_dir, "*.cif")
            ):
        
        # Nucleotide file name 3-letter code
        residue = Path(cif_file).stem
        
        try:
            outfile = os.path.join(
                fig_dir,
                f"{residue}.png")
            
            metadata = process_ccd_file(
                work_dir,
                cif_file,
                fig_dir,
                outfile
                )
            
            # append to metadeta to a list
            metadata_records.append(metadata)

            logger.info(
                f"{residue}: OK"
            )

            success += 1

        except Exception as e:

            logger.error(
                f"{residue}: FAILED ({e})"
            )

            failed += 1
    
    metadata_df = pd.DataFrame(metadata_records)
    # Write a CSV file containing metadeta
    metadata_df.to_csv(
        metadatatocsv,
        index=False
        )

    logger.info(f"Wrote metadata to {metadatatocsv} file\n")

    logger.info(
        f"Completed processing. "
        f"Success={success}, "
        f"Failed={failed}\n"
    )

if __name__ == "__main__":
    main()
