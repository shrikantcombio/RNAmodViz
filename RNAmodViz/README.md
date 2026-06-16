# RNAmodViz

RNAmodViz is a lightweight toolkit for the retrieval, parsing, visualization, and cataloging of modified RNA nucleotides from the Protein Data Bank (PDB) Chemical Component Dictionary (CCD).

The toolkit automates:

* Downloading CCD entries for modified nucleotides
* Parsing CCD CIF files using Gemmi
* Generating publication-quality 2D structures using RDKit
* Generating 3D molecular illustrations using PyMOL
* Building multi-panel figures for reports, manuscripts, and presentations
* Creating metadata catalogs of RNA modifications

---

## Requirements

### Python packages

Install the required packages:

```bash
pip install numpy pandas gemmi rdkit pillow
```

### Additional software

For 3D molecular rendering:

```bash
PyMOL
```

---

## Workflow

### Step 1: Download CCD entries

Download modified nucleotide definitions from the PDB Chemical Component Dictionary.

```bash
python3 fetch_ccd_entries.py
```

Output:

```text
CCD_Modified_Nucleotides/
├── PSU.cif
├── 5MC.cif
├── M2G.cif
└── ...
```

---

### Step 2: Generate 2D molecular structures

Generate 2D chemical structures from CCD files using RDKit.

```bash
python3 generate_2d_structures.py
```

Features:

* Parses CCD CIF files using Gemmi
* Uses canonical SMILES whenever available
* Falls back to CCD atom/bond records when required
* Produces PNG images
* Generates a metadata catalog

Output:

```text
figures/
├── PSU.png
├── 5MC.png
├── M2G.png
└── ...
```

Metadata:

```text
modified_nucleotide_catalog.csv
```

---

### Step 3: Build figure panels

Generate publication-quality multi-panel figures.

```bash
python3 build_structure_panels.py
```

Features:

* Dynamic panel generation
* Configurable rows and columns
* Residue code annotation
* Full modification names
* High-resolution PNG export

Example output:

```text
panel_selected.png
```

---

### Step 4: Generate 3D molecular illustrations

Generate PyMOL-based 3D molecular figures.

```bash
python3 generate_3d_structures.py
```

Output:

```text
figures_PyMol/
├── PSU.png
├── 5MC.png
├── M2G.png
└── ...
```

---

## Directory Structure

```text
RNAmodViz/
├── CCD_Modified_Nucleotides/
├── figures/
├── figures_PyMol/
│
├── fetch_ccd_entries.py
├── generate_2d_structures.py
├── build_structure_panels.py
├── generate_3d_structures.py
│
├── modified_nucleotide_catalog.csv
├── RNAmodViz.log
└── README.md
```

---

## Metadata Catalog

The generated catalog contains:

| Field          | Description                          |
| -------------- | ------------------------------------ |
| id             | CCD identifier                       |
| name           | Modification name                    |
| type           | CCD chemical type                    |
| formula        | Chemical formula                     |
| formula_weight | Molecular weight                     |
| parent_comp_id | Parent nucleotide                    |
| parent_base    | Adenosine/Guanosine/Cytidine/Uridine |
| smiles         | Canonical SMILES                     |
| inchi          | InChI                                |
| inchikey       | InChIKey                             |
| atom_count     | Number of atoms                      |
| bond_count     | Number of bonds                      |
| build_method   | SMILES or CCD                        |
| png_file       | Output image                         |

---

## Example Applications

* RNA modification surveys
* Structural bioinformatics analyses
* RNA–protein interaction studies
* Preparation of publication figures
* Construction of RNA modification databases
* Educational visualization of modified nucleotides

---

## Author

Shri Kant

Computational Approaches to Protein Science Laboratory
National Centre for Biological Sciences (NCBS)
Tata Institute of Fundamental Research (TIFR)
Bengaluru, India
