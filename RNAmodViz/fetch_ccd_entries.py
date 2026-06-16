import os
import requests

# Update list of modified nucleotides as set()
mods = ('4SU', 'H2U', 'QUO', 'G7M', '5MU', 'PSU', '1MG', '2MG', 'M2G', '5MC', 
        '1MA', 'OMG', 'OMC', '6IA', 'IU ', 'DG ', 'N5M', 'GTP', 'A23', 'OMU', 
        'A2M', '5BU', 'FHU', 'P5P', 'FMU', 'FOU', 'CTP', 'U37', 'MIA', 'U34', 
        '2MA', 'CCC', 'IU', 'AVC', 'AET', 'UR3', 'YYG', '7MG', 'PPU', 'T6A', 
        'UD5', 'PYO', 'SUR', '5FU', 'ONE', '6MZ', 'G3A', 'DG', 'SSU', '8AN', 
        'GDP', 'PGP', '70U', '12A', 'N79', '23G', 'GDO', 'RPC', 'A9Z', '1RN', 
        'UBD', 'G5J', 'N7X', 'T39', 'C5L', ' SC', 'PST', ' AS', '6OO', 'RFJ', 
        '6NW', '8AZ', 'GTG', 'CM0', 'ZJS')

# Update list of modified nucleotides as set()
modres = ('G7M', '5MU', 'PSU')

# Path where you would like to store your downloded CIF files
outdir = "CCD_Modified_Nucleotides"
os.makedirs(outdir, exist_ok=True)

# Download CIF MODRES files
for mod in mods:

    url = f"https://files.rcsb.org/ligands/download/{mod}.cif"

    r = requests.get(url)

    if r.status_code == 200:
        outfile = os.path.join(outdir, f"{mod}.cif")

        with open(outfile, "wb") as f:
            f.write(r.content)

        print(f"Downloaded {mod}")

    else:
        print(f"Failed {mod}")