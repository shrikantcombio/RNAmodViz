import os
from pymol import cmd

# Workdir
work_dir = "/home/labuser/Projects/shrikant_libs/modres_helper/"
data_subdir = "CCD_Modified_Nucleotides"

# Output fig
fig_dir = os.path.join(work_dir, "figures_PyMol")
os.makedirs(fig_dir, exist_ok=True)

mods = ["5FU", "5MC", "6IA","2MG"]

for mod in mods:
    ccd_cif_file = os.path.join(work_dir, data_subdir, f"{mod}.cif")
    cmd.reinitialize()
    cmd.load(ccd_cif_file)

    cmd.select("mod", f"resn {mod}")

    if cmd.count_atoms("mod") == 0:
        continue

    cmd.hide("everything")
    cmd.show("sticks", "mod")

    cmd.orient("mod")

    cmd.bg_color("white")

    cmd.ray(1200,1200)

    cmd.png(f"{fig_dir}/{mod}.png", dpi=600)