import os
import math
from PIL import Image, ImageDraw, ImageFont


# Update the dictionary of modified nucleotides
MODIFICATION_NAMES = {
    "PSU": "Pseudouridine",
    "5MC": "5-Methylcytidine",
    "M2G": "N2-Methylguanosine",
    "OMG": "O2'-Methylguanosine",
    "5FU": "5-Fluorouridine",
    "6IA": "6-Iodoadenosine",
    "2MG": "N2-Methylguanosine",
}


def create_modification_panel(
        fig_dir,
        output_file,
        residue_list=None,
        ncols=4,
        padding=40,
        title_height=40,
        subtitle_height=40,
        bg_color="white",
        dpi=600
):
    """
    Create a publication-quality panel of modified nucleotides.

    Parameters
    ----------
    fig_dir : str
        Directory containing PNG files.

    residue_list : list or None
        Example:
            ["PSU","5MC","M2G"]

        If None, all PNG files are used.

    ncols : int
        Number of columns.

    output_file : str
        Output panel PNG.
    """

    # --------------------------------------------------
    # Determine files
    # --------------------------------------------------

    if residue_list is None:

        files = sorted(
            f for f in os.listdir(fig_dir)
            if f.endswith(".png")
        )

    else:

        files = []

        for residue in residue_list:

            png_file = f"{residue}.png"

            if os.path.exists(
                os.path.join(fig_dir, png_file)
            ):
                files.append(png_file)

            else:
                print(
                    f"WARNING: {png_file} not found"
                )

    if len(files) == 0:
        raise ValueError("No PNG files found.")

    # --------------------------------------------------
    # Load images
    # --------------------------------------------------

    images = [
        Image.open(os.path.join(fig_dir, f))
        for f in files
    ]

    img_w, img_h = images[0].size

    nrows = math.ceil(len(images) / ncols)

    cell_h = (
        img_h +
        title_height +
        subtitle_height
    )

    panel_w = (
        ncols * img_w +
        (ncols + 1) * padding
    )

    panel_h = (
        nrows * cell_h +
        (nrows + 1) * padding
    )

    panel = Image.new(
        "RGB",
        (panel_w, panel_h),
        bg_color
    )

    draw = ImageDraw.Draw(panel)

    # --------------------------------------------------
    # Fonts
    # --------------------------------------------------

    try:

        title_font = ImageFont.truetype(
            "DejaVuSans-Bold.ttf",
            40
        )

        subtitle_font = ImageFont.truetype(
            "DejaVuSans.ttf",
            30
        )

    except:

        title_font = ImageFont.load_default()
        subtitle_font = ImageFont.load_default()

    # --------------------------------------------------
    # Draw panel
    # --------------------------------------------------

    for idx, (img, fname) in enumerate(
            zip(images, files)
    ):

        residue = os.path.splitext(fname)[0]

        full_name = MODIFICATION_NAMES.get(
            residue,
            "Unknown Modification"
        )

        row = idx // ncols
        col = idx % ncols

        x = padding + col * (img_w + padding)

        y = padding + row * (cell_h + padding)

        panel.paste(img, (x, y))

        # -----------------------------
        # Residue code
        # -----------------------------

        title_bbox = draw.textbbox(
            (0, 0),
            residue,
            font=title_font
        )

        title_w = (
            title_bbox[2] -
            title_bbox[0]
        )

        draw.text(
            (
                x + (img_w - title_w) / 2,
                y + img_h + 5
            ),
            residue,
            fill="black",
            font=title_font
        )

        # -----------------------------
        # Full modification name
        # -----------------------------

        subtitle_bbox = draw.textbbox(
            (0, 0),
            full_name,
            font=subtitle_font
        )

        subtitle_w = (
            subtitle_bbox[2] -
            subtitle_bbox[0]
        )

        draw.text(
            (
                x + (img_w - subtitle_w) / 2,
                y + img_h + title_height
            ),
            full_name,
            fill="black",
            font=subtitle_font
        )

    panel.save(
        output_file,
        dpi=(dpi, dpi)
    )

    print(
        f"Saved panel: {output_file}"
    )

# Main function
work_dir = "/home/labuser/Projects/shrikant_libs/modres_helper/"
fig_dir = os.path.join(work_dir, "figures")
#fig_dir = os.path.join(work_dir, "figures_PyMol")

# Create panel for your residues of intrests
create_modification_panel(
    fig_dir=fig_dir,
    output_file=os.path.join(work_dir, "panel_selected.png"),
    #output_file=os.path.join(work_dir, "panel_selected_PyMol.png"),
    residue_list=[
        "5FU",
        "5MC",
        "6IA",
        "2MG"
    ],
    ncols=2
)

# If want to build a panel with all the PNG in fig_dir
"""
create_modification_panel(
    fig_dir=fig_dir,
    output_file="panel_all.png",
    ncols=4
)
"""