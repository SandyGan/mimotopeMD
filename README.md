# mimotopeMD

This repository contains GROMACS molecular dynamics (MD) input files and analysis data for four antibody–mimotope systems. The data include RMSD, radius of gyration (Rg), free-energy landscape (FEL), and peptide–antibody contact-frequency results derived from three independent replicate trajectories.

The repository also includes frame-normalized contact-frequency analysis and a heatmap for the ChiLob 7/4–CD40 crystal complex (PDB: 6FAX), together with a description and correction of the XPM row-order error affecting the rendering of the free-energy landscapes, together with robustness improvements to `xpm2png.py`.

## Directory index

| Directory or file | System | Contents |
| --- | --- | --- |
| [`case1/`](case1/README.md) | ChiLob 7/4 | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case2/`](case2/README.md) | aεFab | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case3/`](case3/README.md) | MW01 | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case4/`](case4/README.md) | m31A7 | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case1-6fax/`](case1-6fax/README.md) | ChiLob 7/4–CD40 (PDB: 6FAX) | GROMACS input files, normalized antibody-contact frequencies, and a heatmap for the CD40 `53-EFTE-56` segment. |
| [`special_notes/`](special_notes/README.md) | MW01 / XPM script notes | Additional MW01 FEL data, the erroneous script, the corrected script, and a detailed explanation. |
| [`merge_xtc.py`](merge_xtc.py) | All systems | Trajectory centering, periodic-boundary correction, and XTC concatenation. |

See the README in each directory for file descriptions, simulation parameters, data details, and usage instructions. From any case directory, run `./RUN` to start the complete GROMACS workflow from scratch.
