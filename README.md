# mimotopeMD

This repository contains GROMACS molecular dynamics (MD) input files and analysis data for four antibody–mimotope systems. The data include RMSD, radius of gyration (Rg), free-energy landscape (FEL), and peptide–antibody contact-frequency results. Cases 1, 2, and 4 each contain three independent 1-μs trajectories. MW01 (Case 3) contains six 1-μs trajectories, distributed as two three-trajectory analysis sets: Traj1–3 and Traj4–6.

The repository also includes a frame-normalized contact-frequency matrix for the ChiLob 7/4–CD40 crystal complex (PDB: 6FAX), together with a description of the XPM row-order error affecting FEL rendering and corrected and reproduction versions of the XPM conversion script.

## Directory index

| Directory or file | System | Contents |
| --- | --- | --- |
| [`case1/`](case1/README.md) | ChiLob 7/4 | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case2/`](case2/README.md) | aεFab | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case3/`](case3/README.md) | MW01 | GROMACS input files and separate RMSD, Rg, FEL, and contact-frequency data sets for Traj1–3 and Traj4–6. |
| [`case4/`](case4/README.md) | m31A7 | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case1-6fax/`](case1-6fax/README.md) | ChiLob 7/4–CD40 (PDB: 6FAX) | GROMACS input files and a normalized antibody-contact frequency matrix for the CD40 `53-EFTE-56` segment. |
| [`special_notes/`](special_notes/README.md) | MW01 / XPM script notes | MW01 FEL reference copies, the all-six-trajectory FEL, the erroneous script, the corrected script, and a detailed explanation. |
| [`merge_xtc.py`](merge_xtc.py) | All systems | Trajectory centering, periodic-boundary correction, and XTC concatenation. |

See the README in each directory for file descriptions, simulation parameters, data details, and usage instructions. From any case directory, run `./RUN` to start the complete GROMACS workflow from scratch.
