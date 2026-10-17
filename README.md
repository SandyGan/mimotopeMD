# mimotopeMD

Code and data accompanying: **Z. Sha et al., "Mimotope-Guided Molecular Dynamics Enables Precise Functional Epitope Mapping," *J. Chem. Inf. Model.* 2026, 66 (17), 11259–11274.**
DOI: [10.1021/acs.jcim.6c01562](https://doi.org/10.1021/acs.jcim.6c01562)

This repository contains GROMACS molecular dynamics (MD) input files and analysis data for four antibody–mimotope systems. The data include RMSD, radius of gyration (Rg), free-energy landscape (FEL), and peptide–antibody contact-frequency results derived from three independent 1-μs trajectories for each system.

The repository also includes a frame-normalized contact-frequency matrix for the ChiLob 7/4–CD40 crystal complex (PDB: 6FAX).

## Directory index

| Directory or file | System | Contents |
| --- | --- | --- |
| [`case1/`](case1/README.md) | ChiLob 7/4 | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case2/`](case2/README.md) | aεFab | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case3/`](case3/README.md) | MW01 | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case4/`](case4/README.md) | m31A7 | GROMACS input files and RMSD, Rg, FEL, and contact-frequency data. |
| [`case1-6fax/`](case1-6fax/README.md) | ChiLob 7/4–CD40 (PDB: 6FAX) | GROMACS input files and a normalized antibody-contact frequency matrix for the CD40 `53-EFTE-56` segment. |
| [`merge_xtc.py`](merge_xtc.py) | All systems | Trajectory centering, periodic-boundary correction, and XTC concatenation. |

See the README in each directory for file descriptions, simulation parameters, data details, and usage instructions. From any case directory, run `./RUN` to start the complete GROMACS workflow from scratch.
