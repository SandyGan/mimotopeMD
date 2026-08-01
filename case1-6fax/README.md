# Case 1–6FAX — ChiLob 7/4–CD40 complex

This directory contains the input files for molecular dynamics simulations of the 6FAX crystal structure (ChiLob 7/4–CD40 complex), together with frame-normalized antibody-contact data for the CD40 `53-EFTE-56` segment. The CHARMM-GUI-generated inputs cover the minimization, equilibration, and production stages. This simulation is not part of the epitope-mapping workflow. It was performed as a reference for Case 1, in order to compare the interaction interface formed by the PWVP mimotope with that formed by the corresponding 53-EFTE-56 segment in the native ChiLob 7/4–CD40 complex.

## Running the simulation

Run `./RUN` directly to start GROMACS from scratch:

```bash
./RUN
```

`RUN` performs energy minimization, equilibration, and five consecutive production segments. In `step5_production.mdp`, `dt = 0.002 ps` and each segment contains `100,000,000` steps. Each segment therefore represents 200 ns, for a total production time of 1,000 ns (1 μs).

The current script uses GPU 0 and 12 CPU threads. Adjust these settings for the available hardware and ensure that `csh` and `gmx` are available before starting the simulation.

## Contents

| File or directory | Description |
| --- | --- |
| `RUN` | Complete C shell workflow from minimization through five production segments. |
| `step3_input.gro/.pdb/.psf` | Initial structure of the 6FAX complex. |
| `step4.0_minimization.mdp` | Energy-minimization parameters. |
| `step4.1_equilibration.mdp` | Equilibration parameters. |
| `step5_production.mdp` | Production parameters for one 200 ns segment. |
| `index.ndx` | GROMACS atom-group index. |
| `topol.top`, `toppar/` | System topology and force-field parameters. |
| `data/efte_contact_map.xlsx` | Normalized contact-frequency matrix between CD40 residues `53E`, `54F`, `55T`, and `56E` and ChiLob 7/4 antibody residues. |

Rows in `efte_contact_map.xlsx` represent the four CD40 EFTE residues, and columns represent antibody residues. Each cell contains a contact frequency obtained by normalizing the original contact count by the total number of analyzed frames. Values range from 0 to 1 and are stored with six decimal places.
