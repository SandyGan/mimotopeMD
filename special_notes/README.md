# special_notes — MW01 FEL data and XPM conversion scripts

This directory provides companion free-energy landscape (FEL) data for MW01 (Case 3) and documents, reproduces, and corrects the XPM row-order error in the original conversion script. The Traj1–3 and Traj4–6 XPM files mirror the corresponding files in [`../case3/data/`](../case3/data/); the combined Traj1–6 XPM file is provided only here.

## Contents

| File or directory | Description |
| --- | --- |
| `Traj1-3/gibbs.xpm` | RMSD–Rg FEL data from the combined MW01 Trajectories 1–3; identical to [`../case3/data/gibbs1-3.xpm`](../case3/data/gibbs1-3.xpm). |
| `Traj4-6/gibbs.xpm` | RMSD–Rg FEL data from the combined MW01 Trajectories 4–6; identical to [`../case3/data/gibbs4-6.xpm`](../case3/data/gibbs4-6.xpm). |
| `Traj1-6/gibbs.xpm` | RMSD–Rg FEL data from all six combined MW01 trajectories. |
| `xpm2png(wrong).py` | Original erroneous implementation, retained only for reproduction and comparison. |
| `xpm2png(correct).py` | Corrected implementation for XPM row order and energy-sign handling; it also retains the FEL point-sorting and local-minimum utilities. |

## Why the original script produces an incorrect FEL

XPM row order is opposite to the numerical coordinate direction. XPM stores pixels in image order: row 0 is the top row and corresponds to the largest y value. The `y-axis` values written by GROMACS are ordered from smallest to largest. The original script constructs `energy_matrix` in file order and pairs it directly with the ascending `y_axis`, without reversing the rows. The result is an FEL mirrored along the y-axis.

```python
energy_matrix = energy_matrix[::-1, :]
```

Omitting this operation places energy basins at incorrect Rg coordinates. This is the primary data-order error.

The default energy convention follows `gmx sham`: the lowest energy is `ΔG = 0`, and all other values are positive.

## Usage

Python 3, NumPy, Matplotlib, and SciPy are required. The script reads a file named `gibbs.xpm` from the current working directory by default. Enter one of the trajectory directories and run the corrected script from the parent directory:

```bash
cd Traj1-3
python '../xpm2png(correct).py'
```

Alternatively, copy `xpm2png(correct).py` into any directory containing a `gibbs.xpm` file and run it there. With its default settings, the script:

- converts XPM coordinates from nm to Å using `k_c = 10`;
- saves the two-dimensional FEL as `fel_from_xmp.svg`;
- reports the lowest-energy FEL points in the terminal.

The plotting range, Gaussian smoothing, contour lines, colormap, output path, and energy convention can be changed in the user-configurable parameters section at the beginning of the script. Use `xpm2png(correct).py` for analysis; do not use `xpm2png(wrong).py` for production results.
