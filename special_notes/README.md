# special_notes — MW01 FEL data and `xpm2png.py`

This directory provides additional free-energy landscape (FEL) data for MW01 (case 3) and documents, reproduces, and corrects the XPM parsing errors in the original `xpm2png.py` implementation.

## Contents

| File or directory | Description |
| --- | --- |
| `Traj1-3/gibbs.xpm` | RMSD–Rg FEL data from the combined MW01 Trajectories 1–3. |
| `Traj4-6/gibbs.xpm` | RMSD–Rg FEL data from the combined MW01 Trajectories 4–6. |
| `Traj1-6/gibbs.xpm` | RMSD–Rg FEL data from all six combined MW01 trajectories. |
| `xpm2png(error).py` | Original erroneous implementation, retained only for reproduction and comparison. |
| `xpm2png.py` | Corrected implementation for XPM row order and coordinate handling, with point-extraction utilities removed. |

## Why the original script produces an incorrect FEL

XPM row order is opposite to the numerical coordinate direction. XPM stores pixels in image order: row 0 is the top row and corresponds to the largest y value. The `y-axis` values written by GROMACS are ordered from smallest to largest. The original script constructs `energy_matrix` in file order and pairs it directly with the ascending `y_axis`, without reversing the rows. The result is an FEL mirrored along the y-axis.

```python
energy_matrix = energy_matrix[::-1, :]
```

Omitting this operation places energy basins at incorrect Rg coordinates. This is the primary data-order error.

The corrected script also includes adaptive fixes for several potential visualization issues involving multi-line axis comments, grid dimensions, characters per pixel, energy-sign conventions, and Matplotlib compatibility. None of these issues occurred with the data analyzed in this study, and the changes were made solely to improve portability and robustness across different input files and runtime environments.

Because point extraction is not needed in this workflow, all minimum-point detection, sorting, reporting, and annotation features have been removed from the revised script. The default energy convention follows `gmx sham`: the lowest energy is `ΔG = 0`, and all other values are positive.

## Usage

Python 3, NumPy, Matplotlib, and SciPy are required. The script reads a file named `gibbs.xpm` from the current working directory by default. Enter one of the trajectory directories and run the corrected script from the parent directory:

```bash
cd Traj1-3
python ../xpm2png.py
```

Alternatively, copy `xpm2png.py` into any directory containing a `gibbs.xpm` file and run it there. With its default settings, the script:

- converts XPM coordinates from nm to Å using `k_c = 10`;
- saves the two-dimensional FEL as `fel_from_xmp.svg`;
- generates and displays a three-dimensional iceberg plot, but does not save it because its default `save_path` is `None`.

The plotting range, Gaussian smoothing, contour lines, colormap, output path, and energy convention can be changed in the user-configurable parameters section at the beginning of the script. Use `xpm2png.py` for analysis; do not use `xpm2png(error).py` for production results.
