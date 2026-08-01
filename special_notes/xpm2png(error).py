#!/usr/bin/env python3
"""
Free Energy Landscape (FEL) Generator
Based on XPM data format from GROMACS
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap, LightSource, Normalize
from scipy.ndimage import gaussian_filter
import re
# from mpl_toolkits.mplot3d import Axes3D  # Enables 3D support; the name is not used directly.

###################### User-configurable parameters #############################################

# Path to the XPM file.
xpm_path = 'gibbs.xpm'

# Output path for the 2D figure; use None to disable saving.
out_fig_path = './fel_from_xmp.svg'

# Axis limits: ax.axis([xmin, xmax, ymin, ymax]), or None for automatic limits.
limit_axis = None
# limit_axis = [1.10, 1.80, 11.52, 11.85]

# Show axis grid lines.
need_grids = False

# Gaussian smoothing; 0 disables smoothing.
gauss = 0.0

# Unit-conversion factor; use 10 for nm to angstrom or 0.1 for angstrom to nm.
k_c = 10

# Contour-line settings; use None to disable contour lines.
# contour_line = None
contour_line = {
    'level': 15,
    'alpha': 0.9,
    'label': True, # True labels the contour lines with their values.
}

# Number of minimum-energy points to report.
top_values = 10

# Colormap.
# whatcmap = 'custom'
whatcmap = 'GnBu_r'

##################################################################################

def make_cmap(wcmap='viridis'):
    if wcmap == 'custom1':
    # Unevenly spaced color anchors matching the 2D palette.
        stops  = [0.0, 0.11, 0.22, 0.33, 0.44, 0.55, 0.66, 0.77, 0.88, 1.0]
        colors = ["#FCFC00", "#FCCF30", "#DCBC28", "#81CC59", "#31C59F",
                "#08B4D1", "#2895EC", "#3F6EFF", "#4646EB", "#3E27AC"][::-1]
        return LinearSegmentedColormap.from_list("my_cmap", list(zip(stops, colors)), N=256)
    elif wcmap == 'custom2':
        stops  = [0.0, 0.14, 0.28, 0.42, 0.56, 0.70, 0.84, 1.0]
        colors = ['#BFDFD2','#51999F','#4198AC','#7BC0CD',
                  '#DBCB92','#ECB66C','#EA9E58','#ED8D5A']
        return LinearSegmentedColormap.from_list("my_cmap", list(zip(stops, colors)), N=256)
    elif wcmap == 'custom3':
        stops  = [0.0, 0.14, 0.28, 0.42, 0.56, 0.70, 0.84, 1.0]
        colors = ['#B6B3D6','#CFCCE3','#D5D3DE','#D5D1D1',
                  '#F6DFD6','#F8B2A2','#F1837A','#E9687A']
        return LinearSegmentedColormap.from_list("my_cmap", list(zip(stops, colors)), N=256)
    elif wcmap == 'custom4':
        stops  = [0.0, 0.20, 0.40, 0.60, 0.80, 1.0]
        colors = ['#61AACF','#98CADD','#EAEFF6',
                  '#F9EFEF','#E9C6C6','#DA9599']
        return LinearSegmentedColormap.from_list("my_cmap", list(zip(stops, colors)), N=256)
    elif wcmap == 'custom5':
        stops  = [0.0, 0.14, 0.28, 0.42, 0.56, 0.70, 0.84, 1.0]
        colors = ['#3C9BC9','#65BDBA','#B0D6A9','#FEE199',
                  '#FDCD94','#FAA26F','#F97F5F','#FC757B']
        return LinearSegmentedColormap.from_list("my_cmap", list(zip(stops, colors)), N=256)
    else:
        return cm.get_cmap(wcmap)


def parse_xpm_file(filename):
    """
    Parse an XPM file and extract its free-energy data.
    """
    with open(filename, 'r') as f:
        content = f.read()
    
    # Extract metadata.
    title_match = re.search(r'/\*\s*title:\s*"(.*?)"\s*\*/', content)
    legend_match = re.search(r'/\*\s*legend:\s*"(.*?)"\s*\*/', content)
    xlabel_match = re.search(r'/\*\s*x-label:\s*"(.*?)"\s*\*/', content)
    ylabel_match = re.search(r'/\*\s*y-label:\s*"(.*?)"\s*\*/', content)
    x_axis_match = re.search(r'/\*\s*x-axis:\s*(.*?)\s*\*/', content)
    y_axis_match = re.search(r'/\*\s*y-axis:\s*(.*?)\s*\*/', content)
    
    title = title_match.group(1)
    legend = legend_match.group(1)
    xlabel = xlabel_match.group(1)
    ylabel = ylabel_match.group(1)
    x_axis = [d*k_c for d in list(map(float, x_axis_match.group(1).split(' ')))]
    y_axis = [d*k_c for d in list(map(float, y_axis_match.group(1).split(' ')))]
    
    # Extract the color map and its numeric values.
    lines = content.split('\n')
    
    # Locate the first grid-data row.
    data_start = -1
    color_map = {}
    
    for i, line in enumerate(lines):
        if line.startswith('"') and ' c #' in line:
            # Parse a color-map entry such as "AA c #000000 " /* "0" */.
            parts = line.split(' c #')
            if len(parts) >= 2:
                symbol = parts[0].strip('"')
                color_part = parts[1].split(' ')[0]
                # Extract the numeric value.
                value_match = re.search(r'/\* "([^"]*)" \*/', line)
                if value_match:
                    try:
                        value = float(value_match.group(1))
                        color_map[symbol] = {'color': '#' + color_part, 'value': value}
                    except ValueError:
                        pass
        elif line.startswith('"') and not ' c #' in line and len(color_map) > 0:
            data_start = i
            break
    
    # Extract grid rows.
    grid_data = []
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if line.startswith('"') and line.endswith('",') or line.endswith('"'):
            # Remove quotes and a trailing comma.
            clean_line = line.replace('"','').rstrip(',')
            if clean_line:
                # Split the row into two-character symbols if explicit decoding is needed.
                # row = [clean_line[j:j+2] for j in range(0, len(clean_line), 2)]
                grid_data.append(clean_line)

    return {
        'title': title,
        'legend': legend,
        'xlabel': xlabel,
        'ylabel': ylabel,
        'color_map': color_map,
        'grid_data': grid_data,
        'x_axis': np.linspace(min(x_axis), max(x_axis), len(grid_data)),
        'y_axis': np.linspace(min(y_axis), max(y_axis), len(grid_data)),
    }

def create_fel_from_data(grid_strings, color_map, gauss=0.5):
    """
    Create an FEL matrix from decoded XPM data.
    """
    # Convert grid symbols into a numeric matrix.
    energy_matrix = np.zeros((len(grid_strings), len(grid_strings)))
    
    for i, row_string in enumerate(grid_strings):
        for j in range(0, len(row_string), 2):
            symbol = row_string[j:j+2]
            if symbol in color_map:
                energy_matrix[i, j//2] = color_map[symbol]['value']
    
    # Shift the energy range so the maximum is zero and all other values are negative.
    max_energy = np.max(energy_matrix)
    energy_matrix = energy_matrix - max_energy
    
    # Apply Gaussian smoothing to reduce jagged edges.
    if gauss > 0:
        energy_matrix = gaussian_filter(energy_matrix, sigma=gauss)
    
    return energy_matrix


def plot_fel(
        energy_matrix, x_axis, y_axis, 
        title, x_label, y_label, bar_label,
        limit_axis=None, 
        save_path=None,
        need_grids=None,
        contour_line=None,
    ):
    """
    Plot a free-energy landscape (FEL).
    - Draw contour lines as solid lines.
    - Use a built-in continuous colormap.
    - Format colorbar ticks as integers.
    """
    # Create the figure.
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # A built-in continuous colormap such as viridis can also be used.
    # cmap = cm.Blues
    
    cmap = make_cmap(wcmap=whatcmap)

    # Draw the filled-contour plot.
    X, Y = np.meshgrid(x_axis, y_axis)
    contour_filled = ax.contourf(X, Y, energy_matrix, levels=100, cmap=cmap)
    
    # Add solid contour lines.
    if contour_line:
        contour_level = contour_line['level']
        contour_alpha = contour_line['alpha']
        contour_label = contour_line['label']
        contour_lines = ax.contour(
            X, Y, energy_matrix, 
            levels=contour_level, 
            alpha=contour_alpha, 
            colors='black', linewidths=1.0, linestyles="solid",
        )
        if contour_label:
            ax.clabel(contour_lines, inline=True, fontsize=12, fmt='%d')  # Use integer contour labels.
    
    # Set labels and title.
    ax.set_xlabel(x_label, fontsize=14, fontweight='bold')
    ax.set_ylabel(y_label, fontsize=14, fontweight='bold')
    ax.set_title(title, fontsize=16, fontweight='bold')
    
    # Add the colorbar.
    cbar = plt.colorbar(contour_filled, ax=ax, shrink=0.5)
    cbar.set_label(bar_label, fontsize=14, fontweight='bold')
    
    # Format the colorbar with no more than eight integer ticks.
    cbar.locator = MaxNLocator(integer=True, nbins=8)
    cbar.update_ticks()
    
    # Set axis ticks and optional limits.
    if limit_axis:
        ax.axis(limit_axis)

    ax.xaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
    ax.tick_params(axis='both', which='major', labelsize=14)
    
    # Add axis grid lines.
    if need_grids:
        ax.grid(True, alpha=0.3, linestyle='--')
    
    # Adjust the layout.
    plt.tight_layout()
    
    # Save the figure.
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"FEL saved to: {save_path}")
    
    return fig, ax

def plot_fel_iceberg(
    energy_matrix, x_axis, y_axis,
    title, x_label, y_label, bar_label,
    save_path=None,
    show_minima=None,
    seabed=dict(
        show=True,
        mode="contourf",   # "contourf" or "surface"
        levels=100,
        alpha=1.0,
        offset="min",      # "zero", "min", or a numeric value such as -5.0
        epsilon_ratio=1e-3 # Relative offset that prevents overlap with the surface.
    ),
    light_azimuth=315,
    light_altitude=45,
    z_exaggeration=1.0,
    zlim_margin=0.02,
):
    """
    Plot an inverted 3D iceberg view.
    - energy_matrix is expected to have a maximum of zero and all other values <= 0.
    - seabed.offset:
        "zero" -> place the plane at z=0
        "min"  -> place it just below the global minimum
        float  -> use the supplied z value directly
    """
    cmap = make_cmap(wcmap=whatcmap)
    X, Y = np.meshgrid(x_axis, y_axis)

    Z = energy_matrix * float(z_exaggeration)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Color normalization and lighting.
    norm = Normalize(vmin=np.min(energy_matrix), vmax=np.max(energy_matrix))
    ls = LightSource(azdeg=light_azimuth, altdeg=light_altitude)
    face_rgb = ls.shade(energy_matrix, cmap=cmap, vert_exag=0.5, blend_mode="soft")

    # Iceberg surface.
    ax.plot_surface(
        X, Y, Z,
        facecolors=face_rgb,
        rstride=1, cstride=1,
        linewidth=0, antialiased=True, shade=False
    )

    # Place the 2D FEL plane at the minimum or a user-specified offset.
    if seabed and seabed.get("show", True):
        mode   = seabed.get("mode", "contourf")
        levels = seabed.get("levels", 100)
        alpha  = seabed.get("alpha", 1.0)
        offset_cfg = seabed.get("offset", "min")
        eps_ratio  = float(seabed.get("epsilon_ratio", 1e-3))

        if isinstance(offset_cfg, (int, float)):
            offset = float(offset_cfg)
        else:
            offset_cfg = str(offset_cfg).lower()
            if offset_cfg == "min":
                zmin = float(Z.min())
                # Apply a relative offset to prevent z-fighting with the surface.
                epsilon = max(1e-9, abs(zmin) * eps_ratio)
                offset = zmin - epsilon
            elif offset_cfg == "zero":
                offset = 0.0
            else:
                raise ValueError(f"Unknown seabed offset: {offset_cfg}")

        if mode == "contourf":
            ax.contourf(
                X, Y, energy_matrix,
                zdir='z', offset=offset,
                levels=levels, cmap=cmap, alpha=alpha
            )
        elif mode == "surface":
            colors_rgba = cmap(norm(energy_matrix))
            Z0 = np.full_like(Z, fill_value=offset)
            ax.plot_surface(
                X, Y, Z0,
                facecolors=colors_rgba,
                rstride=1, cstride=1,
                linewidth=0, antialiased=False, shade=False, alpha=alpha
            )
        else:
            raise ValueError(f"Unknown seabed mode: {mode}")

    # Minimum-energy points.
    if show_minima:
        xs = [p["x"] for p in show_minima]
        ys = [p["y"] for p in show_minima]
        zs = [p["z"] * float(z_exaggeration) for p in show_minima]
        ax.scatter(xs, ys, zs, s=35, c='crimson', depthshade=True, label="Local minima")
        for k, (x, y, z) in enumerate(zip(xs, ys, zs), 1):
            ax.text(x, y, z, f"{k}", color='crimson', fontsize=9)

    # Axes and title.
    z_label = bar_label
    ax.set_xlabel(x_label, labelpad=8, fontweight='bold')
    ax.set_ylabel(y_label, labelpad=8, fontweight='bold')
    ax.set_zlabel(z_label, labelpad=8, fontweight='bold')
    ax.set_title(title, pad=12, fontweight='bold')

    # Fit the axes tightly without extra padding.
    ax.set_xlim(x_axis.min(), x_axis.max())
    ax.set_ylim(y_axis.min(), y_axis.max())
    zmin = Z.min()
    zmax = 0.0
    if seabed and seabed.get("show", True):
        zmax = max(zmax, offset)
        zmin = min(zmin, offset)
    ax.set_zlim(zmin - abs(zmin) * zlim_margin, zmax)

    ax.set_box_aspect([1,1,1])  # Preserve equal aspect ratios.
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')

    # Colorbar.
    mappable = cm.ScalarMappable(norm=norm, cmap=cmap)
    mappable.set_array([])
    cbar = fig.colorbar(mappable, shrink=0.6, pad=0.08)
    cbar.set_label(z_label, fontsize=12, fontweight='bold')
    cbar.locator = MaxNLocator(integer=True, nbins=8)
    cbar.update_ticks()

    ax.view_init(elev=28, azim=-60)
    ax.dist = 9
    ax.grid(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"3D FEL (iceberg) saved to: {save_path}")

    return fig, ax


def sort_fel_points(energy_matrix, x_axis, y_axis, top_n=None):
    """
    Sort all points in energy_matrix from lowest to highest energy.
    Return dictionaries containing matrix indices, coordinates, and energy values.
    
    Parameters:
    - energy_matrix: 2D numpy array, shape = (len(y_axis), len(x_axis))
    - x_axis, y_axis: corresponding coordinate arrays
    - top_n: optional number of lowest-energy points to return
    """
    E = np.array(energy_matrix, dtype=float)
    nrows, ncols = E.shape

    results = []
    for i in range(nrows):
        for j in range(ncols):
            results.append({
                'i': i,
                'j': j,
                'x': float(x_axis[j]),
                'y': float(y_axis[i]),
                'z': float(E[i, j]),
            })

    results.sort(key=lambda d: d['z'])  # Sort from lowest to highest energy.

    if top_n:
        results = results[:top_n]

    return results

def find_fel_minima(energy_matrix, x_axis, y_axis,
                    top_n=5,
                    exclude_edges=True,
                    nan_policy='omit'):
    """
    Find all strict local minima in a 2D FEL using the eight-neighbor rule.
    Return dictionaries containing matrix indices, coordinates, and energy values.
    Assumption: energy_matrix.shape = (len(y_axis), len(x_axis)).
    """
    E = np.array(energy_matrix, dtype=float, copy=True)

    # Handle NaN values.
    if nan_policy == 'omit':
        # Replace NaN with +inf so it cannot become a minimum.
        E = np.where(np.isnan(E), np.inf, E)
    elif nan_policy == 'propagate':
        pass  # Preserve NaN; comparisons prevent it from becoming a minimum.

    # Strictly compare each point with its eight neighbors.
    neighbors = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)]
    mask = np.ones_like(E, dtype=bool)
    for dy, dx in neighbors:
        shifted = np.roll(E, shift=(dy, dx), axis=(0,1))
        mask &= (E < shifted)

    # Remove false edge minima introduced by np.roll wrapping.
    if exclude_edges:
        mask[0, :] = False
        mask[-1, :] = False
        mask[:, 0] = False
        mask[:, -1] = False

    # Collect the indices of all minima.
    idxs = np.argwhere(mask)
    minima = []
    for (i, j) in idxs:
        minima.append({
            'i': int(i),
            'j': int(j),
            'x': float(x_axis[j]),
            'y': float(y_axis[i]),
            'z': float(E[i, j]),
        })
    # Sort minima from lowest to highest energy.
    minima.sort(key=lambda d: d['z'])
    return minima[:top_n]

###################################

def main():
    """
    Run the FEL generation workflow.
    """
    costom_file = xpm_path

    print("Free Energy Landscape (FEL) Generator")
    print("=====================================")
    
    # Generate the FEL from the data in the XPM file.
    try:
        print("Generating FEL from XPM data...")
        read_data = parse_xpm_file(costom_file)
        energy_matrix = create_fel_from_data(
            grid_strings=read_data['grid_data'],
            color_map=read_data['color_map'],
            gauss=gauss,
        )
        x_axis = read_data['x_axis']
        y_axis = read_data['y_axis']

        print('############ FEL_points ############')
        # value_list = sort_fel_points(energy_matrix, x_axis, y_axis, top_n=top_values)
        value_list = find_fel_minima(energy_matrix, x_axis, y_axis, top_n=top_values)
        for k, p in enumerate(value_list, 1):
            print(f"{k:2d}. i={p['i']:>4d}, j={p['j']:>4d}, "
                f"x={p['x']:.6g}, y={p['y']:.6g}, energy={p['z']:.6g}")
        print('####################################')

        fig, ax = plot_fel(
            energy_matrix, x_axis, y_axis, 
            title=read_data['title'],
            x_label=read_data['xlabel'],
            y_label=read_data['ylabel'],
            bar_label='Gibbs (kJ/mol)',
            save_path=out_fig_path,
            limit_axis=limit_axis,
            need_grids=need_grids,
            contour_line=contour_line,
        )
        print("✓ FEL generated successfully from XMP data")

        # Optionally generate a 3D iceberg plot.
        fig3d, ax3d = plot_fel_iceberg(
            energy_matrix, x_axis, y_axis,
            title=read_data['title'] + " (3D Iceberg)",
            x_label=read_data['xlabel'],
            y_label=read_data['ylabel'],
            bar_label='Gibbs (kJ/mol)',
            save_path=None,
            # show_minima=value_list,
            show_minima=None,
            seabed=dict(
                show=True,
                mode="contourf",
                levels=100,
                alpha=1.0,
                offset="min",       # Place the plane at the minimum (the default).
                epsilon_ratio=1e-3, # Increase to 1e-2 for more visual separation.
            ),
            light_azimuth=315,
            light_altitude=45,
            z_exaggeration=1.0,
        )

    except Exception as e:
        print(f"✗ Error generating FEL from XMP data: {e}")
    
    plt.show()
    
    print("\nProgram completed!")

if __name__ == "__main__":
    main()
