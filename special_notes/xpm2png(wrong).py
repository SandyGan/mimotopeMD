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
# from mpl_toolkits.mplot3d import Axes3D  # Enables 3D support even though the name is not used directly.

############################ User-configurable parameters ############################

# Path to the XPM file.
xpm_path = 'gibbs.xpm'

# Output image path; set to None to disable saving.
out_fig_path = './fel_from_xmp.svg'

# Axis limits: ax.axis([xmin, xmax, ymin, ymax]) or None.
limit_axis = None
# limit_axis = [0.70, 2.10, 8.55, 9.25]

# Display grid lines on the axes.
need_grids = False

# Gaussian smoothing; 0 disables smoothing.
gauss = 0.0

# Unit conversion factor, e.g. multiply nm by 10 to obtain angstroms.
k_c = 10

# Contour-line settings; set to None to disable contour lines.
# contour_line = None
contour_line = {
    'level': 15,
    'alpha': 0.2,
    'label': False, # Set to True to label the contour lines.
}

# Number of lowest-energy points to report.
top_values = 10

# Colormap selection.
# whatcmap = 'custom'
whatcmap = 'GnBu_r'
# whatcmap = 'custom5'

##################################################################################

def make_cmap(wcmap='viridis'):
    if wcmap == 'custom1':
    # Nonuniform color anchors matching the 2D palette.
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
    Parse an XPM file and extract the free-energy data.
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

    # Extract the color mapping and associated values.
    lines = content.split('\n')

    # Locate the start of the grid data.
    data_start = -1
    color_map = {}

    for i, line in enumerate(lines):
        if line.startswith('"') and ' c #' in line:
            # Parse a color entry such as "AA c #000000 " /* "0" */.
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
    
    # Extract the grid rows.
    grid_data = []
    for i in range(data_start, len(lines)):
        line = lines[i].strip()
        if line.startswith('"') and line.endswith('",') or line.endswith('"'):
            # Remove quotation marks and the trailing comma.
            clean_line = line.replace('"','').rstrip(',')
            if clean_line:
                # Split the row into two-character symbols if individual symbols are needed.
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
    Create an FEL matrix from the supplied XPM grid data.
    """
    # Convert the grid data to a numeric matrix.
    energy_matrix = np.zeros((len(grid_strings), len(grid_strings)))

    for i, row_string in enumerate(grid_strings):
        for j in range(0, len(row_string), 2):
            symbol = row_string[j:j+2]
            if symbol in color_map:
                energy_matrix[i, j//2] = color_map[symbol]['value']

    # Shift the energy range so that the maximum is zero and other values are negative.
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
        contour_line=None | dict,
    ):
    """
    Plot a free-energy landscape (FEL).
    - Draw solid contour lines.
    - Use a built-in continuous colormap.
    - Use integer ticks on the colorbar.
    """
    # Create the figure.
    fig, ax = plt.subplots(figsize=(5, 4))
    
    # A built-in continuous colormap can be selected here.
    # cmap = cm.Blues
    
    cmap = make_cmap(wcmap=whatcmap)

    # Draw the filled contour plot.
    X, Y = np.meshgrid(x_axis, y_axis)
    contour_filled = ax.contourf(X, Y, energy_matrix, levels=150, cmap=cmap)
    
    # Add solid contour lines.
    if contour_line:
        contour_level = contour_line['level']
        contour_alpha = contour_line['alpha']
        contour_label = contour_line['label']
        contour_lines = ax.contour(
            X, Y, energy_matrix, 
            levels=contour_level, 
            alpha=contour_alpha, 
            colors='black', linewidths=0.3, linestyles="solid",
        )
        if contour_label:
            ax.clabel(contour_lines, inline=True, fontsize=6, fmt='%d')  # Format contour labels as integers.
    
    # Set the width of all four axes spines.
    for spine in ax.spines.values():
        spine.set_linewidth(0.3)  # Spine width in points.

    # Set axis labels and the title.
    ax.set_xlabel(x_label, fontsize=12, fontweight='normal')
    ax.set_ylabel(y_label, fontsize=12, fontweight='normal')
    # ax.set_title(title, fontsize=16, fontweight='bold')
    # ax.set_xlim(0, 14.9)
    # ax.set_ylim(13.7, 19.5)
    
    # Add the colorbar.
    cbar = plt.colorbar(contour_filled, ax=ax, shrink=0.5)
    cbar.set_label(bar_label, fontsize=12, fontweight='normal')
    cbar.ax.tick_params(labelsize=8)   # Tick-label font size.
    cbar.outline.set_linewidth(0.3)       # Border width.
    cbar.outline.set_edgecolor('black') # Border color.
    # Use integer ticks on the colorbar.
    cbar.locator = MaxNLocator(integer=True, nbins=8)  # Use at most eight integer ticks.
    cbar.update_ticks()
    
    # Apply the requested axis limits.
    if limit_axis:
        ax.axis(limit_axis)

    # ax.xaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
    # ax.yaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))

    ax.tick_params(axis='both', which='major', labelsize=10)
    # ax.tick_params(axis='x', labelsize=12)  # Change only the x-axis labels.
    # ax.tick_params(axis='y', labelsize=16)  # Change only the y-axis labels.
    
    # Add axis grid lines.
    if need_grids:
        ax.grid(True, alpha=0.3, linestyle='--')
    
    # Adjust the layout.
    plt.tight_layout()
    
    # Save the image.
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        print(f"FEL saved to: {save_path}")
    
    return fig, ax

def plot_fel_iceberg(
    energy_matrix, x_axis, y_axis,
    title, x_label, y_label, bar_label,
    save_path=None,
    show_minima=None,
    seabed=dict(
        show=True,
        mode="contourf",   # Either "contourf" or "surface".
        levels=100,
        alpha=1.0,
        offset="min",      # Accepts "zero", "min", or a numeric value such as -5.0.
        epsilon_ratio=1e-3 # Relative offset that prevents overlap with the surface.
    ),
    light_azimuth=315,
    light_altitude=45,
    z_exaggeration=1.0,
    zlim_margin=0.02,
):
    """
    Create an inverted 3D iceberg view.
    - energy_matrix is expected to have a maximum of zero and other values at or below zero.
    - seabed.offset:
        "zero" -> place the plane at z=0.
        "min"  -> place the plane just below the global minimum.
        float  -> use the specified z value directly.
    """
    cmap = make_cmap(wcmap=whatcmap)
    X, Y = np.meshgrid(x_axis, y_axis)

    Z = energy_matrix * float(z_exaggeration)

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Normalize colors and configure lighting.
    norm = Normalize(vmin=np.min(energy_matrix), vmax=np.max(energy_matrix))
    ls = LightSource(azdeg=light_azimuth, altdeg=light_altitude)
    face_rgb = ls.shade(energy_matrix, cmap=cmap, vert_exag=0.5, blend_mode="soft")

    # Draw the iceberg surface.
    ax.plot_surface(
        X, Y, Z,
        facecolors=face_rgb,
        rstride=1, cstride=1,
        linewidth=0, antialiased=True, shade=False
    )

    # Place the 2D FEL plane at the minimum or at a user-defined offset.
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

    # Plot local minima.
    if show_minima:
        xs = [p["x"] for p in show_minima]
        ys = [p["y"] for p in show_minima]
        zs = [p["z"] * float(z_exaggeration) for p in show_minima]
        ax.scatter(xs, ys, zs, s=35, c='crimson', depthshade=True, label="Local minima")
        for k, (x, y, z) in enumerate(zip(xs, ys, zs), 1):
            ax.text(x, y, z, f"{k}", color='crimson', fontsize=9)

    # Configure axes and the title.
    z_label = bar_label
    ax.set_xlabel(x_label, labelpad=8, fontweight='bold')
    ax.set_ylabel(y_label, labelpad=8, fontweight='bold')
    ax.set_zlabel(z_label, labelpad=8, fontweight='bold')
    ax.set_title(title, pad=12, fontweight='bold')

    # Fit the axes tightly to the data range.
    ax.set_xlim(x_axis.min(), x_axis.max())
    ax.set_ylim(y_axis.min(), y_axis.max())
    zmin = Z.min()
    zmax = 0.0
    if seabed and seabed.get("show", True):
        zmax = max(zmax, offset)
        zmin = min(zmin, offset)
    ax.set_zlim(zmin - abs(zmin) * zlim_margin, zmax)

    ax.set_box_aspect([1,1,1])  # Keep equal aspect ratios.
    ax.xaxis.pane.set_edgecolor('black')
    ax.yaxis.pane.set_edgecolor('black')
    ax.zaxis.pane.set_edgecolor('black')

    # Add the colorbar.
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
    Sort every point in energy_matrix from lowest to highest energy.
    Return dictionaries containing row index, column index, x, y, and energy values.

    Parameters:
    - energy_matrix: 2D NumPy array with shape (len(y_axis), len(x_axis)).
    - x_axis, y_axis: Coordinate arrays.
    - top_n: Optional number of lowest-energy points to return; otherwise return all points.
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
    Find local minima that are strictly below all eight neighbors in a 2D FEL matrix.
    Return dictionaries containing row index, column index, x, y, and energy values.
    energy_matrix.shape is expected to equal (len(y_axis), len(x_axis)).
    """
    E = np.array(energy_matrix, dtype=float, copy=True)

    # Handle NaN values.
    if nan_policy == 'omit':
        # Replace NaN values with positive infinity so they cannot be minima.
        E = np.where(np.isnan(E), np.inf, E)
    elif nan_policy == 'propagate':
        pass  # Preserve NaN values; no minimum can be located at a NaN position.

    # Compare strictly against all eight neighbors.
    neighbors = [(1,0),(-1,0),(0,1),(0,-1),(1,1),(1,-1),(-1,1),(-1,-1)]
    mask = np.ones_like(E, dtype=bool)
    for dy, dx in neighbors:
        shifted = np.roll(E, shift=(dy, dx), axis=(0,1))
        mask &= (E < shifted)

    # Exclude false edge minima introduced by array wrapping.
    if exclude_edges:
        mask[0, :] = False
        mask[-1, :] = False
        mask[:, 0] = False
        mask[:, -1] = False

    # Collect all local-minimum indices.
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
    # Sort from lowest to highest energy for easier inspection.
    minima.sort(key=lambda d: d['z'])
    return minima[:top_n]

###################################

def main():
    """
    Run the FEL conversion and plotting workflow.
    """
    costom_file = xpm_path

    print("Free Energy Landscape (FEL) Generator")
    print("=====================================")
    
    # Generate the FEL from the data stored in the XPM file.
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
            bar_label='-ΔG (kJ/mol)',
            save_path=out_fig_path,
            limit_axis=limit_axis,
            need_grids=need_grids,
            contour_line=contour_line,
        )
        print("✓ FEL generated successfully from XPM data")

        # Optional: generate a 3D iceberg plot.
        fig3d, ax3d = plot_fel_iceberg(
            energy_matrix, x_axis, y_axis,
            title=read_data['title'] + " (3D Iceberg)",
            x_label=read_data['xlabel'],
            y_label=read_data['ylabel'],
            bar_label='-ΔG (kJ/mol)',
            save_path=None,
            # show_minima=value_list,
            show_minima=None,
            seabed=dict(
                show=True,
                mode="contourf",
                levels=100,
                alpha=1.0,
                offset="min",       # Place the plane at the minimum; this is the default.
                epsilon_ratio=1e-3, # Increase to 1e-2 for greater separation.
            ),
            light_azimuth=315,
            light_altitude=45,
            z_exaggeration=1.0,
        )

    except Exception as e:
        print(f"✗ Error generating FEL from XPM data: {e}")
    
    plt.show()
    
    print("\nProgram completed!")

if __name__ == "__main__":
    main()
