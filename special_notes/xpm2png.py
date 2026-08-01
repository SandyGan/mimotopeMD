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
# limit_axis = [0.50, 2.50, 17.90, 18.35]

# Show axis grid lines.
need_grids = False

# Gaussian smoothing; 0 disables smoothing.
gauss = 0.5

# Unit-conversion factor; use 10 for nm to angstrom or 0.1 for angstrom to nm.
k_c = 10

# Contour-line settings; use None to disable contour lines.
# contour_line = None
contour_line = {
    'level': 20,
    'alpha': 0.2,
    'label': False, # True labels the contour lines with their values.
}

# Colormap.
# whatcmap = 'custom'
whatcmap = 'GnBu_r'
# whatcmap = 'custom5'

# Energy convention:
#   'min_zero': minimum ΔG=0 and all other values positive (the native gmx sham convention)
#   'max_zero': maximum=0 and all other values negative (legacy script behavior)
energy_convention = 'min_zero'

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
        stops  = [0.0, 0.10, 0.30, 0.45, 0.60, 0.75, 0.92, 1.0]
        colors = ['#B6B3D6','#CFCCE3','#D5D3DE','#D5D1D1',
                  '#F6DFD6','#F8B2A2','#F1837A','#E9687A']
        return LinearSegmentedColormap.from_list("my_cmap", list(zip(stops, colors)), N=256)
    elif wcmap == 'custom4':
        stops  = [0.0, 0.20, 0.40, 0.60, 0.80, 1.0]
        colors = ['#61AACF','#98CADD','#EAEFF6',
                  '#F9EFEF','#E9C6C6','#DA9599']
        return LinearSegmentedColormap.from_list("my_cmap", list(zip(stops, colors)), N=256)
    elif wcmap == 'custom5':
        stops  = [0.0, 0.20, 0.35, 0.50, 0.75, 0.85, 0.90, 1.0]
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
    # [FIX 1] GROMACS may split long axis arrays across multiple comments.
    # A single re.search captured only the first line, truncating min/max and miscalibrating the axes.
    x_chunks = re.findall(r'/\*\s*x-axis:\s*(.*?)\s*\*/', content)
    y_chunks = re.findall(r'/\*\s*y-axis:\s*(.*?)\s*\*/', content)
    x_axis = [float(v)*k_c for ch in x_chunks for v in ch.split()]
    y_axis = [float(v)*k_c for ch in y_chunks for v in ch.split()]
    
    # [FIX 2] Read nx, ny, and chars-per-pixel from the XPM header instead of assuming a square grid and cpp=2.
    lines = content.split('\n')
    _h = re.search(r'"\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s*"', content)
    nx, ny, ncolors, cpp = (int(v) for v in _h.groups())
    
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
                # Split by chars-per-pixel later when decoding the matrix.
                # row = [clean_line[j:j+2] for j in range(0, len(clean_line), 2)]
                grid_data.append(clean_line)

    print(f"  [xpm] nx={nx} ny={ny} ncolors={ncolors} chars/pixel={cpp} | "
          f"x-axis lines={len(x_chunks)} ({len(x_axis)} vals), "
          f"y-axis lines={len(y_chunks)} ({len(y_axis)} vals)")
    assert len(grid_data) == ny, f'XPM row count {len(grid_data)} != ny {ny}'

    # [FIX 12] gmx sham writes nbin+1 grid boundaries rather than nbin bin centers.
    # The old np.linspace(min, max, nbin) approach changed the interval denominator from nbin
    # to nbin-1, stretching the coordinates and shifting them by roughly half a bin.
    # Convert boundaries to bin centers explicitly.
    def _to_centers(a, n, name):
        a = np.asarray(a, dtype=float)
        if len(a) == n + 1:                       # Boundaries -> centers
            print(f"         {name}: {len(a)} boundaries -> {n} bin centers "
                  f"(bin width = {np.mean(np.diff(a)):.6g})")
            return (a[:-1] + a[1:]) / 2.0
        if len(a) == n:                           # Already bin centers
            return a
        raise ValueError(f'{name} has {len(a)} values for a grid with {n} bins; '
                         f'expected {n} centers or {n+1} boundaries. Check the XPM file.')

    x_axis = _to_centers(x_axis, nx, 'x-axis')
    y_axis = _to_centers(y_axis, ny, 'y-axis')

    return {
        'nx': nx, 'ny': ny, 'cpp': cpp,
        'title': title,
        'legend': legend,
        'xlabel': xlabel,
        'ylabel': ylabel,
        'color_map': color_map,
        'grid_data': grid_data,
        # [FIX 3] The old code sized x_axis from len(grid_data)=ny, which misaligned non-square grids.
        'x_axis': x_axis,
        'y_axis': y_axis,
    }

def create_fel_from_data(grid_strings, color_map, gauss=0.5,
                         cpp=2, convention='min_zero'):
    """
    Create an FEL matrix from decoded XPM data.
    Matrix shape is (len(y_axis), len(x_axis)); row indices correspond to ascending y.
    """
    ny = len(grid_strings)
    nx = len(grid_strings[0]) // cpp
    energy_matrix = np.zeros((ny, nx))

    for i, row_string in enumerate(grid_strings):
        for j in range(0, len(row_string), cpp):
            symbol = row_string[j:j+cpp]
            if symbol in color_map:
                energy_matrix[i, j//cpp] = color_map[symbol]['value']

    # ==================== [FIX 4] Critical correction ====================
    # XPM is an image format: row 0 is at the top and corresponds to the maximum y value.
    # Pairing file rows directly with ascending y_axis mirrors the entire y axis:
    #     y_displayed = y_min + y_max - y_true
    # The resulting basin appears at a coordinate unsupported by the source data.
    energy_matrix = energy_matrix[::-1, :]
    
    # Shift the energy range so the maximum is zero and all other values are negative.
    max_energy = np.max(energy_matrix)
    energy_matrix = energy_matrix - max_energy
    # =========================================================

    if convention == 'min_zero':
        energy_matrix = energy_matrix - np.min(energy_matrix)   # Native gmx sham convention
    else:
        energy_matrix = energy_matrix - np.max(energy_matrix)   # Legacy script behavior

    # Apply Gaussian smoothing for rendering.
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
            colors='black', linewidths=0.3, linestyles="solid",
        )
        if contour_label:
            ax.clabel(contour_lines, inline=True, fontsize=14, fmt='%d')  # Use integer contour labels.
    
    # Set the width of all four axis spines.
    for spine in ax.spines.values():
        spine.set_linewidth(0.3)  # Line width in points.

    # Set labels and title.
    ax.set_xlabel(x_label, fontsize=16, fontweight='normal')
    ax.set_ylabel(y_label, fontsize=16, fontweight='normal')
    # ax.set_title(title, fontsize=16, fontweight='bold')
    # [FIX 6] The old script hard-coded ax.set_xlim(0,14) and ax.set_ylim(13.9,19.5),
    # making the displayed range an arbitrary plotting choice. Derive it from the data instead;
    # use limit_axis above and document it in the figure caption when fixed limits are required.
    ax.set_xlim(x_axis.min(), x_axis.max())
    ax.set_xlim(0, x_axis.max())
    ax.set_ylim(y_axis.min(), y_axis.max())
    
    # Add the colorbar.
    cbar = plt.colorbar(contour_filled, ax=ax, shrink=0.5)
    cbar.set_label(bar_label, fontsize=14, fontweight='normal')
    cbar.ax.tick_params(labelsize=14)   # Tick-label font size.
    cbar.outline.set_linewidth(0.3)     # Border width.
    cbar.outline.set_edgecolor('black') # Border color.
    # Format the colorbar with no more than eight integer ticks.
    cbar.locator = MaxNLocator(integer=True, nbins=8)
    cbar.update_ticks()
    
    # Set axis ticks and optional limits.
    if limit_axis:
        ax.axis(limit_axis)

    # ax.xaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))
    # ax.yaxis.set_major_locator(MaxNLocator(nbins=8, integer=True))

    ax.tick_params(axis='both', which='major', labelsize=16)
    # ax.tick_params(axis='x', labelsize=12)  # Change only the x-axis.
    # ax.tick_params(axis='y', labelsize=16)  # Change only the y-axis.
    
    # Add axis grid lines.
    if need_grids:
        ax.grid(True, alpha=0.3, linestyle='--')
    
    # Adjust the layout.
    plt.tight_layout()
    
    # Save the figure.
    if save_path:
        plt.savefig(save_path, dpi=600, bbox_inches='tight')
        print(f"FEL saved to: {save_path}")
    
    return fig, ax

def plot_fel_iceberg(
    energy_matrix, x_axis, y_axis,
    title, x_label, y_label, bar_label,
    save_path=None,
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
    cbar = fig.colorbar(mappable, ax=ax, shrink=0.6, pad=0.08)  # [FIX 10] New Matplotlib requires an explicit Axes.
    cbar.set_label(z_label, fontsize=12, fontweight='bold')
    cbar.locator = MaxNLocator(integer=True, nbins=8)
    cbar.update_ticks()

    ax.view_init(elev=28, azim=-60)
    try:
        ax.dist = 9      # [FIX 11] Removed in Matplotlib >=3.9.
    except Exception:
        pass
    ax.grid(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"3D FEL (iceberg) saved to: {save_path}")

    return fig, ax


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
            cpp=read_data['cpp'],
            convention=energy_convention,
        )
        x_axis = read_data['x_axis']
        y_axis = read_data['y_axis']

        fig, ax = plot_fel(
            energy_matrix, x_axis, y_axis, 
            title=read_data['title'],
            x_label=read_data['xlabel'],
            y_label=read_data['ylabel'],
            bar_label=('ΔG (kJ/mol)' if energy_convention == 'min_zero'
                       else '-ΔG (kJ/mol)'),   # [FIX 8] Keep the label consistent with the energy convention.
            save_path=out_fig_path,
            limit_axis=limit_axis,
            need_grids=need_grids,
            contour_line=contour_line,
        )
        print("✓ FEL generated successfully from XMP data")

        # Optionally generate a 3D iceberg plot.
        # [FIX 9] The inverted view requires max=0 and all other values negative; convert either convention here.
        fig3d, ax3d = plot_fel_iceberg(
            energy_matrix - np.max(energy_matrix), x_axis, y_axis,
            title=read_data['title'] + " (3D Iceberg)",
            x_label=read_data['xlabel'],
            y_label=read_data['ylabel'],
            bar_label='-ΔG (kJ/mol)',
            save_path=None,
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
