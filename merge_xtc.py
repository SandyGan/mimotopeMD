import subprocess
import os
import argparse

def find_center_atom(file='step5_1.gro', bottom=None):
    '''
    Return the atom index closest to the protein geometric center in a .gro file.
    '''
    print(f"Arguments: file={file}; bottom={bottom}, type={type(bottom)}")
    aminoacids = [
        'ALA','ASN','CYS','GLN','GLY','ILE','LEU','MET',
        'PHE','PRO','SER','THR','TRP','TYR','VAL','ASP',
        'ARG','HIS','LYS','GLU',
        'PYL','SEC',    # Pyrrolysine and selenocysteine
        'HID','HSD',    # Histidine protonated at delta nitrogen (ND1)
        'HIE','HSE',    # Histidine protonated at epsilon nitrogen (NE2)
        'HIP','HSP',    # Doubly protonated, positively charged histidine
        'ASH', # Protonated, neutral aspartic acid
        'GLH', # Protonated glutamic acid
        'CYX', # Cysteine in a disulfide bond (CYS-CYS)
        'CYN', # Deprotonated, negatively charged cysteine
        'LYN', # Deprotonated lysine
        'TYM', # Deprotonated, negatively charged tyrosine
        'ARN', # Deprotonated arginine
        ]
    # Initialize the geometric-center accumulators.
    center_X, center_Y, center_Z = 0, 0, 0 
    i = 0
    position = []
    with open(file, "r", encoding='utf-8') as fo:
        try:
            next(fo)  # Skip the title line.
            next(fo)  # Skip the atom-count line.
        except StopIteration:
            pass  # Continue when the file has fewer than two lines.
        if bottom:
            for line in fo:
                if line[5:8] == bottom:
                    break
                else:
                    # Protein coordinates occupy the fixed-width x, y, and z fields.
                    x, y, z, atom = float(line[20:28]), float(line[28:36]), float(line[36:44]), int(line[15:20])
                    center_X += x
                    center_Y += y
                    center_Z += z
                    i += 1
                    position.append([atom, x, y, z])
        else:
            for line in fo:
                if line[5:8] in aminoacids:
                    # Protein coordinates occupy the fixed-width x, y, and z fields.
                    x, y, z, atom = float(line[20:28]), float(line[28:36]), float(line[36:44]), int(line[15:20])
                    center_X += x
                    center_Y += y
                    center_Z += z
                    i += 1
                    position.append([atom, x, y, z])
                else:
                    break            
    # Calculate the geometric-center coordinates.
    center_X = center_X / i
    center_Y = center_Y / i
    center_Z = center_Z / i
    print(f'Geometric center: {center_X, center_Y, center_Z}')
    dist = 10
    atom_info = ''
    # Search the atoms again for the one closest to the geometric center.
    for entry in position:
        # Update the nearest atom whenever a shorter distance is found.
        atom_center_dist = ((entry[1]-center_X)**2 + (entry[2]-center_Y)**2 + (entry[3]-center_Z)**2 )**0.5
        if atom_center_dist < dist:
            dist = atom_center_dist
            # Record the atom index.
            atom_info = entry[0]
    print(f'Atom index closest to the geometric center: {atom_info}')
    return atom_info


def add_first_line(file_path, new_line):
    '''
    Prepend content to a file.
    '''
    with open(file_path, 'r+', encoding='utf-8') as f:
        content = f.read()
        if new_line in content:
            print(f'[Read {file_path}] The [center] group is already up to date')
        else:
            f.seek(0)  # Move the file pointer to the beginning.
            f.write(new_line + '\n' + content)
            print(f'[Read {file_path}] Updated the [center] group')
        f.close()


def gmx_trjconv_trans(start, end, prod_step, number):
    '''
    Center each trajectory geometrically with gmx trjconv.
    '''
    print('[gmx_trjconv_trans] Centering XTC trajectories')
    for i in range(start, end + 1):
        istee = f"{prod_step}_{i}"
        isteq = f"{istee}_{number}_center"        
        # Remove an existing output file, with exception handling.
        try:
            for file in [f"{isteq}.xtc"]:
                if os.path.exists(file):
                    os.remove(file)
                    print(f'> Removed existing file: {file}')
        except Exception as e:
            print(f"Failed to remove file: {e}")
            continue  # Continue with the next trajectory segment.
        # Run commands sequentially.
        try:
            # subprocess.run keeps execution synchronous and explicit.
            result = subprocess.run(
                [
                    "gmx", "trjconv",
                    "-s", f"{istee}.tpr",
                    "-f", f"{istee}.xtc",
                    "-n", "index.ndx",
                    "-o", f"{isteq}.xtc",
                    "-pbc", "atom", "-center",
                ],
                input="0\n3\n".encode(),  # Send both interactive selections at once.
                stdout=subprocess.PIPE,     # Capture standard output.
                stderr=subprocess.PIPE,     # Capture standard error.
                check=True,                 # Raise an exception on a nonzero exit code.
            )
            print(result.stdout.decode())
            print(f'[gmx_trjconv] Wrote {isteq}.xtc\n'+'-'*10)
        
        except subprocess.CalledProcessError as e:
            print(f"Command failed with exit code {e.returncode}")
            print("Error output:", e.stderr.decode())
            continue  # Continue with the next segment; use raise to stop instead.
        except subprocess.TimeoutExpired:
            print("Command timed out")
            continue
        except Exception as e:
            print(f"Unexpected error: {e}")
            continue


def gmx_trjconv_noline(start, end, prod_step, number):
    '''
    Reconstruct molecules across PBC to prevent line artifacts in VMD.
    '''
    print('[gmx_trjconv_noline] Removing PBC line artifacts for visualization')
    for i in range(start, end + 1):
        istee = f"{prod_step}_{i}"
        isteq = f"{istee}_{number}_center"
        istep = f"{isteq}_noline"
        # Remove an existing output file, with exception handling.
        try:
            for file in [f"{istep}.xtc"]:
                if os.path.exists(file):
                    os.remove(file)
                    print(f'> Removed existing file: {file}')
        except Exception as e:
            print(f"Failed to remove file: {e}")
            continue  # Continue with the next trajectory segment.
        try:
            result2 = subprocess.run(
            [
                "gmx", "trjconv", 
                "-f", f"{isteq}.xtc",
                "-o", f"{istep}.xtc",
                "-s", f"{istee}.tpr",
                "-pbc", "mol",
            ],
            input="0\n".encode(),  # Send the interactive selection at once.
            stdout=subprocess.PIPE,     # Capture standard output.
            stderr=subprocess.PIPE,     # Capture standard error.
            check=True,                 # Raise an exception on a nonzero exit code.
            )
            print(result2.stdout.decode())
            print(f'[gmx_trjconv] Wrote {istep}.xtc\n'+'-'*10)
        
        except subprocess.CalledProcessError as e:
            print(f"Command failed with exit code {e.returncode}")
            print("Error output:", e.stderr.decode())
            continue  # Continue with the next segment; use raise to stop instead.
        except subprocess.TimeoutExpired:
            print("Command timed out")
            continue
        except Exception as e:
            print(f"Unexpected error: {e}")
            continue

def gmx_trjcat(start, end, prod_step, number):
    '''
    Concatenate the processed XTC trajectory segments.
    '''
    out_xtc = f"{prod_step}_all.xtc"
    if os.path.exists(out_xtc):
        os.remove(out_xtc)
        print(f'> Removed existing file: {out_xtc}')
    entf = [f"{prod_step}_{j}_{number}_center_noline.xtc" 
            for j in range(start, end + 1)]
    trajectory_list = '\n - '.join(entf)
    print(f'[gmx_trjcat] Merging \n - {trajectory_list}')
    try:
        result3 = subprocess.run(
        ["gmx", "trjcat", "-f",] + entf + 
        ["-o", out_xtc, "-sort", "yes", "-cat", "-settime"],
        input=str("c\n"*int(1+end-start)).encode(),  # Send all time selections at once.
        stdout=subprocess.PIPE,     # Capture standard output.
        stderr=subprocess.PIPE,     # Capture standard error.
        check=True,
        )
        print(result3.stdout.decode())
        print(f'[gmx_trjcat] Wrote merged trajectory: {out_xtc}')
    
    except subprocess.CalledProcessError as e:
        print(f"Command failed with exit code {e.returncode}")
        print("Error output:", e.stderr.decode())
    except subprocess.TimeoutExpired:
        print("Command timed out")
    except Exception as e:
        print(f"Unexpected error: {e}")

gro_file = 'step5_1.gro'
step_start = 1
step_end = 5
xtc_step = 'step5'


if __name__ == '__main__':
    # Create the command-line argument parser.
    parser = argparse.ArgumentParser(description="Center, reconstruct, and merge GROMACS XTC files")
    # Add command-line options and their defaults.
    parser.add_argument('--gro', type=str, default=gro_file, help=f".gro file (default: {gro_file})")
    parser.add_argument('--gro_until', type=str, default=None, help="Stop reading the GRO file at this residue name (default: read recognized protein residues)")
    parser.add_argument('--start', type=int, default=step_start, help=f"Start step5_*.xtc index (default: {step_start})")
    parser.add_argument('--end', type=int, default=step_end, help=f"End step5_*.xtc index (default: {step_end})")
    parser.add_argument('--prod_step', type=str, default=xtc_step, help=f"Production step (default: {xtc_step})")

    # Parse command-line arguments.
    args = parser.parse_args()
    
    # Find the atom index closest to the geometric center.
    center_id = find_center_atom(file=args.gro, bottom=args.gro_until)
    # Prepend the new group to index.ndx.
    add_first_line("index.ndx", f'[ center ]\n{str(center_id)}\n')
    # Preprocess and concatenate the trajectory files.
    gmx_trjconv_trans(start=args.start, end=args.end, prod_step=args.prod_step, number=center_id)
    gmx_trjconv_noline(start=args.start, end=args.end, prod_step=args.prod_step, number=center_id)
    gmx_trjcat(start=args.start, end=args.end, prod_step=args.prod_step, number=center_id)
